//! BAM records as the bytes on disk.
//!
//! `assign` needs a handful of fixed-position fields, the CIGAR and one or two
//! tags per record, and with `--updated_sam` it copies most records to the
//! output unchanged. Decoding every record into an owned structure and
//! re-encoding it costs more than everything else the program does, so
//! records stay as raw BAM bytes: fields are read in place, unchanged records
//! are written back verbatim, and the few edits Telescope makes (flag, MAPQ,
//! appended tags) are applied to the bytes directly.
//!
//! SAM input is converted to the same encoding as it is read.

use std::fs::File;
use std::io::{self, BufReader, Read, Write};

use noodles::{
    bam, bgzf,
    sam::{self, Header, alignment::RecordBuf, alignment::io::Write as AlignmentWrite},
};

pub use bgzf::io::writer::CompressionLevel;

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

fn u16_at(b: &[u8], i: usize) -> u16 {
    u16::from_le_bytes([b[i], b[i + 1]])
}

fn u32_at(b: &[u8], i: usize) -> u32 {
    u32::from_le_bytes([b[i], b[i + 1], b[i + 2], b[i + 3]])
}

fn i32_at(b: &[u8], i: usize) -> i32 {
    u32_at(b, i) as i32
}

/// A view of one BAM record (the bytes after the 4-byte length prefix).
#[derive(Clone, Copy)]
pub struct Raw<'a>(pub &'a [u8]);

/// One auxiliary field located in a record.
struct Aux {
    /// start of the whole field (its two-letter tag)
    start: usize,
    kind: u8,
    /// the value's bytes
    value: std::ops::Range<usize>,
    /// one past the whole field
    end: usize,
}

impl<'a> Raw<'a> {
    /// Rejects records too short to hold the fields their own header declares.
    pub fn check(self) -> io::Result<()> {
        let b = self.0;
        if b.len() < 32 || self.aux_start() > b.len() {
            return Err(invalid("truncated BAM record"));
        }
        Ok(())
    }

    pub fn ref_id(self) -> i32 {
        i32_at(self.0, 0)
    }
    /// 0-based leftmost position, -1 if unset (pysam's `reference_start`).
    pub fn pos(self) -> i32 {
        i32_at(self.0, 4)
    }
    fn name_len(self) -> usize {
        self.0[8] as usize
    }
    fn n_cigar(self) -> usize {
        u16_at(self.0, 12) as usize
    }
    pub fn flag(self) -> u16 {
        u16_at(self.0, 14)
    }
    fn seq_len(self) -> usize {
        u32_at(self.0, 16) as usize
    }
    pub fn mate_ref_id(self) -> i32 {
        i32_at(self.0, 20)
    }
    pub fn mate_pos(self) -> i32 {
        i32_at(self.0, 24)
    }
    pub fn template_len(self) -> i32 {
        i32_at(self.0, 28)
    }
    /// Read name without its terminating NUL.
    pub fn name(self) -> &'a [u8] {
        &self.0[32..32 + self.name_len().saturating_sub(1)]
    }
    fn cigar_start(self) -> usize {
        32 + self.name_len()
    }
    fn aux_start(self) -> usize {
        self.cigar_start() + 4 * self.n_cigar() + self.seq_len().div_ceil(2) + self.seq_len()
    }

    /// Calls `f(op, len)` for each CIGAR operation (op codes as in the BAM
    /// spec: 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X). A CIGAR too long
    /// for the fixed field lives in the `CG` tag behind a placeholder, which
    /// htslib resolves transparently; so does this.
    pub fn for_each_cigar_op(self, mut f: impl FnMut(u8, u32)) {
        let b = self.0;
        let c = self.cigar_start();
        let op_at = |i: usize| u32_at(b, c + 4 * i);
        if self.n_cigar() == 2
            && op_at(0) & 0xf == 4
            && (op_at(0) >> 4) as usize == self.seq_len()
            && op_at(1) & 0xf == 3
            && let Some(a) = self.find_aux(*b"CG")
            && a.kind == b'B'
            && b[a.value.start] == b'I'
        {
            let n = u32_at(b, a.value.start + 1) as usize;
            for i in 0..n {
                let v = u32_at(b, a.value.start + 5 + 4 * i);
                f((v & 0xf) as u8, v >> 4);
            }
            return;
        }
        for i in 0..self.n_cigar() {
            let v = op_at(i);
            f((v & 0xf) as u8, v >> 4);
        }
    }

    fn find_aux(self, tag: [u8; 2]) -> Option<Aux> {
        let b = self.0;
        let mut i = self.aux_start();
        while i + 3 <= b.len() {
            let kind = b[i + 2];
            let v = i + 3;
            let len = match kind {
                b'A' | b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                b'Z' | b'H' => b[v..].iter().position(|&x| x == 0)? + 1,
                b'B' => {
                    let size = match *b.get(v)? {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        _ => 4,
                    };
                    5 + size * u32_at(b, v + 1) as usize
                }
                _ => return None,
            };
            if v + len > b.len() {
                return None;
            }
            if b[i] == tag[0] && b[i + 1] == tag[1] {
                return Some(Aux { start: i, kind, value: v..v + len, end: v + len });
            }
            i = v + len;
        }
        None
    }

    /// An integer tag's value, whatever width it was stored in.
    pub fn aux_int(self, tag: [u8; 2]) -> Option<i64> {
        let a = self.find_aux(tag)?;
        let v = &self.0[a.value];
        Some(match a.kind {
            b'c' => v[0] as i8 as i64,
            b'C' => v[0] as i64,
            b's' => u16_at(v, 0) as i16 as i64,
            b'S' => u16_at(v, 0) as i64,
            b'i' => i32_at(v, 0) as i64,
            b'I' => u32_at(v, 0) as i64,
            _ => return None,
        })
    }

    /// A string tag's value without its terminating NUL.
    pub fn aux_str(self, tag: [u8; 2]) -> Option<&'a [u8]> {
        let a = self.find_aux(tag)?;
        (a.kind == b'Z').then(|| &self.0[a.value.start..a.value.end - 1])
    }
}

/// pysam's `set_tag(..., replace=True)`: drop any existing field, append.
fn replace_aux(rec: &mut Vec<u8>, tag: [u8; 2], kind: u8, value: &[u8]) {
    if let Some(a) = Raw(rec).find_aux(tag) {
        rec.drain(a.start..a.end);
    }
    rec.extend_from_slice(&tag);
    rec.push(kind);
    rec.extend_from_slice(value);
}

pub fn set_tag_str(rec: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) {
    replace_aux(rec, tag, b'Z', value);
    rec.push(0);
}

pub fn set_tag_u8(rec: &mut Vec<u8>, tag: [u8; 2], value: u8) {
    replace_aux(rec, tag, b'C', &[value]);
}

pub fn set_flag(rec: &mut [u8], flag: u16) {
    rec[14..16].copy_from_slice(&flag.to_le_bytes());
}

pub fn set_mapq(rec: &mut [u8], mapq: u8) {
    rec[9] = mapq;
}

enum Source {
    Bam(bam::io::Reader<bgzf::io::Reader<File>>),
    /// SAM records are re-encoded as BAM through an in-memory writer.
    Sam {
        reader: sam::io::Reader<BufReader<File>>,
        record: RecordBuf,
        encoder: bam::io::Writer<Vec<u8>>,
    },
}

/// Reads SAM or BAM (picked by content, as pysam does) one raw record at a time.
pub struct RawReader {
    pub header: Header,
    source: Source,
}

impl RawReader {
    pub fn open(path: &str) -> io::Result<Self> {
        let mut magic = [0u8; 2];
        let n = File::open(path)?.read(&mut magic)?;
        if n == 2 && magic == [0x1f, 0x8b] {
            let mut reader = bam::io::Reader::new(File::open(path)?);
            let header = reader.read_header()?;
            Ok(RawReader { header, source: Source::Bam(reader) })
        } else {
            let mut reader = sam::io::Reader::new(BufReader::new(File::open(path)?));
            let header = reader.read_header()?;
            let source =
                Source::Sam { reader, record: RecordBuf::default(), encoder: bam::io::Writer::from(Vec::new()) };
            Ok(RawReader { header, source })
        }
    }

    /// Fills `buf` with the next record; `false` at end of input.
    pub fn read(&mut self, buf: &mut Vec<u8>) -> io::Result<bool> {
        match &mut self.source {
            Source::Bam(reader) => {
                let inner = reader.get_mut();
                let mut len = [0u8; 4];
                // A clean end of file falls between records.
                let mut got = 0;
                while got < 4 {
                    match inner.read(&mut len[got..])? {
                        0 if got == 0 => return Ok(false),
                        0 => return Err(invalid("truncated BAM record length")),
                        n => got += n,
                    }
                }
                buf.resize(u32::from_le_bytes(len) as usize, 0);
                inner.read_exact(buf)?;
            }
            Source::Sam { reader, record, encoder } => {
                if reader.read_record_buf(&self.header, record)? == 0 {
                    return Ok(false);
                }
                encoder.get_mut().clear();
                encoder.write_alignment_record(&self.header, record)?;
                buf.clear();
                buf.extend_from_slice(&encoder.get_ref()[4..]);
            }
        }
        Raw(buf).check()?;
        Ok(true)
    }
}

/// A BAM being written with the input's header.
pub struct RawWriter {
    writer: bam::io::Writer<bgzf::io::Writer<File>>,
}

impl RawWriter {
    pub fn create(path: &str, header: &Header, level: CompressionLevel) -> io::Result<Self> {
        let inner =
            bgzf::io::writer::Builder::default().set_compression_level(level).build_from_writer(File::create(path)?);
        let mut writer = bam::io::Writer::from(inner);
        writer.write_header(header)?;
        Ok(RawWriter { writer })
    }

    pub fn write(&mut self, rec: &[u8]) -> io::Result<()> {
        let out = self.writer.get_mut();
        out.write_all(&(rec.len() as u32).to_le_bytes())?;
        out.write_all(rec)
    }

    pub fn finish(mut self) -> io::Result<()> {
        self.writer.try_finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// name "r1", 2 CIGAR ops (5M 3S), 8 bases, tags AS:i:-7 (as int8) and XS:Z:ab
    fn sample() -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(&3i32.to_le_bytes()); // ref
        b.extend_from_slice(&99i32.to_le_bytes()); // pos
        b.push(3); // name length incl. NUL
        b.push(42); // mapq
        b.extend_from_slice(&0u16.to_le_bytes()); // bin
        b.extend_from_slice(&2u16.to_le_bytes()); // n_cigar
        b.extend_from_slice(&0x53u16.to_le_bytes()); // flag
        b.extend_from_slice(&8u32.to_le_bytes()); // l_seq
        b.extend_from_slice(&(-1i32).to_le_bytes()); // mate ref
        b.extend_from_slice(&(-1i32).to_le_bytes()); // mate pos
        b.extend_from_slice(&(-120i32).to_le_bytes()); // tlen
        b.extend_from_slice(b"r1\0");
        b.extend_from_slice(&(5u32 << 4).to_le_bytes());
        b.extend_from_slice(&((3u32 << 4) | 4).to_le_bytes());
        b.extend_from_slice(&[0x12, 0x48, 0x12, 0x48]); // 8 bases
        b.extend_from_slice(&[30; 8]); // quals
        b.extend_from_slice(b"ASc");
        b.push((-7i8) as u8);
        b.extend_from_slice(b"XSZab\0");
        b
    }

    #[test]
    fn reads_fixed_fields_cigar_and_tags() {
        let b = sample();
        let r = Raw(&b);
        r.check().unwrap();
        assert_eq!((r.ref_id(), r.pos(), r.flag(), r.template_len()), (3, 99, 0x53, -120));
        assert_eq!(r.name(), b"r1");
        let mut ops = Vec::new();
        r.for_each_cigar_op(|op, len| ops.push((op, len)));
        assert_eq!(ops, vec![(0, 5), (4, 3)]);
        assert_eq!(r.aux_int(*b"AS"), Some(-7));
        assert_eq!(r.aux_str(*b"XS"), Some(&b"ab"[..]));
        assert_eq!(r.aux_int(*b"NM"), None);
    }

    #[test]
    fn set_tag_replaces_then_appends() {
        let mut b = sample();
        set_tag_str(&mut b, *b"XS", b"xyz");
        set_tag_u8(&mut b, *b"XP", 87);
        set_flag(&mut b, 0x153);
        set_mapq(&mut b, 0);
        let r = Raw(&b);
        assert_eq!(r.aux_str(*b"XS"), Some(&b"xyz"[..]));
        assert_eq!(r.aux_int(*b"XP"), Some(87));
        assert_eq!(r.aux_int(*b"AS"), Some(-7));
        assert_eq!((r.flag(), b[9]), (0x153, 0));
        // the replaced tag moved to the end, ahead of the newer one
        assert!(b.ends_with(b"XSZxyz\0XPC\x57"));
    }
}
