"""Build the small alignment/annotation fixtures under tests/data/legacy.

This is a record of how the fixtures were made; the tests do not run it.
It needs the local testbed (alignments of public CCLE and ENCODE reads, see
README.md in this directory) and `samtools` on PATH:

    TESTBED=~/herv_working_data/telescope_ccle python build_fixtures.py

Each fixture is a subset of whole fragments (every record sharing a read
name) from a real alignment file, in the original order, with SEQ and QUAL
blanked: Telescope never reads them and they are most of the file size.
Annotations are cut down to the rows a fixture can touch and to the three
attributes that matter.
"""
import hashlib
import os
import re
import subprocess
import sys
from collections import defaultdict

TESTBED = os.path.expanduser(os.environ.get('TESTBED', '~/herv_working_data/telescope_ccle'))
OUT = os.path.dirname(os.path.abspath(__file__))
# A gene-dense, HERV-containing stretch of hg38 shared by most fixtures.
WINDOW = ('chr19', int(os.environ.get('WINDOW_START', 0)), int(os.environ.get('WINDOW_END', 0)))
S1 = 'SRR8618300'


def records(path):
    """Yield (header_lines, None) once, then (None, bundle) per read name."""
    p = subprocess.Popen(['samtools', 'view', '-h', '--no-PG', path], stdout=subprocess.PIPE, text=True)
    header, bundle, name = [], [], None
    for line in p.stdout:
        if line.startswith('@'):
            header.append(line)
            continue
        if header is not None:
            yield header, None
            header = None
        f = line.rstrip('\n').split('\t')
        if f[0] != name and bundle:
            yield None, bundle
            bundle = []
        name = f[0]
        bundle.append(f)
    if bundle:
        yield None, bundle


def keep_fraction(name, one_in):
    """Deterministic 1-in-N sample by read name."""
    return int(hashlib.md5(name.encode()).hexdigest()[:8], 16) % one_in == 0


def blocks(rec):
    """Reference spans of the aligned (M/=/X) runs of a SAM record, 0-based half-open."""
    pos, out = int(rec[3]) - 1, []
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', rec[5]):
        n = int(n)
        if op in 'M=X':
            out.append((pos, pos + n))
            pos += n
        elif op in 'DN':
            pos += n
    return out


def subset(src, dst, choose, limit=None, as_sam=False):
    """Write the fragments `choose(bundle)` accepts; return the records written."""
    kept, n = [], 0
    tmp = dst + '.sam'
    with open(tmp, 'w') as out:
        for header, bundle in records(src):
            if header is not None:
                out.writelines(header)
                continue
            if limit is not None and n >= limit:
                break
            if not choose(bundle):
                continue
            n += 1
            for f in bundle:
                f[9] = f[10] = '*'
                out.write('\t'.join(f) + '\n')
                kept.append(f)
    if as_sam:
        os.replace(tmp, dst)
    else:
        subprocess.check_call(['samtools', 'view', '-b', '--no-PG', '-o', dst, tmp])
        os.remove(tmp)
    flags = [int(f[1]) for f in kept]
    print('%-34s %6d fragments %7d records (%d supplementary, %d secondary, %d unmapped)' % (
        os.path.relpath(dst, OUT), n, len(kept), sum(1 for x in flags if x & 0x800),
        sum(1 for x in flags if x & 0x100), sum(1 for x in flags if x & 0x4)))
    return kept


def in_window(bundle):
    c, s, e = WINDOW
    return any(f[2] == c and not int(f[1]) & 0x4 and s <= int(f[3]) < e for f in bundle)


def chimeric_or_sampled(bundle):
    """Long reads: every in-window read with a supplementary (chimeric) alignment, and one in ten of the rest."""
    return in_window(bundle) and (any(int(f[1]) & 0x800 for f in bundle) or keep_fraction(bundle[0][0], 10))


def names_in(bam):
    return {b[0][0] for h, b in records(bam) if b is not None}


def slim(line):
    """A GTF row with only the attributes the tests use."""
    f = line.rstrip('\n').split('\t')
    attrs = dict(re.findall(r'(\w+)\s+"(.+?)";', f[8]))
    f[8] = ' '.join('%s "%s";' % (k, attrs[k]) for k in ('gene_id', 'transcript_id', 'locus') if k in attrs)
    return '\t'.join(f) + '\n'


def gtf_for_window(src, dst):
    c, s, e = WINDOW
    n = 0
    with open(dst, 'w') as out:
        for line in open(src):
            if line.startswith('#'):
                continue
            f = line.split('\t', 5)
            if f[0] == c and int(f[3]) < e and int(f[4]) > s:
                out.write(slim(line))
                n += 1
    print('%-34s %6d rows' % (os.path.relpath(dst, OUT), n))


def gtf_for_records(src, dst, recs):
    """Keep every row of every locus that any kept alignment comes within a base of."""
    rows = [l for l in open(src) if not l.startswith('#')]
    bins = defaultdict(list)
    for i, l in enumerate(rows):
        f = l.split('\t', 5)
        for b in range(int(f[3]) // 10000, int(f[4]) // 10000 + 1):
            bins[(f[0], b)].append((int(f[3]) - 2, int(f[4]) + 2, i))
    locus = [dict(re.findall(r'(\w+)\s+"(.+?)";', l.rstrip('\n').split('\t')[8]))['locus'] for l in rows]
    touched = set()
    for f in recs:
        if int(f[1]) & 0x4:
            continue
        for bs, be in blocks(f):
            for b in range(bs // 10000, be // 10000 + 1):
                for s, e, i in bins.get((f[2], b), ()):
                    if s < be and bs < e:
                        touched.add(locus[i])
    n = 0
    with open(dst, 'w') as out:
        for l, name in zip(rows, locus):
            if name in touched:
                out.write(slim(l))
                n += 1
    print('%-34s %6d rows, %d loci' % (os.path.relpath(dst, OUT), n, len(touched)))


def main():
    T = TESTBED
    if not WINDOW[2]:
        sys.exit('set WINDOW_START and WINDOW_END (hg38 chr19 coordinates)')

    # 1. Paired-end, unstranded, T2T-CHM13, HERV annotation: every fragment Telescope found
    #    overlapping the annotation in sample 1, plus a sprinkling of the rest. Large enough
    #    (over 8,192 fragments, rows of 9+ candidates, a weight that underflows) to make
    #    numpy's summation order and scipy's zero-dropping matter.
    if os.environ.get('ONLY'):
        return window_fixtures(T)
    hits = names_in('%s/telescope_updated_sam/hs1_herv/%s/%s-tmp_tele.bam' % (T, S1, S1))
    recs = subset('%s/bam/hs1/%s.bam' % (T, S1), OUT + '/pe_hs1_herv.bam',
                  lambda b: b[0][0] in hits or keep_fraction(b[0][0], 3000 if not int(b[0][1]) & 0x4 else 1500))
    gtf_for_records('%s/annot/hs1_herv.genome.gtf' % T, OUT + '/hs1_herv_subset.gtf', recs)

    # 2. The same library against the HML2-only GTF, whose attribute quoting is mangled
    #    (locus """"HML2.LTR1"""";) and whose chromosomes are named 1, 2, ...
    hits = names_in('%s/telescope_updated_sam/hg38_hml2/%s/%s-tmp_tele.bam' % (T, S1, S1))
    subset('%s/bam/hg38_ensembl/%s.bam' % (T, S1), OUT + '/pe_hg38_hml2.bam',
           lambda b: b[0][0] in hits or keep_fraction(b[0][0], 4000))

    # 3. One hg38 window with genes and HERVs, annotation with overlapping loci (ties and
    #    interval-tree merges), seen through each aligner and library type.
    gtf_for_window('%s/annot/hg38_herv_genes.genome.gtf' % T, OUT + '/hg38_window_herv_genes.gtf')
    window_fixtures(T)


def window_fixtures(T):
    for name, src, limit, sam, choose in [
        ('pe_hisat2', 'bam/hg38/%s.bam' % S1, 4000, False, in_window),
        ('pe_bowtie2', 'aligners/bowtie2.bam', 2500, False, in_window),
        ('pe_star', 'aligners/star_x86.bam', 2500, False, in_window),
        ('se_hisat2', 'se_lr/se_hisat2.bam', 2500, False, in_window),
        ('se_bowtie2', 'se_lr/se_bowtie2.bam', 2500, False, in_window),
        ('pe_stranded_hisat2', 'encode/bam/A549_hg38.bam', 2500, False, in_window),
        ('se_stranded_hisat2', 'se_lr/se_stranded_hisat2.bam', 2500, False, in_window),
        # kept as SAM so SAM input is exercised
        ('long_minimap2', 'se_lr/lr_minimap2.sam', None, True, chimeric_or_sampled),
    ]:
        if os.environ.get('ONLY') and os.environ['ONLY'] != name:
            continue
        subset('%s/%s' % (T, src), '%s/%s.%s' % (OUT, name, 'sam' if sam else 'bam'), choose, limit, sam)


if __name__ == '__main__':
    main()
