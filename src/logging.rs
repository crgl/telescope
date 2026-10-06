use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verbosity {
    Quiet,
    Normal,
    Verbose,
}

pub struct Logger {
    start: Instant,
    verbosity: Verbosity,
}

impl Logger {
    pub fn new(verbosity: Verbosity) -> Self {
        Logger {
            start: Instant::now(),
            verbosity,
        }
    }

    /// Log a pipeline stage checkpoint (Normal and Verbose).
    pub fn stage(&self, msg: &str) {
        if self.verbosity == Verbosity::Quiet {
            return;
        }
        let elapsed = self.start.elapsed().as_secs_f64();
        eprintln!("[{elapsed:7.2}s] {msg}");
    }

    /// Log extra detail (Verbose only).
    pub fn detail(&self, msg: &str) {
        if self.verbosity != Verbosity::Verbose {
            return;
        }
        let elapsed = self.start.elapsed().as_secs_f64();
        eprintln!("[{elapsed:7.2}s]   {msg}");
    }

    pub fn elapsed_secs(&self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }
}

/// Return the peak resident set size of this process so far, in bytes.
/// `None` on platforms where we cannot read it (e.g. Windows builds, where
/// `libc` is not pulled in via the target-conditional dependency).
#[cfg(unix)]
pub fn peak_rss_bytes() -> Option<u64> {
    // SAFETY: getrusage with RUSAGE_SELF and a writable, properly-aligned
    // libc::rusage buffer is sound. Zero-initialization is valid for rusage
    // (all fields are integer types).
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) };
    if rc != 0 {
        return None;
    }
    let raw = usage.ru_maxrss as u64;
    // macOS reports ru_maxrss in bytes; Linux/BSD report it in kilobytes.
    #[cfg(target_os = "macos")]
    let bytes = raw;
    #[cfg(not(target_os = "macos"))]
    let bytes = raw * 1024;
    Some(bytes)
}

#[cfg(not(unix))]
pub fn peak_rss_bytes() -> Option<u64> {
    None
}

/// Format a byte count as a human-readable string (GB / MB / KB).
pub fn format_bytes(n: u64) -> String {
    const KB: u64 = 1024;
    const MB: u64 = KB * 1024;
    const GB: u64 = MB * 1024;
    if n >= GB {
        format!("{:.2} GB", n as f64 / GB as f64)
    } else if n >= MB {
        format!("{:.2} MB", n as f64 / MB as f64)
    } else {
        format!("{} KB", n / KB)
    }
}
