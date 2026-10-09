//! Bit-exact stand-ins for the numpy/scipy primitives Telescope leans on.
//!
//! Telescope's numbers depend on *how* numpy adds floats (pairwise blocks, not
//! left-to-right), on numpy's legacy Mersenne Twister, and on Python's float
//! formatting. Each helper here mirrors one of those behaviours and was checked
//! against numpy 1.26 / scipy 1.15.

use std::sync::atomic::{AtomicBool, Ordering};

/// When set, every sum below is a plain left-to-right loop instead of
/// numpy's blocked pairwise scheme (`--float_sums sequential`).
static SEQUENTIAL_SUMS: AtomicBool = AtomicBool::new(false);

pub fn set_sequential_sums(on: bool) {
    SEQUENTIAL_SUMS.store(on, Ordering::Relaxed);
}

fn sequential() -> bool {
    SEQUENTIAL_SUMS.load(Ordering::Relaxed)
}

/// expm1, log1p and log10 are not correctly rounded in system math libraries,
/// and their last digit differs between platforms (macOS and Linux disagree
/// on about one expm1 result in ten over the model's range). numpy calls the
/// system's, so Telescope's own output varies by machine. By default these
/// come from CORE-MATH instead, which is correctly rounded and therefore the
/// same everywhere; `--math system` uses the platform's, to match a Python
/// Telescope run on the same machine.
static SYSTEM_MATH: AtomicBool = AtomicBool::new(false);

pub fn set_system_math(on: bool) {
    SYSTEM_MATH.store(on, Ordering::Relaxed);
}

fn system_math() -> bool {
    SYSTEM_MATH.load(Ordering::Relaxed)
}

pub fn expm1(x: f64) -> f64 {
    if system_math() { x.exp_m1() } else { core_math::expm1(x) }
}

pub fn log1p(x: f64) -> f64 {
    if system_math() { x.ln_1p() } else { core_math::log1p(x) }
}

pub fn log10(x: f64) -> f64 {
    if system_math() { x.log10() } else { core_math::log10(x) }
}

/// `ndarray.sum()` of a 1-D f64 array. numpy feeds the reduction through its
/// 8192-element buffer, so the pairwise sum restarts every 8192 values and
/// the block totals are added left to right.
pub fn np_sum(a: &[f64]) -> f64 {
    const BUFSIZE: usize = 8192;
    if sequential() {
        return a.iter().sum();
    }
    let mut chunks = a.chunks(BUFSIZE);
    let mut acc = chunks.next().map_or(0.0, pairwise_sum);
    for c in chunks {
        acc += pairwise_sum(c);
    }
    acc
}

/// numpy's `pairwise_sum` over one contiguous block.
pub fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8;
        let stop = n - (n % 8);
        while i < stop {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

/// One segment of `np.add.reduceat`, which is how scipy sums a CSR row: the
/// first element seeds the accumulator and the rest are pairwise-summed.
pub fn reduceat_sum(a: &[f64]) -> f64 {
    if sequential() {
        return a.iter().sum();
    }
    match a.len() {
        0 => 0.0,
        1 => a[0],
        _ => a[0] + pairwise_sum(&a[1..]),
    }
}

/// `sparse_plus._recip0`: reciprocal with infinities replaced by zero.
pub fn recip0(v: f64) -> f64 {
    let r = 1.0 / v;
    if r.is_infinite() { 0.0 } else { r }
}

/// numpy's legacy `RandomState` generator (MT19937 seeded by `init_genrand`).
pub struct Mt19937 {
    mt: [u32; 624],
    idx: usize,
}

impl Mt19937 {
    /// `np.random.seed(seed)` for an integer seed.
    pub fn new(seed: u32) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = seed;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Mt19937 { mt, idx: 624 }
    }

    fn generate(&mut self) {
        for k in 0..624 {
            let y = (self.mt[k] & 0x8000_0000) | (self.mt[(k + 1) % 624] & 0x7fff_ffff);
            let mag = if y & 1 == 1 { 0x9908_b0df } else { 0 };
            self.mt[k] = self.mt[(k + 397) % 624] ^ (y >> 1) ^ mag;
        }
        self.idx = 0;
    }

    fn next_u32(&mut self) -> u32 {
        if self.idx >= 624 {
            self.generate();
        }
        let mut y = self.mt[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// `np.random.random_sample()`: a double in [0, 1) from two 32-bit draws.
    pub fn random_f64(&mut self) -> f64 {
        let a = (self.next_u32() >> 5) as f64;
        let b = (self.next_u32() >> 6) as f64;
        (a * 67108864.0 + b) / 9007199254740992.0
    }

    /// Index drawn by `np.random.choice(range(n))`: masked rejection sampling
    /// on 32-bit draws. Consumes nothing when `n == 1`.
    pub fn choice(&mut self, n: u32) -> u32 {
        let rng = n - 1;
        if rng == 0 {
            return 0;
        }
        let mut mask = rng;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        loop {
            let v = self.next_u32() & mask;
            if v <= rng {
                return v;
            }
        }
    }
}

fn nonfinite(x: f64) -> Option<String> {
    if x.is_nan() {
        Some("nan".to_string())
    } else if x.is_infinite() {
        Some(if x > 0.0 { "inf" } else { "-inf" }.to_string())
    } else {
        None
    }
}

/// Python `'{:.2f}'.format(x)`.
pub fn fmt_f2(x: f64) -> String {
    nonfinite(x).unwrap_or_else(|| format!("{x:.2}"))
}

/// Python `'{:.3g}'.format(x)`.
pub fn fmt_g3(x: f64) -> String {
    if let Some(s) = nonfinite(x) {
        return s;
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    // Round to 3 significant digits first; the exponent of the *rounded*
    // value picks fixed vs scientific notation.
    let sci = format!("{x:.2e}");
    let (mant, exp) = sci.split_once('e').expect("scientific format");
    let exp: i32 = exp.parse().expect("exponent");
    let strip = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if !(-4..3).contains(&exp) {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", strip(mant), sign, exp.abs())
    } else {
        let decimals = (2 - exp).max(0) as usize;
        strip(&format!("{x:.decimals$}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g3_matches_python() {
        // Expected strings produced by CPython's '{:.3g}'.format.
        let cases: &[(f64, &str)] = &[
            (0.998, "0.998"),
            (7.28e-15, "7.28e-15"),
            (0.00161, "0.00161"),
            (1.0, "1"),
            (0.5, "0.5"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (0.99999, "1"),
            (999.5, "1e+03"),
            (123.456, "123"),
            (0.000123456, "0.000123"),
            (1.5e-300, "1.5e-300"),
            (0.0, "0"),
            (0.0799, "0.0799"),
            (0.1005, "0.101"),
        ];
        for &(x, want) in cases {
            assert_eq!(fmt_g3(x), want, "x={x}");
        }
    }

    #[test]
    fn random_f64_matches_numpy() {
        // np.random.seed(12345); np.random.random_sample()
        let mut rng = Mt19937::new(12345);
        assert_eq!(rng.random_f64(), 0.9296160928171479);
    }

    #[test]
    fn choice_matches_numpy() {
        // np.random.seed(12345); [np.random.choice(range(n)) for n in (2,3,2,5,7,100,2,2,9,33)]
        let mut rng = Mt19937::new(12345);
        let got: Vec<u32> = [2, 3, 2, 5, 7, 100, 2, 2, 9, 33]
            .iter()
            .map(|&n| rng.choice(n))
            .collect();
        assert_eq!(got, vec![0, 1, 1, 1, 4, 41, 0, 0, 5, 29]);
    }

    #[test]
    fn np_sum_restarts_every_8192() {
        // Chosen so one-pass pairwise and chunked sums round differently.
        let a: Vec<f64> = (0..20000).map(|i| ((i * 37 % 101) as f64).exp()).collect();
        let chunked = pairwise_sum(&a[..8192]) + pairwise_sum(&a[8192..16384]) + pairwise_sum(&a[16384..]);
        assert_eq!(np_sum(&a), chunked);
        assert_eq!(np_sum(&a[..8192]), pairwise_sum(&a[..8192]));
    }

    #[test]
    fn pairwise_small_is_sequential() {
        let a = [1e16, 1.0, -1e16, 1.0];
        assert_eq!(pairwise_sum(&a), ((0.0 + 1e16 + 1.0) - 1e16) + 1.0);
    }
}
