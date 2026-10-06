//! CPython 3.10's `set`, reduced to what decides iteration order.
//!
//! Telescope picks the feature for an alignment from
//! `Counter.most_common()[0]`, and the Counter is filled by iterating a Python
//! `set` of intervals. When two loci overlap an alignment equally, the winner
//! is whichever the set yields first, which depends on CPython's hash-table
//! layout: slot = hash & mask, linear probing, perturbation, resize points and
//! dummy entries. This mirrors `Objects/setobject.c` closely enough to yield
//! the same order.
//!
//! Members are interval ids; two members are equal only if they are the same
//! id, though different ids may share a hash (same begin and end).

const EMPTY: u32 = u32::MAX;
const DUMMY: u32 = u32::MAX - 1;
const DUMMY_HASH: u64 = u64::MAX; // CPython stores hash -1 in dummy entries
const MINSIZE: usize = 8;
const LINEAR_PROBES: usize = 9;
const PERTURB_SHIFT: u32 = 5;

#[derive(Clone, Copy)]
struct Entry {
    key: u32,
    hash: u64,
}

const VACANT: Entry = Entry { key: EMPTY, hash: 0 };

#[derive(Clone)]
pub struct PySet {
    table: Vec<Entry>,
    mask: usize,
    /// active + dummy entries
    fill: usize,
    /// active entries
    used: usize,
}

impl Default for PySet {
    fn default() -> Self {
        Self::new()
    }
}

/// `hash((begin, end))` for two non-negative Python ints below 2**61 - 1
/// (CPython's xxHash-style tuple hash; an int hashes to itself).
pub fn hash_interval(begin: u64, end: u64) -> u64 {
    const P1: u64 = 11400714785074694791;
    const P2: u64 = 14029467366897019727;
    const P5: u64 = 2870177450012600261;
    let mut acc = P5;
    for lane in [begin, end] {
        acc = acc.wrapping_add(lane.wrapping_mul(P2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(P1);
    }
    acc = acc.wrapping_add(2 ^ (P5 ^ 3527539));
    if acc == u64::MAX { 1546275796 } else { acc }
}

fn insert_clean(table: &mut [Entry], mask: usize, key: u32, hash: u64) {
    let mut perturb = hash as usize;
    let mut i = hash as usize & mask;
    loop {
        if table[i].key == EMPTY {
            table[i] = Entry { key, hash };
            return;
        }
        if i + LINEAR_PROBES <= mask {
            for j in 1..=LINEAR_PROBES {
                if table[i + j].key == EMPTY {
                    table[i + j] = Entry { key, hash };
                    return;
                }
            }
        }
        perturb >>= PERTURB_SHIFT;
        i = i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb) & mask;
    }
}

fn growth_target(used: usize) -> usize {
    if used > 50000 { used * 2 } else { used * 4 }
}

impl PySet {
    pub fn new() -> Self {
        PySet { table: vec![VACANT; MINSIZE], mask: MINSIZE - 1, fill: 0, used: 0 }
    }

    /// Back to a freshly constructed `set()`, keeping the allocation.
    pub fn reset(&mut self) {
        self.table.clear();
        self.table.resize(MINSIZE, VACANT);
        self.mask = MINSIZE - 1;
        self.fill = 0;
        self.used = 0;
    }

    pub fn is_empty(&self) -> bool {
        self.used == 0
    }

    /// Members in iteration order.
    pub fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        self.table.iter().filter(|e| e.key != EMPTY && e.key != DUMMY).map(|e| e.key)
    }

    /// `set_table_resize`
    fn resize(&mut self, minused: usize) {
        let mut newsize = MINSIZE;
        while newsize <= minused {
            newsize <<= 1;
        }
        if newsize == MINSIZE && self.table.len() == MINSIZE && self.fill == self.used {
            return;
        }
        let old = std::mem::replace(&mut self.table, vec![VACANT; newsize]);
        self.mask = newsize - 1;
        self.fill = self.used;
        for e in old {
            if e.key != EMPTY && e.key != DUMMY {
                insert_clean(&mut self.table, self.mask, e.key, e.hash);
            }
        }
    }

    /// `set_add_entry`
    pub fn add(&mut self, key: u32, hash: u64) {
        let mask = self.mask;
        let mut i = hash as usize & mask;
        let mut perturb = hash as usize;
        let mut freeslot: Option<usize> = None;
        let slot = 'search: loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            for j in i..=i + probes {
                let e = self.table[j];
                if e.key == EMPTY {
                    break 'search j;
                }
                if e.hash == hash {
                    if e.key == key {
                        return;
                    }
                } else if e.key == DUMMY {
                    freeslot = Some(j);
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb) & mask;
        };
        if let Some(f) = freeslot {
            self.used += 1;
            self.table[f] = Entry { key, hash };
            return;
        }
        self.fill += 1;
        self.used += 1;
        self.table[slot] = Entry { key, hash };
        if self.fill * 5 < mask * 3 {
            return;
        }
        self.resize(growth_target(self.used));
    }

    /// `set_lookkey`: index of `key`'s entry, if present.
    fn lookup(&self, key: u32, hash: u64) -> Option<usize> {
        let mask = self.mask;
        let mut i = hash as usize & mask;
        let mut perturb = hash as usize;
        loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            for j in i..=i + probes {
                let e = self.table[j];
                if e.key == EMPTY {
                    return None;
                }
                if e.hash == hash && e.key == key {
                    return Some(j);
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb) & mask;
        }
    }

    /// `set.remove` / `set.discard`: leaves a dummy entry behind.
    pub fn discard(&mut self, key: u32, hash: u64) -> bool {
        match self.lookup(key, hash) {
            Some(j) => {
                self.table[j] = Entry { key: DUMMY, hash: DUMMY_HASH };
                self.used -= 1;
                true
            }
            None => false,
        }
    }

    /// `set.update(other_set)` (`set_merge`).
    pub fn update(&mut self, other: &PySet) {
        if other.used == 0 {
            return;
        }
        if (self.fill + other.used) * 5 >= self.mask * 3 {
            self.resize((self.used + other.used) * 2);
        }
        if self.fill == 0 && self.mask == other.mask && other.fill == other.used {
            self.table.copy_from_slice(&other.table);
            self.fill = other.fill;
            self.used = other.used;
            return;
        }
        if self.fill == 0 {
            self.fill = other.used;
            self.used = other.used;
            for e in &other.table {
                if e.key != EMPTY && e.key != DUMMY {
                    insert_clean(&mut self.table, self.mask, e.key, e.hash);
                }
            }
            return;
        }
        for e in &other.table {
            if e.key != EMPTY && e.key != DUMMY {
                self.add(e.key, e.hash);
            }
        }
    }

    /// `set(other_set)`
    pub fn copy_of(other: &PySet) -> PySet {
        let mut s = PySet::new();
        s.update(other);
        s
    }

    /// `self -= other`
    pub fn difference_update(&mut self, other: &PySet) {
        for e in &other.table {
            if e.key != EMPTY && e.key != DUMMY {
                self.discard(e.key, e.hash);
            }
        }
        if self.fill - self.used <= self.mask / 4 {
            return;
        }
        self.resize(growth_target(self.used));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_hash_matches_cpython() {
        // hash((1, 2)) etc. on 64-bit CPython 3.10, as unsigned 64-bit.
        assert_eq!(hash_interval(1, 2) as i64, -3550055125485641917);
        assert_eq!(hash_interval(1410684, 1410774) as i64, -7529434180462733536);
    }

    #[test]
    fn small_ints_iterate_in_slot_order() {
        // {5, 1, 3} iterates 1, 3, 5: small hashes land in their own slots.
        let mut s = PySet::new();
        for k in [5u32, 1, 3] {
            s.add(k, k as u64);
        }
        assert_eq!(s.iter().collect::<Vec<_>>(), vec![1, 3, 5]);
    }

    #[test]
    fn colliding_hashes_keep_insertion_order_and_resize() {
        // {8, 0, 16}: all hash to slot 0 of an 8-slot table; iteration is 8, 0, 16.
        let mut s = PySet::new();
        for k in [8u32, 0, 16] {
            s.add(k, k as u64);
        }
        assert_eq!(s.iter().collect::<Vec<_>>(), vec![8, 0, 16]);
        // The fifth insert grows the table to 32 slots and re-homes everything.
        for k in [24u32, 32] {
            s.add(k, k as u64);
        }
        assert_eq!(s.iter().collect::<Vec<_>>(), vec![0, 32, 8, 16, 24]);
        // A removed member leaves a dummy that the next colliding insert reuses.
        s.discard(8, 8);
        s.add(40, 40);
        assert_eq!(s.iter().collect::<Vec<_>>(), vec![0, 32, 40, 16, 24]);
    }
}
