//! The `intervaltree` 3.1.0 package, ported operation for operation.
//!
//! Only the *order* in which `overlap()` yields intervals matters here (see
//! [`super::pyset`]), but that order depends on which tree node each interval
//! sits in, and so on every rotation and prune since the tree was built. The
//! node logic below therefore follows `intervaltree/node.py` line by line,
//! including where it is unusual, rather than being a fresh interval tree.

use std::collections::BTreeMap;

use super::pyset::{PySet, hash_interval};

const NIL: u32 = u32::MAX;

#[derive(Clone, Copy)]
pub struct Iv {
    pub begin: u64,
    pub end: u64,
    pub locus: u32,
    pub strand: u8,
    hash: u64,
}

impl Iv {
    pub fn new(begin: u64, end: u64, locus: u32, strand: u8) -> Self {
        Iv { begin, end, locus, strand, hash: hash_interval(begin, end) }
    }
    fn contains_point(&self, p: u64) -> bool {
        self.begin <= p && p < self.end
    }
}

struct Node {
    x_center: u64,
    s_center: PySet,
    child: [u32; 2],
    depth: i32,
    balance: i32,
}

/// Interval storage shared by the per-chromosome trees. Slots of removed
/// intervals are reused, so the table stays the size of the live annotation.
#[derive(Default)]
pub struct IvTable {
    ivs: Vec<Iv>,
    free: Vec<u32>,
}

impl IvTable {
    pub fn push(&mut self, iv: Iv) -> u32 {
        match self.free.pop() {
            Some(id) => {
                self.ivs[id as usize] = iv;
                id
            }
            None => {
                self.ivs.push(iv);
                (self.ivs.len() - 1) as u32
            }
        }
    }
    /// Call once the interval is out of every tree.
    pub fn release(&mut self, id: u32) {
        self.free.push(id);
    }
    pub fn get(&self, id: u32) -> &Iv {
        &self.ivs[id as usize]
    }
    pub fn shrink(&mut self) {
        self.ivs.shrink_to_fit();
        self.free = Vec::new();
    }
}

pub struct PyIntervalTree {
    nodes: Vec<Node>,
    /// Slots of nodes that dropped out of the tree, for reuse.
    free: Vec<u32>,
    top: u32,
    boundary: BTreeMap<u64, u32>,
}

impl Default for PyIntervalTree {
    fn default() -> Self {
        PyIntervalTree { nodes: Vec::new(), free: Vec::new(), top: NIL, boundary: BTreeMap::new() }
    }
}

impl PyIntervalTree {
    /// `Node(x_center, s_center)`
    fn new_node(&mut self, x_center: u64, s_center: PySet, ivs: &IvTable) -> u32 {
        let node = Node { x_center, s_center, child: [NIL, NIL], depth: 0, balance: 0 };
        let n = match self.free.pop() {
            Some(n) => {
                self.nodes[n as usize] = node;
                n
            }
            None => {
                self.nodes.push(node);
                (self.nodes.len() - 1) as u32
            }
        };
        self.rotate(n, ivs)
    }

    /// The node has left the tree (Python would garbage-collect it).
    fn discard_node(&mut self, n: u32) {
        self.nodes[n as usize].s_center = PySet::new();
        self.free.push(n);
    }

    pub fn shrink(&mut self) {
        self.nodes.shrink_to_fit();
        self.free = Vec::new();
    }

    fn refresh_balance(&mut self, n: u32) {
        let depth_of = |c: u32| if c == NIL { 0 } else { self.nodes[c as usize].depth };
        let [l, r] = self.nodes[n as usize].child;
        let (ld, rd) = (depth_of(l), depth_of(r));
        let node = &mut self.nodes[n as usize];
        node.depth = 1 + ld.max(rd);
        node.balance = rd - ld;
    }

    fn rotate(&mut self, n: u32, ivs: &IvTable) -> u32 {
        self.refresh_balance(n);
        let bal = self.nodes[n as usize].balance;
        if bal.abs() < 2 {
            return n;
        }
        let my_heavy = bal > 0;
        let child_bal = self.nodes[self.nodes[n as usize].child[my_heavy as usize] as usize].balance;
        let child_heavy = child_bal > 0;
        if my_heavy == child_heavy || child_bal == 0 { self.srotate(n, ivs) } else { self.drotate(n, ivs) }
    }

    fn srotate(&mut self, n: u32, ivs: &IvTable) -> u32 {
        let heavy = (self.nodes[n as usize].balance > 0) as usize;
        let light = 1 - heavy;
        let save = self.nodes[n as usize].child[heavy];
        self.nodes[n as usize].child[heavy] = self.nodes[save as usize].child[light];
        let rotated = self.rotate(n, ivs);
        self.nodes[save as usize].child[light] = rotated;
        let save_x = self.nodes[save as usize].x_center;
        let promotees: Vec<u32> = self.nodes[rotated as usize]
            .s_center
            .iter()
            .filter(|&id| ivs.get(id).contains_point(save_x))
            .collect();
        if !promotees.is_empty() {
            for &id in &promotees {
                let sub = self.nodes[save as usize].child[light];
                let sub = self.remove(sub, id, ivs);
                self.nodes[save as usize].child[light] = sub;
            }
            for &id in &promotees {
                self.nodes[save as usize].s_center.add(id, ivs.get(id).hash);
            }
        }
        self.refresh_balance(save);
        save
    }

    fn drotate(&mut self, n: u32, ivs: &IvTable) -> u32 {
        let my_heavy = (self.nodes[n as usize].balance > 0) as usize;
        let c = self.nodes[n as usize].child[my_heavy];
        self.nodes[n as usize].child[my_heavy] = self.srotate(c, ivs);
        self.refresh_balance(n);
        self.srotate(n, ivs)
    }

    fn add_at(&mut self, n: u32, id: u32, ivs: &IvTable) -> u32 {
        let iv = *ivs.get(id);
        let x = self.nodes[n as usize].x_center;
        if iv.contains_point(x) {
            self.nodes[n as usize].s_center.add(id, iv.hash);
            return n;
        }
        let dir = (iv.begin > x) as usize;
        let c = self.nodes[n as usize].child[dir];
        if c == NIL {
            let mut s = PySet::new();
            s.add(id, iv.hash);
            let leaf = self.new_node(iv.begin, s, ivs);
            self.nodes[n as usize].child[dir] = leaf;
            self.refresh_balance(n);
            n
        } else {
            let sub = self.add_at(c, id, ivs);
            self.nodes[n as usize].child[dir] = sub;
            self.rotate(n, ivs)
        }
    }

    /// `Node.remove`: returns the subtree's new root (NIL if it emptied).
    fn remove(&mut self, n: u32, id: u32, ivs: &IvTable) -> u32 {
        let mut done = false;
        self.remove_helper(n, id, &mut done, ivs)
    }

    fn remove_helper(&mut self, n: u32, id: u32, done: &mut bool, ivs: &IvTable) -> u32 {
        let iv = *ivs.get(id);
        let x = self.nodes[n as usize].x_center;
        if iv.contains_point(x) {
            let found = self.nodes[n as usize].s_center.discard(id, iv.hash);
            assert!(found, "interval missing from tree node");
            if !self.nodes[n as usize].s_center.is_empty() {
                *done = true;
                return n;
            }
            return self.prune(n, ivs);
        }
        let dir = (iv.begin > x) as usize;
        let c = self.nodes[n as usize].child[dir];
        assert!(c != NIL, "interval missing from tree");
        let sub = self.remove_helper(c, id, done, ivs);
        self.nodes[n as usize].child[dir] = sub;
        if !*done { self.rotate(n, ivs) } else { n }
    }

    fn prune(&mut self, n: u32, ivs: &IvTable) -> u32 {
        let [l, r] = self.nodes[n as usize].child;
        self.discard_node(n);
        if l == NIL || r == NIL {
            return if l == NIL { r } else { l };
        }
        let (heir, rest) = self.pop_greatest_child(l, ivs);
        self.nodes[heir as usize].child = [rest, r];
        self.refresh_balance(heir);
        self.rotate(heir, ivs)
    }

    fn pop_greatest_child(&mut self, n: u32, ivs: &IvTable) -> (u32, u32) {
        let right = self.nodes[n as usize].child[1];
        if right == NIL {
            // This node is the greatest child.
            let mut sorted: Vec<u32> = self.nodes[n as usize].s_center.iter().collect();
            sorted.sort_by_key(|&id| (ivs.get(id).end, ivs.get(id).begin));
            let max_end = ivs.get(sorted.pop().expect("non-empty s_center")).end;
            let mut new_x = self.nodes[n as usize].x_center;
            while let Some(id) = sorted.pop() {
                let end = ivs.get(id).end;
                if end == max_end {
                    continue;
                }
                new_x = new_x.max(end);
            }
            let mut moved = PySet::new();
            for id in self.nodes[n as usize].s_center.iter() {
                if ivs.get(id).contains_point(new_x) {
                    moved.add(id, ivs.get(id).hash);
                }
            }
            let remaining = {
                let s = &mut self.nodes[n as usize].s_center;
                s.difference_update(&moved);
                !s.is_empty()
            };
            if remaining {
                (self.new_node(new_x, moved, ivs), n)
            } else {
                let left = self.nodes[n as usize].child[0];
                self.discard_node(n);
                (self.new_node(new_x, moved, ivs), left)
            }
        } else {
            let (greatest, rest) = self.pop_greatest_child(right, ivs);
            self.nodes[n as usize].child[1] = rest;
            let gx = self.nodes[greatest as usize].x_center;
            let snapshot = PySet::copy_of(&self.nodes[n as usize].s_center);
            for id in snapshot.iter() {
                let iv = *ivs.get(id);
                if iv.contains_point(gx) {
                    self.nodes[n as usize].s_center.discard(id, iv.hash);
                    // greatest.add(iv): always a centre hit, result unused
                    self.nodes[greatest as usize].s_center.add(id, iv.hash);
                }
            }
            if !self.nodes[n as usize].s_center.is_empty() {
                self.refresh_balance(n);
                (greatest, self.rotate(n, ivs))
            } else {
                (greatest, self.prune(n, ivs))
            }
        }
    }

    fn search_point(&self, mut n: u32, point: u64, result: &mut PySet, ivs: &IvTable) {
        loop {
            let node = &self.nodes[n as usize];
            for id in node.s_center.iter() {
                let iv = ivs.get(id);
                if iv.contains_point(point) {
                    result.add(id, iv.hash);
                }
            }
            let next = if point < node.x_center {
                node.child[0]
            } else if point > node.x_center {
                node.child[1]
            } else {
                NIL
            };
            if next == NIL {
                return;
            }
            n = next;
        }
    }

    /// Ids of every interval in the tree (in no particular order).
    pub fn interval_ids(&self) -> Vec<u32> {
        let mut out = Vec::new();
        let mut stack = if self.top == NIL { Vec::new() } else { vec![self.top] };
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n as usize];
            out.extend(node.s_center.iter());
            stack.extend(node.child.iter().copied().filter(|&c| c != NIL));
        }
        out
    }

    /// `IntervalTree.add`
    pub fn add(&mut self, id: u32, ivs: &IvTable) {
        let iv = *ivs.get(id);
        if self.top == NIL {
            let mut s = PySet::new();
            s.add(id, iv.hash);
            self.top = self.new_node(iv.begin, s, ivs);
        } else {
            self.top = self.add_at(self.top, id, ivs);
        }
        *self.boundary.entry(iv.begin).or_insert(0) += 1;
        *self.boundary.entry(iv.end).or_insert(0) += 1;
    }

    /// `IntervalTree.remove`
    pub fn remove_interval(&mut self, id: u32, ivs: &IvTable) {
        let iv = *ivs.get(id);
        self.top = self.remove(self.top, id, ivs);
        for p in [iv.begin, iv.end] {
            let c = self.boundary.get_mut(&p).expect("boundary present");
            *c -= 1;
            if *c == 0 {
                self.boundary.remove(&p);
            }
        }
    }

    /// `IntervalTree.overlap(begin, end)`, leaving the result in `result` in
    /// Python's iteration order. `scratch` is a second reusable set.
    pub fn overlap(&self, begin: u64, end: u64, result: &mut PySet, scratch: &mut PySet, ivs: &IvTable) {
        result.reset();
        if self.top == NIL || begin >= end {
            return;
        }
        self.search_point(self.top, begin, result, ivs);
        scratch.reset();
        for (&p, _) in self.boundary.range(begin..end) {
            self.search_point(self.top, p, scratch, ivs);
        }
        result.update(scratch);
    }
}
