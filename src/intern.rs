use std::collections::HashMap;

/// A compact, Copy-able handle to an interned string.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Symbol(u32);

/// String interner: maps unique strings to compact u32 symbols.
/// Built during GTF loading (single-threaded), then used read-only during parallel phases.
pub struct Interner {
    map: HashMap<String, Symbol>,
    strings: Vec<String>,
}

impl Interner {
    pub fn new() -> Self {
        Interner {
            map: HashMap::new(),
            strings: Vec::new(),
        }
    }

    /// Intern a string, returning its symbol. Returns existing symbol if already interned.
    pub fn intern(&mut self, s: &str) -> Symbol {
        if let Some(&sym) = self.map.get(s) {
            return sym;
        }
        let sym = Symbol(self.strings.len() as u32);
        self.strings.push(s.to_string());
        self.map.insert(s.to_string(), sym);
        sym
    }

    /// Look up the string for a symbol.
    pub fn resolve(&self, sym: Symbol) -> &str {
        &self.strings[sym.0 as usize]
    }

    /// Number of interned symbols. Equal to one past the highest valid `Symbol::idx()`.
    pub fn len(&self) -> usize {
        self.strings.len()
    }
}

impl Symbol {
    /// Dense `usize` index into a vector keyed by symbol (e.g. EM's theta).
    pub fn idx(self) -> usize {
        self.0 as usize
    }
}
