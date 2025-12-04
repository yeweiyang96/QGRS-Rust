use std::collections::HashMap;
use std::collections::hash_map::Entry;

use crate::qgrs::data::SequenceSlice;
use crate::qgrs::search::G4;

#[derive(Eq, PartialEq, Hash)]
struct DedupKey {
    start: usize,
    end: usize,
    slice: SequenceSlice,
}

impl DedupKey {
    fn new(g4: &G4) -> Self {
        Self {
            start: g4.start,
            end: g4.end,
            slice: g4.sequence_view.clone(),
        }
    }
}

fn overlapped(a: &G4, b: &G4) -> bool {
    let a_start = a.start as isize;
    let a_end = (a.start + a.length) as isize;
    let b_start = b.start as isize;
    let b_end = (b.start + b.length) as isize;
    (a_start >= b_start && a_start <= b_end)
        || (a_end >= b_start && a_end <= b_end)
        || (b_start >= a_start && b_start <= a_end)
        || (b_end >= a_start && b_end <= a_end)
}

fn belongs_in(g4: &G4, family: &[G4]) -> bool {
    family.iter().any(|member| overlapped(g4, member))
}

pub fn consolidate_g4s(mut raw_g4s: Vec<G4>) -> Vec<G4> {
    raw_g4s.sort_by(|a, b| a.start.cmp(&b.start));

    let mut best_by_key: HashMap<DedupKey, G4> = HashMap::new();
    for g in raw_g4s.into_iter() {
        let key = DedupKey::new(&g);
        match best_by_key.entry(key) {
            Entry::Vacant(slot) => {
                slot.insert(g);
            }
            Entry::Occupied(mut slot) => {
                if g.gscore > slot.get().gscore {
                    slot.insert(g);
                }
            }
        }
    }
    let mut deduped: Vec<G4> = best_by_key.into_values().collect();
    deduped.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));

    let mut families: Vec<Vec<G4>> = Vec::new();

    for g4 in deduped.into_iter() {
        let mut inserted = false;
        for family in &mut families {
            if belongs_in(&g4, family) {
                family.push(g4.clone());
                inserted = true;
                break;
            }
        }
        if !inserted {
            families.push(vec![g4]);
        }
    }

    let mut results = Vec::new();
    for family in families.into_iter() {
        if family.is_empty() {
            continue;
        }
        let mut best = family[0].clone();
        for member in &family {
            if member.gscore > best.gscore {
                best = member.clone();
            }
        }
        results.push(best);
    }

    results
}
