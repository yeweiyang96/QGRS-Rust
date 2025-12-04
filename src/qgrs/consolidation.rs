use crate::qgrs::search::G4;

fn is_better_candidate(current: &G4, candidate: &G4) -> bool {
    candidate.gscore > current.gscore
    // || (candidate.gscore == current.gscore && candidate.length < current.length)
}

pub fn consolidate_g4s(raw_g4s: Vec<G4>) -> Vec<G4> {
    if raw_g4s.is_empty() {
        return Vec::new();
    }

    debug_assert!(
        raw_g4s
            .windows(2)
            .all(|pair| pair[0].start <= pair[1].start),
        "consolidate_g4s expects raw hits sorted by start"
    );

    let mut consolidated = Vec::with_capacity(raw_g4s.len());
    let mut iter = raw_g4s.into_iter();
    let mut current_best = iter.next().expect("iterator is non-empty");
    let mut family_end = current_best.end;

    for candidate in iter {
        if candidate.start <= family_end {
            family_end = family_end.max(candidate.end);
            if is_better_candidate(&current_best, &candidate) {
                current_best = candidate;
            }
        } else {
            consolidated.push(current_best);
            current_best = candidate;
            family_end = current_best.end;
        }
    }

    consolidated.push(current_best);
    consolidated
}
