use std::collections::BTreeMap;

use crate::qgrs::data::SequenceTopology;
use crate::qgrs::search::G4;

fn is_better_candidate(current: &G4, candidate: &G4) -> bool {
    candidate.gscore > current.gscore
        || (candidate.gscore == current.gscore && candidate.length < current.length)
}

pub fn consolidate_g4s(raw_g4s: Vec<G4>) -> (Vec<G4>, Vec<(usize, usize)>) {
    consolidate_linear(raw_g4s)
}

pub fn consolidate_g4s_with_topology(
    raw_g4s: Vec<G4>,
    topology: SequenceTopology,
    sequence_len: usize,
) -> (Vec<G4>, Vec<(usize, usize)>) {
    if topology.is_circular() {
        return consolidate_circular(raw_g4s, sequence_len);
    }
    consolidate_linear(raw_g4s)
}

fn consolidate_linear(raw_g4s: Vec<G4>) -> (Vec<G4>, Vec<(usize, usize)>) {
    if raw_g4s.is_empty() {
        return (Vec::new(), Vec::new());
    }

    debug_assert!(
        raw_g4s
            .windows(2)
            .all(|pair| pair[0].start <= pair[1].start),
        "consolidate_g4s expects raw hits sorted by start"
    );

    let mut consolidated = Vec::with_capacity(raw_g4s.len());
    let mut family_ranges: Vec<(usize, usize)> = Vec::new();
    let mut iter = raw_g4s.into_iter();
    let mut current_best = iter.next().expect("iterator is non-empty");
    let mut family_start = current_best.start;
    let mut family_end = current_best.end;

    for candidate in iter {
        if candidate.start <= family_end {
            family_end = family_end.max(candidate.end);
            if is_better_candidate(&current_best, &candidate) {
                current_best = candidate;
            }
        } else {
            family_ranges.push((family_start, family_end));
            consolidated.push(current_best);
            current_best = candidate;
            family_start = current_best.start;
            family_end = current_best.end;
        }
    }

    family_ranges.push((family_start, family_end));
    consolidated.push(current_best);
    (consolidated, family_ranges)
}

fn consolidate_circular(raw_g4s: Vec<G4>, sequence_len: usize) -> (Vec<G4>, Vec<(usize, usize)>) {
    if raw_g4s.is_empty() || sequence_len == 0 {
        return (Vec::new(), Vec::new());
    }
    debug_assert!(
        raw_g4s.iter().all(|g4| g4.start <= sequence_len),
        "circular consolidation expects start coordinates within sequence length"
    );

    let mut dsu = DisjointSet::new(raw_g4s.len());
    let mut segments = Vec::with_capacity(raw_g4s.len() * 2);
    for (owner, g4) in raw_g4s.iter().enumerate() {
        segments.push(Segment {
            start: g4.start,
            end: g4.end,
            owner,
        });
        segments.push(Segment {
            start: g4.start + sequence_len,
            end: g4.end + sequence_len,
            owner,
        });
    }
    segments.sort_by_key(|segment| (segment.start, segment.end, segment.owner));

    let mut current_owner: Option<usize> = None;
    let mut current_end = 0usize;
    for segment in segments {
        match current_owner {
            Some(owner) if segment.start <= current_end => {
                dsu.union(segment.owner, owner);
                current_end = current_end.max(segment.end);
            }
            _ => {
                current_owner = Some(segment.owner);
                current_end = segment.end;
            }
        }
    }

    let mut members_by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for index in 0..raw_g4s.len() {
        let root = dsu.find(index);
        members_by_root.entry(root).or_default().push(index);
    }

    let mut grouped: Vec<(usize, usize, G4)> = Vec::with_capacity(members_by_root.len());
    for members in members_by_root.values() {
        let mut best_index = members[0];
        for &candidate_index in members.iter().skip(1) {
            if is_better_candidate(&raw_g4s[best_index], &raw_g4s[candidate_index]) {
                best_index = candidate_index;
            }
        }
        let family_range = circular_family_range(&raw_g4s, members, sequence_len);
        grouped.push((family_range.0, family_range.1, raw_g4s[best_index].clone()));
    }

    grouped.sort_by_key(|(start, end, best)| (*start, *end, best.start, best.end));
    let mut consolidated = Vec::with_capacity(grouped.len());
    let mut family_ranges = Vec::with_capacity(grouped.len());
    for (start, end, best) in grouped {
        family_ranges.push((start, end));
        consolidated.push(best);
    }
    (consolidated, family_ranges)
}

fn circular_family_range(raw_g4s: &[G4], members: &[usize], sequence_len: usize) -> (usize, usize) {
    let mut entries: Vec<(usize, usize)> = members
        .iter()
        .map(|&index| {
            let g4 = &raw_g4s[index];
            (g4.start, g4.end)
        })
        .collect();
    entries.sort_unstable_by_key(|(start, end)| (*start, *end));

    let mut prefix_max_end = vec![0usize; entries.len()];
    let mut suffix_max_end = vec![0usize; entries.len()];
    let mut running_prefix = 0usize;
    for (idx, &(_, end)) in entries.iter().enumerate() {
        running_prefix = running_prefix.max(end);
        prefix_max_end[idx] = running_prefix;
    }
    let mut running_suffix = 0usize;
    for idx in (0..entries.len()).rev() {
        running_suffix = running_suffix.max(entries[idx].1);
        suffix_max_end[idx] = running_suffix;
    }

    let mut best: Option<(usize, usize, usize)> = None;
    let mut group_start = 0usize;
    while group_start < entries.len() {
        let cut = entries[group_start].0;
        let prefix_end = if group_start == 0 {
            0
        } else {
            prefix_max_end[group_start - 1].saturating_add(sequence_len)
        };
        let end = suffix_max_end[group_start].max(prefix_end);
        let span = end.saturating_sub(cut);
        let replace = match best {
            None => true,
            Some((best_span, best_start, best_end)) => {
                span < best_span
                    || (span == best_span
                        && (cut < best_start || (cut == best_start && end < best_end)))
            }
        };
        if replace {
            best = Some((span, cut, end));
        }

        let mut next = group_start + 1;
        while next < entries.len() && entries[next].0 == cut {
            next += 1;
        }
        group_start = next;
    }

    let (_, start, end) = best.expect("members is non-empty");
    (start, end)
}

#[derive(Clone, Copy)]
struct Segment {
    start: usize,
    end: usize,
    owner: usize,
}

struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, node: usize) -> usize {
        if self.parent[node] != node {
            let root = self.find(self.parent[node]);
            self.parent[node] = root;
        }
        self.parent[node]
    }

    fn union(&mut self, lhs: usize, rhs: usize) {
        let mut lhs_root = self.find(lhs);
        let mut rhs_root = self.find(rhs);
        if lhs_root == rhs_root {
            return;
        }
        if self.rank[lhs_root] < self.rank[rhs_root] {
            std::mem::swap(&mut lhs_root, &mut rhs_root);
        }
        self.parent[rhs_root] = lhs_root;
        if self.rank[lhs_root] == self.rank[rhs_root] {
            self.rank[lhs_root] += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::{circular_family_range, consolidate_circular};
    use crate::qgrs::{ScanLimits, SequenceTopology, find_owned_bytes_with_topology};

    fn arc_from_sequence(seq: &str) -> Arc<Vec<u8>> {
        Arc::new(seq.bytes().map(|b| b.to_ascii_lowercase()).collect())
    }

    #[test]
    fn circular_family_ranges_have_at_most_one_wraparound_interval() {
        let sequence_len = 20;
        let raw_g4s = find_owned_bytes_with_topology(
            arc_from_sequence("GAGGGGAGGGGAGGGGGGG"),
            4,
            17,
            ScanLimits::default(),
            SequenceTopology::Circular,
        );

        let (_hits, family_ranges) = consolidate_circular(raw_g4s, sequence_len);
        let wraparound_count = family_ranges
            .iter()
            .filter(|(_, end)| *end > sequence_len)
            .count();

        assert!(
            wraparound_count <= 1,
            "expected at most one wraparound family, got {wraparound_count:?} from {family_ranges:?}"
        );
    }

    #[test]
    fn circular_family_range_is_a_continuous_expanded_interval() {
        let sequence_len = 20;
        let raw_g4s = find_owned_bytes_with_topology(
            arc_from_sequence("GAGGGGAGGGGAGGGGGGG"),
            4,
            17,
            ScanLimits::default(),
            SequenceTopology::Circular,
        );
        let members: Vec<usize> = (0..raw_g4s.len()).collect();

        let family_range = circular_family_range(&raw_g4s, &members, sequence_len);
        assert!(family_range.0 <= sequence_len);
        assert!(family_range.1 > sequence_len);

        for &index in &members {
            let g4 = &raw_g4s[index];
            let (start, end) = if g4.start < family_range.0 {
                (g4.start + sequence_len, g4.end + sequence_len)
            } else {
                (g4.start, g4.end)
            };
            assert!(
                start >= family_range.0 && end <= family_range.1,
                "member {:?} does not fit continuously inside family range {:?} after expansion to ({start}, {end})",
                (g4.start, g4.end),
                family_range
            );
        }
    }

    #[test]
    fn circular_dense_overlaps_collapse_into_one_family() {
        let sequence = "GGGGAGGGGAGGGGAGGGGAGGGG";
        let raw_g4s = find_owned_bytes_with_topology(
            arc_from_sequence(sequence),
            4,
            17,
            ScanLimits::default(),
            SequenceTopology::Circular,
        );

        let (hits, family_ranges) = consolidate_circular(raw_g4s.clone(), sequence.len());

        eprintln!(
            "sequence_len={} raw_hits={} consolidated_hits={} family_ranges={:?}",
            sequence.len(),
            raw_g4s.len(),
            hits.len(),
            family_ranges
        );

        assert!(
            raw_g4s.len() > 10,
            "expected many overlapping raw hits, got {}",
            raw_g4s.len()
        );
        assert_eq!(hits.len(), 1, "expected one best hit after consolidation");
        assert_eq!(
            family_ranges.len(),
            1,
            "expected one family for densely overlapping circular hits"
        );
    }
}
