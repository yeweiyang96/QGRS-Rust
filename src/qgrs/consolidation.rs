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

    let mut active: Vec<(usize, usize)> = Vec::new();
    for segment in segments {
        active.retain(|(end, _)| *end >= segment.start);
        for &(_, owner) in &active {
            dsu.union(segment.owner, owner);
        }
        active.push((segment.end, segment.owner));
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
    let mut cut_points: Vec<usize> = members.iter().map(|&index| raw_g4s[index].start).collect();
    cut_points.sort_unstable();
    cut_points.dedup();

    let mut best: Option<(usize, usize, usize)> = None;
    for cut in cut_points {
        let mut min_start = usize::MAX;
        let mut max_end = 0usize;
        for &index in members {
            let g4 = &raw_g4s[index];
            let mut start = g4.start;
            let mut end = g4.end;
            if start < cut {
                start += sequence_len;
                end += sequence_len;
            }
            min_start = min_start.min(start);
            max_end = max_end.max(end);
        }
        let span = max_end.saturating_sub(min_start);
        let replace = match best {
            None => true,
            Some((best_span, best_start, best_end)) => {
                span < best_span
                    || (span == best_span
                        && (min_start < best_start
                            || (min_start == best_start && max_end < best_end)))
            }
        };
        if replace {
            best = Some((span, min_start, max_end));
        }
    }

    let (_, mut start, mut end) = best.expect("members is non-empty");
    while start > sequence_len {
        start -= sequence_len;
        end -= sequence_len;
    }
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
