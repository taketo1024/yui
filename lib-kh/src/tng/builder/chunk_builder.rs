//! Divide-and-conquer chunked build, shared by [`TngComplexBuilder`](super::TngComplexBuilder)
//! and [`SymTngBuilder`](super::SymTngBuilder): plan k chunks up front at thin cutwidth
//! interfaces, build each via a child builder, and merge the reduced chunks into the parent —
//! so the parent never materializes the full dense slice.

use itertools::Itertools;
use log::info;
use yui_link::Node;

use super::{boundary_edges, CutOption};

/// The planning view of a boundary-cutwidth profile: node-index units in build order (1 node
/// plain / on-axis, 2 for a τ-pair) and the boundary cutwidth after each unit. Cuts land
/// between units, so units stay whole.
pub(crate) struct PlanProfile {
    pub units: Vec<Vec<usize>>,
    pub widths: Vec<usize>,
}

/// A builder that can be chunk-built: planning inputs plus the child build.
pub(crate) trait Chunkable: Sized {
    type Elem;

    fn nodes(&self) -> &[Node];
    fn cut_option(&self) -> &CutOption;
    fn profile(&self) -> PlanProfile;
    fn stat(&self) -> String;

    /// A child builder over `chunk` (a sub-tangle).
    fn init_child(&self, chunk: &[Node]) -> Self;
    fn run(self) -> Self;
    fn take_elements(&mut self) -> Vec<Self::Elem>;
}

pub(crate) struct ChunkBuilder<'a, B> {
    builder: &'a B,
}

impl<'a, B: Chunkable> ChunkBuilder<'a, B> {
    pub(crate) fn new(builder: &'a B) -> Self {
        Self { builder }
    }

    // Plan the chunks and build each into a reduced sub-complex, paired with the crossings it
    // covers. Building is independent of the parent, so the caller merges them afterwards.
    pub(crate) fn build_chunks(&self) -> Vec<(Vec<Node>, B, Vec<B::Elem>)> {
        let plan = self.plan();
        info!("chunk plan: {} pieces {:?}", plan.len(),
            plan.iter().map(|c| c.len()).collect_vec());

        plan.into_iter().map(|chunk| {
            let (child, elems) = self.build_chunk(&chunk);
            (chunk, child, elems)
        }).collect()
    }

    // Partition the crossings into pieces: `Auto(k)` at the deepest cutwidth valleys of the
    // MinCut order (contiguous, thin interface), `AtCrossings` at the requested crossing counts.
    fn plan(&self) -> Vec<Vec<Node>> {
        let k = match self.builder.cut_option() {
            CutOption::AtCrossings(counts) => return self.at_crossings_plan(counts),
            CutOption::Auto(k) => (*k).max(1),
            CutOption::None => 1,
        };
        let prof = self.builder.profile();
        let cuts = select_cuts(&prof.widths, k - 1);
        self.segment_plan(&prof, &cuts)
    }

    // Cut after the unit positions whose cumulative crossing count is closest to each requested
    // count — direct control over chunk balance (units stay whole, so counts land within ±1).
    fn at_crossings_plan(&self, counts: &[usize]) -> Vec<Vec<Node>> {
        let prof = self.builder.profile();
        let cum: Vec<usize> = prof.units.iter()
            .scan(0, |acc, unit| {
                *acc += unit.len();
                Some(*acc)
            })
            .collect();

        // the final position is excluded — cutting after everything leaves an empty tail chunk.
        let cuts = counts.iter()
            .filter_map(|&c| (0..cum.len().saturating_sub(1)).min_by_key(|&p| cum[p].abs_diff(c)))
            .sorted().dedup().collect_vec();
        for (&c, &p) in counts.iter().sorted().zip(cuts.iter()) {
            info!("cut requested at {c} crossings -> position {p} ({} crossings, width {})", cum[p], prof.widths[p]);
        }
        self.segment_plan(&prof, &cuts)
    }

    // Segment the profile order after each cut position; expand each unit to its nodes.
    fn segment_plan(&self, prof: &PlanProfile, cuts: &[usize]) -> Vec<Vec<Node>> {
        let nodes = self.builder.nodes();
        let starts = std::iter::once(0).chain(cuts.iter().map(|&v| v + 1));
        let ends = cuts.iter().map(|&v| v + 1).chain(std::iter::once(prof.units.len()));
        starts.zip(ends)
            .map(|(s, e)| prof.units[s..e].iter().flatten().map(|&i| nodes[i].clone()).collect())
            .filter(|c: &Vec<Node>| !c.is_empty())
            .collect()
    }

    // Build `chunk` into a reduced sub-complex via a child builder.
    fn build_chunk(&self, chunk: &[Node]) -> (B, Vec<B::Elem>) {
        let ends = boundary_edges(&chunk.iter().collect::<Vec<_>>()).into_iter().sorted().collect_vec();
        info!("build chunk (n: {}, nb: {} {:?}): {}", chunk.len(), ends.len(), ends, chunk.iter().join(", "));

        let mut child = self.builder.init_child(chunk).run();
        info!("chunk built: {}", child.stat());
        let elems = child.take_elements();
        (child, elems)
    }
}

// ---- Auto chunk planning: pick the cut positions in a cutwidth profile ----

// Valley positions (a descent that turns back up) in `widths`, ascending. The `scan` carries a
// `descended` flag so flats are ignored — a mid-descent plateau isn't mistaken for the bottom.
fn valleys(widths: &[usize]) -> Vec<usize> {
    widths.windows(2).enumerate()
        .scan(false, |descended, (i, w)| Some(
            if w[1] < w[0] {
                *descended = true;
                None
            } else if w[1] > w[0] && *descended {
                *descended = false;
                Some(i)
            } else {
                None
            }
        ))
        .flatten()
        .collect()
}

// `n_cuts` cut positions in `widths`: take the deepest valleys first; if more cuts than valleys
// are needed, keep every valley and fill the rest with evenly spaced positions.
fn select_cuts(widths: &[usize], n_cuts: usize) -> Vec<usize> {
    let n = widths.len();
    let n_cuts = n_cuts.min(n.saturating_sub(1));
    if n_cuts == 0 {
        return vec![];
    }

    let vs = valleys(widths); // ascending positions
    if vs.len() >= n_cuts {
        // enough valleys: keep the `n_cuts` deepest, back in position order.
        return vs.into_iter()
            .sorted_by_key(|&p| widths[p])
            .take(n_cuts)
            .sorted()
            .collect();
    }

    let extra = n_cuts - vs.len();
    let even = (1..=extra).map(|j| j * (n - 1) / (extra + 1));
    vs.into_iter().chain(even).sorted().dedup().collect()
}

