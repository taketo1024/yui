//! Build planning shared by [`TngComplexBuilder`](super::TngComplexBuilder) and
//! [`SymTngBuilder`](super::SymTngBuilder): given the crossing units (single nodes, or τ-pairs
//! for the sym builder), decide the build order and the chunk decomposition. The builders then
//! build each chunk via a child builder and merge the reduced results.

use itertools::Itertools;
use log::info;
use rustc_hash::FxHashSet;
use yui_link::{Edge, Node};

use super::{cutwidth_after, toggle_boundary, CutOption, NodeOrder};

/// Plans the build of a set of crossing `units` (1 node plain / on-axis, 2 for a τ-pair):
/// the build order (cutwidth-greedy, or as given), cut into chunks per the [`CutOption`].
pub(crate) struct BuildPlanner<'a> {
    nodes: &'a [Node],
    units: Vec<Vec<usize>>, // node indices per unit
    cut: &'a CutOption,
    node_order: NodeOrder,
    open: FxHashSet<Edge>, // boundary of the already-built part (usually empty)
}

impl<'a> BuildPlanner<'a> {
    pub(crate) fn new(nodes: &'a [Node], units: Vec<Vec<usize>>, cut: &'a CutOption, node_order: NodeOrder, open: impl IntoIterator<Item = Edge>) -> Self {
        let open = open.into_iter().collect();
        Self { nodes, units, cut, node_order, open }
    }

    /// The planned chunks, each a node list in build order.
    pub(crate) fn plan(&self) -> Vec<Vec<Node>> {
        let k = match self.cut {
            CutOption::At(counts) => return self.at_plan(counts),
            CutOption::Auto(k) => (*k).max(1),
            CutOption::None => 1,
        };
        let (order, widths) = self.unit_order();
        let cuts = select_cuts(&widths, k - 1);
        self.segment_plan(&order, &cuts)
    }

    // Cut after the unit positions whose cumulative crossing count is closest to each requested
    // count — direct control over chunk balance (units stay whole, so counts land within ±1).
    fn at_plan(&self, counts: &[usize]) -> Vec<Vec<Node>> {
        let (order, widths) = self.unit_order();
        let cum: Vec<usize> = order.iter()
            .scan(0, |acc, &u| {
                *acc += self.units[u].len();
                Some(*acc)
            })
            .collect();

        // the final position is excluded — cutting after everything leaves an empty tail chunk.
        let cuts = counts.iter()
            .filter_map(|&c| (0..cum.len().saturating_sub(1)).min_by_key(|&p| cum[p].abs_diff(c)))
            .sorted().dedup().collect_vec();
        for (&c, &p) in counts.iter().sorted().zip(cuts.iter()) {
            info!("cut requested at {c} crossings -> position {p} ({} crossings, width {})", cum[p], widths[p]);
        }
        self.segment_plan(&order, &cuts)
    }

    // Segment the unit order after each cut position; expand each unit to its nodes.
    fn segment_plan(&self, order: &[usize], cuts: &[usize]) -> Vec<Vec<Node>> {
        let starts = std::iter::once(0).chain(cuts.iter().map(|&v| v + 1));
        let ends = cuts.iter().map(|&v| v + 1).chain(std::iter::once(order.len()));
        starts.zip(ends)
            .map(|(s, e)| order[s..e].iter()
                .flat_map(|&u| &self.units[u])
                .map(|&i| self.nodes[i].clone())
                .collect())
            .filter(|c: &Vec<Node>| !c.is_empty())
            .collect()
    }

    // The build order of the units — cutwidth-greedy or as given — with its width profile.
    fn unit_order(&self) -> (Vec<usize>, Vec<usize>) {
        match self.node_order {
            NodeOrder::MinCut => self.greedy_order(),
            NodeOrder::Given => self.given_order(),
        }
    }

    // Greedy MinCut order over the units (ties by first node index), with its width profile.
    fn greedy_order(&self) -> (Vec<usize>, Vec<usize>) {
        let mut open = self.open.clone();
        let mut remaining: Vec<usize> = (0..self.units.len()).collect();

        std::iter::from_fn(|| {
            let pos = remaining.iter()
                .position_min_by_key(|&&u| (cutwidth_after(&open, &self.unit_nodes(u)), self.units[u][0]))?;
            let u = remaining.swap_remove(pos);
            toggle_boundary(&mut open, &self.unit_nodes(u));
            Some((u, open.len()))
        }).unzip()
    }

    // The units as given, with the replayed width profile.
    fn given_order(&self) -> (Vec<usize>, Vec<usize>) {
        let mut open = self.open.clone();
        let widths = (0..self.units.len()).map(|u| {
            toggle_boundary(&mut open, &self.unit_nodes(u));
            open.len()
        }).collect();
        ((0..self.units.len()).collect(), widths)
    }

    fn unit_nodes(&self, u: usize) -> Vec<&Node> {
        self.units[u].iter().map(|&i| &self.nodes[i]).collect()
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

