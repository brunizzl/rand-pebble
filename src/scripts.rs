use itertools::{Itertools, izip};

use crate::bool_csr::BoolCSR;

type PlotRes = Result<(), Box<dyn std::error::Error>>;

pub struct Graph {
    /// `edges[v]` lists all vertex neighboring the vertex `v`.
    edges: BoolCSR,
    /// either the number of vertices or the squareroot thereof.
    /// this value is only used for prettier printing.
    len: usize,
    /// only used for logging
    name: &'static str,
}

impl Graph {
    #[allow(dead_code)]
    pub fn new_grid(len: usize) -> Self {
        let mut edges = BoolCSR::new();
        for y in 0..len {
            for x in 0..len {
                edges.start_new_row();
                for (nx, ny) in [
                    (x, y.wrapping_sub(1)),
                    (x.wrapping_sub(1), y),
                    (x + 1, y),
                    (x, y + 1),
                ] {
                    if nx < len && ny < len {
                        edges.add_entry_in_last_row(ny * len + nx);
                    }
                }
            }
        }
        let name = "grid";
        Self { edges, len, name }
    }

    #[allow(dead_code)]
    pub fn new_torus_grid(len: usize) -> Self {
        let mut edges = BoolCSR::new();
        for y in len..(2 * len) {
            for x in len..(2 * len) {
                edges.start_new_row();
                for (nx, ny) in [(x, y - 1), (x - 1, y), (x + 1, y), (x, y + 1)] {
                    let nx = nx % len;
                    let ny = ny % len;
                    edges.add_entry_in_last_row(ny * len + nx);
                }
            }
        }
        let name = "torus";
        Self { edges, len, name }
    }

    /// in order to allow printing as a square, we construct a complete graph over
    /// `len` squared many vertices.
    #[allow(dead_code)]
    pub fn new_complete_len_times_len(len: usize) -> Self {
        let nr_vertices = len * len;
        let mut edges = BoolCSR::new();
        for v in 0..nr_vertices {
            edges.add_row((0..nr_vertices).filter(|&u| u != v));
        }
        let name = "complete";
        Self { edges, len, name }
    }

    /// unlike most other graphs, this one has only `len` many vertices.
    /// for simplicity, a minimum `len` of `2` is required.
    #[allow(dead_code)]
    pub fn new_path(len: usize) -> Self {
        assert!(len >= 2);
        let mut edges = BoolCSR::new();
        edges.add_row([1].into_iter());
        for v in 1..(len - 1) {
            edges.add_row([v - 1, v + 1].into_iter());
        }
        edges.add_row([len - 2].into_iter());

        let name = "path";
        Self { edges, len, name }
    }

    /// unlike most other graphs, this one has only `len` many vertices.
    /// there are no circles with fewer than `3` vertices.
    #[allow(dead_code)]
    pub fn new_circle(len: usize) -> Self {
        assert!(len >= 3);
        let mut edges = BoolCSR::new();
        edges.add_row([len - 1, 1].into_iter());
        for v in 1..(len - 1) {
            edges.add_row([v - 1, v + 1].into_iter());
        }
        edges.add_row([len - 2, 0].into_iter());

        let name = "circle";
        Self { edges, len, name }
    }

    pub fn nr_vertices(&self) -> usize {
        self.edges.nr_rows()
    }

    #[allow(dead_code)]
    pub fn debug_print(&self) {
        for (v, neighs) in izip!(0.., self.edges.iter_rows()) {
            if self.len < self.nr_vertices() && v % self.len == 0 {
                println!("--------------------------------");
            }
            println!("{v:>4} | {neighs:?}");
        }
    }
}

/// experiment to determine stationary distribution:
/// walk a random walk and keep track of how often each vertex was visited.
/// thus, not exactly the stationary distribution, but this log is returned.
/// because the number of steps is not predetermined, the number of walked steps is retuned alongside the graph:
/// ```
/// use rand_pebble::scripts::*;
/// let g = Graph::new_grid(4);
/// let (steps_per, steps_total) = find_stationary_distribution(&g, 0.1);
///
/// assert_eq!(steps_per.len(), g.nr_vertices());
/// assert_eq!(steps_per.iter().sum::<isize>() as usize, steps_total);
/// ```
pub fn find_stationary_distribution(graph: &Graph, eps: f64) -> (Vec<isize>, usize) {
    let mut curr_visits = vec![0isize; graph.nr_vertices()];
    let mut all_visits = vec![0isize; graph.nr_vertices()];
    let mut nr_steps = 1;
    let mut curr_pos = 0;
    let mut lcg = crate::rand_lcg::Lcg::new(0);
    loop {
        for _ in 0..nr_steps {
            curr_visits[curr_pos] += 1;
            let curr_neighs = &graph.edges[curr_pos];
            curr_pos = curr_neighs[(lcg.gen_u32() as usize) % curr_neighs.len()];
        }
        let cum_relative_diff = izip!(&curr_visits, &all_visits)
            .map(|(&a, &b)| ((a - b).abs() as f64) / (isize::max(a, b) as f64))
            .sum::<f64>();
        let avg_relative_diff = cum_relative_diff / (curr_visits.len() as f64);
        println!("steps: {nr_steps:>20} diff: {avg_relative_diff}");
        if avg_relative_diff < eps {
            break;
        }
        for (cv, av) in izip!(&mut curr_visits, &mut all_visits) {
            *av += *cv;
            *cv = 0;
        }
        nr_steps *= 2;
    }
    for (cv, av) in izip!(&curr_visits, &mut all_visits) {
        *av += cv;
    }

    (all_visits, nr_steps * 2 - 1)
}

/// prints for each vertex how often it is visited, relative to the average.
/// => above average iff shown value is greater than 1.0
#[allow(dead_code)]
fn print_stationary_distribution(eps: f64) {
    let len = FIXED_LEN;
    let graph = build_graph(len);
    let (visits, nr_steps) = find_stationary_distribution(&graph, eps);
    let avg = graph.nr_vertices() as f64 / nr_steps as f64;
    for line in visits.chunks(len) {
        for &entry in line {
            print!(" {:>.2}", entry as f64 * avg);
        }
        println!();
    }
}

fn initial_print(graph: &Graph) {
    if PRINT_PROGRESS {
        let nr_lines = graph.nr_vertices().div_ceil(graph.len);
        for _ in 0..nr_lines {
            println!();
        }
    }
}

/// print progress in process.
/// free vertices are shown as empty space, filled vertices as `'.'`
/// and vertices with markers show how many markers are held (if > 9 an `'X'` is shown.)
fn print_free_slots<'a>(
    graph: &Graph,
    free_slots: &[bool],
    markers: impl Iterator<Item = &'a usize>,
) {
    if !PRINT_PROGRESS {
        return;
    }
    use crossterm::{ExecutableCommand, cursor::MoveUp};
    let nr_lines = graph.nr_vertices().div_ceil(graph.len);
    std::io::stdout().execute(MoveUp(nr_lines as u16)).ok();

    let mut output = free_slots
        .iter()
        .map(|&free| if free { ' ' } else { '.' })
        .collect_vec();
    for &marker in markers {
        output[marker] = match output[marker] {
            ' ' | '.' => '1',
            '1' => '2',
            '2' => '3',
            '3' => '4',
            '4' => '5',
            '5' => '6',
            '6' => '7',
            '7' => '8',
            '8' => '9',
            _ => 'X',
        };
    }
    let mut line_str = String::with_capacity(2 * graph.len);
    for line in output.chunks(graph.len) {
        line_str.clear();
        line_str.extend(
            line.iter()
                .copied()
                .interleave(std::iter::repeat_n(' ', graph.len)),
        );
        println!("| {line_str}|");
    }
    std::thread::sleep(PRINT_PAUSE_TIME);
}

/// simulate parallel IDLA, returned are the lengths of each walk.
fn simulate_walk_lengths_parallel(graph: &Graph, start_vertex: usize) -> Vec<usize> {
    use rand::Rng;
    let mut rng = rand::rng();
    initial_print(graph);
    let mut active_pebbles = vec![start_vertex; graph.nr_vertices()];
    let mut new_active_pebbles = active_pebbles.clone();
    let mut walk_lengths = Vec::with_capacity(graph.nr_vertices());
    let mut free_slots = vec![true; graph.nr_vertices()];
    let mut nr_free_when_last_printed = usize::MAX;
    let mut settled_in_this_step = Vec::new();
    for step in 0.. {
        settled_in_this_step.clear();
        if active_pebbles.is_empty() {
            break;
        }
        new_active_pebbles.clear();
        for &pebble in &active_pebbles {
            if free_slots[pebble] {
                free_slots[pebble] = false;
                walk_lengths.push(step);
                settled_in_this_step.push(pebble);
            } else {
                let new_pebble = if WALK_LAZILY && rng.random_range(0..16) == 0 {
                    pebble
                } else {
                    let neighs = &graph.edges[pebble];
                    neighs[rng.random_range(0..neighs.len())]
                };
                new_active_pebbles.push(new_pebble);
            }
        }
        std::mem::swap(&mut active_pebbles, &mut new_active_pebbles);

        if PRINT_PROGRESS
            && (PRINT_EVERY_FRAME || nr_free_when_last_printed != active_pebbles.len())
        {
            nr_free_when_last_printed = active_pebbles.len();

            print_free_slots(graph, &free_slots, active_pebbles.iter());
        }
    }

    walk_lengths
}

/// a random walk is a vec of `Step`'s
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Step {
    /// (original) index of pebble.
    /// this is mostly metadata.
    pebble: usize,
    /// vertex index
    pos: usize,
}

/// simulate parallel IDLA, returned are the walks of all pebbles.
fn simulate_walks_parallel(graph: &Graph, start_vertex: usize) -> Vec<Vec<Step>> {
    use rand::Rng;
    let mut rng = rand::rng();
    initial_print(graph);
    let mut nr_active_pebbles = graph.nr_vertices();
    let mut active_pebbles = vec![true; graph.nr_vertices()];
    let mut walks = (0..graph.nr_vertices())
        .map(|pebble| {
            vec![Step {
                pebble,
                pos: start_vertex,
            }]
        })
        .collect_vec();
    let mut free_slots = vec![true; graph.nr_vertices()];

    let mut nr_free_when_last_printed = nr_active_pebbles + 1;
    let mut curr_pebble_positions = Vec::new();
    while nr_active_pebbles > 0 {
        for (i, active, walk) in izip!(0.., &mut active_pebbles, &mut walks) {
            if !*active {
                continue;
            }
            let Step {
                pebble,
                pos: pebble_pos,
            } = *walk.last().unwrap();
            debug_assert_eq!(pebble, i);
            curr_pebble_positions.push(pebble_pos);
            if free_slots[pebble_pos] {
                free_slots[pebble_pos] = false;
                *active = false;
                nr_active_pebbles -= 1;
            } else {
                let new_pos = if WALK_LAZILY && rng.random_range(0..16) == 0 {
                    pebble_pos
                } else {
                    let neighs = &graph.edges[pebble_pos];
                    neighs[rng.random_range(0..neighs.len())]
                };
                walk.push(Step {
                    pebble,
                    pos: new_pos,
                });
            }
        }

        if PRINT_PROGRESS && (PRINT_EVERY_FRAME || nr_free_when_last_printed != nr_active_pebbles) {
            nr_free_when_last_printed = nr_active_pebbles;

            print_free_slots(graph, &free_slots, curr_pebble_positions.iter());
        }
        curr_pebble_positions.clear();
    }

    walks
}

/// simulate sequential IDLA, returned are the lengths of each walk.
fn simulate_walk_lengths_serial(graph: &Graph, start_vertex: usize) -> Vec<usize> {
    use rand::Rng;
    let mut rng = rand::rng();
    initial_print(graph);
    let mut walk_lengths = Vec::with_capacity(graph.nr_vertices());
    let mut free_slots = vec![true; graph.nr_vertices()];
    let mut path = Vec::new();
    for _ in 0..graph.nr_vertices() {
        let mut curr = start_vertex;
        path.push(curr);
        let mut nr_steps = 0;
        while !free_slots[curr] {
            curr = if !WALK_LAZILY || rng.random_range(0..16) == 0 {
                let curr_neighs = &graph.edges[curr];
                curr_neighs[rng.random_range(0..curr_neighs.len())]
            } else {
                curr
            };
            path.push(curr);
            nr_steps += 1;
        }
        free_slots[curr] = false;
        walk_lengths.push(nr_steps);

        print_free_slots(graph, &free_slots, path.iter());
        path.clear();
    }

    walk_lengths
}

/// simulate sequential IDLA, returned are the walks of each pebble.
fn simulate_walks_serial(graph: &Graph, start_vertex: usize) -> Vec<Vec<Step>> {
    use rand::Rng;
    let mut rng = rand::rng();
    initial_print(graph);
    let mut walks = Vec::with_capacity(graph.nr_vertices());
    let mut free_slots = vec![true; graph.nr_vertices()];
    for pebble in 0..graph.nr_vertices() {
        let mut path = Vec::new();
        let mut pos = start_vertex;
        path.push(Step { pebble, pos });
        while !free_slots[pos] {
            pos = if !WALK_LAZILY || rng.random_range(0..16) == 0 {
                let curr_neighs = &graph.edges[pos];
                curr_neighs[rng.random_range(0..curr_neighs.len())]
            } else {
                pos
            };
            path.push(Step { pebble, pos });
        }
        free_slots[pos] = false;

        print_free_slots(graph, &free_slots, path.iter().map(|s| &s.pos));
        walks.push(path);
    }

    walks
}

#[allow(dead_code)]
fn test_walk_lengths_once() {
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = build_graph(len);

    let mut serial_walk_lengths = simulate_walk_lengths_serial(&graph, start);
    println!("serial walk lengths (unsorted): {serial_walk_lengths:?}");
    serial_walk_lengths.sort();
    println!("serial walk lengths (sorted): {serial_walk_lengths:?}");

    let parallel_walk_lengths = simulate_walk_lengths_parallel(&graph, start);
    println!("parallel walk lengths: {parallel_walk_lengths:?}");
}

/// mostly made obsolete by [`plot_walks`], but will actually show which vertices where taken.
#[allow(dead_code)]
fn print_walks(walks: &[Vec<Step>]) {
    if !PRINT_WALKS {
        return;
    }
    for (i, walk) in izip!(0.., walks) {
        print!("{i:>4} |");
        for Step { pebble: _, pos } in walk {
            print!(" {pos:>3}");
        }
        println!();
    }
}

#[allow(dead_code)]
fn plot_walks(plot_name: &str, walks: &[Vec<Step>]) -> PlotRes {
    use plotters::prelude::*;
    if !PRINT_WALKS {
        return Ok(());
    }
    let nr_vertices = walks.len();
    let longest = max_walk_len(walks);

    let file_name = format!("{plot_name}.png");
    let root = BitMapBackend::new(&file_name, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(plot_name, ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(30)
        .build_cartesian_2d(0..longest, 0..nr_vertices)?;
    chart.configure_mesh().draw()?;

    for (y, walk) in izip!(1.., walks) {
        let draw_y = nr_vertices - y;
        chart.draw_series(izip!(0.., walk).map(|(x, Step { pebble, pos: _ })| {
            let pebble_val = *pebble as f64 / nr_vertices as f64;
            //let color = plotters::style::HSLColor(0.4, 1.0, pebble_val * 0.7 + 0.15);
            let color = plotters::style::HSLColor(pebble_val * 0.8 + 0.1, 1.0, 0.5);
            Rectangle::new([(x, draw_y), (x + 1, draw_y + 1)], color.filled())
        }))?;
    }

    root.present()?;

    Ok(())
}

/// returns how many steps each pebble took and which pebble took most.
/// this always considers the original pebble index, thus may be different from
/// `walks.iter().map(Vec::len)`
fn pebble_frequencies(walks: &[Vec<Step>]) -> (Vec<usize>, usize) {
    let mut frequencies = vec![0; walks.len()];
    for walk in walks {
        for Step { pebble, pos: _ } in walk {
            frequencies[*pebble] += 1;
        }
    }
    let max_pebble = {
        let mut best_i = 0;
        let mut best_val = 0;
        for (i, &val) in izip!(0.., &frequencies) {
            if val > best_val {
                best_i = i;
                best_val = val;
            }
        }
        best_i
    };
    (frequencies, max_pebble)
}

fn max_walk_len(walks: &[Vec<Step>]) -> usize {
    walks.iter().map(Vec::len).max().unwrap_or(0)
}

fn plot_walks_variations(plot_name: &str, walks: &[Vec<Step>]) -> PlotRes {
    plot_walks(plot_name, walks)?;
    let mut walks_clone = Vec::from(walks);
    walks_clone.sort_by_key(Vec::len);
    let sorted_name = format!("{plot_name}-sorted");
    plot_walks(&sorted_name, &walks_clone)?;

    let nr_vertices = walks.len();
    let (_, max_pebble) = pebble_frequencies(walks);
    walks_clone = Vec::from(walks);
    for walk in &mut walks_clone {
        for Step { pebble, pos: _ } in walk {
            // abuse of pebble meaning: we print different pebbles as different colors,
            // if only the max_pebble should have a special color, all others are the same.
            *pebble = if *pebble == max_pebble {
                nr_vertices
            } else {
                (nr_vertices * 55) / 100
            };
        }
    }
    let max_name = format!("{plot_name}-max");
    plot_walks(&max_name, &walks_clone)?;

    Ok(())
}

#[allow(dead_code)]
fn test_walks_once() {
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = build_graph(len);

    println!("simulate serial...");
    let serial_walks = simulate_walks_serial(&graph, start);
    plot_walks_variations("serial", &serial_walks).ok();
    //print_walks(&serial_walks);

    println!("simulate parallel...");
    let parallel_walks = simulate_walks_parallel(&graph, start);
    plot_walks_variations("parallel", &parallel_walks).ok();
    //print_walks(&parallel_walks);
}

/// algorithm from paper
fn transform_walks_serial_to_parallel(walks: &mut [Vec<Step>]) {
    let nr_vertices = walks.len();
    // if `v` is a vertex, then `walks[ends_in_vertex[v]].last() == Some(&v)` must hold.
    // this cashes information, that we would otherwise compute in O(nr_vertices),
    // whenever a cut-and-paste occurs.
    let mut ends_in_vertex = vec![usize::MAX; nr_vertices];
    for (i, walk) in izip!(0.., walks.iter()) {
        let last_vertex = walk.last().unwrap().pos;
        debug_assert_eq!(ends_in_vertex[last_vertex], usize::MAX);
        ends_in_vertex[last_vertex] = i;
    }
    debug_assert!(!ends_in_vertex.contains(&usize::MAX));

    let mut nr_vertices_visited = 0;
    let mut vertices_visited = vec![false; nr_vertices];
    let mut t = 0;
    while nr_vertices_visited < nr_vertices {
        for i_donor in 0..walks.len() {
            if walks[i_donor].len() > t {
                let Step { pebble: _, pos } = walks[i_donor][t];
                if !vertices_visited[pos] {
                    vertices_visited[pos] = true;
                    nr_vertices_visited += 1;
                    let i_reciever = ends_in_vertex[pos];
                    debug_assert!(walks[i_reciever].last().is_some_and(|s| s.pos == pos));
                    if i_reciever != i_donor {
                        let mut donor_walk = std::mem::take(&mut walks[i_donor]);
                        ends_in_vertex.swap(pos, donor_walk.last().unwrap().pos);
                        let to_transfer = &donor_walk[(t + 1)..];
                        walks[i_reciever].extend_from_slice(to_transfer);
                        donor_walk.truncate(t + 1);
                        walks[i_donor] = donor_walk;
                    }
                }
            }
        }
        t += 1;
    }
}

/// algorithm from paper
fn transform_walks_parallel_to_serial(walks: &mut [Vec<Step>]) {
    let nr_vertices = walks.len();
    // if `v` is a vertex, then `walks[ends_in_vertex[v]].last() == Some(&v)` must hold.
    // this cashes information, that we would otherwise compute in O(nr_vertices),
    // whenever a cut-and-paste occurs.
    let mut ends_in_vertex = vec![usize::MAX; nr_vertices];
    for (i, walk) in izip!(0.., walks.iter()) {
        let last_vertex = walk.last().unwrap().pos;
        debug_assert_eq!(ends_in_vertex[last_vertex], usize::MAX);
        ends_in_vertex[last_vertex] = i;
    }
    debug_assert!(!ends_in_vertex.contains(&usize::MAX));

    let mut vertices_visited = vec![false; nr_vertices];
    for i_donor in 0..walks.len() {
        let mut t = 0;
        while t < walks[i_donor].len() {
            let Step { pebble: _, pos } = walks[i_donor][t];
            if !vertices_visited[pos] {
                vertices_visited[pos] = true;
                let i_reciever = ends_in_vertex[pos];
                debug_assert!(walks[i_reciever].last().is_some_and(|s| s.pos == pos));
                if i_reciever != i_donor {
                    let mut donor_walk = std::mem::take(&mut walks[i_donor]);
                    ends_in_vertex.swap(pos, donor_walk.last().unwrap().pos);
                    let to_transfer = &donor_walk[(t + 1)..];
                    walks[i_reciever].extend_from_slice(to_transfer);
                    donor_walk.truncate(t + 1);
                    walks[i_donor] = donor_walk;
                }
            }
            t += 1;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn transforms_are_inverses() {
        for len in 10..20 {
            let graph = build_graph(len);
            {
                let serial = simulate_walks_serial(&graph, 0);
                let mut round_trip = serial.clone();
                transform_walks_serial_to_parallel(&mut round_trip);
                transform_walks_parallel_to_serial(&mut round_trip);
                assert_eq!(serial, round_trip);
            }
            {
                let parallel = simulate_walks_parallel(&graph, 0);
                let mut round_trip = parallel.clone();
                transform_walks_parallel_to_serial(&mut round_trip);
                transform_walks_serial_to_parallel(&mut round_trip);
                assert_eq!(parallel, round_trip);
            }
        }
    }
}

#[allow(dead_code)]
fn test_transformed_walk_lengths_serial() {
    println!("compute serial, then transform to parallel");
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = build_graph(len);

    let mut walks = simulate_walks_serial(&graph, start);
    //print_walks(&walks);
    plot_walks_variations("serial", &walks).ok();
    let max_serial = max_walk_len(&walks);
    println!("max serial len = {max_serial}");

    transform_walks_serial_to_parallel(&mut walks);
    //print_walks(&walks);
    plot_walks_variations("parallel", &walks).ok();
    let max_parallel = max_walk_len(&walks);
    println!("max parallel len = {max_parallel}");

    println!("ratio: {}", max_serial as f64 / max_parallel as f64);
}

#[allow(dead_code)]
fn test_transformed_walk_lengths_parallel() {
    println!("compute parallel, then transform to serial");
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = build_graph(len);

    let mut walks = simulate_walks_parallel(&graph, start);
    //print_walks(&walks);
    plot_walks_variations("parallel", &walks).ok();
    let max_parallel = max_walk_len(&walks);
    println!("max parallel len = {max_parallel}");

    transform_walks_parallel_to_serial(&mut walks);
    //print_walks(&walks);
    plot_walks_variations("serial", &walks).ok();
    let max_serial = max_walk_len(&walks);
    println!("max serial len = {max_serial}");

    println!("ratio: {}", max_serial as f64 / max_parallel as f64);
}

/// examine for the graph family passed as argument, which values t_par and t_ser aproxximately have.
#[allow(dead_code)]
fn find_max_expected_serial_parallel() {
    let mut len_f = 10.0;

    while len_f < 200.0 {
        use rayon::prelude::*;
        let len = len_f as usize;
        len_f *= 1.1;

        let range = [(); SAMPLES_PER_EXPERIMENT];
        let graph = build_graph(len);

        let last_v_0 = if START_AT_0 {
            0
        } else {
            println!("-------------------------------------------------");
            graph.nr_vertices().div_ceil(2)
        };
        for v_0 in 0..=last_v_0 {
            let serial_sum: usize = range
                .par_iter()
                .map(|_| {
                    let mut serial = simulate_walk_lengths_serial(&graph, v_0);
                    serial.sort();
                    *serial.last().unwrap()
                })
                .sum();

            let parallel_sum: usize = range
                .par_iter()
                .map(|_| {
                    let parallel = simulate_walk_lengths_parallel(&graph, v_0);
                    *parallel.last().unwrap()
                })
                .sum();

            let exp_serial = serial_sum as f64 / SAMPLES_PER_EXPERIMENT as f64;
            let exp_parallel = parallel_sum as f64 / SAMPLES_PER_EXPERIMENT as f64;
            let ratio = serial_sum as f64 / parallel_sum as f64;
            println!(
                "len ={len:>4}, \
                t_par ={exp_parallel:>11.2}, \
                t_ser ={exp_serial:>11.2}, \
                ratio ={ratio:>6.3}, \
                v_0 = {v_0}"
            );
        }
    }
}

/// simulate the parallel process, transform the walks to serial processes.
/// this function should be examined together with [`find_max_expected_serial_parallel`],
/// where both processes are simulated.
#[allow(dead_code)]
fn compare_expected_transformed_serial_parallel() {
    let mut len_f = 10.0;

    while len_f < 200.0 {
        use rayon::prelude::*;
        let len = len_f as usize;
        len_f *= 1.1;

        let range = [(); SAMPLES_PER_EXPERIMENT];
        let graph = build_graph(len);

        let last_v_0 = if START_AT_0 {
            0
        } else {
            println!("-------------------------------------------------");
            graph.nr_vertices().div_ceil(2)
        };
        for v_0 in 0..=last_v_0 {
            let res = range
                .par_iter()
                .map(|_| {
                    let mut walks = simulate_walks_parallel(&graph, v_0);
                    let max_parallel = max_walk_len(&walks);
                    let slowest_pebble = walks.iter().map(Vec::len).position_max().unwrap();
                    transform_walks_parallel_to_serial(&mut walks);
                    let max_serial = max_walk_len(&walks);
                    let nr_cuts = walks
                        .iter()
                        .map(|w| {
                            // if possible, it should be extremly uncommon for the walk of single new pebble
                            // to be made up of multiple seperated pieces of the walk of the same old pebble.
                            // yet, we also handle this case here.
                            w.iter()
                                .map(|s| s.pebble)
                                .dedup()
                                .filter(|&p| p == slowest_pebble)
                                .count()
                        })
                        .sum::<usize>()
                        - 1; // nr cuts is nr segments minus one
                    (max_parallel, max_serial, nr_cuts)
                })
                .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
            let avg_max_parallel = res.0 as f64 / SAMPLES_PER_EXPERIMENT as f64;
            let avg_max_serial = res.1 as f64 / SAMPLES_PER_EXPERIMENT as f64;
            let ratio = avg_max_serial / avg_max_parallel;
            let avg_nr_cuts = res.2 as f64 / SAMPLES_PER_EXPERIMENT as f64;
            println!(
                "len ={len:>4}, \
                t_par ={avg_max_parallel:>11.2}, \
                (transformed) t_ser ={avg_max_serial:>11.2}, \
                cuts ={avg_nr_cuts:>6.3}, \
                ratio ={ratio:>6.3}, \
                v_0 = {v_0}"
            );
        }
    }
}

fn plot_functions(name: &str, xs: &[usize], yss: &[(Vec<f64>, String)]) -> PlotRes {
    use plotters::prelude::*;
    let float_ord = |a: &&f64, b: &&f64| {
        if let Some(ord) = a.partial_cmp(b) {
            return ord;
        }
        std::cmp::Ordering::Equal
    };
    let x_min = *xs.iter().min().unwrap_or(&0);
    let x_max = *xs.iter().max().unwrap_or(&0);
    let y_min = *yss
        .iter()
        .flat_map(|ys| ys.0.iter())
        .min_by(float_ord)
        .unwrap_or(&0.0);
    let y_max = *yss
        .iter()
        .flat_map(|ys| ys.0.iter())
        .max_by(float_ord)
        .unwrap_or(&0.0);

    let plot_name = format!("{name}.png");
    let root = BitMapBackend::new(&plot_name, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE).ok();
    let mut chart = ChartBuilder::on(&root)
        .caption(name, ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;

    let mut draw_legend = false;
    const COLORS: [RGBColor; 6] = [RED, BLUE, YELLOW, CYAN, GREEN, MAGENTA];
    for ((ys, name), color) in izip!(yss, COLORS.into_iter().cycle()) {
        let series =
            chart.draw_series(LineSeries::new(izip!(xs, ys).map(|(&x, &y)| (x, y)), color))?;
        if !name.is_empty() {
            series.label(name);
            draw_legend = true;
        }
    }

    if draw_legend {
        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.8))
            .border_style(BLACK)
            .draw()?;
    }

    root.present()?;

    Ok(())
}

/// when simulating parallel IDLA and transforming the walks to serial,
/// if one finds a non-zero length segment of the originally longest walk in a transformed row,
/// what is the expected length of such a segment?
#[allow(dead_code)]
fn longest_walk_transformed_distribution() {
    let len = FIXED_LEN;
    let graph = build_graph(len);

    const NR_BUCKETS: usize = 400;
    // sum lengths of segments in `.0` and count how often a segment was found in `.1`
    let mut histogram = [(0, 0); NR_BUCKETS];

    for experiment_nr in 0..SAMPLES_PER_EXPERIMENT {
        if experiment_nr % (SAMPLES_PER_EXPERIMENT / 20) == 0 {
            println!("compute experiment {experiment_nr}");
        }
        let mut walks = simulate_walks_parallel(&graph, 0);
        let slowest_pebble = walks.iter().map(Vec::len).position_max().unwrap();
        transform_walks_parallel_to_serial(&mut walks);
        for (row, walk) in izip!(0.., walks) {
            let nr_longest = walk.iter().filter(|s| s.pebble == slowest_pebble).count();
            if nr_longest > 0 {
                histogram[(row * NR_BUCKETS) / graph.nr_vertices()].0 += nr_longest;
                histogram[(row * NR_BUCKETS) / graph.nr_vertices()].1 += 1;
            }
        }
    }

    let bucket_size = graph.nr_vertices() as f64 / NR_BUCKETS as f64;
    let corrected_histogram = histogram
        .iter()
        .map(|&(length, frequency)| length as f64 / (frequency as f64 * bucket_size + 1e-50))
        .collect_vec();
    println!("corrected histogram: {corrected_histogram:?}");
    let xs = (0..NR_BUCKETS)
        .map(|x| x * bucket_size as usize)
        .collect_vec();
    let yss = [(corrected_histogram, String::new())];
    plot_functions("histogram-longest-walk", &xs, &yss).ok();
}

#[allow(dead_code)]
fn walk_lengths_distributions() {
    use rayon::prelude::*;

    let len = FIXED_LEN;
    let range = [(); SAMPLES_PER_EXPERIMENT];
    let graph = build_graph(len);
    let init_vals = || [vec![0; graph.nr_vertices()], vec![0; graph.nr_vertices()]];
    let [serial_res, parallel_res] = range
        .par_iter()
        .map(|_| {
            let serial = simulate_walk_lengths_serial(&graph, 0);
            let parallel = simulate_walk_lengths_parallel(&graph, 0);
            [serial, parallel]
        })
        .reduce(init_vals, |mut acc, new| {
            for (acc_vec, new_vec) in izip!(&mut acc, new) {
                for (acc_val, new_val) in izip!(acc_vec, new_vec) {
                    *acc_val += new_val;
                }
            }
            acc
        });
    let serial_f64 = serial_res.into_iter().map(|y| y as f64).collect_vec();
    let parallel_f64 = parallel_res.into_iter().map(|y| y as f64).collect_vec();

    let xs = (0..graph.nr_vertices()).collect_vec();
    let yss = [
        (serial_f64, "serial".to_string()),
        (parallel_f64, "parallel".to_string()),
    ];
    plot_functions("walk-lengths-distribution", &xs, &yss).ok();
}

fn print_config_info() {
    let graph = build_graph(FIXED_LEN);
    println!("graph family: {}", graph.name);
    println!("samples per experiment: {SAMPLES_PER_EXPERIMENT}");
    println!("walk lazily: {WALK_LAZILY}");
}

// ------------------------------------------------------------------
//                            config
// ------------------------------------------------------------------

/// which family of graphs should be examined
fn build_graph(len: usize) -> Graph {
    //Graph::new_grid(len)
    //Graph::new_torus_grid(len)
    //Graph::new_complete_len_times_len(len)
    Graph::new_path(len)
    //Graph::new_circle(len)
}
/// if an expected value is estimated by averaging over random samples, this is the number of samples taken.
const SAMPLES_PER_EXPERIMENT: usize = 1_000_000;
/// fix vertex `0` as starting position, not test (half of) all vertices.
/// note: testing half of the vertices suffice in the other case, because all examined graphs have enough symmetry.
const START_AT_0: bool = false;
/// if some function assembles only a single graph, this is the used len.
const FIXED_LEN: usize = 40;
/// after finishing the simulation, should all walks be printed / plotted
const PRINT_WALKS: bool = true;
/// during a simulation, should the screen show the progress (this is a VERY significant slowdown)
const PRINT_PROGRESS: bool = !cfg!(test) && false;
/// currently set up with 1/16 chance to not move.
const WALK_LAZILY: bool = false;
/// for parallel IDLA and with PRINT_PROGRESS enabled, should every time step be printed
/// or should ony steps where a pebble finds it's hole be printed?
const PRINT_EVERY_FRAME: bool = true;
/// pause computation after PRINT_PROGRESS happened, so one can actually admire the picture.
const PRINT_PAUSE_TIME: std::time::Duration = std::time::Duration::from_millis(100);

pub fn scripts_main() {
    print_config_info();

    print_stationary_distribution(0.00000001);
    //test_walk_lengths_once();
    //test_walks_once();
    //test_transformed_walk_lengths_serial();
    //test_transformed_walk_lengths_parallel();
    //longest_walk_transformed_distribution();
    //walk_lengths_distributions();

    //find_max_expected_serial_parallel();
    //compare_expected_transformed_serial_parallel();
}
