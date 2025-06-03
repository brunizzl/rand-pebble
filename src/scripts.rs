use itertools::{Itertools, izip};
use rand::Rng;

use crate::bool_csr::BoolCSR;

pub struct Graph {
    edges: BoolCSR,
}

impl Graph {
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
        Self { edges }
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
        Self { edges }
    }

    pub fn nr_vertices(&self) -> usize {
        self.edges.nr_rows()
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
            curr_pos = curr_neighs[(lcg.next() as usize) % curr_neighs.len()];
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
fn print_stationary_distribution() {
    let len = FIXED_LEN;
    let graph = Graph::new_grid(len);
    let (visits, nr_steps) = find_stationary_distribution(&graph, 0.0002);
    for y in 0..len {
        for x in 0..len {
            let val = visits[y * len + x];
            print!(
                " {:>.2}",
                (val as usize * len * len) as f64 / (nr_steps as f64)
            );
        }
        println!();
    }
}

/// print progress in process.
/// free vertices are shown as empty space, filled vertices as `'.'`
/// and vertices with markers show how many markers are held (if > 9 an `'X'` is shown.)
fn print_free_slots<'a>(free_slots: &[bool], len: usize, markers: impl Iterator<Item = &'a usize>) {
    use crossterm::{ExecutableCommand, cursor::MoveUp};
    std::io::stdout().execute(MoveUp(len as u16)).ok();

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
    let mut line_str = String::with_capacity(2 * len);
    for line in output.chunks(len) {
        line_str.clear();
        line_str.extend(
            line.iter()
                .copied()
                .interleave(std::iter::repeat_n(' ', len)),
        );
        println!("| {line_str}|");
    }
    std::thread::sleep(PRINT_PAUSE_TIME);
}

fn simulate_walk_lengths_parallel(graph: &Graph, start_vertex: usize) -> Vec<usize> {
    let mut rng = rand::rng();
    let len = f64::sqrt(graph.nr_vertices() as f64) as usize;
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
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

            print_free_slots(&free_slots, len, active_pebbles.iter());
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

fn simulate_walks_parallel(graph: &Graph, start_vertex: usize) -> Vec<Vec<Step>> {
    let mut rng = rand::rng();
    let len = f64::sqrt(graph.nr_vertices() as f64) as usize;
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
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
    for _ in 0.. {
        if nr_active_pebbles == 0 {
            break;
        }
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

        if PRINT_PROGRESS
            && (PRINT_EVERY_FRAME || nr_free_when_last_printed != active_pebbles.len())
        {
            nr_free_when_last_printed = active_pebbles.len();

            print_free_slots(&free_slots, len, curr_pebble_positions.iter());
        }
        curr_pebble_positions.clear();
    }

    walks
}

fn simulate_walk_lengths_serial(graph: &Graph, start_vertex: usize) -> Vec<usize> {
    let mut rng = rand::rng();
    let len = f64::sqrt(graph.nr_vertices() as f64) as usize;
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
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

        if PRINT_PROGRESS {
            print_free_slots(&free_slots, len, path.iter());
        }
        path.clear();
    }

    walk_lengths
}

fn simulate_walks_serial(graph: &Graph, start_vertex: usize) -> Vec<Vec<Step>> {
    let mut rng = rand::rng();
    let len = f64::sqrt(graph.nr_vertices() as f64) as usize;
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
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

        if PRINT_PROGRESS {
            print_free_slots(&free_slots, len, path.iter().map(|s| &s.pos));
        }
        walks.push(path);
    }

    walks
}

#[allow(dead_code)]
fn test_walk_lengths_once() {
    let len = FIXED_LEN;
    let start = len * (len / 2) + (len / 2);
    let graph = Graph::new_grid(len);

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
fn plot_walks(plot_name: &str, walks: &[Vec<Step>]) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;
    if !PRINT_WALKS {
        return Ok(());
    }
    let nr_vertices = walks.len();
    let longest = walks.iter().map(|w| w.len()).max().unwrap();
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

/// returns many steps each pebble took and which pebble took most.
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

fn plot_walks_variations(
    plot_name: &str,
    walks: &[Vec<Step>],
) -> Result<(), Box<dyn std::error::Error>> {
    plot_walks(plot_name, walks)?;
    let mut walks_clone = Vec::from(walks);
    walks_clone.sort_by_key(Vec::len);
    let sorted_name = format!("{plot_name}-sorted");
    plot_walks(&sorted_name, &walks_clone)?;

    let nr_vertices = walks.len();
    let (_, max_pebble) = pebble_frequencies(walks);
    for walk in &mut walks_clone {
        for Step { pebble, pos: _ } in walk {
            // abuse of pebble meaning: we print different pebbles as different colors,
            // if only the max_pebble should have a special color, all others are the same.
            *pebble = if *pebble == max_pebble {
                nr_vertices
            } else {
                nr_vertices / 10
            };
        }
    }
    let max_name = format!("{plot_name}-sorted-max");
    plot_walks(&max_name, &walks_clone)?;

    Ok(())
}

#[allow(dead_code)]
fn test_walks_once() {
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = Graph::new_grid(len);

    let serial_walks = simulate_walks_serial(&graph, start);
    plot_walks("serial", &serial_walks).ok();
    //print_walks(&serial_walks);

    let parallel_walks = simulate_walks_parallel(&graph, start);
    plot_walks("parallel", &parallel_walks).ok();
    //print_walks(&parallel_walks);
}

/// algorithm from paper
fn transform_walks_serial_to_parallel(graph: &Graph, walks: &mut [Vec<Step>]) {
    let mut nr_vertices_visited = 0;
    let mut vertices_visited = vec![false; graph.nr_vertices()];
    let mut t = 0;
    while nr_vertices_visited < graph.nr_vertices() {
        for i in 0..graph.nr_vertices() {
            if walks[i].len() > t {
                let Step { pebble: _, pos } = walks[i][t];
                if !vertices_visited[pos] {
                    vertices_visited[pos] = true;
                    nr_vertices_visited += 1;
                    let paste_pos = walks
                        .iter()
                        .position(|w| w.last().is_some_and(|s| s.pos == pos))
                        .unwrap();
                    if paste_pos != i {
                        let mut walk_i = std::mem::take(&mut walks[i]);
                        let to_copy = &walk_i[(t + 1)..];
                        walks[paste_pos].extend_from_slice(to_copy);
                        walk_i.truncate(t + 1);
                        walks[i] = walk_i;
                    }
                }
            }
        }
        t += 1;
    }
}

/// algorithm from paper
fn transform_walks_parallel_to_serial(graph: &Graph, walks: &mut [Vec<Step>]) {
    let mut vertices_visited = vec![false; graph.nr_vertices()];
    for i in 0..graph.nr_vertices() {
        let mut t = 0;
        while t < walks[i].len() {
            let Step { pebble: _, pos } = walks[i][t];
            if !vertices_visited[pos] {
                vertices_visited[pos] = true;
                let paste_pos = walks
                    .iter()
                    .position(|w| w.last().is_some_and(|s| s.pos == pos))
                    .unwrap();
                if paste_pos != i {
                    let mut walk_i = std::mem::take(&mut walks[i]);
                    let to_copy = &walk_i[(t + 1)..];
                    walks[paste_pos].extend_from_slice(to_copy);
                    walk_i.truncate(t + 1);
                    walks[i] = walk_i;
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
            let graph = Graph::new_grid(len);
            {
                let serial = simulate_walks_serial(&graph, 0);
                let mut round_trip = serial.clone();
                transform_walks_serial_to_parallel(&graph, &mut round_trip);
                transform_walks_parallel_to_serial(&graph, &mut round_trip);
                assert_eq!(serial, round_trip);
            }
            {
                let parallel = simulate_walks_parallel(&graph, 0);
                let mut round_trip = parallel.clone();
                transform_walks_parallel_to_serial(&graph, &mut round_trip);
                transform_walks_serial_to_parallel(&graph, &mut round_trip);
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
    let graph = Graph::new_torus_grid(len);

    let mut walks = simulate_walks_serial(&graph, start);
    //print_walks(&walks);
    plot_walks_variations("serial", &walks).ok();
    let max_serial = max_walk_len(&walks);
    println!("max serial len = {max_serial}");

    transform_walks_serial_to_parallel(&graph, &mut walks);
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
    let graph = Graph::new_grid(len);

    let mut walks = simulate_walks_parallel(&graph, start);
    //print_walks(&walks);
    plot_walks_variations("parallel", &walks).ok();
    let max_parallel = max_walk_len(&walks);
    println!("max parallel len = {max_parallel}");

    transform_walks_parallel_to_serial(&graph, &mut walks);
    //print_walks(&walks);
    plot_walks_variations("serial", &walks).ok();
    let max_serial = max_walk_len(&walks);
    println!("max serial len = {max_serial}");

    println!("ratio: {}", max_serial as f64 / max_parallel as f64);
}

#[allow(dead_code)]
fn find_max_expected_serial_parallel() {
    let mut len_f = 10.0;

    while len_f < 200.0 {
        use rayon::prelude::*;
        let len = len_f as usize;
        len_f *= 1.1;

        const EXPERIMENTS_PER_LEN: usize = 10_000;
        let range = [(); EXPERIMENTS_PER_LEN];
        let graph = Graph::new_grid(len);

        let serial_sum: usize = range
            .par_iter()
            .map(|_| {
                let mut serial = simulate_walk_lengths_serial(&graph, 0);
                serial.sort();
                *serial.last().unwrap()
            })
            .sum();

        let parallel_sum: usize = range
            .par_iter()
            .map(|_| {
                let parallel = simulate_walk_lengths_parallel(&graph, 0);
                *parallel.last().unwrap()
            })
            .sum();

        let exp_serial = serial_sum as f64 / EXPERIMENTS_PER_LEN as f64;
        let exp_parallel = parallel_sum as f64 / EXPERIMENTS_PER_LEN as f64;
        let ratio = serial_sum as f64 / parallel_sum as f64;
        println!(
            "len = {len:>3}, t_ser = {exp_serial:>10.2}, t_par = {exp_parallel:>10.2}, ratio = {ratio}"
        );
    }
}

#[allow(dead_code)]
fn compare_expected_transformed_serial_parallel() {
    let mut len_f = 10.0;

    while len_f < 200.0 {
        use rayon::prelude::*;
        let len = len_f as usize;
        len_f *= 1.1;

        const EXPERIMENTS_PER_LEN: usize = 10_000;
        let range = [(); EXPERIMENTS_PER_LEN];
        let graph = Graph::new_grid(len);

        let res = range
            .par_iter()
            .map(|_| {
                let mut walks = simulate_walks_parallel(&graph, 0);
                let max_parallel = max_walk_len(&walks);
                let slowest_pebble = walks.iter().map(Vec::len).position_max().unwrap();
                transform_walks_parallel_to_serial(&graph, &mut walks);
                let max_serial = max_walk_len(&walks);
                let nr_cuts = walks
                    .iter()
                    .filter(|w| w.iter().any(|s| s.pebble == slowest_pebble))
                    .count()
                    - 1;
                (max_parallel, max_serial, nr_cuts)
            })
            .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
        let avg_max_parallel = res.0 as f64 / EXPERIMENTS_PER_LEN as f64;
        let avg_max_serial = res.1 as f64 / EXPERIMENTS_PER_LEN as f64;
        let ratio = avg_max_serial / avg_max_parallel;
        let avg_nr_cuts = res.2 as f64 / EXPERIMENTS_PER_LEN as f64;
        println!(
            "len = {len:>3}, \
            t_par = {avg_max_parallel:>10.2}, \
            (transformed) t_ser = {avg_max_serial:>10.2}, \
            cuts = {avg_nr_cuts:>3.2}, \
            ratio = {ratio}"
        );
    }
}

const FIXED_LEN: usize = 15;
const PRINT_WALKS: bool = true;
const PRINT_PROGRESS: bool = false;
const WALK_LAZILY: bool = false;
const PRINT_EVERY_FRAME: bool = true;
const PRINT_PAUSE_TIME: std::time::Duration = std::time::Duration::from_millis(10);

pub fn scrips_main() {
    //print_static_distribution();
    //find_max_expected_serial_parallel();
    //test_walk_lengths_once();
    //test_walks_once();
    //test_transformed_walk_lengths_serial();
    //test_transformed_walk_lengths_parallel();
    compare_expected_transformed_serial_parallel();
}
