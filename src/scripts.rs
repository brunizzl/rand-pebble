use itertools::{Itertools, izip};
use rand::Rng;

use crate::bool_csr::BoolCSR;

struct Graph {
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

fn find_stationary_distribution(len: usize) -> (Vec<isize>, usize) {
    let grid = Graph::new_grid(len);
    let mut curr_visits = vec![0isize; len * len];
    let mut all_visits = vec![0isize; len * len];
    let mut nr_steps = 1;
    let mut curr_pos = len * (len / 2) + len / 2;
    let mut lcg = crate::rand_lcg::Lcg::new(0);
    loop {
        for _ in 0..nr_steps {
            curr_visits[curr_pos] += 1;
            let curr_neighs = &grid.edges[curr_pos];
            curr_pos = curr_neighs[(lcg.next() as usize) % curr_neighs.len()];
        }
        let cum_relative_diff = izip!(&curr_visits, &all_visits)
            .map(|(&a, &b)| ((a - b).abs() as f64) / (isize::max(a, b) as f64))
            .sum::<f64>();
        let avg_relative_diff = cum_relative_diff / (curr_visits.len() as f64);
        println!("steps: {nr_steps:>20} diff: {avg_relative_diff}");
        if avg_relative_diff < 0.0002 {
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

#[allow(dead_code)]
fn print_stationary_distribution() {
    let len = FIXED_LEN;
    let (visits, nr_steps) = find_stationary_distribution(len);
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

fn print_free_slots(free_slots: &[bool], len: usize, markers: &[usize]) {
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

            print_free_slots(&free_slots, len, &active_pebbles);
        }
    }

    walk_lengths
}

fn simulate_walks_parallel(graph: &Graph, start_vertex: usize) -> Vec<Vec<usize>> {
    let mut rng = rand::rng();
    let len = f64::sqrt(graph.nr_vertices() as f64) as usize;
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
    let mut nr_active_pebbles = graph.nr_vertices();
    let mut active_pebbles = vec![true; graph.nr_vertices()];
    let mut walks = vec![vec![start_vertex]; graph.nr_vertices()];
    let mut free_slots = vec![true; graph.nr_vertices()];

    let mut nr_free_when_last_printed = nr_active_pebbles + 1;
    let mut curr_pebble_positions = Vec::new();
    for _ in 0.. {
        if nr_active_pebbles == 0 {
            break;
        }
        for (active, walk) in izip!(&mut active_pebbles, &mut walks) {
            if !*active {
                continue;
            }
            let pebble = *walk.last().unwrap();
            curr_pebble_positions.push(pebble);
            if free_slots[pebble] {
                free_slots[pebble] = false;
                *active = false;
                nr_active_pebbles -= 1;
            } else {
                let new_pebble = if WALK_LAZILY && rng.random_range(0..16) == 0 {
                    pebble
                } else {
                    let neighs = &graph.edges[pebble];
                    neighs[rng.random_range(0..neighs.len())]
                };
                walk.push(new_pebble);
            }
        }

        if PRINT_PROGRESS
            && (PRINT_EVERY_FRAME || nr_free_when_last_printed != active_pebbles.len())
        {
            nr_free_when_last_printed = active_pebbles.len();

            print_free_slots(&free_slots, len, &curr_pebble_positions);
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
            print_free_slots(&free_slots, len, &path);
        }
        path.clear();
    }

    walk_lengths
}

fn simulate_walks_serial(graph: &Graph, start_vertex: usize) -> Vec<Vec<usize>> {
    let mut rng = rand::rng();
    let len = f64::sqrt(graph.nr_vertices() as f64) as usize;
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
    let mut walks = Vec::with_capacity(graph.nr_vertices());
    let mut free_slots = vec![true; graph.nr_vertices()];
    for _ in 0..graph.nr_vertices() {
        let mut path = Vec::new();
        let mut curr = start_vertex;
        path.push(curr);
        while !free_slots[curr] {
            curr = if !WALK_LAZILY || rng.random_range(0..16) == 0 {
                let curr_neighs = &graph.edges[curr];
                curr_neighs[rng.random_range(0..curr_neighs.len())]
            } else {
                curr
            };
            path.push(curr);
        }
        free_slots[curr] = false;

        if PRINT_PROGRESS {
            print_free_slots(&free_slots, len, &path);
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

fn print_walks(walks: &[Vec<usize>]) {
    if !PRINT_WALKS {
        return;
    }
    for (i, walk) in izip!(0.., walks) {
        print!("{i:>4} |");
        for step in walk {
            print!(" {step:>3}");
        }
        println!();
    }
}

#[allow(dead_code)]
fn test_walks_once() {
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = Graph::new_grid(len);

    let serial_walks = simulate_walks_serial(&graph, start);
    print_walks(&serial_walks);

    let parallel_walks = simulate_walks_parallel(&graph, start);
    print_walks(&parallel_walks);
}

fn transform_walks_serial_to_parallel(graph: &Graph, walks: &mut [Vec<usize>]) {
    let mut nr_vertices_visited = 0;
    let mut vertices_visited = vec![false; graph.nr_vertices()];
    let mut t = 0;
    while nr_vertices_visited < graph.nr_vertices() {
        for i in 0..graph.nr_vertices() {
            if walks[i].len() > t {
                let pebble = walks[i][t];
                if !vertices_visited[pebble] {
                    vertices_visited[pebble] = true;
                    nr_vertices_visited += 1;
                    let paste_pos = walks
                        .iter()
                        .position(|w| w.last() == Some(&pebble))
                        .unwrap();
                    if paste_pos != i {
                        let mut walk_i = std::mem::take(&mut walks[i]);
                        let to_copy = &walk_i[(t + 1)..];
                        walks[paste_pos].extend_from_slice(to_copy);
                        walk_i.resize(t + 1, 0);
                        walks[i] = walk_i;
                    }
                }
            }
        }
        t += 1;
    }
}

#[allow(dead_code)]
fn transform_walks_parallel_to_serial(graph: &Graph, walks: &mut [Vec<usize>]) {
    let mut vertices_visited = vec![false; graph.nr_vertices()];
    for i in 0..graph.nr_vertices() {
        let mut t = 0;
        while t < walks[i].len() {
            let pebble = walks[i][t];
            if !vertices_visited[pebble] {
                vertices_visited[pebble] = true;
                let paste_pos = walks
                    .iter()
                    .position(|w| w.last() == Some(&pebble))
                    .unwrap();
                if paste_pos != i {
                    let mut walk_i = std::mem::take(&mut walks[i]);
                    let to_copy = &walk_i[(t + 1)..];
                    walks[paste_pos].extend_from_slice(to_copy);
                    walk_i.resize(t + 1, 0);
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
    let start = 0;//len * (len / 2) + (len / 2);
    let graph = Graph::new_torus_grid(len);
    let mut ratio_sum = 0.0;
    const NR_SAMPLES: usize = 1000;
    for _ in 0..NR_SAMPLES {
        let mut walks = simulate_walks_serial(&graph, start);
        print_walks(&walks);
        let max_serial = walks.iter().map(|w| w.len()).max().unwrap_or(0);
        println!("max serial len = {max_serial}");

        transform_walks_serial_to_parallel(&graph, &mut walks);
        print_walks(&walks);
        let max_parallel = walks.iter().map(|w| w.len()).max().unwrap_or(0);
        println!("max parallel len = {max_parallel}");
        let ratio = max_serial as f64 / max_parallel as f64;
        println!("ratio: {ratio}");
        ratio_sum += ratio;
    }
    println!("=> avg ratio: {}", ratio_sum / NR_SAMPLES as f64);
}

#[allow(dead_code)]
fn test_transformed_walk_lengths_parallel() {
    println!("compute parallel, then transform to serial");
    let len = FIXED_LEN;
    let start = 0; //len * (len / 2) + (len / 2);
    let graph = Graph::new_grid(len);

    let mut walks = simulate_walks_parallel(&graph, start);
    print_walks(&walks);
    let max_parallel = walks.iter().map(|w| w.len()).max().unwrap_or(0);
    println!("max parallel len = {max_parallel}");

    transform_walks_parallel_to_serial(&graph, &mut walks);
    print_walks(&walks);
    let max_serial = walks.iter().map(|w| w.len()).max().unwrap_or(0);
    println!("max serial len = {max_serial}");
    println!("ratio: {}", max_serial as f64 / max_parallel as f64);
}

#[allow(dead_code)]
fn find_max_expected_serial_parallel() {
    let mut expectation_serial = Vec::new();
    let mut expectation_parallel = Vec::new();

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
        expectation_serial.push((len * len, exp_serial));
        expectation_parallel.push((len * len, exp_parallel));
        let ratio = serial_sum as f64 / parallel_sum as f64;
        println!(
            "len = {len:>3}, t_ser = {exp_serial:>10.2}, t_par = {exp_parallel:>10.2}, ratio = {ratio}"
        );
    }
}

const FIXED_LEN: usize = 100;
const PRINT_WALKS: bool = false;
const PRINT_PROGRESS: bool = false;
const WALK_LAZILY: bool = false;
const PRINT_EVERY_FRAME: bool = true;
const PRINT_PAUSE_TIME: std::time::Duration = std::time::Duration::from_millis(10);

pub fn scrips_main() {
    //print_static_distribution();
    //find_max_expected_serial_parallel();
    //test_walk_lengths_once();
    //test_walks_once();
    test_transformed_walk_lengths_serial();
    //test_transformed_walk_lengths_parallel();
}
