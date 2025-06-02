use itertools::{Itertools, izip};

struct Graph {
    edges: rand_pebble::bool_csr::BoolCSR,
}

impl Graph {
    pub fn new_grid(len: usize) -> Self {
        let mut edges = rand_pebble::bool_csr::BoolCSR::new();
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
    let mut lcg = rand_pebble::rand_lcg::Lcg::new(0);
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
    let len = 15;
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
        println!("|{line_str}|");
    }
    std::thread::sleep(PRINT_PAUSE_TIME);
}

fn simulate_pebbles_parallel(mut rng: impl FnMut() -> usize, len: usize, start_vertex: usize) -> Vec<usize> {
    let grid = Graph::new_grid(len);
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
    let mut active_pebbles = vec![start_vertex; grid.nr_vertices()];
    let mut new_active_pebbles = active_pebbles.clone();
    let mut walk_lengths = Vec::with_capacity(grid.nr_vertices());
    let mut free_slots = vec![true; grid.nr_vertices()];
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
                let new_pebble = if WALK_LAZILY && rng() % 16 == 0 {
                    pebble
                } else {
                    let neighs = &grid.edges[pebble];
                    neighs[rng() % neighs.len()]
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

fn simulate_pebbles_serial(mut rng: impl FnMut() -> usize, len: usize, start_vertex: usize) -> Vec<usize> {
    let grid = Graph::new_grid(len);
    if PRINT_PROGRESS {
        for _ in 0..len {
            println!();
        }
    }
    let mut walk_lengths = Vec::with_capacity(grid.nr_vertices());
    let mut free_slots = vec![true; grid.nr_vertices()];
    let mut path = Vec::new();
    for _ in 0..grid.nr_vertices() {
        let mut curr = start_vertex;
        path.push(curr);
        let mut nr_steps = 0;
        while !free_slots[curr] {
            curr = if !WALK_LAZILY || rng() % 16 == 0 {
                let curr_neighs = &grid.edges[curr];
                curr_neighs[rng() % curr_neighs.len()]
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

#[allow(dead_code)]
fn test_serial_parallel_once() {
    let len = 30;
    let start = len * (len / 2) + (len / 2);

    let mut lcg = rand_pebble::rand_lcg::Lcg::new(0);

    let rng = || lcg.next_usize();
    let mut serial_walk_lengths = simulate_pebbles_serial(rng, len, start);

    let rng = || lcg.next_usize();
    let parallel_walk_lengths = simulate_pebbles_parallel(rng, len, start);
    
    println!("{parallel_walk_lengths:?}");
    println!("{serial_walk_lengths:?}");
    serial_walk_lengths.sort();
    println!("{serial_walk_lengths:?}");
}

#[allow(dead_code)]
fn find_max_expected_serial_parallel() {
    let mut expectation_serial = Vec::new();
    let mut expectation_parallel = Vec::new();

    let mut len_f = 10.0;
    
    while len_f < 200.0 {
        use rayon::prelude::*;
        use rand::prelude::*;
        let len = len_f as usize;
        len_f *= 1.1;

        const EXPERIMENTS_PER_LEN: usize = 10_000;
        let range = [(); EXPERIMENTS_PER_LEN];

        let serial_sum: usize = range.par_iter().map(|_| {
            let mut rng = rand::rng();
            let rng_f = || rng.random_range(1..10000);
            let mut serial = simulate_pebbles_serial(rng_f, len, 0);
            serial.sort();
            *serial.last().unwrap()
        }).sum();

        
        let parallel_sum: usize = range.par_iter().map(|_| {
            let mut rng = rand::rng();
            let rng_f = || rng.random_range(1..10000);
            let parallel = simulate_pebbles_parallel(rng_f, len, 0);
            *parallel.last().unwrap()
        }).sum();

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

const PRINT_PROGRESS: bool = true;
const WALK_LAZILY: bool = false;
const PRINT_EVERY_FRAME: bool = true;
const PRINT_PAUSE_TIME: std::time::Duration = std::time::Duration::from_millis(100);

fn main() {
    //print_static_distribution();
    //find_max_expected_serial_parallel();
    test_serial_parallel_once();
}
