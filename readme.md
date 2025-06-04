
note: this is a script, not a program.

## how to use
- install rust via rustup
- if on linux: install dependencies required by plotters crate
- depending on what should be examined, please scroll to the end of `src/scripts.rs` and change things
- run `cargo run --release` in project folder
- enjoy

## assumptions in proof of t_par = O(t_ser * log(n))
- every time the remnants of the longest row are encountered, it is again cut, except when it has reached the last row.
- every time the remnants of the longest row are pasted somwhere, they are pasted onto an otherwise empty row.
- the expected length of largest row remnant, one walks until the next cut happens, is constant throughut the process.
    (it actually increases with each finished row, as each row removes a vertex from the set of those causing a cut.)

in parallel IDLA, each row has the same length distribution.
in serial IDLA, the expected length increases with the row index.
when transforming parallel to serial, _every_ row is potentially partially pasted further back.
-> sections of all rows accumulate at the end.

