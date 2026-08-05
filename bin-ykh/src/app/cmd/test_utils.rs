//! Shared helpers for the per-command tests.

use std::error::Error;
use std::fmt::Debug;
use yui_link::{InvLink, Link};

// PD code of a `lib-link` test knot, in the JSON form the `link` argument takes.
pub fn pd(name: &str) -> String {
    format!("{:?}", Link::test_data(name).pd_code())
}

// PD code of a `lib-link` symmetric test diagram, for the involutive commands.
pub fn inv_pd(name: &str) -> String {
    format!("{:?}", InvLink::test_data(name).pd_code())
}

// Compare a command's output cell by cell, so column padding and blank lines don't matter
// but every printed entry does.
pub fn assert_out(res: Result<String, Box<dyn Error>>, expected: &str) {
    let cells = |s: &str| s.lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect::<Vec<_>>())
        .filter(|l: &Vec<String>| !l.is_empty())
        .collect::<Vec<_>>();

    assert_eq!(cells(&res.expect("dispatch failed")), cells(expected));
}

// `Args::default()` must agree with clap's `default_value`s — the command tests build `Args` by
// struct literal, so a mismatch would have them testing something the CLI never runs.
pub fn assert_cli_default<T>(parsed: &T, default: &T)
where T: PartialEq + Debug {
    assert_eq!(parsed, default);
}
