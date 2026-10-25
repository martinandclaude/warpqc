//! Shared harness: runs a fixture in-process, once single-threaded with one
//! block and once with 4 workers and tiny blocks, and requires identical
//! JSON (every fixture is built so that its numbers are deterministic).

#![allow(dead_code)]

use std::path::PathBuf;

use serde_json::Value;
use warpqc::pipeline::{RunConfig, run};

pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn reference() -> PathBuf {
    root().join("tests/fixtures/ref/mini.fa")
}

pub fn fixture(name: &str) -> PathBuf {
    root().join("tests/fixtures/sam").join(name)
}

pub type Depth = Vec<(String, Vec<u16>, Option<Vec<u16>>)>;

pub fn config(name: &str) -> RunConfig {
    let mut c = RunConfig::new("T", reference(), Some(fixture(name)));
    c.contigs = Some("chr1,chr2,chrM".into());
    c.keep_depth = true;
    c
}

pub fn run_cfg(c: &RunConfig) -> (Value, Depth) {
    let out = run(c).unwrap_or_else(|e| panic!("{e}"));
    (serde_json::to_value(&out.report).unwrap(), out.depth)
}

pub fn run_fixture(name: &str) -> (Value, Depth) {
    let mut one = config(name);
    one.threads = 1;
    one.block_size = 1 << 20;
    let mut many = config(name);
    many.threads = 4;
    many.block_size = 64;
    let (a, da) = run_cfg(&one);
    let (b, db) = run_cfg(&many);
    assert_eq!(a, b, "{name}: 1 thread and 4 threads disagree");
    assert_eq!(da, db, "{name}: depth arrays disagree");
    (a, da)
}

pub fn get<'a>(v: &'a Value, ptr: &str) -> &'a Value {
    v.pointer(ptr).unwrap_or_else(|| panic!("missing {ptr}"))
}

pub fn u(v: &Value, ptr: &str) -> u64 {
    get(v, ptr)
        .as_u64()
        .unwrap_or_else(|| panic!("{ptr} is not an integer: {}", get(v, ptr)))
}

#[track_caller]
pub fn close(v: &Value, ptr: &str, want: f64) {
    let got = get(v, ptr)
        .as_f64()
        .unwrap_or_else(|| panic!("{ptr} is not a number"));
    assert!((got - want).abs() < 1e-9, "{ptr}: got {got}, want {want}");
}

/// Depth at a 1-based position.
pub fn depth_at(d: &Depth, contig: &str, pos1: usize) -> (u16, Option<u16>) {
    let (_, f, r) = d.iter().find(|(n, _, _)| n == contig).unwrap();
    (f[pos1 - 1], r.as_ref().map(|r| r[pos1 - 1]))
}
