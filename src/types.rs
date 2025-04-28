use std::{collections::HashSet, path::PathBuf};

#[derive(Clone, Hash, PartialEq, Eq)]
pub struct ContigStats {
    pub id: String,
    pub length: u64,
}

pub struct Bin {
    pub name: String,
    pub file: PathBuf,
    pub contig_ids: HashSet<ContigStats>,
    pub completeness: Option<f64>,
    pub contamination: Option<f64>,
}

pub struct BinIntersection {
    pub binner_a: String,
    pub bin_a: String,
    pub binner_b: String,
    pub bin_b: String,
    pub intersection_count: usize,
    pub union_count: usize,
    pub jaccard_index: f64,
    pub intersection_size: u64,
    pub union_size: u64,
    pub weighted_jaccard_index: f64,
    pub ani: ANI,
}

#[derive(Default)]
pub struct ANI {
    pub ani: f64,
    pub alignment_fraction_ref: f64,
    pub alignment_fraction_query: f64,
}
