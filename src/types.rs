use std::collections::HashSet;

#[derive(Clone, Hash, PartialEq, Eq)]
pub struct ContigStats {
    pub id: String,
    pub length: u64,
}

pub struct Bin {
    pub name: String,
    // pub contigs: Vec<Contig>,
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
}
