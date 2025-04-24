use log::info;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::types::{Bin, BinIntersection};

pub fn all_pairwise(input: HashMap<String, Vec<Bin>>) -> Vec<BinIntersection> {
    let binner_labels: Vec<_> = input.keys().cloned().collect();

    let binner_pairs = binner_labels
        .iter()
        .enumerate()
        .flat_map(|(i, label1)| {
            binner_labels
                .iter()
                .skip(i + 1)
                .map(move |label2| (label1, label2))
        })
        .collect::<Vec<(_, _)>>();

    let intersections: Vec<BinIntersection> = binner_pairs
        .into_par_iter()
        .map(|(label1, label2)| {
            let mut local_vec = Vec::new();
            info!("Comparing: {} vs {}", label1, label2);

            if let (Some(bins1), Some(bins2)) = (input.get(label1), input.get(label2)) {
                for bin1 in bins1 {
                    let set1 = bin1.contig_ids.clone();

                    for bin2 in bins2 {
                        let set2 = bin2.contig_ids.clone();

                        let intersection_count = set1.intersection(&set2).count();
                        let union_count = set1.len() + set2.len() - intersection_count;
                        let jaccard_index = if union_count > 0 {
                            intersection_count as f64 / union_count as f64
                        } else {
                            0.0
                        };

                        local_vec.push(BinIntersection {
                            binner_a: label1.clone(),
                            bin_a: bin1.name.clone(),
                            binner_b: label2.clone(),
                            bin_b: bin2.name.clone(),
                            intersection_count,
                            union_count,
                            jaccard_index,
                        });
                    }
                }
            }
            local_vec
        })
        .reduce_with(|mut acc, mut next| {
            acc.append(&mut next);
            acc
        })
        .unwrap_or_default();

    intersections
}
