use anyhow::anyhow;
use log::info;
use rayon::prelude::*;
use std::collections::HashMap;
use std::process::Command;
use std::str;

use crate::types::{Bin, BinIntersection, ANI};

pub fn compute_ani(genome1: &Bin, genome2: &Bin) -> anyhow::Result<ANI> {
    let output = Command::new("skani")
        .args(&[
            "dist",
            genome1.file.as_os_str().to_str().unwrap(),
            genome2.file.as_os_str().to_str().unwrap(),
        ])
        .output()?;

    if !output.status.success() {
        return Err(anyhow!("Skani failed: {}", String::from_utf8_lossy(&output.stderr)).into());
    }

    let stdout = str::from_utf8(&output.stdout)?;

    let mut lines = stdout.lines();

    let _header = lines.next();

    let maybe_result_line = lines.next();
    if let Some(result_line) = maybe_result_line {
        let fields: Vec<&str> = result_line.split('\t').collect();
        if fields.len() != 7 {
            return Err(anyhow!("Unexpected output from skani"));
        }

        let ani: f64 = fields[2].parse()?;
        let afr: f64 = fields[3].parse()?;
        let afq: f64 = fields[4].parse()?;

        let ani_result = ANI {
            ani,
            alignment_fraction_ref: afr,
            alignment_fraction_query: afq,
        };

        Ok(ani_result)
    } else {
        Ok(ANI::default())
    }
}

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

                        let intersection = set1.intersection(&set2);

                        let intersection_count = intersection.clone().count();
                        let union_count = set1.len() + set2.len() - intersection_count;
                        let jaccard_index = if intersection_count > 0 {
                            intersection_count as f64 / union_count as f64
                        } else {
                            0.0
                        };

                        let intersection_size: u64 =
                            intersection.into_iter().map(|c| c.length).sum();

                        let bin1_size: u64 = set1.clone().into_iter().map(|c| c.length).sum();
                        let bin2_size: u64 = set2.clone().into_iter().map(|c| c.length).sum();

                        let union_size = bin1_size + bin2_size - intersection_size;

                        let weighted_jaccard_index = if intersection_count > 0 {
                            intersection_size as f64 / union_size as f64
                        } else {
                            0.0
                        };

                        let ani = if intersection_count > 0 {
                            compute_ani(&bin1, &bin2).unwrap()
                        } else {
                            ANI::default()
                        };

                        local_vec.push(BinIntersection {
                            binner_a: label1.clone(),
                            bin_a: bin1.name.clone(),
                            binner_b: label2.clone(),
                            bin_b: bin2.name.clone(),
                            intersection_count,
                            union_count,
                            jaccard_index,
                            intersection_size,
                            union_size,
                            weighted_jaccard_index,
                            ani,
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
