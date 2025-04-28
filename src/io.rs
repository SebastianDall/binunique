use anyhow::{anyhow, Context};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::{
    fasta::read_fasta,
    types::{Bin, ContigStats},
};

pub fn load_bins(
    bin_dirs: &Vec<String>,
    labels: &Vec<String>,
    exts: &Vec<String>,
) -> anyhow::Result<HashMap<String, Vec<Bin>>> {
    let mut binner_inputs: HashMap<String, Vec<Bin>> = HashMap::new();
    let style =
        ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("##-");

    for bl in bin_dirs.clone().iter().zip(labels.clone().iter_mut()) {
        let (dir, label) = bl;

        let entries =
            std::fs::read_dir(dir).with_context(|| format!("Failed to list dir: {}", dir))?;

        let bin_paths = entries
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| {
                if let Some(filename) = path.file_name().and_then(|n| n.to_str()) {
                    exts.iter()
                        .any(|ext| filename.to_lowercase().ends_with(ext))
                } else {
                    false
                }
            })
            .collect::<Vec<_>>();

        if bin_paths.is_empty() {
            return Err(anyhow!("No bins found in path: {}", dir));
        }
        let nbins = bin_paths.len() as u64;

        let pb = ProgressBar::new(nbins);
        pb.set_style(style.clone());
        let msg = format!("{}: Loading bins...", &label);
        pb.set_message(msg);

        let bins: Vec<Bin> = bin_paths
            .into_par_iter()
            .progress_with(pb)
            .map(|bp| {
                let name = bp
                    .file_stem()
                    .and_then(|n| n.to_str())
                    .ok_or_else(|| anyhow!("Invalid UTF-8 in filename: {:#?}", bp))?;

                let contigs = read_fasta(&bp)?;

                let contig_ids = contigs
                    .iter()
                    .map(|c| c.clone())
                    .collect::<HashSet<ContigStats>>();

                Ok(Bin {
                    name: name.to_string(),
                    file: bp,
                    contig_ids,
                    completeness: None,
                    contamination: None,
                })
            })
            .collect::<Result<Vec<Bin>, anyhow::Error>>()?;

        binner_inputs.insert(label.to_owned(), bins);
    }
    Ok(binner_inputs)
}
