use anyhow::{anyhow, Context};
use log::info;
use rayon::ThreadPoolBuilder;
use std::collections::HashSet;

use binunique::{cli, compare, io, output};

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = cli::parse_args();
    info!("Using extensions to search for bins {:?}", &args.extensions);

    if args.bin_dirs.len() != args.labels.len() {
        return Err(anyhow!(
            "The provided labels does not match the provided dirs"
        ));
    }

    if args
        .labels
        .clone()
        .into_iter()
        .collect::<HashSet<_>>()
        .len()
        != args.labels.clone().len()
    {
        return Err(anyhow!("Duplicate entries in labels"));
    }

    // Create output dir
    std::fs::create_dir_all(&args.output)
        .with_context(|| anyhow!("Could not create output directory: {}", args.output))?;

    ThreadPoolBuilder::new()
        .num_threads(args.threads as usize)
        .build_global()
        .unwrap();

    let binner_inputs = io::load_bins(&args.bin_dirs, &args.labels, &args.extensions)?;

    // Compare contig overlap between bins for different binners.
    info!("Comparing bin overlaps");
    let intersections = compare::all_pairwise(binner_inputs);

    output::write_tsv(&args.output, &intersections)?;

    Ok(())
}
