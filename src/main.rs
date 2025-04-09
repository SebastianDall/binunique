use std::{collections::HashMap, path::Path};

use anyhow::{anyhow, Context};
use clap::Parser;
use seq_io::fasta::{Reader, Record};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    #[arg(short, long, num_args(2..), required=true, help="directories with bins")]
    bin_dirs: Vec<String>,

    #[arg(short, long, num_args(2..), required=true, help="Labels for binners in the same order as above")]
    labels: Vec<String>,

    #[arg(short, long, num_args(2..), help="Path to qualities provided by CheckM2")]
    qualities: Option<Vec<String>>,

    #[arg(short, long, num_args(1..), default_values=vec![".fa", ".fna", ".fasta"], help="Extensions to look for in directories")]
    extensions: Option<Vec<String>>,

    #[arg(
        short,
        long,
        default_value = "output",
        help = "Path to output directory"
    )]
    output: String,
}

struct Contig {
    pub id: String,
    pub sequence: String,
}

impl Contig {
    pub fn new(id: String, sequence: String) -> Self {
        Self { id, sequence }
    }
}

pub struct Bin {
    name: String,
    contigs: Vec<Contig>,
    completeness: Option<f64>,
    contamination: Option<f64>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let allowed_exts: Vec<String> = args
        .extensions
        .unwrap_or_else(|| vec![".fa".into(), ".fna".into(), ".fasta".into()])
        .into_iter()
        .map(|e| e.to_lowercase())
        .collect();

    if args.bin_dirs.len() != args.labels.len() {
        return Err(anyhow!(
            "The provided labels does not match the provided dirs"
        ));
    }

    let mut binner_inputs: HashMap<String, Vec<Bin>> = HashMap::new();

    for bl in args
        .bin_dirs
        .clone()
        .iter()
        .zip(args.labels.clone().iter_mut())
    {
        let (dir, label) = bl;

        let entries =
            std::fs::read_dir(dir).with_context(|| format!("Failed to list dir: {}", dir))?;

        let bin_paths = entries
            .filter_map(Result::ok)
            .filter(|entry| {
                let path = entry.path();
                if let Some(ext_os) = path.extension() {
                    let ext = ext_os.to_string_lossy().to_lowercase();
                    allowed_exts.iter().any(|allowed| {
                        let allowed_norm = allowed.trim_start_matches(".");
                        ext == allowed_norm
                    })
                } else {
                    false
                }
            })
            .collect::<Vec<_>>();

        if bin_paths.is_empty() {
            return Err(anyhow!("No bins found in path: {}", dir));
        }

        for bp 
    }

    Ok(())
}
