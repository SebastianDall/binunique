use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

use anyhow::{anyhow, Context};
use clap::Parser;
use flate2::read::GzDecoder;
use seq_io::fasta::{Reader, Record};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    #[arg(short, long, num_args(2..), required=true, help="directories with bins")]
    bin_dirs: Vec<String>,

    #[arg(short, long, num_args(2..), required=true, help="Labels for binners in the same order as above")]
    labels: Vec<String>,

    #[arg(long, num_args(2..), help="Path to qualities provided by CheckM2")]
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

pub struct Contig {
    pub id: String,
    // pub sequence: String,
}

pub struct Bin {
    pub name: String,
    pub contigs: Vec<Contig>,
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

fn read_fasta<P>(path: P) -> anyhow::Result<Vec<Contig>>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();

    let file = File::open(path)?;
    let buf_reader: Box<dyn std::io::BufRead> = if path.extension().map_or(false, |ext| ext == "gz")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut reader = Reader::new(buf_reader);
    let mut contig_ids = Vec::new();

    while let Some(record) = reader.next() {
        let record = record?;

        let rid = record.id()?;
        let id = rid.to_string();
        // let sequence = String::from_utf8(record.owned_seq())?;
        let c = Contig { id }; //, sequence
        contig_ids.push(c);
    }

    Ok(contig_ids)
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
            .map(|entry| entry.path())
            .filter(|path| {
                if let Some(filename) = path.file_name().and_then(|n| n.to_str()) {
                    allowed_exts
                        .iter()
                        .any(|ext| filename.to_lowercase().ends_with(ext))
                } else {
                    false
                }
            })
            .collect::<Vec<_>>();

        if bin_paths.is_empty() {
            return Err(anyhow!("No bins found in path: {}", dir));
        }

        let mut bins = Vec::new();
        for bp in bin_paths.iter() {
            let name = bp
                .file_stem()
                .and_then(|n| n.to_str())
                .ok_or_else(|| anyhow!("Invalid UTF-8 in filename: {:#?}", bp))?;

            let contigs = read_fasta(bp)?;

            let bin = Bin {
                name: name.to_string(),
                contigs,
                completeness: None,
                contamination: None,
            };

            bins.push(bin);
        }
        binner_inputs.insert(label.to_owned(), bins);
    }

    // Compare contig overlap between bins for different binners.
    let mut intersections: Vec<BinIntersection> = Vec::new();

    let binner_labels: Vec<_> = binner_inputs.keys().cloned().collect();
    for (i, label1) in binner_labels.iter().enumerate() {
        for label2 in binner_labels.iter().skip(i + 1) {
            if let (Some(bins1), Some(bins2)) =
                (binner_inputs.get(label1), binner_inputs.get(label2))
            {
                for bin1 in bins1 {
                    let set1: HashSet<String> = bin1.contigs.iter().map(|c| c.id.clone()).collect();

                    for bin2 in bins2 {
                        let set2: HashSet<String> =
                            bin2.contigs.iter().map(|c| c.id.clone()).collect();

                        let intersection_count = set1.intersection(&set2).count();
                        let union_count = set1.len() + set2.len() - intersection_count;
                        let jaccard_index = if union_count > 0 {
                            intersection_count as f64 / union_count as f64
                        } else {
                            0.0
                        };

                        intersections.push(BinIntersection {
                            binner_a: label1.clone(),
                            bin_a: bin1.name.clone(),
                            binner_b: label2.clone(),
                            bin_b: bin2.name.clone(),
                            intersection_count,
                            union_count,
                            jaccard_index,
                        })
                    }
                }
            }
        }
    }

    let outpath = Path::new(&args.output).join("bin-comparisons.tsv");
    let outfile = std::fs::File::create(outpath.clone())
        .with_context(|| anyhow!("Could not create file: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(
        writer,
        "binner_a\tbin_a\tbinner_b\tbin_b\tintersection\tunion\tjaccard_index"
    )?;

    for i in &intersections {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i.binner_a,
            i.bin_a,
            i.binner_b,
            i.bin_b,
            i.intersection_count,
            i.union_count,
            i.jaccard_index
        )?;
    }

    writer.flush()?;
    Ok(())
}
