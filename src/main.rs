use anyhow::{anyhow, Context};
use clap::Parser;
use flate2::read::GzDecoder;
use log::info;
use rayon::{prelude::*, ThreadPoolBuilder};
use seq_io::fasta::{Reader, Record};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

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

    #[arg(short, long, default_value_t = 1, help = "Num processes")]
    threads: usize,

    #[arg(
        short,
        long,
        default_value = "output",
        help = "Path to output directory"
    )]
    output: String,
}

#[derive(Clone, Hash, PartialEq, Eq)]
pub struct Contig {
    pub id: String,
    // pub sequence: String,
}

pub struct Bin {
    pub name: String,
    // pub contigs: Vec<Contig>,
    pub contig_ids: HashSet<Contig>,
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
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    let allowed_exts: Vec<String> = args
        .extensions
        .unwrap_or_else(|| vec![".fa".into(), ".fna".into(), ".fasta".into()])
        .into_iter()
        .map(|e| e.to_lowercase())
        .collect();
    info!("Using extensions to search for bins {:?}", &allowed_exts);

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
        info!("{}: {} bins found.", &label, bin_paths.len());

        let bins: Vec<Bin> = bin_paths
            .into_par_iter()
            .map(|bp| {
                let name = bp
                    .file_stem()
                    .and_then(|n| n.to_str())
                    .ok_or_else(|| anyhow!("Invalid UTF-8 in filename: {:#?}", bp))?;

                let contigs = read_fasta(&bp)?;

                let contig_ids = contigs
                    .iter()
                    .map(|c| c.clone())
                    .collect::<HashSet<Contig>>();

                Ok(Bin {
                    name: name.to_string(),
                    contig_ids,
                    completeness: None,
                    contamination: None,
                })
            })
            .collect::<Result<Vec<Bin>, anyhow::Error>>()?;

        binner_inputs.insert(label.to_owned(), bins);
    }

    // Compare contig overlap between bins for different binners.
    info!("Comparing bin overlaps");
    // let mut intersections: Vec<BinIntersection> = Vec::new();

    let binner_labels: Vec<_> = binner_inputs.keys().cloned().collect();

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

            if let (Some(bins1), Some(bins2)) =
                (binner_inputs.get(label1), binner_inputs.get(label2))
            {
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

    // for (i, label1) in binner_labels.iter().enumerate() {
    //     for label2 in binner_labels.iter().skip(i + 1) {
    //         if let (Some(bins1), Some(bins2)) =
    //             (binner_inputs.get(label1), binner_inputs.get(label2))
    //         {
    //             for bin1 in bins1 {
    //                 let set1 = bin1.contig_ids.clone();

    //                 for bin2 in bins2 {
    //                     let set2 = bin2.contig_ids.clone();

    //                     let intersection_count = set1.intersection(&set2).count();
    //                     let union_count = set1.len() + set2.len() - intersection_count;
    //                     let jaccard_index = if union_count > 0 {
    //                         intersection_count as f64 / union_count as f64
    //                     } else {
    //                         0.0
    //                     };

    //                     intersections.push(BinIntersection {
    //                         binner_a: label1.clone(),
    //                         bin_a: bin1.name.clone(),
    //                         binner_b: label2.clone(),
    //                         bin_b: bin2.name.clone(),
    //                         intersection_count,
    //                         union_count,
    //                         jaccard_index,
    //                     })
    //                 }
    //             }
    //         }
    //     }
    // }

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
