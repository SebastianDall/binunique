use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub struct Args {
    #[arg(short, long, num_args(2..), required=true, help="directories with bins")]
    pub bin_dirs: Vec<String>,

    #[arg(short, long, num_args(2..), required=true, help="Labels for binners in the same order as above")]
    pub labels: Vec<String>,

    #[arg(long, num_args(2..), help="Path to qualities provided by CheckM2")]
    pub qualities: Option<Vec<String>>,

    #[arg(short, long, num_args(1..), default_values=vec![".fa", ".fna", ".fasta"], help="Extensions to look for in directories. Can be .gz")]
    pub extensions: Vec<String>,

    #[arg(short, long, default_value_t = 1, help = "Num processes")]
    pub threads: usize,

    #[arg(
        short,
        long,
        default_value = "output",
        help = "Path to output directory"
    )]
    pub output: String,
}

pub fn parse_args() -> Args {
    let mut args = Args::parse();

    args.extensions = args
        .extensions
        .into_iter()
        .map(|e| e.to_lowercase())
        .collect();
    args
}
