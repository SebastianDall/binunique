use std::{fs::File, io::BufReader, path::Path};

use flate2::read::GzDecoder;
use seq_io::fasta::{Reader, Record};

use crate::types::ContigStats;

pub fn read_fasta<P>(path: P) -> anyhow::Result<Vec<ContigStats>>
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
        let sequence = String::from_utf8(record.owned_seq())?;
        let seq_length = sequence.len() as u64;
        let c = ContigStats {
            id,
            length: seq_length,
        }; //, sequence
        contig_ids.push(c);
    }

    Ok(contig_ids)
}
