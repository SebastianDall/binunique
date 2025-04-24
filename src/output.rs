use std::{
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{anyhow, Context, Result};

use crate::types::BinIntersection;

pub fn write_tsv<P: AsRef<Path>>(outdir: P, data: &[BinIntersection]) -> Result<()> {
    let outpath = outdir.as_ref().join("bin-comparisons.tsv");
    let outfile = std::fs::File::create(outpath.clone())
        .with_context(|| anyhow!("Could not create file: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(
        writer,
        "binner_a\tbin_a\tbinner_b\tbin_b\tintersection\tunion\tjaccard_index\tintersection_size\tunion_size\tweighted_jaccard_index"
    )?;

    for i in data {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i.binner_a,
            i.bin_a,
            i.binner_b,
            i.bin_b,
            i.intersection_count,
            i.union_count,
            i.jaccard_index,
            i.intersection_size,
            i.union_size,
            i.weighted_jaccard_index
        )?;
    }

    writer.flush()?;
    Ok(())
}
