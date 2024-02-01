use clap::{Arg, Command};
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::error::Error;
use std::sync::{Arc, Mutex};
use streaming_algorithms::HyperLogLog;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("K-mer Counter")
        .version("0.1.0")
        .author("Your Name")
        .about("Counting unique k-mers in sequence files via HLL")
        .arg(Arg::new("fasta_file1")
            .help("The first FASTA file to read")
            .required(true)
            .index(1))
        .arg(Arg::new("fasta_file2")
            .help("The second FASTA file to read")
            .required(true)
            .index(2))
        .arg(Arg::new("kmer_length")
            .help("The length of the k-mers")
            .required(true)
            .index(3))
        .get_matches();

    let fasta_file1 = matches.get_one::<String>("fasta_file1").unwrap().as_str();
    let fasta_file2 = matches.get_one::<String>("fasta_file2").unwrap().as_str();
    let kmer_length: usize = matches.get_one::<String>("kmer_length").unwrap().parse()?;

    let hll1 = process_file(fasta_file1, kmer_length)?;
    let hll2 = process_file(fasta_file2, kmer_length)?;

    // Combine the HyperLogLogs from both files
    let mut combined_hll = hll1.clone();
    combined_hll.union(&hll2);

    let cardinality1 = hll1.len();
    let cardinality2 = hll2.len();
    let combined_cardinality = combined_hll.len();
    let jaccard_index = (cardinality1 + cardinality2 - combined_cardinality) / combined_cardinality;

    println!("Number of unique k-mers in first file: {}", cardinality1);
    println!("Number of unique k-mers in second file: {}", cardinality2);
    println!("Combined number of unique k-mers: {}", combined_cardinality);
    println!("Jaccard Index: {}", jaccard_index);
    Ok(())
}

fn process_file(file_name: &str, kmer_length: usize) -> Result<HyperLogLog<Vec<u8>>, Box<dyn Error>> {
    let mut reader = parse_fastx_file(file_name)?;
    let hll = Arc::new(Mutex::new(HyperLogLog::<Vec<u8>>::new(0.00408)));

    let mut sequences = Vec::new();
    while let Some(result) = reader.next() {
        let record = result?;
        sequences.push(record.seq().to_vec());
    }

    // Determine chunk size based on the number of available CPUs
    let mut chunk_size = sequences.len() / rayon::current_num_threads();
    //just incase there is not enough number of sequences to run on multiple threads
    if chunk_size == 0 {
        chunk_size += 1;
    }

    sequences
        .chunks(chunk_size)
        .collect::<Vec<_>>()
        .into_par_iter()
        .for_each(|chunk| {
            let mut local_hll = HyperLogLog::<Vec<u8>>::new(0.00408);
            for seq in chunk {
                let kmer_length_u8 = kmer_length.try_into().unwrap(); // Handle error appropriately

                for kmer in needletail::kmer::Kmers::new(seq, kmer_length_u8) {
                    local_hll.push(&kmer.to_vec());
                }
            }

            let mut global_hll = hll.lock().unwrap(); // Lock the global HyperLogLog for update
            global_hll.union(&local_hll); // Merge local HyperLogLog into global one
        });

    let final_hll = Arc::try_unwrap(hll).expect("Lock still has multiple owners");
    final_hll.into_inner().map_err(Into::into)
}
