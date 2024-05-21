use clap::{Arg, Command};
use needletail::{parse_fastx_file, Sequence};
use needletail::kmer::Kmers;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::thread;
use crossbeam_channel::{bounded, unbounded};
use streaming_algorithms::HyperLogLog;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("K-mer Counter")
        .version("0.1.0")
        .author("Jianshu Zhao")
        .about("Counting unique k-mers in sequence files via HyperLogLog and HyperMinHash")
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
    let (sender, receiver) = bounded(20);

    // Spawn a thread to read sequences and batch them
    let reader_thread = thread::spawn(move || {
        let mut batch = Vec::new();
        while let Some(result) = reader.next() {
            if let Ok(record) = result {
                batch.push(record.seq().to_vec());
                if batch.len() == 40000 {
                    sender.send(batch.clone()).expect("Failed to send batch");
                    batch = Vec::new(); // Reset the batch
                }
            }
        }
        if !batch.is_empty() {
            sender.send(batch).expect("Failed to send the final batch");
        }
    });

    // Process batches of sequences and update the global hll once per batch
    let processor_thread = thread::spawn({
        let hll_clone = Arc::clone(&hll);
        move || {
            receiver.into_iter().for_each(|batch| {
                // Use fold to create a local HyperLogLog for each thread
                let local_hll = batch.par_iter().fold(
                    || HyperLogLog::<Vec<u8>>::new(0.00408),
                    |mut acc, seq| {
                        let kmer_length_u8 = kmer_length.try_into().unwrap();
                        for kmer in Kmers::new(seq, kmer_length_u8) {
                            acc.push(&kmer.to_vec());
                        }
                        acc
                    }
                );

                // Use reduce to merge all local HyperLogLogs into one
                let reduced_hll = local_hll.reduce(
                    || HyperLogLog::<Vec<u8>>::new(0.00408),
                    |mut a, b| {
                        a.union(&b);
                        a
                    }
                );

                let mut global_hll = hll_clone.lock().unwrap();
                global_hll.union(&reduced_hll);  // Perform the union operation once per batch
            });
        }
    });

    reader_thread.join().unwrap();
    processor_thread.join().unwrap();

    Arc::try_unwrap(hll)
        .map_err(|_arc| Box::new(std::io::Error::new(std::io::ErrorKind::Other, "Lock still has multiple owners")) as Box<dyn Error>)
        .and_then(|mutex| mutex.into_inner().map_err(|e| Box::new(e) as Box<dyn Error>))
}