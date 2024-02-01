# kmerHLL
This crate provides implementation of probabilistic algorithms to count unique kmers in sequence files via HyperLogLog(HLL)-like algorithms.

It is developed for daily sequence processing. I did not find a stand-alone Rust program for counting kmers in sequence files. It is also very useful for MinHash algorithms like FracMinHash to estimate Jaccard index from containment index, in which the carginality of each genome files must be known. Exact algorithms for real-world genome/metagenomes files are not practical. Raw fasta and zipped fasta file format are supported.

It also estiamtes the Jaccard index via the inclusion-exclusion rule after knowing the cardinality of two genomes and their union.

## Install (nightly Rust must be installed)

```bash
git clone https://github.com/jianshu93/kmerHLL
cd kmerHLL
cargo build --release
./target/release/kmerHLL -h
```

```bash
Counting unique k-mers in sequence files via HLL

Usage: kmerHLL <fasta_file1> <fasta_file2> <kmer_length>

Arguments:
  <fasta_file1>  The first FASTA file to read
  <fasta_file2>  The second FASTA file to read
  <kmer_length>  The length of the k-mers

Options:
  -h, --help     Print help
  -V, --version  Print version
```

## Details of implementation
We use the streaming algorithms crate, which utilize SIMD instructions (nightly Rust for all platforms) for speedup certain operations. 

Since HyperLogLog are mergeable sketches, parallelism is very easy. Merging sketches from each thread require lock of the global sketches. For large files like metagenomes, the sequence file will be split into chunks so that each thread can processing part of the sequence file.

## Performance and accuracy

It took about 1.5 minutes for processing a 15G metagenome using 24 threads on a Intel(R) Xeon(R) Gold 6226 CPU@2.70GHz. IO is the main limit now.

You can adjust the accuracy via initializing HLL with a given error rate, which controls the number of counters (registers):

$$MSE=\frac{3*log2 -1}{m}\approx \frac{1.07944}{m}$$

However, the accuracy can never go beyond the limit below no matter how many counters:

$$MSE=\frac{log2/(\frac{\pi ^{2}}{6}-1)}{m}\approx \frac{1.07475}{m}$$

which is the Cramér Rao lower bound for such counting algorithms.

## References
Flajolet, Philippe, Éric Fusy, Olivier Gandouet, and Frédéric Meunier. "Hyperloglog: the analysis of a near-optimal cardinality estimation algorithm." Discrete mathematics & theoretical computer science Proceedings (2007).

Heule, Stefan, Marc Nunkesser, and Alexander Hall. "Hyperloglog in practice: Algorithmic engineering of a state of the art cardinality estimation algorithm." In Proceedings of the 16th International Conference on Extending Database Technology, pp. 683-692. 2013.


