[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://www.zenodo.org/badge/643454906.svg)](https://www.zenodo.org/badge/latestdoi/643454906)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/fast-adapter-trimming)
![GitHub Continuous Integration Status](https://github.com/linsalrob/fast-adapter-trimming/actions/workflows/c.yml/badge.svg)

# Fast Adapter Trimming

Sequencing adapters occassionally end up in the sequences, and we need to trim them out. There are several tools for doing this, but `fast-adapter-trimming` has been designed to be _fast_ and to be _memory efficient_. 


# Should I use _fast-adapter-trimming_?

Several other tools offer adapter trimming options, but none met our specifications, so we wrote this. In particular, if you are processing lots of `fastq` files, this is a solution for you to try.

Advantages:
 - Written in C
    - Portable
    - Fast
    - Memory efficient
 - Threaded to process R1 and R2 files in parallel
 - Reconcile paired-end reads
 - Adapters provided in fasta file
 - Checks for mismatches between adapter and sequence
 - Summarizes all adapters found in the R1 and R2 files
 - Trims sequences and writes gzip compressed files


Disadvantages:
 - Truncates all adapters to a maximum length of 31 nucleotides
 - Only performs adapter trimming, does not perform sequence QC (we recommend [prinseq++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) for that!


# Basic installation and use

We are working on a conda installation, and so at the moment you need to clone this repository:

```
git clone https://github.com/linsalrob/fast-adapter-trimming.git
cd fast-adapter-trimming
make all
make install PREFIX=$HOME/bin
fast-adapter-trimming
```

Once you have it installed, you can trim the adapters from sequence reads. Note that we provide a few files of [adapter sequences](https://github.com/linsalrob/fast-adapter-trimming/tree/main/adapters) but if you know of some we are missing, please [post an issue](https://github.com/linsalrob/fast-adapter-trimming/issues)

For example, to find the adapters in the [test sequences](https://github.com/linsalrob/fast-adapter-trimming/tree/main/fastq) you can use:

```
fast-adapter-trimming -1 fastq/788707_20180313_S_R1.small.fastq.gz -2 fastq/788707_20180313_S_R2.small.fastq.gz -f adapters/truseq.fa  --primeroccurrences 50 --matchesR1 matches.R1.tsv --matchesR2 matches.R2.tsv --outputR1 R1.trimmed.fastq.gz --outputR2 R2.trimmed.fastq.gz
```

# Two and a half ways of removing adapters

We have three different ways of removing adapters, but two of those are related!

You can either remove the adapters separately, or as a paired end library. The fastest way to use `fast-adapter-removal` is to remove the adapters from R1 and R2 reads independently from each other. In this approach, we look at each read and trim it starting with the adapter that lies closest to the start of the sequence (most 5'). We have two ways of doing this: we can either use seperate threads for the R1 and R2 files, which is the fastest and default mode for `fast-adapter-trimming`, or process them sequentially by using the `--nothreads` option.

Alternatively, you can remove the adapters in `--paired_end` mode. In this case, we identify all the adapters that match to all the reads, and then reconcile the R1 and R2 reads so that we remove the same amount of sequence from both. Currently, we just trim to whichever is shorter, because (a) that is by far quicker than aligning the two reads, and (b) in most of our test cases that is what was required. Other tools, like [fastp](https://github.com/OpenGene/fastp#base-correction-for-pe-data) will align the reads and attempt base correction.


# Options

The only _required_ options are a sequencing file and a file of adapter sequences. Of course, if you want to run with the `--paired_end` flag, you need both an R1 and an R2 file!

If you run `fast-adapter-trimming` with no options you will be provided with a reminder:

```
USAGE: search-paired-snp -1 -2 --primers -outputR1 --outputR2 --matchesR1 --matchesR2

Search for primers listed in --primers, allowing for 1-bp mismatches, against all the reads in --R1 and --R2
-1 --R1 R1 file (required)
-2 --R2 R2 file (required)
-f --primers fasta file of primers (required)
-p --outputR1 R1 output fastq file (will be gzip compressed)
-q --outputR2 R2 output fastq file (will be gzip compressed)
-j --matchesR1 Write the R1 matches to this file. Default: stdout
-k --matchesR2 Write the R2 matches to this file. Default: stdout
-n --noreverse Do not reverse the sequences
-a --adjustments Write the trimming adjustments here
--paired_end use a paired end (slower) search.
--primeroccurrences minimum number of times a primer was matched to include in the report
--nothreads use a single thread only. We typically want upto 4 threads to read and write R1 and R2 files
--verbose more output (but less than --debug)
--debug more more output
-v --version print the version and exit
```



Short option | Long option | Required? | Meaning
---|---|---|---
`-1` | `--R1` | Optional | The R1 (left) reads file. This can be gzip compressed or not compressed. Note that one R1 or R2 file is required, or else there is nothing to do.
`-2` | `--R2` | Optional | The R2 (right) reads file. This can be gzip compressed or not compressed.
`-f` | `--primers` | Required | A (typically) fasta file with adapters sequences. This can also be gzip compressed. For examples, see the [adapter](https://github.com/linsalrob/fast-adapter-trimming/tree/main/adapters) directory.
`-p` | `--outputR1` | Optional | Where to write the trimmed fastq reads from R1. This will be gzip compressed.
`-q` | `--outputR2` | Optional | Where to write the trimmed fastq reads from R2. This will be gzip compressed.
`-j` | `--matchesR1` |  Optional | Where to write a list of the adapters that match the R1 reads. This is a tab separated output of `adapter name`, `R1 sequence ID`, `matched position`, `offset from the right end`.
`-k` | `--matchesR2` | Optional | Where to write a list of the adapters that match the R2 reads. Same format as above.
`-n` | `--noreverse` | Optional | Only consider the forward direction of the adapers. By default we look for both the adapter sequences as specified in `--primers` and their reverse complement.
`-a` | `--adjustments` | Optional | Only valid with `--paired_end`. Where to write a summary of the adjustments to the R1 or R2 read trimming locations. If we find adapters in different locations in the R1 and R2 mate pairs, this file summarises the changes we made to accomodate both primers.
  &nbsp; | `--paired_end` | Optional | Use a paired end search which is slower and requires slightly more RAM.
  &nbsp; | `--primeroccurrences` | Optional | At the end we summarise the adapters that we found. This limits that output to those adapters found _n_ times or more. We often find one read that matches a single adapter (e.g. because there is a sequencing error), and so this just limits that output.
  &nbsp; | `--nothreads` | Optional | Only use a single thread for searching for the adapters.
  &nbsp; | `--versbose` | Optional | Write a lot more output
  &nbsp; | `--debug` | Optional | Write a lot, lot more output
`-v` | `--version` | Optional | Print the version and exit.




# How does it work?

We use 2-bit encoding of the DNA sequences to convert the sequences to a number:


Base | Two-bit encoding | Number
--- | --- | ---
A | 00 | 0
C | 01 | 1
G | 10 | 2
T | 11 | 3


We use a 64-bit integer and We can encode any sequence, upto 31 bases (because the encoding is 0-63 bits), using two bit encoding. By default, we start with the low bits (eg.0000000000000000000000000000000000000000000000000000000000000011 represents a single T). We have provided the vanity function `char* int_to_binary(uint64_t);` that will take any encoding and return a string representation of the binary code for that sequence.

We start by calculating the value for a 31bp sequence, using a technique called [bit-shifting](https://en.wikipedia.org/wiki/Bitwise_operation#Bit_shifts), essentially multiplying the number by 2 each time. So we start wirh our T, like about, and then add a C, so we get 0000000000000000000000000000000000000000000000000000000000001101.  Notice how the two 11's that represent the T have moved to the left by two spots, and then we added a 01 for the C.

So we calculate the numbers for all the adapter sequences and remember them. Then, we look through a file, and at every position we calculate the encoding for 31bp at a time. If that number is the same as one of the adapters, we have found the adapter, and can remove it. However, there is a cool CS trick here: we don't calculate that number for all 31 positions for every base in the sequence, We start at position 0 and calculate our number that represents from bases 0 to 30 (31 in total). Now for bases from 1 to 31 we remove the number for the base at position 0, move everything to the left two spots, and add the number for the base at position 31. Similarly, for the next string - positions 2-32 (inclusive) we remove the base at position 1 from the left of the number, move everything over, and add the number for the base at position 32. Thus, instead of doing 31 operations for each base in the sequence, we actually only do three (subtraction, bit-shifting, and addition).

The second step that we had to implement was to check for integers representing each of the variants in the adapter sequence, so we check for:


```
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
AATCGGAAGAGCACACGTCTGAACTCCAGTCAC
CATCGGAAGAGCACACGTCTGAACTCCAGTCAC
TATCGGAAGAGCACACGTCTGAACTCCAGTCAC
GCTCGGAAGAGCACACGTCTGAACTCCAGTCAC
GGTCGGAAGAGCACACGTCTGAACTCCAGTCAC
GTTCGGAAGAGCACACGTCTGAACTCCAGTCAC
...
```

To do this, we use a binary search tree. If you want more information about that, watch [Rob's YouTube channel](https://www.youtube.com/watch?v=lhTCSGRAlXI).

Now that we have found the adapters, we trim the sequences at those positions and print out the trimmed sequences. 

*Adjustments*. There is one more adjustment that we make. Refer to the figure above and it becomes apparent that the sequence on the left from the forward read and the sequence on the right from the reverse read that are complimentary to each other should be the same length. We compare those sequences and use this logic:

- If we find one adapter but not the other, we assume there is an additional error, so we trim both reads to the same point. This is a bit of a belt and braces option.
- If we find adapters in both reads, but at different positions, we trim to the shorter sequence. We haven't (yet) figured out why this happens, but it is only a really small subset of reads. We're working on that!


## Wow! How do I cite this amazing piece of work

Check our out DOI code and please cite it as:

```
Edwards, J.A. and Edwards, R.A. 2023 Fast Adapter Trimming: Efficient Removal of Sequencing Adapters. DOI: 
```



# Contributors

Most of the code was written by Rob Edwards with some help from John Edwards.

