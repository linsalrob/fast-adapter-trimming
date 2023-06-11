
/*
 * A fast search that does not compare R1 and R2 but only trims the most 5' adapter sequence for each read.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#include "colours.h"
#include "create-snps.h"
#include "definitions.h"
#include "hash.h"
#include "kseq.h"
#include "primer-match-counts.h"
#include "primers.h"
#include "print-sequences.h"
#include "rob_dna.h"
#include "search.h"
#include "seqs_to_ints.h"
#include "structs.h"
#include "version.h"


KSEQ_INIT(gzFile, gzread);


void fast_search(struct options *opt) {
	/*
	 * opt contains our variables for this search:
	 * 	opt->primers = name of a file of primers
	 * 	opt->reverse = include the reverse complement of the primers
	 * 	opt->verbose = sometimes more output!
	 */	

	typedef struct COUNTS {
		int R1_seqs;
		int R2_seqs;
		int R1_found;
		int R2_found;
		int R1_adjusted;
		int R2_adjusted;
		int R1_trimmed;
		int R2_trimmed;
		int same;
	} COUNTS;

	COUNTS counts = {};

	// create an array of kmer_bsts and malloc them
	kmer_bst_t *all_primers[MAXKMER+1];

	for (int i = 0; i<=MAXKMER; i++) {
		all_primers[i] = malloc(sizeof(kmer_bst_t));

		if (all_primers[i] == NULL) {
			fprintf(stderr, "Can't malloc memory for all_primers %d\n", i);
			exit(1);
		}
		all_primers[i]->bigger = NULL;
		all_primers[i]->smaller = NULL;
		all_primers[i]->value = 0;
		all_primers[i]->id = "";
	}

	// read the primer file
	read_primers_create_snps(opt->primers, all_primers, opt->reverse, opt->verbose);

	if (opt->debug) {	
		fprintf(stderr, "%sWe have read the primers%s\n", GREEN, ENDC);
		for (int i = 0; i<=MAXKMER; i++)
			print_all_primers(all_primers[i], i);
	}

	// A place to store our encdoded kmers
	// 
	// First, we make an array with just the lengths of the kmers we need to encode
	// and then we encode those kmers. This reduces our search from ~30 to ~3
	// (depending on how many primer lengths we have).
	// We also do this longer primers to shorter, so that we initially trim off the
	// longest possible primers

	int kmer_lengths[MAXKMER];
	int unique_kmer_count = 0;
	for (int i=MAXKMER; i>=0; i--) 
		if (all_primers[i]->value > 0) 
			kmer_lengths[unique_kmer_count++] = i; 	// we need to remember this kmer length

	// now we just add those kmers to this array and then we test them each time
	uint64_t encoded_kmers[unique_kmer_count];

	// Initialize a primer count structure
	primer_counts_t *pc;
	pc = malloc(sizeof(primer_counts_t));
	pc->id = NULL;
	pc->count = 0;
	pc->next_primer = NULL;
	for (int i=0; i<5; i++) {
		pc->before[i] = 0;
		pc->after[i] = 0;
	}

	// Step 1. Read the R1 file and find the matches to any primer
	if (opt->R1_file) {
		if( access( opt->R1_file, R_OK ) == -1 ) {
			// file doesn't exist
			fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, opt->R1_file, ENDC);
			return;
		}

		if (opt->verbose)
			fprintf(stderr, "%sReading %s%s\n", GREEN, opt->R1_file, ENDC);

		gzFile fp1 = gzopen(opt->R1_file, "r");
		if (fp1 == NULL) {
			fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R1_file, ENDC);
			exit(3);
		}
		kseq_t *seq = kseq_init(fp1);
		int l;

		FILE *match_out = NULL;
		if (opt->R1_matches)
			match_out = fopen(opt->R1_matches, "w");

		// do we need to write to R1
		FILE *pipe = NULL;
		if (opt->R1_output) {
			char* pipe_file = malloc(sizeof(char) * (strlen(opt->R1_output) + 10));
			strcpy(pipe_file, "gzip - > ");
			strcat(pipe_file, opt->R1_output);
			pipe = popen(pipe_file, "w");
			free(pipe_file);
		}

		bool warning_printed = false;

		while ((l = kseq_read(seq)) >= 0) {
			counts.R1_seqs++;
			if (opt->debug)
				fprintf(stderr, "Reading %s\n", seq->name.s);

			// housekeeping warnings and definitions
			if (opt->verbose && !warning_printed && has_n(seq->seq.s)) {
				fprintf(stderr, "%sWARNING: sequences have an N but we don't deal with them. They are encoded as A%s\n", BLUE, ENDC);
				warning_printed = true;
			}

			int trim = -1;
			char *primerid;
			primerid = malloc(sizeof(char) * MAXLINELEN);
			char before;
			char after;

			// end housekeeping warnings and definitions

			//  encode the first kmers in the sequence
			for (int i = 0; i<unique_kmer_count; i++)
				encoded_kmers[i] = kmer_encoding(seq->seq.s, 0, kmer_lengths[i]);

			// test for the first kmers in our data structure. We need to iterate the 
			// array of encoded primers and see if we find an encoding in all_primers[k]
			// for k being the length of the encoding
			for (int i=0; i<unique_kmer_count; i++) {
				uint64_t enc  = encoded_kmers[i];
				kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);
				if (ks) {
					strcpy(primerid, ks->id);
					before = '^';
					after = seq->seq.s[kmer_lengths[i]+1];
					trim = 0;
				} 
			}


			for (int posn=1; posn<seq->seq.l - MAXKMER + 1; posn++) {
				for (int i=0; i<unique_kmer_count; i++) {
					// calculate the next encoding for this kmer length
					uint64_t enc  = next_kmer_encoding(seq->seq.s, posn, kmer_lengths[i], encoded_kmers[i]);
					encoded_kmers[i] = enc; // remember it for next time!
					kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);
					if (ks) {
						if (trim == -1 || posn < trim) {
							strcpy(primerid, ks->id);
							before = seq->seq.s[posn-1];
							after = seq->seq.s[kmer_lengths[i]+1];
							trim = posn;
						}
					}
				}
			}

			if (trim > -1) {
				counts.R1_found++;
				count_primer_occurrence(pc, primerid, before, after);
				if (opt->R1_matches)
					fprintf(match_out, "R1\t%s\t%s\t%d\t-%ld\n", primerid, seq->name.s, trim, strlen(seq->seq.s)-trim);
				seq->seq.s[trim] = '\0';
				seq->qual.s[trim] = '\0';
				counts.R1_trimmed++;
			}
			if (pipe)
				fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);

		}

		// I am going to reset kseq so we have to initiate it again later
		kseq_destroy(seq);
		gzclose(fp1);

		if (opt->R1_matches)
			fclose(match_out);

		if (pipe)
			pclose(pipe);
	}
	// Step 2. Read the R2 file and find the locations of any of the primers.

	if (opt->R2_file) {
		if( access( opt->R2_file, R_OK ) == -1 ) {
			// file doesn't exist
			fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, opt->R2_file, ENDC);
			return;
		}

		if (opt->verbose)
			fprintf(stderr, "%sReading %s%s\n", GREEN, opt->R2_file, ENDC);

		// Open R2 for reading
		gzFile fp2 = gzopen(opt->R2_file, "r");
		if (fp2 == NULL) {
			fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R2_file, ENDC);
			exit(3);
		}
		kseq_t *seq = kseq_init(fp2);

		// if we want to write the files, we open a pipe
		// otherwise it is null. so we just need to check before writing
		// Open R2 for writing
		FILE *pipe = NULL;
		if (opt->R2_output) {
			char* pipe_file = malloc(sizeof(char) * (strlen(opt->R2_output) + 10));
			strcpy(pipe_file, "gzip - > ");
			strcat(pipe_file, opt->R2_output);
			pipe = popen(pipe_file, "w");
			free(pipe_file);
		}

		// open our log files
		FILE *match_out = NULL;
		if (opt->R2_matches)
			match_out = fopen(opt->R2_matches, "w");

		int l;
		while ((l = kseq_read(seq)) >= 0) {
			counts.R2_seqs++;
			//  encode the first kmers in the sequence
			for (int i = 0; i<unique_kmer_count; i++)
				encoded_kmers[i] = kmer_encoding(seq->seq.s, 0, kmer_lengths[i]);

			int trim = -1;
			char *primerid;
			primerid = malloc(sizeof(char) * MAXLINELEN);
			char before;
			char after;

			for (int i=0; i<unique_kmer_count; i++) {
				uint64_t enc  = encoded_kmers[i];
				kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);

				if (ks) {
					counts.R2_found++;
					strcpy(primerid, ks->id);
					before = '^';
					after = seq->seq.s[kmer_lengths[i]+1];
					trim = 0;
				} 
			}

			for (int posn=1; posn<seq->seq.l - MAXKMER + 1; posn++) {
				for (int i=0; i<unique_kmer_count; i++) {
					// calculate the next encoding for this kmer length
					uint64_t enc  = next_kmer_encoding(seq->seq.s, posn, kmer_lengths[i], encoded_kmers[i]);
					encoded_kmers[i] = enc; // remember it for next time!
					kmer_bst_t *ks = find_primer(enc, all_primers[kmer_lengths[i]]);
					if (ks) {
						if (trim == -1 || trim < posn)	{
							strcpy(primerid, ks->id);
							before = seq->seq.s[posn-1];
							after = seq->seq.s[kmer_lengths[i]+1];
							trim = posn;
						}
					}
				}
			}
			if (trim > -1) {
				counts.R2_found++;
				count_primer_occurrence(pc, primerid, before, after);
				if (opt->R2_matches)
					fprintf(match_out, "R2\t%s\t%s\t%d\t-%ld\n", primerid, seq->name.s, trim, strlen(seq->seq.s)-trim);
				seq->seq.s[trim] = '\0';
				seq->qual.s[trim] = '\0';
				counts.R2_trimmed++;
			}
			if (pipe)
				fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);


		}
		if (pipe)
			pclose(pipe);


		if (opt->R2_matches)
			fclose(match_out);
		kseq_destroy(seq);
		gzclose(fp2);
	}
	for (int i=0; i<=MAXKMER; i++)
		free(all_primers[i]);

	printf("Total sequences: R1 %d R2 %d\n", counts.R1_seqs, counts.R2_seqs);
	printf("Primer found: R1 %d R2 %d\n", counts.R1_found, counts.R2_found);
	printf("Sequences trimmed: R1 %d R2 %d\n", counts.R1_trimmed, counts.R2_trimmed);


	printf("\nAdapter occurrences:\n");
	print_primers(pc, opt->primer_occurrences);
}

