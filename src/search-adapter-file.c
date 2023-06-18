/*
 * Search for all adapters using a fasta file of adapters and paired end files
 *
 * (c) Rob. 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include "structs.h"
#include "search.h"
#include "colours.h"
#include "version.h"

void help() {
	printf("USAGE: search-paired-snp -1 -2 --primers -outputR1 --outputR2 --matchesR1 --matchesR2\n");
	printf("\nSearch for primers listed in %s--primers%s, allowing for 1-bp mismatches, against all the reads in %s--R1%s and %s--R2%s\n", 
			GREEN, ENDC, GREEN, ENDC, GREEN, ENDC);
	printf("-1 --R1 R%s1%s file (%srequired%s)\n", GREEN, ENDC, RED, ENDC);
	printf("-2 --R2 R%s2%s file (%srequired%s)\n", GREEN, ENDC, RED, ENDC);
	printf("-f --primers fasta file of primers (%srequired%s)\n", RED, ENDC);
	printf("-p --outputR1 R1 output fastq file (will be gzip compressed)\n");
	printf("-q --outputR2 R2 output fastq file (will be gzip compressed)\n");
	printf("-j --matchesR1 Write the R1 matches to this file. Default: stdout\n");
	printf("-k --matchesR2 Write the R2 matches to this file. Default: stdout\n");
	printf("-l --length Minimum sequence length (bp). Sequences shorter than this will be filtered out (Default 100)\n");
	printf("--noreverse Do not reverse the sequences\n");
	printf("--adjustments Write the trimming adjustments here\n");
	printf("--paired_end use a paired end (slower) search.\n");
	printf("--primeroccurrences minimum number of times a primer was matched to include in the report\n");
	printf("--nothreads use a single thread only. We typically want upto 4 threads to read and write R1 and R2 files\n");
	printf("--verbose more output (but less than --debug)\n");
	printf("--debug more more output\n");
	printf("-v --version print the version and exit\n");
}


int main(int argc, char* argv[]) {
	if (argc < 2) {
		help();
		exit(0);
	}

	if (argc == 2 && ((strcmp(argv[1], "-v") == 0) || (strcmp(argv[1], "--version") == 0))) {
		printf("%s version: %f\n", argv[0], __version__);
		exit(0);
	}

	struct options *opt;
	opt = malloc(sizeof(struct options));

	opt->tablesize = 17976787; // this is the next prime number after my largest sequence :)

	opt->R1_file = NULL;
	opt->R2_file = NULL;
	opt->R1_output = NULL;
	opt->R2_output = NULL;
	opt->R1_matches = NULL;
	opt->R2_matches = NULL;
	opt->min_sequence_length = 100;
	opt->primer_occurrences = 50;
	opt->reverse = true;
	opt->primers = NULL;
	opt->debug = false;
	opt->verbose = false;
	opt->adjustments = NULL;

	bool nothreads = false;
	bool paired_end = false;

	int gopt = 0;
	static struct option long_options[] = {
		{"R1",  required_argument, 0, '1'},
		{"R2",  required_argument, 0, '2'},
		{"primers",  required_argument, 0, 'f'},
		{"outputR1",  required_argument, 0, 'p'},
		{"outputR2",  required_argument, 0, 'q'},
		{"matchesR1",  required_argument, 0, 'j'},
		{"matchesR2",  required_argument, 0, 'k'},
		{"length", required_argument, 0, 'l'},
		{"primeroccurrences", required_argument, 0, 3},
		{"paired_end", no_argument, 0, 4},
		{"nothreads", no_argument, 0, 5},
		{"adjustments", required_argument, 0, 6},
		{"noreverse", required_argument, 0, 7},
		{"debug", no_argument, 0, 'd'},
		{"version", no_argument, 0, 'v'},
		{"verbose", no_argument, 0, 'b'},
		{0, 0, 0, 0}
	};
	int option_index = 0;
	while ((gopt = getopt_long(argc, argv, "1:2:p:q:f:j:k:l:bndv", long_options, &option_index )) != -1) {
		switch (gopt) {
			case '1' :
				opt->R1_file = strdup(optarg);
				break;
			case '2':
				opt->R2_file = strdup(optarg);
				break;
			case 'p':
				opt->R1_output = strdup(optarg);
				break;
			case 'q':
				opt->R2_output = strdup(optarg);
				break;
			case 'j' :
				opt->R1_matches = strdup(optarg);
				break;
			case 'k':
				opt->R2_matches = strdup(optarg);
				break;
			case 'f':
				opt->primers = strdup(optarg);
				break;
			case 'l':
				opt->min_sequence_length = atoi(optarg);
				break;
			case 'd': 
				opt->debug = true;
				break;
			case 'b':
				opt->verbose = true;
				break;
			case 'v':
				printf("Version: %f\n", __version__);
				return 0;
			case 3:
				opt->primer_occurrences = atoi(optarg);
				break;
			case 4:
				paired_end = true;
				break;
			case 5:
				nothreads = true;
				break;
			case 6:
				opt->adjustments = strdup(optarg);
				break;
			case 7:
				opt->reverse = false;
				break;
			default: help();
				 exit(EXIT_FAILURE);
		}
	}

	if (opt->R1_file == NULL && opt->R2_file == NULL) {
		fprintf(stderr, "Please provide at least one R1 or R2 read file\n");
		help();
		exit(EXIT_FAILURE);
	}

	if (opt->primers == NULL) {
		fprintf(stderr, "Please provide a primer file\n");
		help();
		exit(EXIT_FAILURE);
	}


	if (nothreads)
		fast_search(opt);
	else if (paired_end)
		paired_end_search(opt);
	else {
		pthread_t threads[2];
		thread_args_t *thread0_args;
		thread0_args = malloc(sizeof(thread_args_t));
		thread0_args->opt = opt;
		thread_args_t *thread1_args;
		thread1_args = malloc(sizeof(thread_args_t));
		thread1_args->opt = opt;
		// process R1
		if (opt->R1_file) {
			thread0_args->fqfile = strdup(opt->R1_file);
			if (opt->R1_matches)
				thread0_args->matches_file = strdup(opt->R1_matches);
			if (opt->R1_output)
				thread0_args->output_file = strdup(opt->R1_output);
			int result_code = pthread_create(&threads[0], NULL, &fast_search_one_file, (void *)thread0_args);
			if (result_code)
				fprintf(stderr, "%sERROR: Starting thread 0 returned the error code %d%s\n", RED, result_code, ENDC);
		}
		// process R2
		if (opt->R2_file) {
			thread1_args->fqfile = strdup(opt->R2_file);
			if (opt->R2_matches)
				thread1_args->matches_file = strdup(opt->R2_matches);
			if (opt->R2_output)
				thread1_args->output_file = strdup(opt->R2_output);
			int result_code = pthread_create(&threads[1], NULL, &fast_search_one_file, (void *)thread1_args);
			if (result_code)
				fprintf(stderr, "%sERROR: Starting thread 1 returned the error code %d%s\n", RED, result_code, ENDC);
		}
		if (opt->R1_file) {
			int result_code = pthread_join(threads[0], NULL);
			if (result_code)
				fprintf(stderr, "%sERROR: Joining thread 0 for it to finish returned the error code %d%s\n", RED, result_code, ENDC);
		}
		if (opt->R2_file) {
			int result_code = pthread_join(threads[1], NULL);
			if (result_code)
				fprintf(stderr, "%sERROR: Joining thread 1 for it to finish returned the error code %d%s\n", RED, result_code, ENDC);
		}
		free(thread0_args);
		free(thread1_args);
	}

	free(opt);
}

