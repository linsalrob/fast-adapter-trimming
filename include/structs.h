#ifndef FAST_SEARCH_STRUCTS_H
#define FAST_SEARCH_STRUCTS_H

#include <stdbool.h>
#include <stdint.h>

/*
 * Structs that are used in searching the sequences
 *
 * The options that we need
 */

struct options {
	char* R1_file;
	char* R2_file;
	char* R1_output;
	char* R2_output;
	char* R1_matches;
	char* R2_matches;
	char* adjustments;
	int min_sequence_length;
	int min_adapter_length;
	char* primers;
	int primer_occurrences;
	bool reverse;
	int tablesize;
	bool verbose;
	bool debug;
};

/*
 * R1_read is a struct with the R1 name, the id string for the sequence, and whether there was a match to I7left
 * next is a pointer to the next R1_read element in the hash.
 */
struct R1_read {
	int trim;
	char *id;
	struct R1_read *next;
};

/*
 * This is an unbalanced binary search tree, and so could devolve into O(n)
 * performance, however with random ints it should be ~O(log n)
 */
typedef struct kmer_bst {
    uint64_t value;
    char* id;
    struct kmer_bst *bigger;
    struct kmer_bst *smaller;
} kmer_bst_t;

/*
 * Some counts and information about primer matches. Currently this is done
 * as an O(n) search each time we find a match because (a) we don't have a 
 * lot of primers (by name), and (b) we don't have a lot of matches.
 *
 * id: primer name
 * count: occurrence
 * before: array of counts of A, C, G, T as the preceeding base
 * after:  array of counts of A, C, G, T as the following base
 * next_primer: the next one in the list
 */

typedef struct primer_counts {
	char* id;
	int count;
	int before[5];
	int after[5];
	struct primer_counts *next_primer;
} primer_counts_t;

/*
 * Count of things that we find
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


/*
 * A struct to pass the data specifically to a pthread_create thread
 * for multi-threaded reading
 */
typedef struct thread_arg_struct {
	struct options *opt;
	char* fqfile;
	char* matches_file;
	char* output_file;
} thread_args_t;


#endif
