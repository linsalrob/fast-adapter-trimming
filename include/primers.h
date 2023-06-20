
#ifndef FAST_SEARCH_PRIMERS_H
#define FAST_SEARCH_PRIMERS_H


#include "structs.h"

/* 
 * Add a sequence encoding (a long long int) to a binary
 * search tree.
 */
void add_primer(uint64_t, char*, kmer_bst_t*);

/*
 * Find an encoding in a bst
 */
kmer_bst_t* find_primer(uint64_t, kmer_bst_t*);

/*
 * recursively print all primers
 */

void print_all_primers(kmer_bst_t*, int);

// read the primers and populate the kmer_bst_t
void read_primers(char*, kmer_bst_t**, int, bool, int);

// read the priemrs and populate the kmer_bst_t with a snp in every position
void read_primers_create_snps(char*, kmer_bst_t**, int, bool, int);

// read the primers but only store a substring of each
void read_trunc_primers(char*, int, kmer_bst_t*, bool, int);

#endif
