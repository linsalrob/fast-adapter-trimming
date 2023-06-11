
#ifndef FAST_SEARCH_PRIMER_MATCH_COUNTS_H
#define FAST_SEARCH_PRIMER_MATCH_COUNTS_H

#include "structs.h"

/*
 * Add a primer match
 * id: primer id
 * before: the base preceeding (or NULL)
 * after: the base following (or NULL)
 */

void count_primer_occurrence(primer_counts_t *pc, char * id, char before, char after);


/* 
 * recursively print all the primer sequences
 */

void print_primers(primer_counts_t *pc, int);

#endif
