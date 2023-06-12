#ifndef FAST_SEARCH_DEFS_H
#define FAST_SEARCH_DEFS_H


//  paired end search
void paired_end_search(struct options *opt);

// fast search without pairing
void fast_search(struct options *opt);

// search one file, for multi-threading
void *fast_search_one_file(void *);


#endif
