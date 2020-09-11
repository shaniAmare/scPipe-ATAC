#include <zlib.h> // for reading compressed .fq file
#include <string>
#include <stdio.h>
#include <iostream>
#include <Rcpp.h>
#include "config_hts.h"
#include "utils.h"

#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif

#ifndef TRIMBARCODE_H
#define TRIMBARCODE_H
#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

static const char empty_header[] = "@HD\tVN:1.4\tSO:unknown\n";

// Conversion functions
void paired_fastq_to_fastq(char *fq1_fn, std::vector<std::string> fq2_fn_list, char *fq3_fn,char *fq_out, const bool write_gz);
void paired_fastq_to_csv(char *fq1_fn, char *fq3_fn,char *fq_out, char *barcode_file, int start, int length, int umi_start, int umi_length, char* umi_in, const bool write_gz);

#endif
