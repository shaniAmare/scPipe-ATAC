//trim_barcode
#include "trimbarcode.h"
#include<string.h>

using namespace Rcpp;

void fq_write(std::ofstream& o_stream, kseq_t *seq, int trim_n)
{
  o_stream << "@" << seq->name.s << " " << seq->comment.s << "\n" << 
    (seq->seq.s+trim_n) << "\n" << 
      "+" << "\n" << 
        (seq->qual.s+trim_n) << "\n";
}

void fq_gz_write(gzFile out_file, kseq_t *seq, int trim_n) {
  std::stringstream stream;
  stream << "@" << seq->name.s << " " << seq->comment.s <<"\n" << 
    (seq->seq.s+trim_n) << "\n" << 
      "+" << "\n" << 
        (seq->qual.s+trim_n) << "\n";
  gzputs(out_file, stream.str().c_str());
}


// paired_fastq_to_fastq ------------------

void paired_fastq_to_fastq(
    char *fq1_fn,
    std::vector<std::string> fq2_fn_list,
    char *fq3_fn,
    char *fq_out,
    const bool write_gz
)
{
  int passed_reads = 0;
  int removed_have_N = 0;
  int removed_low_qual = 0;
  
  bool R3 = false;
  
  if(!strcmp(fq3_fn,"")){
    R3 = false;
  }else{
    R3 = true;
  }
  
  int l1 = 0;
  int l2 = 0;
  int l3 = 0;
  gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
  if (!fq1) {
    file_error(fq1_fn);
  }
  
  std::vector<kseq_t*> seq2_list;
  std::vector<gzFile> fq2_list;
  for(int i=0;i<fq2_fn_list.size();i++){
    char* fq2_fn = (char *)fq2_fn_list[i].c_str();
    gzFile fq2 = gzopen(fq2_fn, "r");
    if (!fq2) {
      file_error(fq2_fn);
    }
    fq2_list.push_back(fq2);
    seq2_list.push_back(kseq_init(fq2));
  }
  
  gzFile o_stream_gz_R1;
  std::ofstream o_stream_R1;
  
  gzFile fq3;
  gzFile o_stream_gz_R3;
  std::ofstream o_stream_R3;
  kseq_t *seq3;
  
  if(R3){
    fq3 = gzopen(fq3_fn, "r");
    if (!fq3) {
      file_error(fq3_fn);
    }
    
    const char* appendR3 = "/demux_" ;
    char *fqoutR3 = (char*)malloc(strlen(fq_out)+ strlen(appendR3) + strlen(getFileName(fq3_fn)) +  1 );
    strcpy(fqoutR3, fq_out);
    strcat(fqoutR3, appendR3);
    strcat(fqoutR3, getFileName(fq3_fn));
    seq3 =  kseq_init(fq3);
    if (write_gz){
      o_stream_gz_R3 = gzopen(fqoutR3, "wb2"); // open gz file
      if (!o_stream_gz_R3) {
        file_error(fqoutR3);
      }
    }else{
      o_stream_R3.open(fqoutR3); // output file
      if (!o_stream_R3.is_open()) {
        file_error(fqoutR3);
      }
    }
  }
  
  kseq_t *seq1;
  seq1 =  kseq_init(fq1);
  
  
  const char* appendR1 = "/demux_";
  char *fqoutR1 = (char*)malloc(strlen(fq_out) + strlen(appendR1) + strlen(getFileName(fq1_fn)) + 1);
  strcpy(fqoutR1, fq_out);
  strcat(fqoutR1,appendR1);
  strcat(fqoutR1, getFileName(fq1_fn));
  
  
  if (write_gz) {
    o_stream_gz_R1 = gzopen(fqoutR1, "wb2"); // open gz file
    if (!o_stream_gz_R1) {
      file_error(fqoutR1);
    }
    
  } else {
    o_stream_R1.open(fqoutR1); // output file
    if (!o_stream_R1.is_open()) {
      file_error(fqoutR1);
    }
  }
  
  size_t _interrupt_ind = 0;
  // Assuming R1, R2, R3 all are of equal lengths.
  while (((l1 = kseq_read(seq1)) >= 0))
  {
    if (++_interrupt_ind % 4096 == 0) checkUserInterrupt();
    passed_reads++;
    
    
    char * const seq1_name = seq1->name.s;
    char * const seq1_seq = seq1->seq.s;
    int seq1_namelen = seq1->name.l;
    int seq1_seqlen = seq1->seq.l;
    
    char * seq3_name;
    char * seq3_seq;
    int seq3_namelen;
    int seq3_seqlen;
    if (R3){
      if(l3 = kseq_read(seq3) >= 0){
        seq3_name = seq3->name.s;
        seq3_seq = seq3->seq.s;
        seq3_namelen = seq3->name.l;
        seq3_seqlen = seq3->seq.l;
      }
      else{
        Rcpp::Rcout << "read2 file is not of the same length as the barcode fastq file: " << "\n";
      }
    } 
    
    for(int i=0;i<seq2_list.size();i++){
      kseq_t* seq2 = seq2_list[i];
      if (l2 = kseq_read(seq2) >= 0){
        char * const seq2_name = seq2->name.s;
        char * const seq2_seq = seq2->seq.s;
        int seq2_namelen = seq2->name.l;
        int seq2_seqlen = seq2->seq.l;
        
        const int new_name_length1 = seq1_namelen + seq2_seqlen+1;
        seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1); // allocate additional memory
        memmove(seq1_name + seq2_seqlen+1, seq1_name, seq1_namelen * sizeof(char));// move original read name
        memcpy(seq1_name, seq2_seq, seq2_seqlen * sizeof(char)); // copy index one
        seq1_name[seq2_seqlen] = '#'; // add separator
        seq1_name[new_name_length1] = '\0';
        
        if(R3){
          const int new_name_length2 = seq3_namelen + seq2_seqlen+1;
          seq3->name.s = (char*)realloc(seq3->name.s, new_name_length2); // allocate additional memory
          memmove(seq3_name + seq2_seqlen+1, seq3_name, seq3_namelen * sizeof(char));// move original read name
          memcpy(seq3_name, seq2_seq, seq2_seqlen * sizeof(char)); // copy index one
          seq3_name[seq2_seqlen] = '#'; // add separator  
          seq3_name[new_name_length1] = '\0';
          
        }
      }else{
        Rcpp::Rcout << "read1 file is not the same length as the barcode fastq file: " << "\n";
      }
    }
    
    if (write_gz) {
      fq_gz_write(o_stream_gz_R1, seq1, 0); // write to gzipped fastq file
      if(R3){
        fq_gz_write(o_stream_gz_R3, seq3, 0); // write to gzipped fastq file
      }
    } else {
      fq_write(o_stream_R1, seq1, 0); // write to fastq file
      if(R3){
        fq_write(o_stream_R3, seq3, 0); // write to fastq file
      }
    }
    
  }
  
  // free subStr1
  kseq_destroy(seq1); 
  for(int i=0;i<seq2_list.size();i++){
    kseq_t* seq2 = seq2_list[i];
    kseq_destroy(seq2);// free seq
    gzclose(fq2_list[i]);
  }
  gzclose(fq1); // close fastq file
  if(R3){
    kseq_destroy(seq3);
    gzclose(fq3);
  }
  if (write_gz){
    gzclose(o_stream_gz_R1);
    if(R3){
      gzclose(o_stream_gz_R3);
    }
  }
  Rcpp::Rcout << "Total Reads: " << passed_reads << "\n";
  //Rcpp::Rcout << "removed_have_N: " << removed_have_N << "\n";
  //Rcpp::Rcout << "removed_low_qual: " << removed_low_qual << "\n";
}

// paired_fastq_to_csv ------------------

void paired_fastq_to_csv(
    char *fq1_fn,
    char *fq3_fn,
    char *fq_out, 
    char *bc_fn, 
    int start,
    int length, 
    int umi_start,
    int umi_length,
    char *umi_in,
    const bool write_gz
)
{
  int passed_reads = 0;
  int removed_have_N = 0;
  int removed_low_qual = 0;
  int exact_match = 0;
  int approx_match = 0;
  
  bool R3 = false;
  bool isUMIR1 = (umi_length > 0 && (strcmp(umi_in,"both") == 0 || strcmp(umi_in,"R1") == 0));
  bool isUMIR2 = (umi_length > 0 && (strcmp(umi_in,"both") == 0 || strcmp(umi_in,"R2") == 0));
    
  if(!strcmp(fq3_fn,"")){
    R3 = false;
  }else{
    R3 = true;
  }
  
  int l1 = 0;
  int l2 = 0;
  int l3 = 0;
  gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
  if (!fq1) {
    file_error(fq1_fn);
  }
  
  
  std::map<std::string, int> barcode_map; 
  std::ifstream bc(bc_fn);
  std::string line;
  if(!bc.is_open()) throw std::runtime_error("Could not open file");
  while(std::getline(bc, line))
  {
    std::stringstream s_stream(line);
    std::string bcode;
    if(s_stream.good()) {
      getline(s_stream, bcode, ','); //get first string delimited by comma
    }
    
    if(s_stream.good()) {
      getline(s_stream, bcode, ','); //get second string delimited by comma
      std::string substr = bcode.substr(0,length);
      barcode_map.insert(std::pair<std::string, int>(substr, 1)); 
    }
    
  }
  
  if(barcode_map.empty()){
    std::stringstream err_msg;
    err_msg << "Error in retrieving barcodes from the barcode File. Please check the barcode file format. " << bc_fn << "\n";
    Rcpp::stop(err_msg.str());
  }
  
  
  
  gzFile fq3;
  gzFile o_stream_gz_R3;
  std::ofstream o_stream_R3;
  gzFile o_stream_gz_R3_Partial;
  std::ofstream o_stream_R3_Partial;
  gzFile o_stream_gz_R3_No;
  std::ofstream o_stream_R3_No;
  
  kseq_t *seq3;
  const char* appendCompleteMatch = "/demultiplexed_completematch_";
  const char* appendPartialMatch = "/demultiplexed_partialmatch_";
  const char* appendNoMatch = "/demultiplexed_nomatch_";
  
  if(R3){
    fq3 = gzopen(fq3_fn, "r");
    if (!fq3) {
      file_error(fq3_fn);
    }
    
    char *fqoutR3 = createFileWithAppend(fq_out,appendCompleteMatch,fq3_fn);
    openFile(o_stream_gz_R3,o_stream_R3,fqoutR3, write_gz);
    
    char *fqoutR3Partial = createFileWithAppend(fq_out,appendPartialMatch,fq3_fn);
    openFile(o_stream_gz_R3_Partial,o_stream_R3_Partial,fqoutR3Partial, write_gz);
    
    char *fqoutR3No = createFileWithAppend(fq_out,appendNoMatch,fq3_fn);
    openFile(o_stream_gz_R3_No,o_stream_R3_No,fqoutR3No, write_gz);
    
    seq3 =  kseq_init(fq3);
    
  }
  
  kseq_t *seq1;
  seq1 =  kseq_init(fq1);
  
  
  gzFile o_stream_gz_R1;
  std::ofstream o_stream_R1;
  char *fqoutR1 = createFileWithAppend(fq_out,appendCompleteMatch,fq1_fn);
  openFile(o_stream_gz_R1,o_stream_R1,fqoutR1, write_gz);
  
  gzFile o_stream_gz_R1_Partial;
  std::ofstream o_stream_R1_Partial;
  char *fqoutR1Partial = createFileWithAppend(fq_out,appendPartialMatch,fq1_fn);
  openFile(o_stream_gz_R1_Partial,o_stream_R1_Partial,fqoutR1Partial, write_gz);
  
  
  gzFile o_stream_gz_R1_No;
  std::ofstream o_stream_R1_No;
  char *fqoutR1No = createFileWithAppend(fq_out,appendNoMatch,fq1_fn);
  openFile(o_stream_gz_R1_No,o_stream_R1_No,fqoutR1No, write_gz);
  
  
  size_t _interrupt_ind = 0;
  // Assuming R1, R2, R3 all are of equal lengths.
  while (((l1 = kseq_read(seq1)) >= 0))
  {
    if (++_interrupt_ind % 4096 == 0) checkUserInterrupt();
    passed_reads++;
    
    
    char *seq1_name = seq1->name.s;
    char *seq1_seq = seq1->seq.s;
    int seq1_namelen = seq1->name.l;
    int seq1_seqlen = seq1->seq.l;
    
    char * seq3_name;
    char * seq3_seq;
    int seq3_namelen;
    int seq3_seqlen;
    if (R3){
      if(l3 = kseq_read(seq3) >= 0){
        seq3_name = seq3->name.s;
        seq3_seq = seq3->seq.s;
        seq3_namelen = seq3->name.l;
        seq3_seqlen = seq3->seq.l;
      }
      else{
        Rcpp::Rcout << "R3 file is not of same length as R2: " << "\n";
      }
    } 
    
    
    char* subStr1;
    int bcUMIlen1 = 0;
    if(isUMIR1){
      subStr1 = (char*)malloc(length+umi_length+1);
      memcpy(subStr1, &seq1_seq[start], length );
      subStr1 [length] = '_';
      memcpy(subStr1+length+1, &seq1_seq[umi_start], umi_length);
      subStr1[length + umi_length +1] = '\0';
      bcUMIlen1 = length + umi_length +1;
      
    }else{
      subStr1 = (char*)malloc(length);
      memcpy( subStr1, &seq1_seq[start], length );
      subStr1[length] = '\0';
      bcUMIlen1 =length;
    }
    
    char barcode[length];
    memcpy( barcode, &seq1_seq[start], length );
    barcode[length] = '\0';
    
    
    
    if ( barcode_map.find(barcode) == barcode_map.end() ) {
      for (std::map<std::string,int>::iterator it=barcode_map.begin(); it!=barcode_map.end(); ++it){
        if(hamming_distance(it->first, barcode) <2){
          approx_match++;
          const int new_name_length1 = seq1_namelen + bcUMIlen1 + 1;
          seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1); // allocate additional memory
          memmove(seq1_name + bcUMIlen1+1, seq1_name, seq1_namelen );// move original read name
          memcpy(seq1_name, subStr1, bcUMIlen1 * sizeof(char)); // copy index one
          seq1_name[bcUMIlen1] = '#'; // add separator
          seq1_name[new_name_length1] = '\0';
          
          if ( (umi_start + umi_length) > (start + length) && isUMIR1){
            memmove (seq1_seq, seq1_seq + umi_start + umi_length, seq1_seqlen-umi_start-umi_length+1); 
          }else{
            memmove (seq1_seq, seq1_seq + start + length, seq1_seqlen-start-length+1); 
          }
          
          if (write_gz) {
            fq_gz_write(o_stream_gz_R1_Partial, seq1, 0); // write to gzipped fastq file
          } else {
            fq_write(o_stream_R1_Partial, seq1, 0); // write to fastq file
          }
          
        }else{
          if (write_gz) {
            fq_gz_write(o_stream_gz_R1_No, seq1, 0); // write to gzipped fastq file
          } else {
            fq_write(o_stream_R1_No, seq1, 0); // write to fastq file
          }
        }
      }
    } else {
      
      exact_match ++;   
      const int new_name_length1 = seq1_namelen + bcUMIlen1 + 1;
      seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1); // allocate additional memory
      memmove(seq1_name + bcUMIlen1+1, seq1_name, seq1_namelen );// move original read name
      memcpy(seq1_name, subStr1, bcUMIlen1 * sizeof(char)); // copy index one
      seq1_name[bcUMIlen1] = '#'; // add separator
      seq1_name[new_name_length1] = '\0';
      
      if ( (umi_start + umi_length) > (start + length) && isUMIR1){
        memmove (seq1_seq, seq1_seq + umi_start + umi_length, seq1_seqlen-umi_start-umi_length+1); 
      }else{
        memmove (seq1_seq, seq1_seq + start + length, seq1_seqlen-start-length+1); 
      }
      
      if (write_gz) {
        fq_gz_write(o_stream_gz_R1, seq1, 0); // write to gzipped fastq file
      } else {
        fq_write(o_stream_R1, seq1, 0); // write to fastq file
      }
      
    }
    
    
    if(R3){
      char* subStr3;
      int bcUMIlen3 = 0;
      if(isUMIR2){
        subStr3 = (char*)malloc(length+umi_length+1);
        memcpy(subStr3, &seq3_seq[start], length );
        subStr3 [length] = '_';
        memcpy(subStr3+length+1, &seq3_seq[umi_start], umi_length);
        subStr3[length + umi_length +1] = '\0';
        bcUMIlen3 = length + umi_length +1;
        
      }else{
        subStr3 = (char*)malloc(length);
        memcpy( subStr3, &seq3_seq[start], length );
        subStr3[length] = '\0';
        bcUMIlen3 =length;
      }
      
      char barcode[length];
      memcpy( barcode, &seq3_seq[start], length );
      barcode[length] = '\0';
      
      if ( barcode_map.find(subStr3) == barcode_map.end() ) {
        for (std::map<std::string,int>::iterator it=barcode_map.begin(); it!=barcode_map.end(); ++it){
          if(hamming_distance(it->first, barcode) < 2){
            const int new_name_length1 = seq3_namelen + bcUMIlen3+1;
            seq3->name.s = (char*)realloc(seq3->name.s, new_name_length1); // allocate additional memory
            memmove(seq3_name + bcUMIlen3+1, seq3_name, seq3_namelen * sizeof(char) );// move original read name
            memcpy(seq3_name, subStr3, bcUMIlen3 * sizeof(char)); // copy index one
            seq3_name[bcUMIlen3] = '#'; // add separator
            seq3_name[new_name_length1] = '\0';
            
            if ( (umi_start + umi_length) > (start + length) && isUMIR2){
              memmove (seq3_seq, seq3_seq + umi_start + umi_length, seq3_seqlen-umi_start-umi_length+1); 
            }else{
              memmove (seq3_seq, seq3_seq + start + length, seq3_seqlen-start-length+1); 
            }
            
            if (write_gz) {
              fq_gz_write(o_stream_gz_R3_Partial, seq3, 0); // write to gzipped fastq file
            } else {
              fq_write(o_stream_R3_Partial, seq3, 0); // write to fastq file
            }
          }else{
            if (write_gz) {
              fq_gz_write(o_stream_gz_R3_No, seq3, 0); // write to gzipped fastq file
            } else {
              fq_write(o_stream_R3_No, seq3, 0); // write to fastq file
            }
          }
        }
      } else {
        const int new_name_length1 = seq3_namelen + bcUMIlen3+1;
        seq3->name.s = (char*)realloc(seq3->name.s, new_name_length1); // allocate additional memory
        memmove(seq3_name + bcUMIlen3+1, seq3_name, seq3_namelen * sizeof(char) );// move original read name
        memcpy(seq3_name, subStr3, bcUMIlen3 * sizeof(char)); // copy index one
        seq3_name[bcUMIlen3] = '#'; // add separator
        seq3_name[new_name_length1] = '\0';
        
        if ( (umi_start + umi_length) > (start + length) && isUMIR2){
          memmove (seq3_seq, seq3_seq + umi_start + umi_length, seq3_seqlen-umi_start-umi_length+1); 
        }else{
          memmove (seq3_seq, seq3_seq + start + length, seq3_seqlen-start-length+1); 
        }
        
        if (write_gz) {
          fq_gz_write(o_stream_gz_R3, seq3, 0); // write to gzipped fastq file
        } else {
          fq_write(o_stream_R3, seq3, 0); // write to fastq file
        }
        
      }
      
      free(subStr3);
      
    }
    free(subStr1);
    
  }
  
  kseq_destroy(seq1); 
  
  bc.close();
  gzclose(fq1); // close fastq file
  if(R3){
    kseq_destroy(seq3);
    gzclose(fq3);
  }
  if (write_gz){
    gzclose(o_stream_gz_R1);
    gzclose(o_stream_gz_R1_Partial);
    gzclose(o_stream_gz_R1_No);
    
    if(R3){
      gzclose(o_stream_gz_R3);
      gzclose(o_stream_gz_R3_Partial);
      gzclose(o_stream_gz_R3_No);
      
    }
  }
  Rcpp::Rcout << "Total Reads: " << passed_reads << "\n";
  Rcpp::Rcout << "Exact match Reads: " << exact_match << "\n";
  Rcpp::Rcout << "Approx Match Reads: " << approx_match << "\n";
  
  //Rcpp::Rcout << "removed_have_N: " << removed_have_N << "\n";
  //Rcpp::Rcout << "removed_low_qual: " << removed_low_qual << "\n";
}

