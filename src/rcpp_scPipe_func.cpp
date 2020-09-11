
#include <Rcpp.h>
#include "trimbarcode.h"
#include "Timer.h"
#include "transcriptmapping.h"


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_atac_trim_barcode(Rcpp::CharacterVector outfq,
                                 Rcpp::CharacterVector r1,
                                 Rcpp::CharacterVector r3,
                                 Rcpp::StringVector barcode_file,
                                 Rcpp::NumericVector start,
                                 Rcpp::NumericVector len,
                                 Rcpp::NumericVector umi_start,
                                 Rcpp::NumericVector umi_len,
                                 Rcpp::CharacterVector umi_in,
                                 Rcpp::LogicalVector write_gz) {
  
  std::string c_outfq = Rcpp::as<std::string>(outfq);
  std::string c_r1 = Rcpp::as<std::string>(r1);
  std::string c_r3 = Rcpp::as<std::string>(r3);
  std::string c_barcode_file = Rcpp::as<std::string>(barcode_file(0));
  //std::string c_barcode_file = Rcpp::as<std::string>(barcode_file);
  std::string c_umi_in = Rcpp::as<std::string>(umi_in);
  
  
  int c_start = Rcpp::as<int>(start);
  int c_len = Rcpp::as<int>(len);
  int c_umi_start = Rcpp::as<int>(umi_start);
  int c_umi_len = Rcpp::as<int>(umi_len);
  
  bool c_write_gz = Rcpp::as<bool>(write_gz);

  Timer timer;
  timer.start();
  
  paired_fastq_to_csv((char *)c_r1.c_str(), (char *)c_r3.c_str(), (char *)c_outfq.c_str(), (char *)c_barcode_file.c_str(), c_start, c_len, c_umi_start, c_umi_len, (char*)c_umi_in.c_str(), c_write_gz);
  
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_atac_trim_barcode_paired(Rcpp::CharacterVector outfq,
                                 Rcpp::CharacterVector r1,
                                 Rcpp::StringVector r2_list,
                                 Rcpp::CharacterVector r3,
                                 Rcpp::LogicalVector write_gz) {

  std::string c_outfq = Rcpp::as<std::string>(outfq);
  std::string c_r1 = Rcpp::as<std::string>(r1);
  std::vector<std::string> c_r2_list; 
  
  for( int i=0; i < r2_list.size(); i++ ){
    c_r2_list.push_back(Rcpp::as<std::string>(r2_list(i)));
  }
  
  std::string c_r3 = Rcpp::as<std::string>(r3);
  bool c_write_gz = Rcpp::as<bool>(write_gz);

  Timer timer;
  timer.start();

  paired_fastq_to_fastq((char *)c_r1.c_str(), c_r2_list, (char *)c_r3.c_str(), (char *)c_outfq.c_str(), c_write_gz);

  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
void rcpp_sc_atac_bam_tagging(Rcpp::CharacterVector inbam,
                          Rcpp::CharacterVector outbam,
                          Rcpp::CharacterVector bc,
                          Rcpp::CharacterVector mb,
                          Rcpp::NumericVector nthreads)
{
  std::string c_outbam = Rcpp::as<std::string>(outbam);

  std::string c_bc = Rcpp::as<std::string>(bc);
  std::string c_mb = Rcpp::as<std::string>(mb);
  std::vector<std::string> c_inbam_vec = Rcpp::as<std::vector<std::string>>(inbam);
  int c_nthreads = Rcpp::as<int>(nthreads);

  Mapping a = Mapping();
  Timer timer;
  timer.start();
  a.parse_align_warpper(c_inbam_vec, c_outbam, c_bc, c_mb, c_nthreads);
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}
