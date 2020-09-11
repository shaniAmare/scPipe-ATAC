// transcriptmapping.cpp
#include "transcriptmapping.h"

using std::atoi;
using std::atomic;
using std::endl;
using std::fixed;
using std::getline;
using std::ifstream;
using std::ostream;
using std::setprecision;
using std::sort;
using std::string;
using std::stringstream;
using std::thread;
using std::unordered_map;
using std::vector;

using namespace std::this_thread;
using namespace std::chrono;
using namespace Rcpp;

namespace {
    void report_every_3_mins(
        atomic<unsigned long long> &cnt,
        atomic<bool> &running,
        atomic<bool> &report_message
    )
    {
        do {
            // sleep thread for a total of 3 minutes (180 seconds)
            // wake up at shorter intervals to check if process has stopped running
            for (int i = 0; i < 180; i++)
            {
                sleep_for(seconds(1));
                if (!running)
                {
                    break;
                }
            }

            report_message = true;
        } while (running);
    }
}

void Mapping::parse_align_warpper(vector<string> fn_vec, string fn_out, string cellular_tag, string molecular_tag, int nthreads)
{
  if (fn_vec.size()>1)
  {
    parse_align(fn_vec[0], fn_out, cellular_tag, molecular_tag, nthreads);
      for (int i=1;i<fn_vec.size();i++)
      {
        parse_align(fn_vec[i], fn_out, cellular_tag, molecular_tag, nthreads);
      }
  }
  else
  {
    parse_align(fn_vec[0], fn_out, cellular_tag, molecular_tag, nthreads);
  }
}

namespace {
    std::pair<int, int> get_bc_umi_lengths(string bam_fn) {
        BGZF *fp = bgzf_open(bam_fn.c_str(), "r"); // input file
        bam_hdr_t *bam_hdr = bam_hdr_read(fp);

        bam1_t *bam_record = bam_init1();

        if (bam_read1(fp, bam_record) >= 0) {
            string read_header = bam_get_qname(bam_record);
            int break_pos = read_header.find("#");
            // start from 1 to exclude @
            string first_section = read_header.substr(1, break_pos);
            int bc_len = first_section.length();
            int umi_len = 0;
            // header has structure {BARCODE}_{UMI}
            if(first_section.find("_") != -1){
              bc_len = first_section.find("_") + 1;
              umi_len = first_section.length() - bc_len - 1;
            }
            return std::make_pair(bc_len, umi_len);
        }
        else
        {
            throw std::runtime_error("BAM file reading failed.");
        }
    }
}

void Mapping::parse_align(string bam_fn, string fn_out, string cellular_tag, string molecular_tag, int nthreads)
{
    check_file_exists(bam_fn); // htslib does not check if file exist so we do it manually

    int bc_len = 0;
    int UMI_len = 0;
    std::tie(bc_len, UMI_len) = get_bc_umi_lengths(bam_fn);

    Rcout << "Detected bc_len: " << bc_len << "  Detected UMI len:  " << UMI_len  << "\n";
    
    string write_mode = "wb";
    const char * c_write_mode = write_mode.c_str();
    // open files
    bam1_t *b = bam_init1();
    BGZF *fp = bgzf_open(bam_fn.c_str(), "r"); // input file
    samFile *of = sam_open(fn_out.c_str(), c_write_mode); // output file

    // set up htslib threadpool for output
    int out_threads = std::max(nthreads - 1, 1);
    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(out_threads);
    hts_set_opt(of, HTS_OPT_THREAD_POOL, &p);

    int hts_retcode;

    bam_hdr_t *header = bam_hdr_read(fp);
    hts_retcode = sam_hdr_write(of, header);

    int tmp_c[4] = {0,0,0,0};

    const char * c_ptr = cellular_tag.c_str();
    const char * m_ptr = molecular_tag.c_str();
    char buf[999] = ""; // assume the length of barcode or UMI is less than 999

    atomic<unsigned long long> cnt{0};
    atomic<bool> running{true};
    atomic<bool> report_message{false};

    Rcout << "updating progress every 3 minutes..." << "\n";
    // spawn thread to set report_message to true every 3 minutes
    thread reporter_thread(
        [&]() {
            report_every_3_mins(cnt, running, report_message);
        }
    );
    Timer timer;
    timer.start();

    while (bam_read1(fp, b) >= 0)
    {
        string gene_id;

            if (cnt % 1000000 == 0)
            {
                Rcout << "number of read processed:" << cnt << "\n";
                Rcout << tmp_c[0] <<"\t"<< tmp_c[1] <<"\t"<<tmp_c[2] <<"\t"<<tmp_c[3] <<"\t" << "\n";
            }
        
        cnt++;
        if (cnt % 32768 == 0) checkUserInterrupt();

        // The Rcout would be conceptually cleaner if it lived inside the spawned thread
        // but only the master thread can interact with R without error so this code CANNOT
        // be run inside a child thread
        if (report_message)
        {
            Rcout
                << cnt << " reads processed" << ", "
                << cnt / timer.seconds_elapsed() / 1000 << "k reads/sec" << endl;
            report_message = false;
        }
        
        if (bc_len > 0)
        {
            memcpy(buf, bam_get_qname(b), bc_len * sizeof(char));
            buf[bc_len] = '\0';
            bam_aux_append(b, c_ptr, 'Z', bc_len+1, (uint8_t*)buf);
        }
        
        if (UMI_len > 0)
        {
            memcpy(buf, bam_get_qname(b)+bc_len+1, UMI_len * sizeof(char)); // `+1` to add separator
            buf[UMI_len] = '\0';
            bam_aux_append(b, m_ptr, 'Z', UMI_len+1, (uint8_t*)buf);
        }
        
        int re = sam_write1(of, header, b);
        if (re < 0)
        {
            stringstream err_msg;
            err_msg << "fail to write the bam file: " << bam_get_qname(b) << "\n";
            err_msg << "return code: " << re << "\n";
            stop(err_msg.str());
        }
    }

    running = false;
    reporter_thread.join();

    // final report of processing speed
    Rcout
        << cnt << " reads processed" << ", "
        << cnt / timer.seconds_elapsed() / 1000 << "k reads/sec" << endl;

    Rcout << "number of read processed: " << cnt << "\n";
    
    sam_close(of);
    bgzf_close(fp);
    if (p.pool) hts_tpool_destroy(p.pool);
}
