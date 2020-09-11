// transcriptmapping.h
//

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <regex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "config_hts.h"
#include "utils.h"
#include "Timer.h"

#ifndef TRANSCRIPTMAPPING_H
#define TRANSCRIPTMAPPING_H


class Mapping
{
public:

    void parse_align_warpper(std::vector<std::string> fn_vec, std::string fn_out,  std::string cellular_tag, std::string molecular_tag, int nthreads);
    void parse_align(std::string fn, std::string fn_out, std::string cellular_tag, std::string molecular_tag, int nthreads);

    /* data */
};

#endif
