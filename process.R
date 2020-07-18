#'
#' Objective : generate a fully-fledged scPipe module for scATAC-Seq data preprocessing
#' Author: S. Aarasinghe
#' Date: 30/06/2020
#' 
  
# Pre-requisites -------------
# install.packages(c("devtools", "roxygen2", "testthat", "knitr", "tidyr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("RSubread");
# BiocManager::install("Rhtslib")

# TO DO: these dependencies should be handled within the package
# TO DO: a log file needs to record what is happening + basic stats
# TO DO: proper error handling

devtools::load_all()

# Demultiplexing -----------------

R1_path      <- system.file(
  "data","testfastq_S1_L001_R1_001.fastq.gz",package="seqATAC") 

R2_path      <- system.file(
  "data","testfastq_S1_L001_R2_001.fastq.gz",package="seqATAC") 

R3_path      <- system.file(
  "data","testfastq_S1_L001_R3_001.fastq.gz",package="seqATAC") 

Barcode_path <- system.file("data","barcode.csv",package="seqATAC")

dir.create("result") 
# TO DO: if !existing create file
# ask user (optionally) to enter the preferred output location, else mention that results are located in <default location>.
output_path  <- system.file("result",package="seqATAC") 

# execute the modified function; TO DO: change function name
# to run with barcode fastq files
sc_trim_barcode(output_path, R1_path, c(R2_path), R3_path)

# to run with barcode.csv file (barcode in column 2)
# TO DO: only extract the length from barcode too (sometimes the barcode contains some more characters e.g. "-1" after the actual barcode sequence)
sc_trim_barcode(outfq = output_path, r1 = R1_path, r3 = R3_path, barcode_file = Barcode_path, start = 3,length = 16)
# currently r terminates if R2 is not added.
# TO DO: order of input has to change + the third file has to be checked to see whether it is a fastq file or a .csv file, 
# else add a commands that takes whether the barcode/s is/are in fastq or csv format
# TO DO: if .csv file needs to check for the next commands


# TO DO: warning to say the files are existing and options to remove them before re-running

# TO DO: replace the R1 and R3 with read1 and read2 for it to be more intuitive
# TO DO: how is using multiple barcode.csv files are handled?
# TO DO: add the functionality to take the UMI in
# TO DO: add the functionality use the hamming distance for matching the barcode if not coming form a fastq file
# The above hamming distance is expected to generate two more demultiplexed fastq files (e.g. adjusted_match.fastq.gz, ambiguous.fastq.gz)

# TO DO: load the libraries that are required

# TO DO: remove the "non-cells" based on EmptyDrops (https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/EmptyDrop.R)
# or filter barodes (https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/filter_barcodes.R) from scATAC-Pro

# Aligning -------------------

library(Rsubread) # TO DO: this needs to be loaded within package, write a wrapper

# TO DO: Need to take the genome.fa in the wrapper
ref    <- system.file("data","genome.fa",package="seqATAC") 

# build aligning index (wrapper within package)

# TO DO: need to save in a different directory not in the directory where the demultiplexed data is stored, 
# a prefix name (optional) should be taken in the wrapper function, or create one
buildindex(basename="result/genome_index",reference=ref)

reads  <- system.file("result","demux_testfastq_S1_L001_R1_001.fastq.gz",package="seqATAC") 

# run Rsubread
# TO DO: need to write a wrapper within package, need to define a different location for the BAM files
# minor point: BAM file extensions are .bam and not.BAM
# remove Rsubread current version; thi is buggy. Ask user to align outside the package if they want to use Rsubread
# Include the bwa-mem2 (https://github.com/bwa-mem2) and star (https://github.com/alexdobin/STAR) aligners instead witht eh ability for user to pick and define paramters

align(index="result/genome_index",readfile1=reads,output_file="result/testfastq_S1_L001_R1_001.BAM",phredOffset=64)


# TO DO: Need to include several QC steps here from scPIpe (MT, convert to chr, MAPQ level, etc.)

# BAM Tagging ----------------

BAM_file        <- system.file(
  "result","testfastq_S1_L001_R1_001.BAM",package="seqATAC") 

file.create("result/testfastq_S1_L001_R1_001_tagged.BAM")
# TO DO: need to create the BAM file if not existing, currently R terminates if the file does not exist
# AND warn if the file is existing and give options to remove it before re-running

BAM_output_file <- system.file(
  "result","testfastq_S1_L001_R1_001_tagged.BAM",package="seqATAC") 

# modify the sc_exon_mapping to only tag the BAM file; TO DO: change function name
sc_exon_mapping(BAM_file, BAM_output_file)
# TO DO: BAM file with windows endings need to be fixed; revert back to previous subread package and see if the problem persists
# TO DO: index the BAM file when created (samtools needs to be a dependency)

# PEAK FILE Generation using MACS2 ------------------

system2("macs2", c("callpeak", "-t", BAM_output_file, "--outdir", "result")) # MACS2 is a dependency, 
# TO DO: wrapper for this function, take the macs2 path in

# Also need to ask for several parameters that can be adjusted by the user before running MACS2: 
# --nomodel, -f, --shift,  --extsize , --gsize, --outdir, --name, --bdg , --SPMR, --cutoff-analysis, --call-summits

# TO DO: remove peaks overlapped with a blacklisted region
# user can produce their own blacklisted region OR have an inbuild one | in https://github.com/tanlabcode/scATAC-pro/tree/master/annotation
# can use bedtools intersected for this

# TO DO: use csaw for sliding window method (upgrade)


# Feature counting ----------------------------------

BAM_output_file         <- system.file("result","testfastq_S1_L001_R1_001_tagged.BAM",package="seqATAC")
BAM_output_sorted_file  <- system.file("result","testfastq_S1_L001_R1_001_tagged_sorted.BAM",package="seqATAC") 
peak_file               <- system.file("result","NA_peaks.narrowPeak",package="seqATAC")
intersect_file          <- system.file("result","intersect.bed",package="seqATAC")

# Sort the bam file 
system2("/Users/300002291/bin/samtools",c("sort",BAM_output_file,"-o",BAM_output_sorted_file))
# samtools sort "result/testfastq_S1_L001_R1_001_tagged.BAM" -o "result/testfastq_S1_L001_R1_001_tagged_sorted.BAM"

#TO DO: Use bedr package and rsamtools package

# Find the intersection between Sorted bam File and the peak File 
system2("bedtools", c("intersect", "-abam", BAM_output_sorted_file, "-b", peak_file, "-bed", "-wa", "-wb", "-sorted" ), stdout = "result/intersect.bed")
# intersectBed -abam "result/testfastq_S1_L001_R1_001_tagged_sorted.BAM" -b "result/NA_peaks.narrowPeak" -bed -wa -wb -sorted > "result/intersect.bed"

# Seprate the barcode and tag using the awk command.
system2("awk", c("-F'#'","'$1=$1'","OFS='\t'","result/intersect.bed"),stdout = "result/intersectSplit.bed")
# awk -F'#' '$1=$1' OFS='\t' "result/intersect.bed" > "result/intersectSplit.bed"

# GroupBy the Bed File 
# TO DO: IMP. NOTE; below approach is not working!
# idea from scATAC-pro: https://github.com/wbaopaul/scATAC-pro/blob/c85326061109e913a25a3f7ca3a9c62c2d4f8a75/scripts/src/get_mtx.R
system2("bedtools", c("groupby","-i","result/intersectSplit.bed","-g","14,15,16,4","-c", "9","-o","count"), stdout= "result/matrix.bed")
# groupBy -i "result/intersectSplit.bed" -g 14,15,16,4 -c 9 -o count > "result/matrix.bed"

#TO DO: Need to load dplyr and tidyr here

# TO DO:Matrix.bed need headers

matrixData <- as.data.frame(read.table("result/matrix.bed", sep="\t",header=FALSE))
names(matrixData) <- c("chromosome", "start", "end", "barcode", "count")

# TO DO: "Feature" is a better word for rownames
matrixData <- matrixData %>%
  unite("chrS", chromosome:start, sep=":") %>%
  unite("feature", chrS:end, sep="-")

# TO DO: 4th column of matrixData is not seq, but barcode
matrixData <- matrixData %>% 
  group_by(feature,barcode) %>% 
  mutate(grouped_id = row_number())

matrixData <- matrixData %>% 
  spread(barcode, count) %>% 
  select(-grouped_id)

# check the value in matrix  
# matrixData[1,'GAACTTGGTTCAGAAA']

# TO DO: load sce library
# TO DO: row.names of the matrix should be features and column names should be barcodes - add features as row.names and remove features from matrixData

# When I tried following; I get the error below, please see what do to do about it:
row.names(matrixData) <- matrixData$feature
#Error in `.rowNamesDF<-`(x, value = value) : 
#  duplicate 'row.names' are not allowed
# In addition: Warning messages:
#  1: Setting row names on a tibble is deprecated. 
#  2: non-unique values when setting 'row.names': ‘chr21:10085049-10085962’, ‘chr21:10086090-1008640

# TO DO: Get all the info properly into the sce (use scPipe last step : https://github.com/LuyiTian/scPipe/blob/master/R/sc_workflow.R)
# Use the ideas in this step tp generate the relevent reports too

# TO DO: the script of generatign the sce needs to change to something like following:
sce <- SingleCellExperiment(assays = list(counts=as.matrix(matrixData)))

# Ignore Commands 
# Sample_BAM_file  <- system.file("result","Mapped_testfastq_S1_L001_R1_001.BAM",package="seqATAC")
# peak_file        <- system.file("result","NA_peaks.narrowPeak",package="seqATAC")
# sc_demultiplex(Sample_BAM_file,output_path, peak_file)
