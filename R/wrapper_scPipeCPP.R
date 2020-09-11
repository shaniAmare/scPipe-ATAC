library("tools")
library("Rsubread")
library("dplyr")
library("tidyr")
library("locStra")
library("Matrix")

set_paths = function(macs2path,bedToolsPath){
  
  old_path <- Sys.getenv("PATH")
  new_path = paste(old_path, macs2Path, sep = ":")
  new_path = paste(old_path, bedToolsPath, sep = ":")
  
  Sys.setenv(PATH = new_path)
  
}

sc_atac_trim_barcode = function(r1, r2, bc_file, output_folder = "", bc_start=-1, bc_length=-1, umi_start=0, umi_length=0, umi_in = "both") {
  
  if(output_folder == ''){
    output_folder = file.path(getwd(), "output")
  }
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  if (substr(r1, nchar(r1) - 2, nchar(r1)) == ".gz") {
    write_gz = TRUE
  }
  else {
    write_gz = FALSE
  }
  
  if (!is.null(bc_file)) {
    if (!file.exists(r1)) {stop("read1 fastq file does not exists.")}
    i=1;
    for (bc in bc_file) {
      if (!file.exists(bc)) {stop("Barcode file does not exists.")}
      bc_file[i] = path.expand(bc)
      i = i+1;
    }
    
    # expand tilde to home path for downstream gzopen() call
    r1 = path.expand(r1)
    
    if(!is.null(r2)){
      if (!file.exists(r2)) {stop("read2 file does not exists.")}
      r2 = path.expand(r2)
    }else{
      r2=""
    }
    cat("Saving the output at location: ")
    cat(output_folder)
    cat("\n")
    
    if(file_ext(bc_file) != 'csv'){
      rcpp_sc_trim_barcode_paired(output_folder, r1, bc_file,r2,write_gz)
    }
    else {
      cat("Using barcode CSV file, since Barcode FastQ file not passed \n ")
      if(bc_start == -1 || bc_length == -1 ){
        stop("Please pass bc_start and bc_length values")
      }
      rcpp_sc_trim_barcode(output_folder, r1, r2, bc_file, bc_start, bc_length, umi_start, umi_length, umi_in, write_gz)
    }
  }else{
    stop("Barcode file is mandatory")
  }
}

sc_atac_bam_tagging = function(inbam, outbam="",
                               bam_tags = list(bc="CB", mb="OX"),
                               nthreads=1) {
  
  if (any(!file.exists(inbam))) {
    stop("At least one input bam file does not exist")
  } else {
    inbam = path.expand(inbam)
  }
  
  if(outbam == ""){
    fileNameWithoutExtension = strsplit(inbam, "\\.")[[1]][1]
    outbam = paste(fileNameWithoutExtension, "_tagged.bam", sep = "")
    outsortedbam = paste(fileNameWithoutExtension, "_sorted", sep = "")
  }
  
  outbam = path.expand(outbam)
  if(!file.exists(outbam)){
    file.create(outbam)
  }
  
  cat("Output Tagged Bam File Location:")
  cat(outbam)
  cat("\n")
  
  rcpp_sc_exon_mapping(inbam, outbam, bam_tags$bc, bam_tags$mb,nthreads)
  
  cat("Output Sorted & Indexed Bam File Location:")
  cat(outsortedbam)
  cat("\n")
  
  sf <- sortBam(outbam,outsortedbam, indexDestination=TRUE)
}

sc_atac_call_peaks = function(inbam, output_folder = ""){
  
  if(output_folder == ''){
    output_folder = file.path(getwd(), "output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  system2("macs2", c("callpeak", "-t", inbam , "--outdir", output_folder)) 
  
}

sc_atac_feature_counting = function(inbam, feature_input, approach = "peak", bin_size = NULL, output_folder = ""){
  
  if(output_folder == ''){
    output_folder = file.path(getwd(), "output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  if(approach == 'genome_bin'){ 
    cat("Approach Selected Genome Bin", "\n")
    if(is.null(bin_size)){
      stop("Bin Size is mandatory for genome_bin approach")
    }
    system2("awk", c("'/^>/{if (l!=\"\") print l; print; l=0; next}{l+=length($0)}END{print l}'", feature_input, " | paste -d , - -") , stdout = paste(output_folder,"/genome_bin.csv",sep = ""))
    genome_length <- read.csv(file = paste(output_folder,"/genome_bin.csv",sep = ""), header=FALSE)
    x <- data.frame(stringsAsFactors = FALSE)
    
    for( i in rownames(genome_length) ){
      barcode = substr(genome_length[i, 1], 2, nchar(genome_length[i, 1]))
      length = genome_length[i, 2]
      j = 0
      while(j < as.integer(length)){
        if (j+bin_size < as.integer(length)){
          x <- rbind(x,list(barcode, as.integer(j), as.integer(j+bin_size)))
        }else{
          x <- rbind(x,list(barcode, as.integer(j), as.integer(length)))
        }
        j = j + bin_size
      }
    }
    
    write.table(x,paste(output_folder,"/genome_bin_feature.csv",sep = ""), row.names = FALSE, sep='\t', quote = FALSE, col.names = FALSE)
    feature_input = paste(output_folder,"/genome_bin_feature.csv",sep = "")
    cat("Generated the genome bin csv: ", feature_input , "\n")
  }
  
  #cat("Approach Selected Peak File")
  
  system2("bedtools", c("intersect", "-abam", inbam, "-b", feature_input, "-bed", "-wa", "-wb", "-sorted" ), stdout = paste(output_folder,"/intersect.bed",sep = ""))
  
  system2("awk", c("-F'#'","'$1=$1'","OFS='\t'", paste(output_folder,"/intersect.bed",sep = "")),stdout = paste(output_folder,"/intersectSplit.bed",sep = ""))
  
  system2("bedtools", c("groupby","-i", paste(output_folder,"/intersectSplit.bed",sep = ""),"-g","14,15,16,4","-c", "9","-o","count"), stdout= paste(output_folder,"/matrix.bed",sep = ""))
  
  matrixData <- as.data.frame(read.table(paste(output_folder,"/matrix.bed",sep = ""), sep="\t",header=TRUE))
  
  cat("Generated the Matrix Bed ", paste(output_folder,"/matrix.bed",sep = "") , "\n")
  
  names(matrixData) <- c("chromosome","start","end","barcode","count")
  
  matrixData <- matrixData %>%
    unite("chrS", chromosome:start, sep=":") %>%
    unite("Feature", chrS:end, sep="-")
  
  matrixData <- matrixData %>% 
    group_by(Feature,barcode) %>% 
    mutate(grouped_id = row_number())
  
  matrixData <- matrixData %>% 
    spread(barcode, count) %>% 
    select(-grouped_id)
  
  matrixData <- as.matrix(matrixData)
  
  saveRDS(matrixData, file = paste(output_folder,"/feature_matrix.rds",sep = ""))
  cat("Feature matrix generated: ", paste(output_folder,"/feature_matrix.rds",sep = "") , "\n")
  
  sparseM <- Matrix(matrixData,sparse=TRUE)
  cat("Sparse matrix generated: ", "\n")
  
  jaccardM <- jaccardMatrix(sparseM)
  
  saveRDS(jaccardM, file = paste(output_folder,"/jaccard_matrix.rds",sep = ""))
  cat("Jaccard matrix generated: ", paste(output_folder,"/jaccard_matrix.rds",sep = "") , "\n")
  
  matrixData[matrixData>0] = 1
  saveRDS(matrixData, file = paste(output_folder,"/binary_matrix.rds",sep = ""))
  cat("Binary matrix generated: ", paste(output_folder,"/binary_matrix.rds",sep = "") , "\n")
  
}

sc_atac_align = function (ref, readFile1, readFile2=NULL, output_folder = ""){
  
  if(output_folder == ''){
    output_folder = file.path(getwd(), "output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  if(!file.exists(readFile1)){
    stop("Input File readFile1 does not exist")    
  }
  
  if(!is.null(readFile2) && !file.exists(readFile2)){
    stop("Input File readFile2 does not exist")    
  }
  
  indexPath =  file.path(output_folder, "genome_index") 
  buildindex(basename=indexPath,reference=ref)
  
  fileNameWithoutExtension = strsplit(readFile1, "\\.")[[1]][1]
  outputPath = paste(fileNameWithoutExtension, "_aligned.bam", sep = "")
  align(index=indexPath,readfile1=readFile1,readfile2=readFile2, output_file=outputPath)
} 