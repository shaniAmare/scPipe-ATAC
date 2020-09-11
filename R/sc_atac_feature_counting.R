#' sc_atac_feature_counting()
#'
#' @return 
#'
#' @examples
#' \dontrun{
#' 
#' 
#' }
#'
#' @export
#'

library(rtracklayer)
library(GenomicAlignments)
library(Rsamtools)
library(SummarizedExperiment) 
library(DelayedArray)
library(BiocParallel)
library(matrixStats)
library(Biobase)
library(GenomicRanges)
library(GenomeInfoDb)
library(IRanges)            
library(S4Vectors)
library(BiocGenerics)  

sc_atac_feature_counting <- function(insortedbam, 
                                    feature_input, 
                                    bam_tags = list(bc="CB", mb="OX"), 
                                    feature_type = "peak", 
                                    bin_size = NULL, 
                                    yieldsize = 1000000,
                                    mapq = 0,
                                    blacklist = NULL,
                                    output_folder = ""){
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  if(feature_type == 'genome_bin'){ 
    # TODO: test the format of the fasta file here and stop if not proper format
    cat("`genome bin` feature type is selected for feature input. reading the genome fasta file ...", "\n")
    
    # TODO: execute following if the index for fasta is not found...
    Rsamtools::indexFa(feature_input)
    cat("Indes for ", feature_input, " not found. Creating one now... ", "\n")
    
    
    if(is.null(bin_size)){
      bin_size <- 2000
      cat("Default bin size of 2000 is selected", "\n")
    }
    
    # Here I want to execute the C++ code to generate the genome_bin file based on bin_size (above) 
    # that would be later used in this same function to 
    # overlap with another data file type after the genome_bin file is created
    # and fed in as the feature_inpui to create feature.gr below.
    
    #####________________________
    # I attempted to create this genome_bin in R, but below crashes my computer! 
    # and that's why I thought to look for an alternative for this step
    # legacyseq <- seqinr::read.fasta(file = feature_input, as.string = TRUE)
    # l         <- nchar(legacyseq)
    
    # why do I save this temp file here?
    # write.table(l, paste(output_folder,"/genome_bin.csv",sep = ""),row.names=TRUE, col.names = FALSE, sep=",", quote = FALSE)
    
    #system2("awk", c("'/^>/{if (l!=\"\") print l; print; l=0; next}{l+=length($0)}END{print l}'", feature_input, " | paste -d , - -") , 
    #        stdout = paste(output_folder,"/genome_bin.csv",sep = ""))
    

    # genome_length   <- read.csv(file = paste(output_folder,"/genome_bin.csv",sep = ""), header=FALSE)
    # x               <- data.frame(stringsAsFactors = FALSE)
    
    # for( i in rownames(genome_length) ){
    #   barcode <- genome_length[i, 1]
    #   length  <- genome_length[i, 2]
    #   j       <- 0
    #   while(j < as.integer(length)){
    #     if (j+bin_size < as.integer(length)){
    #       x   <- rbind(x,list(barcode, as.integer(j), as.integer(j+bin_size)))
    #     }else{
    #       x   <- rbind(x,list(barcode, as.integer(j), as.integer(length)))
    #     }
    #     j    <- j + bin_size
    #   }
    # }
    
    # write.table(x,paste(output_folder,"/genome_bin_features.csv",sep = ""), row.names = FALSE, sep='\t', quote = FALSE, col.names = FALSE)
    # feature_input <- paste(output_folder,"/genome_bin_feature.csv",sep = "")
    # cat("Generated the genome bin csv: ", feature_input , "\n")
    ####_____________________
    
    # finally what is needed is....
    # feature_input <- genome_bin
  }
  
  
  # generate the GAlignments object for further use
  
  if(feature_type == 'peak'){ 
    cat("`peak` feature_type is selected for feature input", "\n")
  } 
    
  cat("Creating GAlignment object for the sorted BAM file...")
  #tags = "CB"
  param <- Rsamtools::ScanBamParam(tag = as.character(bam_tags),  mapqFilter=mapq)
  bamfl <- Rsamtools::BamFile(insortedbam, yieldSize = yieldsize)
  
  
  open(bamfl)
  
  yld    <- GenomicAlignments::readGAlignments(bamfl,use.names = TRUE, param = param)
  yld.gr <- makeGRangesFromDataFrame(yld,keep.extra.columns=TRUE) 
  #yld.gr <- as(GRanges(yld), "GAlignments")
  
  # ______________ get the average lines per CB (average reads per cell barcode) here for the GAlignment object
  
  saveRDS(yld.gr, file = paste(output_folder,"/BAM_GAlignmentsObject.rds",sep = ""))
  cat("Galignment object is created and saved: ", paste(output_folder,"/BAM_GAlignmentsObject.rds",sep = "") , "\n")
  
  # generate the GAalignment file from the feature_input file
  cat("Creating Galignment object for the feature input...")
 
  feature.gr <- rtracklayer::import(feature_input)
  
  # ______________ need to add the blacklist filtering here , relevent param has been added to the function i.e. blacklist.
  # this param should take either a tab delimited file in the format of feature_input or keywords "hs", "mm" which would load a stored GAlignment file related to them
  
  overlaps <- findOverlaps( query = yld.gr, subject = feature.gr, type = "within", minoverlap = 0, ignore.strand = TRUE)
  
  # generate the matrix using this overlap results above.
  
  #yld.gr[-queryHits(findOverlaps(yld.gr, feature.gr, type="any", ignore.strand = TRUE)),] 
  #feature.gr[-queryHits(findOverlaps(feature.gr, yld.gr, type="any", ignore.strand = TRUE)),] 
  # mergeByOverlaps(feature.gr, yld.gr)
  
  mcols(yld.gr)[queryHits(overlaps), "peakStart"] <- start(ranges(feature.gr)[subjectHits(overlaps)])
  mcols(yld.gr)[queryHits(overlaps), "peakEnd"]   <- end(ranges(feature.gr)[subjectHits(overlaps)])
  
  
  overlap.df <- data.frame(yld.gr) %>% select(seqnames, peakStart, peakEnd, CB)
  matrixData <- overlap.df %>% group_by(seqnames, peakStart, peakEnd, CB) %>% summarise(count = n())
  
  # generate the matrix format
  names(matrixData) <- c("chromosome","start","end","barcode","count")
  
  matrixData <- matrixData %>%
    unite("chrS", chromosome:start, sep=":") %>%
    unite("feature", chrS:end, sep="-")
  
  matrixData <- matrixData %>% 
    group_by(feature,barcode) %>% 
    mutate(grouped_id = row_number())
  
  matrixData <- matrixData %>% 
    spread(barcode, count) %>% 
    select(-grouped_id)
  
  # TODO: sanity check to see if all features are unique here
  #if(!unique(matrixData$feature)) {stop("There are duplicate values in the feature input. Please check and rerun this step again")}
  
  matrixData <- matrixData %>% as_tibble(rownames = "feature")
  matrixData <- matrixData[, -1]
  
  #TODO: Matrix format is not yet implemeted properly. Need to fix this
  #matrixData <- as.matrix(matrixData)
  #matrixData <- matrix(matrixData, dimnames = list(matrixData$feature, colnames(matrixData)))
  
  # call sc_atac_cell_callling.R here ... still ongoing
  
  saveRDS(matrixData, file = paste(output_folder,"/feature_matrix.rds",sep = ""))
  cat("Feature matrix generated: ", paste(output_folder,"/feature_matrix.rds",sep = "") , "\n")
  
  sparseM <- Matrix(matrixData,sparse=TRUE)
  cat("Sparse matrix generated: ", "\n")
  
  jaccardM <- jaccardMatrix(sparseM)
  cat("Jaccard matrix generated: ", "\n")
  
  saveRDS(jaccardM, file = paste(output_folder,"/jaccard_matrix.rds",sep = ""))
  cat("Jaccard matrix generated: ", paste(output_folder,"/jaccard_matrix.rds",sep = "") , "\n")
  
  matrixData[matrixData>0] = 1
  saveRDS(matrixData, file = paste(output_folder,"/binary_matrix.rds",sep = ""))
  cat("Binary matrix generated: ", paste(output_folder,"/binary_matrix.rds",sep = "") , "\n")
  
}