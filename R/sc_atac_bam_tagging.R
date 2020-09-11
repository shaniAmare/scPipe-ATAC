#' sc_atac_bam_tagging()
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

sc_atac_bam_tagging = function(inbam, outbam="",
                               bam_tags = list(bc="CB", mb="OX"),
                               nthreads=1) {
  
  if (any(!file.exists(inbam))) {
    stop("At least one input bam file should be present")
  } else {
    inbam = path.expand(inbam)
  }
  
  if(outbam == ""){
    fileNameWithoutExtension <- strsplit(inbam, "\\.")[[1]][1]
    outbam                   <- paste(fileNameWithoutExtension, "_tagged.bam", sep = "")
    outsortedbam             <- paste(fileNameWithoutExtension, "_sorted", sep = "")
  }
  
  outbam <- path.expand(outbam)
  if(!file.exists(outbam)){
    file.create(outbam)
  }
  
  # if(output_folder == ''){
  #   output_folder <- file.path(getwd(), "scPipe-atac-output")
  # }
  # 
  # if (!dir.exists(output_folder)){
  #   dir.create(output_folder,recursive=TRUE)
  #   cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  # }
  
  
  rcpp_sc_atac_bam_tagging(inbam, outbam, bam_tags$bc, bam_tags$mb,nthreads)
  
  cat("Output Tagged Bam File Location:")
  cat(outbam)
  cat("\n")

  
  sf <- sortBam(outbam,outsortedbam, indexDestination=TRUE)
  
  cat("Output Sorted & Indexed Bam File Location:")
  cat(outsortedbam)
  cat("\n")
  
  # generate the fragment file for the BAM file
  # need bedtools v2.26.0 or later
  #system2("bedtools", c("bamToBed", "i", outsortedbam), "|", "awk", c(-F"#" '{print $1"\t"$2}'), stdout = paste(output_folder,"/fragments.bed",sep = ""))
}
