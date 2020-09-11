#' sc_atac_align()
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

sc_atac_aligning <- function (ref, readFile1, readFile2=NULL, readDir = NULL, output_folder = ""){
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
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
  
  indexPath <-  file.path(output_folder, "genome_index") 
  buildindex(basename=indexPath,reference=ref)
  
  fileNameWithoutExtension <- strsplit(readFile1, "\\.")[[1]][1]
  outbam                   <- paste(fileNameWithoutExtension, "_aligned.bam", sep = "")
  
  #execute Rsubread align()
  Rsubread::align(index=indexPath,readfile1=readFile1,readfile2=readFile2, output_file=outbam)
  
  Rsamtools::indexBam(outbam)
  
  # get the unmapped mapped stats to be output and stored in a log file
  #can use Rsamtools::idxstatsBam()
  
} 