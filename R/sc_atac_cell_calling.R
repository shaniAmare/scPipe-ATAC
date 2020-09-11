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

sc_atac_cell_calling = function (feature_matrix, call_method){
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  mat <- readMM(feature_matrix)
  
  if(call_method == 'emptyDrops'){ 
  
      set.seed(2019)
      cell.out     <- emptyDrops(mat)
      
      filter.out   <- cell.out[complete.cases(cell.out), ]
      
      saveRDS(filter.out, file = paste0(output_folder, '/EmptyDrop_obj.rds'))
      
      filter.out   <- filter.out[filter.out$FDR <= fdr, ]
      
      select.cells <- rownames(filter.out)
      
      out_mat      <- mat[, colnames(mat) %in% select.cells]
      barcodes     <- colnames(out_mat)
      features     <- rownames(out_mat)
      
      writeMM(out_mat, file = paste0(output_folder, '/matrix.mtx'))  
      write.table(barcodes, file = paste0(output_folder, '/non_empty_barcodes.txt'), sep = '\t', 
                row.names = FALSE, quote = FALSE, col.names = FALSE)
      write.table(features, file = paste0(output_folder, '/non_empty_features.txt'), sep = '\t',
              row.names = FALSE, quote = FALSE, col.names = FALSE)
      }
  
  if(call_method == 'callCells'){ 
    
  
    }
  
  
  }

