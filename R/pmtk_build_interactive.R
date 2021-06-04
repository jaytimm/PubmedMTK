#' Build interactive html for topic summary
#'
#' @name pmtk_build_interactive
#' 
#' @param pmtk_lda Object returned from pmtk_summarize_lds 
#' @param pmtk_2d Object returned from pmtk_2d 
#' @param out_dir Directory to output html
#' @param file_name Name to assign html file
#' 
#' @return An html file
#' 
#' @export
#' @rdname pmtk_build_interactive
#' 
#' 
pmtk_build_interactive <- function(pmtk_lda,
                                   pmtk_2d,
                                   out_dir,
                                   file_name){
  
  fx <- system.file("rmd", "pmtk_topic_viz.Rmd", package = "PubmedMTK")
  
  rmarkdown::render(
    input  = fx,
    
    params = list(topic_model = pmtk_lda,
                  tsne = pmtk_2d),
    
    output_file = paste0(out_dir, file_name)
  )
  
  #htmltools::includeHTML('/home/jtimm/Desktop//my.html') ## this is cool -- !! -- 
  
  rstudioapi::viewer(url = paste0(out_dir, file_name))
  
}