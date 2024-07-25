library(rmarkdown)
# set working directory to this file before using the function


# set directory to this file folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# example running with test data
render_de <- function(column, numerator, denominator, subset_value = NULL,
                      params_file = 'params_de-testdata.R'){
                        
  rmarkdown::render(input = "DEG.Rmd",
                    output_dir = ".",
                    output_format = "html_document",
                    output_file = ifelse(!is.null(subset_value),
                                         paste0('DE_', subset_value, '_', numerator, '_vs_', denominator, '.html'),
                                         paste0('DE_', numerator, '_vs_', denominator, '.html')
                                         ),
                    clean = TRUE,
                    envir = new.env(),
                    params = list(
                      column = column,
                      subset_value = subset_value,
                      numerator = numerator,
                      denominator = denominator,
                      params_file = params_file,
                      project_file = '../information.R',
                      functions_file = 'load_data.R'
                    )
  )
}
#Example data
render_de("sample_type","tumor", "normal")
