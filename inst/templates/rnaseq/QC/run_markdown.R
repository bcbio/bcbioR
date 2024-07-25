library(rmarkdown)

# set directory to this file folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# example running with test data
rmarkdown::render("QC_nf-core.Rmd",
                  output_dir = ".",
                  clean = TRUE,
                  output_format = "html_document",
                  params = list(
                    params_file = 'params_qc_nf-core-testdata.R',
                    project_file = '../information.R')
                  )
