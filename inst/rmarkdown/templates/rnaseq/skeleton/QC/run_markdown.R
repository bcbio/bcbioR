library(rmarkdown)

rmarkdown::render("QC_nf-core.Rmd",
                  output_dir = ".",
                  clean = TRUE,
                  output_format = "html_document",
                  params = list(
                    params_file = 'params_qc_nf-core.R',
                    project_file = '../information.R')
                  )
