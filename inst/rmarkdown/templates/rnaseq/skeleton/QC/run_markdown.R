library(rmarkdown)

rmarkdown::render("./inst/rmarkdown/templates/rnaseq/skeleton/QC/QC_nf-core.Rmd",
                  output_dir = "./inst/rmarkdown/templates/rnaseq/skeleton/QC/", 
                  clean = TRUE,
                  output_format = "html_document",
                  params = list(
                    params_file = 'params_qc_nf-core.R',
                    project_file = '../information.R')
                  )
