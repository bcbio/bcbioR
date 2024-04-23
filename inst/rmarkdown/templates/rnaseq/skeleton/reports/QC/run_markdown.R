library(rmarkdown)

rmarkdown::render("./inst/rmarkdown/templates/rnaseq/skeleton/reports/QC/QC.Rmd",
                  output_dir = "./inst/rmarkdown/templates/rnaseq/skeleton/reports/QC/", 
                  clean = TRUE,
                  output_format = "html_document",
                  params = list(
                    params_file = '../../params_qc.R')
                  )
