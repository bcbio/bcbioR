library(rmarkdown)

rmarkdown::render("QC.Rmd",
                  output_dir = "./", 
                  clean = TRUE,
                  output_format = "html_document",
                  params = list(
                    params_file = 'params_qc.R',
                    project_file = '../information.R')
                  )
