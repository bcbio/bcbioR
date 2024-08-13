library(bcbioR)


test_that("rnaseq deg",{
  path <- withr::local_tempdir()
  print(path)
  bcbio_templates(type="rnaseq", outpath=path)
  fs::dir_ls(path,all=T)
  numerator="tumor"
  denominator="normal"
  subset_value=NA
  rmarkdown::render(input = file.path(path,"DE/DEG.Rmd"),
                    output_dir = file.path(path,"DE"),
                    output_format = "html_document",
                    output_file = ifelse(!is.na(subset_value),
                                         paste0('DE_', subset_value, '_', numerator, '_vs_', denominator, '.html'),
                                         paste0('DE_', numerator, '_vs_', denominator, '.html')
                    ),
                    clean = TRUE,
                    # envir = new.env(),
                    params = list(
                      subset_value = subset_value,
                      numerator = numerator,
                      denominator = denominator,
                      params_file = file.path(path,'DE/params_de-example.R'),
                      project_file = file.path(path,'information.R'),
                      functions_file = file.path(path,'libs/load_data.R')
                    )
  )
  # browseURL(file.path(path, "DE/DE_tumor_vs_normal.html"))
  # usethis::proj_activate(path)
})

test_that("rnaseq qc",{
  path <- withr::local_tempdir()
  print(path)
  bcbio_templates(type="rnaseq", outpath=path)
  fs::dir_ls(path,all=T)
  rmarkdown::render(input = file.path(path,"QC/QC_nf-core.Rmd"),
                    output_dir = file.path(path,"QC"),
                    output_format = "html_document",
                    clean = TRUE,
                    params = list(
                      params_file = file.path(path,'QC/params_qc_nf-core-example.R'),
                      project_file = file.path(path,'information.R'),
                      functions_file = file.path(path,'libs/load_data.R')
                    )
  )
  # browseURL(file.path(path, "QC/QC_nf-core.html"))
  # usethis::proj_activate(path)
})
