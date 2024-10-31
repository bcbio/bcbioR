library(bcbioR)


test_that("rnaseq deg",{
  path <- withr::local_tempdir()
  print(path)
  bcbio_templates(type="rnaseq", outpath=path)
  fs::dir_ls(path,all=T)
  rmarkdown::render(input = file.path(path,"DE/DEG.Rmd"),
                    output_dir = file.path(path,"DE"),
                    output_format = "html_document",
                    clean = TRUE
                    # envir = new.env(),
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
