library(bcbioR)

test_that("rnaseq copy",{
  path <- withr::local_tempdir()
  use_bcbio_projects(path, nfcore="nf-core/rnaseq", copy=TRUE, git=FALSE)
})

# test_that("rnaseq testing", {
#   path <- withr::local_tempdir()
#   print(path)
#   bcbio_templates(type="rnaseq", outpath=path)
#   numerator="tumor"
#   denominator="normal"
#   subset_value=NA
#   rmarkdown::render(input = file.path(path,"DE/DEG.Rmd"),
#                     output_dir = file.path(path,"DE"),
#                     output_format = "html_document",
#                     output_file = ifelse(!is.na(subset_value),
#                                          paste0('DE_', subset_value, '_', numerator, '_vs_', denominator, '.html'),
#                                          paste0('DE_', numerator, '_vs_', denominator, '.html')
#                     ),
#                     clean = TRUE,
#                     envir = new.env(),
#                     params = list(
#                       subset_value = subset_value,
#                       numerator = numerator,
#                       denominator = denominator,
#                       params_file = file.path(path,'DE/params_de.R'),
#                       project_file = file.path(path,'information.R'),
#                       functions_file = file.path(path,'DE/load_data.R')
#                     )
#   )
# })
