library(rmarkdown)

render_de <- function(subset_value, numerator, denominator){
  rmarkdown::render(input = "DEG.Rmd",
                    output_dir = "reports",
                    output_format = "html_document",
                    output_file = paste0('DE_', subset_value, '_', numerator, '_vs_', denominator, '.html'),
                    clean = TRUE,
                    envir = new.env(),
                    params = list(
                      subset_value = subset_value,
                      numerator = numerator,
                      denominator = denominator
                    )
  )
}

render_de('HDFn', "TNF", "untreated")
