library(rmarkdown)

render_de <- function(numerator, denominator, subset_value = NA,
                      params_file = 'params_de.R'){

  rmarkdown::render(input = "DEG.Rmd",
                    output_dir = "./",
                    output_format = "html_document",
                    output_file = ifelse(!is.na(subset_value),
                                         paste0('DE_', subset_value, '_', numerator, '_vs_', denominator, '.html'),
                                         paste0('DE_', numerator, '_vs_', denominator, '.html')
                                         ),
                    clean = TRUE,
                    envir = new.env(),
                    params = list(
                      subset_value = subset_value,
                      numerator = numerator,
                      denominator = denominator,
                      params_file = params_file,
                      project_file = '../information.R'
                    )
  )
}

render_de("tumor", "normal")
