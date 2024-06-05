library(rmarkdown)

rmarkdown::render("reports/QC.Rmd",
                  output_dir = "reports", 
                  clean = TRUE,
                  output_format = "html_document",
                  output_file = "QC.html",
                  params = list(
                    seurat_fn = '../data/CTCL_PAT8_DEV3_LVL7.RDS',
                    project_file = './information.R',
                    umap_dim = 'approximateumap_8c6f278e.b9f4.4535.aeca.8955c1dff614_1'
                    )
)