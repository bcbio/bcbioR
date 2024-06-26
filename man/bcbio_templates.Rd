% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{bcbio_templates}
\alias{bcbio_templates}
\title{Function to help deploy analysis folder inside a project folder}
\usage{
bcbio_templates(type = "rnaseq", outpath)
}
\arguments{
\item{type}{string indicating the type of analysis, supported:
\itemize{
\item base
\item rnaseq, scrnaseq,
\item teaseq
\item cosmx
}}

\item{outpath}{string path indicating where to copy all the files to}
}
\description{
This function contains Rmd, R, md, files that help to structure
an analysis following HCBC best-practices.
For rnaseq, it will deploy: QC and DE Rmd with additional files to help
to facilitate the analysis as needed.
}
\details{
Normally these helper files are inside a report folder inside a
project folder.
}
\examples{
 \dontrun{
 bcbio_templates("rnaseq", "path_to_projects/project1/reports")
 }
}
