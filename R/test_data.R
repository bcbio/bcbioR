#' Function to help download test data for QC chipseq report
#'
#' It downloads files from our testdata repository: [bcbio/bcbioR-test-data](https://github.com/bcbio/bcbioR-test-data/tree/main/chipseq)
#'
#' It downloads the `bowtie2/mergedLibrary/macs2/narrowPeak` output
#'
#' @export
bcbio_qc_chipseq_testdata <- function(){
  # if using example data to render report, download peaks from github
  api_url <- "https://api.github.com/repos/bcbio/bcbioR-test-data/contents/chipseq/bowtie2/mergedLibrary/macs2/narrowPeak"
  response <- GET(api_url)

  if (status_code(response) == 200) {
    content <- content(response, as = "text")
    files_info <- fromJSON(content) %>% filter(name != 'consensus')

    # Filter out file paths and construct raw URLs
    file_paths <- files_info$path
    raw_base_url <- "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/"

    raw_file_urls <- paste0(raw_base_url, file_paths)

    # Function to download a file from a URL
    download_file <- function(url) {
      file_name <- basename(url)
      download.file(url, destfile = file_name, mode = "wb")
    }

    # Download all files using the constructed raw URLs
    for (url in raw_file_urls) {
      download_file(url)
    }
    peaks_dir = '.'
  }
}
