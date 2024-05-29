.fix <- function(x){
  x <- tolower(x) %>% str_replace_all(., "[[:punct:]]", "_")
  x <- str_replace_all(x, " ", "_")
  return(x)
}


#' Function to check samplesheet for nf-core
#'
#' @param file path to CSV file for nf-core
#' @examples
#'
#' bcbio_nfcore_check(system.file("extdata", "rnaseq_good.csv", package = "bcbioR") )
#'
#' @export
bcbio_nfcore_check <- function(file){
  required=c("sample","fastq_1","fastq_2","strandedness")
  samplesheet=read_csv(file)

  if (!(all(required %in% colnames(samplesheet)))){
    print(colnames(samplesheet))
    stop("Missing required columns ", paste(required, collapse = " "))
  }else if (any(grepl("^[1-9]", samplesheet[["sample"]]))){
    stop("Avoid samples starting with numbers ")
  }else if (any(is.na(samplesheet))){
    warning("Columns with missing values")
  }else{
    message("All good.")
  }
}

#' Function to help deploy analysis folder inside a project folder
#'
#' This function contains Rmd, R, md, files that help to structure
#' an analysis following HCBC best-practices.
#' For rnaseq, it will deploy: QC and DE Rmd with additional files to help
#'   to facilitate the analysis as needed.
#'
#' Normally these helper files are inside a report folder inside a
#' project folder.
#'
#' @param type string indicating the type of analysis, supported: rnaseq.
#'
#' @param outpath string path indicating where to copy all the files to
#' @examples
#'  \dontrun{
#'  bcbio_templates("rnaseq", "path_to_projects/project1/reports")
#'  }
#' @export
bcbio_templates <- function(type="rnaseq", outpath){
  switch(type,
         rnaseq={

           fpath <- system.file("rmarkdown/templates/rnaseq", "skeleton", package="bcbioR")
           #file.copy(fpath, outpath, recursive = T)
           copyDirectory(fpath, outpath)
         },
         scrnaseq={

           fpath <- system.file("rmarkdown/templates/singlecell", "skeleton", package="bcbioR")
           #file.copy(fpath, outpath, recursive = T)
           copyDirectory(fpath, outpath)
         },
         {
           stop('project type not recognize, please choose: ', 'rnaseq', 'scrnaseq')
         }
  )
}

#' Function to help with project name used for parent folder
#'
#' This function will ask for user input about:
#'   * numeric code
#'   * PI full name
#'   * technology
#'   * tissue
#'   * organism
#'   * project description
#'
#' It removes special character with `_`. The output is a guideline to
#'   what the folder used can be.
#'
#' @returns A string list with hbc_code, and project folder name
#' @export
bcbio_set_project <- function() {
  hbc_code <- readline("What is the hbc code (only numbers):\n")
  hbc_code <- paste0("hbc", hbc_code)
  pi <- readline("What is PI last name:\n")
  technology <- readline("What is the technology:\n")
  tissue <- readline("What is the tissue:\n")
  org <- readline("What is the organism:\n")
  project <- readline("What is the project name:\n")
  #dropbox <- readline("What is the dropbox name:\n")
  #github_org <- readline("What is the github organization:\n")
  #hbc_$technology_of_$pilastname_$intervention_on_$tissue_in_$organism_$hbccode
  project_full <- paste(technology, .fix(pi), .fix(project), tissue, org, hbc_code, sep="_")
  #github <- c(github_org,project_full)
  opts <- list(code=hbc_code, project=project_full)
               #dropbox=file.path(dropbox,project_full),
               #github=github)
  print(opts)
  return(opts)
}


guess_analysis <- function(path){
  if (!fs::dir_exists(path))
    ui_abort("{ui_val(path)} doesn't exist")

  # This file is inside star_salmon/ folder
  counts_fn <- fs::path_join(path, '/star_salmon/salmon.merged.gene_counts.tsv')
  # This folder called "multiqc_report_data" is inside the output directory star_salmon inside multiqc folder
  multiqc_data_dir <- fs::path_join(path, 'star_salmon/multiqc_report_data')
  # This file is inside star_salmon/ folder
  se_object <- fs::path_join(path, 'star_salmon/salmon.merged.gene_counts.rds')

}

read_pipeline_info <- function(path){
  # pipeline_info/params_2024-05-28_12-28-51.json
  config <- fs::path_join(nfcore, "pipeline_info")
  params <- fs::dir_ls(config, regexp = "params")
  metadata <- jsonlite::read_json(params)[["input"]]
  # input
  # tmp_rna/pipeline_info/software_versions.yml
  software <- fs::path_join(nfcore, "pipeline_info", "software_versions.yml")
  software_txt <- yaml::read_yaml(software)
  pipeline <- grep("nf-core", names(software_text$Workflow), value = TRUE)
  # Workflow:
  #   Nextflow: 24.04.1
  # nf-core/rnaseq: 3.14.0
  # check only rnaseq is supported
  if (!(pipeline %in% c("nf-core/rnasew"))){
    iu_abort("Sorry, we don't yet support {.ui_value(pipeline)}")
  }
  list(metadata=metadata, pipeline=pipeline)
}


bcbio_params <-function(path, pipeline, metadata, copy){

  if (pipeline=="nf-core/rnaseq"){
    if (!copy){
      se_object <- fs::path_join(path, "star_salmon/salmon.merged.gene_counts.rds")
      metadata_fn <- metadata
      counts_fn <- fs::path_join(path, "star_salmon/salmon.merged.gene_counts.tsv")
      multiqc_data_dir <- fs::path_join(path, "multiqc/star_salmon/multiqc-report-data/")
      gtf_fn <- fs::path_join(path, "genome/genome.filtered.gtf")
    }

    analysis_template <- fs::path_package("bcbioR", "templates", "rnaseq", "qc")
    fs::dir_copy(analysis_template, fs::path_join(path, "reports"), overwrite = FALSE)
    analysis_template <- fs::path_package("bcbioR", "templates", "rnaseq", "de")
    fs::dir_copy(analysis_template, fs::path_join(path, "reports"), overwrite = FALSE)

    ui_info("Please, to start the analysis, modify these parameter in QC/QC.rmd")
    ui_todo("set genome to hg38, mm10, mm39, or other")
    ui_todo("set factor_of_interest to a column in your metadata")
  }

}

#' @export
use_bcbio_analysis <- function(path, nfcore=NULL, copy=FALSE, metadata=NULL){

  if (copy){
    # deploy files
    ui_info("Rmd templates will be copied but variables path won't be filled automatically.")
  }else{
    if (!fs::dir_exists(nfcore))
      ui_abort("{ui_value(nfcore)} doesn't exist. point to nfcore path or turn on copy mode.")

    #guess analysis from pipeline file
    information <- read_pipeline_info(nfcore)
    fs::dir_create(fs::path_join(path, "meta"))
    meta_path <- fs::path_join(path, "meta", fs::path_file(information$metadata))
    pipeline <- information$pipeline
    if (!is.null(metadata)){
      if (!(fs::file_exists(metadata)))
          ui_abort("{ui_value(metadata)} doesn't exist.")
      fs::file_copy(metadata, meta_path)
    }else{
      if (!fs::file_exists(information$metadata)){
        ui_warn("{ui_value(metadata)} not found. We can only work with local filesytems. For now.")
        ui_todo("Please, copy {ui_value(metadata)} to {ui_value(meta_path)}.")
        ui_warn("If this file is not in the folder, the code will fail.")
      }else{
        fs::file_copy(information$metadata, meta_path)
      }
      metadata <- meta_path
    }

    # set all files from analysis
    bcbio_params <- set_bcbio_params(nfcore, pipeline, metadata, copy=copy)
  }

}

#' @export
#' @examples
#' path <- withr::local_tempdir()
#' # use_bcbio_projects(path,nfcore="nf-core/rnaseq",copy=TRUE)
use_bcbio_projects <- function(path, nfcore=NULL, metadata=NULL, git=TRUE, gh=FALSE, org=NULL, copy=FALSE) {

  ui_info("Creating project at {ui_value(path)}")
  if (!fs::dir_exists(path))
    fs::dir_create(path, mode = "u=xrw,g=xwr,o=r", recurse = TRUE)

  ui_info("Populating base project")
  base_template <- fs::path_package("bcbioR", "templates", "base")
  fs::dir_copy(base_template, path, overwrite = FALSE)

  if (is.null(nfcore)){
    is_nfcore_ready <- ui_yeah("Have you already run nf-core pipeline?",
                               n_yes=1, n_no =1)
    if (is_nfcore_ready){
      nfcore <- readline("? Enter path to nf-core output: ")
    }else{
      ui_warn("Please, turn copy = TRUE to only deploy files or")
      ui_abort("Please use {.run use_bcbio_projects} again when you have the nf-core output.")
    }
    use_bcbio_analysis(path, nfcore, copy, metadata)
  }else{
    if (fs::dir_exists(nfcore)){
      ui_info("Checking {.ui_value(nfcore)} as nf-core output directory")
      use_bcbio_analysis(path, nfcore, copy, metadata)
    }else if (copy){
      # deploy only files
      ui_info("Deploying only templates without pipeline information.")
      use_bcbio_analysis(path, nfcore, metadata=metadata, copy = TRUE)
    }else{
      ui_warn("Please, provide nfcore working directory or")
      ui_warn("turn copy = TRUE to only deploy files.")
    }
  }

  if (git){
    ui_info("Create Git local repo at {ui_value(path)}")
    use_git()
  }
  if (gh){
    ui_info("Create GitHub repo at {ui_value(path)}")
    whoami <- suppressMessages(gh::gh_whoami())
    if (is.null(whoami)) {
      ui_abort(c(
        "x" = "Unable to discover a GitHub personal access token.",
        "i" = "A token is required in order to create and push to a new repo.",
        "_" = "Call {.run usethis::gh_token_help()} for help configuring a token."
      ))
    }
    use_github(organisation=org)
  }

  answer <- ui_yeah("Please, read the README.md file as the session starts.Are you ready?",
                    n_yes=1, n_no =1)
  if (answer)
    proj_activate(path)
  if (!answer)
    ui_info("Please use {.run proj_activate({ui_value(path)})} to start this project.")

}
