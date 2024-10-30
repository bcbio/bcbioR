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
  }else if (any(grep("[^a-zA-Z0-9_]", samplesheet[["sample"]]))){
    stop("Sample names should contain only letters, numbers, and underscores")
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
#' @param org string with the organization name. To deploy specific files.
#' @examples
#'  \dontrun{
#'   path <- withr::local_tempdir()
#'   bcbio_templates(type="base",outpath=path)
#'   fs::dir_ls(path,all=T)
#'  }
#' @export
bcbio_templates <- function(type="rnaseq", outpath=NULL, org=NULL){
  if (type=="all"){
    usethis::ui_info("Showing analysis:")
    msg <- basename(fs::dir_ls(fs::path_package("bcbioR", "templates")))
    return(msg)
  }
  if (is.null(outpath)){
    usethis::ui_stop("outpath needs to be defined.")
  }
  fs::dir_create(outpath)
  switch(type,
         base={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "base", org)
         },
         rnaseq={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "nf-core/rnaseq", org)
         },
         singlecell={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "singlecell", org)
         },
         singlecell_delux={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "singlecell_delux", org)
         },
         spatial={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "spatial", org)
         },
         chipseq={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "chipseq", org)
         },
         multiomics={
           #file.copy(fpath, outpath, recursive = T)
           copy_templates(outpath, "multiomics", org)
         },
         {
           stop(paste('project type not recognize, please choose: ',
                      'rnaseq', 'chipseq',
                      'singlecell','singlecell_delux','spatial'))
         }
  )
}

read_pipeline_info <- function(nfcore){
  # pipeline_info/params_2024-05-28_12-28-51.json
  config <- fs::path_join(c(nfcore, "pipeline_info"))
  params <- fs::dir_ls(config, regexp = "params")
  metadata <- jsonlite::read_json(params)[["input"]]
  # input
  # tmp_rna/pipeline_info/software_versions.yml
  software <- fs::path_join(c(nfcore, "pipeline_info", "software_versions.yml"))
  software_txt <- yaml::read_yaml(software)
  pipeline <- grep("nf-core", names(software_txt$Workflow), value = TRUE)
  # Workflow:
  #   Nextflow: 24.04.1
  # nf-core/rnaseq: 3.14.0
  # check only rnaseq is supported
  if (!(pipeline %in% c("nf-core/rnaseq"))){
    ui_stop("Sorry, we don't yet support {ui_value(pipeline)}")
  }
  list(metadata=metadata, pipeline=pipeline)
}

render_rmd <- function(infile, outfile, ls_data){
  whisker.render(read_file(infile),
                 ls_data) %>%
    write_file(outfile)
}

bcbio_params <-function(nfcore_path, pipeline, metadata, copy){
  ui_info("Reading input files from {ui_value(nfcore_path)}")
  if (pipeline=="nf-core/rnaseq"){
    if (!copy){
      ls_data<-list(
        se_object =fs::path_join(c(nfcore_path, "star_salmon/salmon.merged.gene_counts.rds")),
        metadata_fn = metadata,
        counts_fn = fs::path_join(c(nfcore_path, "star_salmon/salmon.merged.gene_counts.tsv")),
        multiqc_data_dir = fs::path_join(c(nfcore_path, "multiqc/star_salmon/multiqc-report-data/")),
        gtf_fn = fs::path_join(c(nfcore_path, "genome/genome.filtered.gtf")))
      return(ls_data)
    }
  }

}

detect_gitignores <- function(path){
  gits <- fs::dir_ls(path, recurse = TRUE, regexp = 'gitignore')
  sapply(gits, function(fn){
    hidden <- file.path(dirname(fn), paste0(".", basename(fn)))
    fs::file_move(fn, hidden)
  })
}

copy_files_in_folder<- function(origin, remote, is_org=FALSE){
  to_copy <- fs::dir_ls(origin,all = TRUE)
  if (!is_org) {
    to_copy <- grep("org", to_copy,
                    value = TRUE, invert = TRUE)
  }else{
    # don't allow doc files
    to_copy <- grep(".doc.*$", to_copy,
                    value = TRUE, invert = TRUE)
  }
  for (element in to_copy){
    full_new_path <- fs::path_join(c(remote, fs::path_file(element)))

    if (fs::is_dir(element)){
      if (!(fs::dir_exists(full_new_path)) | is_org)
        fs::dir_copy(element, full_new_path, overwrite = is_org)
    }
    if (fs::is_file(element)){
      if (!(fs::file_exists(full_new_path)) | is_org)
        fs::file_copy(element, full_new_path, overwrite = is_org)
    }
  }
  detect_gitignores(remote)
}

deploy_apps <- function(apps, path){
  fs::dir_create(file.path(path, "apps"))
  sapply(names(apps), function(app){
    dest_file=file.path(path, "apps", paste0(app, ".zip"))
    download.file(url = apps[[app]],
                  destfile = dest_file)
    unzip(zipfile = dest_file, exdir = dirname(dest_file))
    fs::file_delete(dest_file)
  })
}

copy_templates <- function(path, pipeline, org=NULL){
  apps=list()
  base = c("bcbioR")
  if (pipeline=="base"){
    parts = c("templates/base")
  }else if(pipeline=="nf-core/rnaseq"){
    parts = c("templates/rnaseq")
  }else if(pipeline=="singlecell"){
    parts = c("templates/singlecell")
    apps=c(apps, scRNAseq_qc="https://github.com/hbc/scRNAseq_qc_app/archive/refs/heads/main.zip")
  }else if(pipeline=="singlecell_delux"){
    parts = c("templates/singlecell_delux")
  }else if(pipeline=="multiomics"){
    parts = c("templates/multiomics")
  }else if(pipeline=="spatial"){
    parts = c("templates/spatial")
  }else if(pipeline=="chipseq"){
    parts = c("templates/chipseq")
  }
  analysis_template <- fs::path_package(base, parts)

  ui_info("Getting templates from {ui_value(analysis_template)}")
  # ls_files <- grep("org", list.files(analysis_template, full.names = TRUE),
  #                  value = TRUE, invert = TRUE)
  # ui_info("{ui_value(length(ls_files))} amount of files to copy")
  copy_files_in_folder(analysis_template, path)
  if (!is.null(org)){
    org_template <- fs::path_package(base, parts, "org", org)
    if (fs::dir_exists(org_template)){
      ui_info("Getting templates from {ui_value(org_template)}")
      copy_files_in_folder(org_template, path, is_org=TRUE)
    }
  }

  # check org folder is in there
  # search for param + _README.md
  # concat file to README.md
  deploy_apps(apps, path)
}

bcbio_render <- function(path, pipeline, data){

  if (pipeline=="nf-core/rnaseq"){
    # analysis_template <- fs::path_package("bcbioR", "templates", "rnaseq", "qc")
    # fs::dir_copy(analysis_template, fs::path_join(c(path, "reports", "qc")), overwrite=TRUE)
    # analysis_template <- fs::path_package("bcbioR", "templates", "rnaseq", "de")
    # fs::dir_copy(analysis_template, fs::path_join(c(path, "reports", "de")), overwrite=TRUE)
    render_rmd(
      fs::path_join(c(path, "reports", "qc", "QC_nf-core.Rmd")),
      fs::path_join(c(path, "reports", "qc", "QC_nf-core.Rmd")),
      data
    )
    render_rmd(
      fs::path_join(c(path, "reports", "de", "DEG.Rmd")),
      fs::path_join(c(path, "reports", "de", "DEG.Rmd")),
      data
    )
    ui_info("Please, to start the analysis, modify these parameter in QC/QC.rmd")
    ui_todo("set genome to hg38, mm10, mm39, or other")
    ui_todo("set factor_of_interest to a column in your metadata")
  }else{
    ui_warn("These are draft templates, are meant to show examples of specific analysis")
    ui_todo("Please, read carefully and adapt to your data and question.")
  }
}

# help with bcbio analysis setup
use_bcbio_analysis <- function(path, pipeline, copy=TRUE, metadata=NULL){

  if (copy){
    # deploy files
    ui_info("Rmd templates will be copied but variables path won't be filled automatically.")
    if (!is.null(metadata)){
      meta_path <- fs::path_join(c(path, "meta", fs::path_file(metadata)))
      if (!(fs::file_exists(metadata)))
        ui_stop("{ui_value(metadata)} doesn't exist.")
      fs::file_copy(metadata, meta_path)
    }
  }
  if (!is.null(pipeline) & fs::dir_exists(pipeline)){
      # ui_stop("{ui_value(nfcore)} doesn't exist. point to nfcore path or turn on copy mode.")
    ui_info("Trying to guess nf-core pipeline at {ui_value(pipeline)}")
    # guess analysis from pipeline file
    information <- read_pipeline_info(pipeline)
    fs::dir_create(fs::path_join(c(path, "meta")))
    meta_path <- fs::path_join(c(path, "meta", fs::path_file(information$metadata)))
    pipeline <- information$pipeline
    if (!is.null(metadata)){
      if (!(fs::file_exists(metadata)))
          ui_stop("{ui_value(metadata)} doesn't exist.")
      fs::file_copy(metadata, meta_path)
    }else{
      if (!fs::file_exists(information$metadata)){
        ui_warn("{ui_value(metadata)} not found. We can only work with local filesytems right now.")
        ui_todo("Please, copy {ui_value(metadata)} to {ui_value(meta_path)}.")
        ui_warn("If this file isn't manually set up, the Rmd code will fail.")
      }else{
        ui_info("Copy metadata to {ui_value(meta_path)}")
        fs::file_copy(information$metadata, meta_path)
      }
      metadata <- meta_path
    }
    path_final <- fs::path_join(c(path, "final"))
    ui_todo("Please, copy nf-core output directory to {ui_value(path_final)}")
  }
  # set all files from analysis
  copy_templates(fs::path_join(c(path, "reports")), pipeline)
  if (fs::dir_exists(pipeline)){
    data <- bcbio_params(nfcore, pipeline, metadata)
    bcbio_render(path, pipeline, data)
  }


}

# Pilot to deploy full projects at once
# path <- withr::local_tempdir()
# use_bcbio_projects(path,pipeline="nf-core/rnaseq",copy=TRUE)
# fs::dir_ls(path)
use_bcbio_projects <- function(path, pipeline=NULL, metadata=NULL,
                               git=TRUE, gh=FALSE, org=NULL, copy=TRUE) {

  ui_info("Creating project at {ui_value(path)}")
  if (!fs::dir_exists(path))
    fs::dir_create(path, mode = "u=xrw,g=xwr,o=r", recurse = TRUE)

  ui_info("Populating base project")
  base_template <- fs::path_package("bcbioR", "templates", "base")
  copy_files_in_folder(base_template, path)

  if (!is.null(pipeline)){
    ui_info("Using this pipeline templates {ui_value(pipeline)}")
    use_bcbio_analysis(path, pipeline, copy = copy, metadata=metadata)
  }
  # is_nfcore_ready <- FALSE
  # if (is.null(pipeline) && rlang::is_interactive()){
  #   is_nfcore_ready <- ui_yeah("Have you already run nf-core pipeline?",
  #                              n_yes=1, n_no =1)
  #   if (is_nfcore_ready && rlang::is_interactive()){
  #     nfcore <- readline("? Enter path to nf-core output: ")
  #   }else{
  #     ui_warn("Please, turn copy = TRUE to only deploy files or,")
  #     ui_stop("Please use {ui_code('use_bcbio_projects')} again when you have the nf-core output.")
  #   }
  #   use_bcbio_analysis(path, nfcore, copy, metadata)
  # }else{
  #   if (fs::dir_exists(nfcore)){
  #     ui_info("Checking {ui_value(nfcore)} as nf-core output directory")
  #     use_bcbio_analysis(path, nfcore, copy, metadata)
  #   }else if (copy){
  #     # deploy only files
  #     ui_info("Deploying only templates without pipeline information.")
  #     use_bcbio_analysis(path, nfcore, copy = TRUE, metadata=metadata)
  #   }else{
  #     ui_warn("Please, provide nfcore working directory or")
  #     ui_warn("turn copy = TRUE to only deploy files.")
  #   }
  # }

  # if (git){
  #   ui_info("Create Git local repo at {ui_value(path)}")
  #   use_git()
  # }
  # if (gh){
  #   ui_info("Create GitHub repo at {ui_value(path)}")
  #   whoami <- suppressMessages(gh::gh_whoami())
  #   if (is.null(whoami)) {
  #     ui_warn(c(
  #       "x" = "Unable to discover a GitHub personal access token.",
  #       "i" = "A token is required in order to create and push to a new repo.",
  #       "_" = "Call {.run usethis::gh_token_help()} for help configuring a token."
  #     ))
  #     ui_todo("Try this later: use_github(organisation=org), private = TRUE")
  #
  #   }
  #   use_github(organisation=org, private = TRUE)
  # }else{
  #   ui_info("You decided not to create a repo, please use this to push when ready")
  #   ui_todo("Try this later: use_github(organisation=org), private = TRUE")
  # }

  answer <- FALSE
  if (rlang::is_interactive())
    answer <- ui_yeah("Please, read the README.md file as the session starts.Are you ready?",
                      n_yes=1, n_no =1, shuffle=FALSE)
  if (answer)
    proj_activate(path)
  if (!answer)
    ui_info("Please use proj_activate({ui_value(path)})} to start this project.")

}
