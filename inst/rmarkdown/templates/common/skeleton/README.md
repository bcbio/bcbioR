# Guidelines

## Set up work-space

-   [ ] Replace the title in this file to match the project's title
-   [ ] Modify `information.R` with the right text for this project, it can be used to source in other `Rmd` files. The main `Rmd` file in this directory can be used to show general information of the project if needed.
-   [ ] Use the same project name to create a folder in *Dropbox* and a repo in *GitHub*
-   [ ] Use the function `bcbio_templates` to create templates inside `reports` for each type of analysis. For instance, for *RNAseq*:
    -   `bcbio_templates(type="rnaseq", outpath="reports")` or
    -   `bcbio_templates(type="rnaseq", outpath="reports/experiment1")`
    -   Then go to that folder and read the `README.md`

## Folders

-   `meta` should contain the CSV/YAML files used by *bcbio* or *nextflow*
-   `scripts` should contain `sbatch` scripts or any custom scripts used in this project
-   `data` contains raw data
-   `reports` contains `Rmd` and `html` together with their files that will be added to *DropBox*. Each type of project have different guidelines.
-   `final` contains the output of *nextflow/bcbio*
-   For any relevant files or papers use the `docs` folder on *DropBox*

## Download

-   [ ] Download data to the `data` directory on O2. Check the md5 checksums if available.

## Analysis

-   [ ] Make sure that final folder is copied from *scratch* or *S3* to `/n/data1/cores/bcbio/PIs/`

## GitHub

-   [ ] Track in *Git* this `README` file
-   [ ] Track in *Git* files in `scripts`, `meta`, and `reports` that belongs to these type:
    -   **Note** Git add `*.Rmd *.R *ipynb *.sh *.yaml`. (feel free use `.gitignore` if you use a GUI for non-tracked files). *DO NOT* use `git add *`. *DO NOT* track `html/csv/figures`
-   [ ] Commit files and push to *Github* as necessary throughout the project, but especially when work is complete

## Dropbox

-   [ ] Add to the *DropBox* folder all files in `reports`
