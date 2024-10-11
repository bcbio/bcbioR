library(stringr)
fsq=list.files("data/GEX_fastq/", pattern = "_S[1-9]")


dir.create("data/GEX_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")
  fnnew=paste0("GEX",ending)
  fs::link_create(path = fs::path_abs(file.path("data/GEX_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/GEX_fastq_fixed",fnnew)))
})


fsq=list.files("data/ATAC_fastq", pattern = "_S[1-9]")
dir.create("data/ATAC_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")
  fnnew=paste0("ATAC",ending)
  fs::link_create(path = fs::path_abs(file.path("data/ATAC_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/ATAC_fastq_fixed",fnnew)))
})


fsq=list.files("data/ADT_fastq/", pattern = "_S[1-9]")
dir.create("data/ADT_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")
  fnnew=paste0("ADT",ending)
  fs::link_create(path = fs::path_abs(file.path("data/ADT_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/ADT_fastq_fixed",fnnew)))
})

fsq=list.files("data/HTO_fastq/", pattern = "_S[1-9]")
dir.create("data/HTO_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")
  fnnew=paste0("HTO",ending)
  fs::link_create(path = fs::path_abs(file.path("data/HTO_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/HTO_fastq_fixed",fnnew)))
})
