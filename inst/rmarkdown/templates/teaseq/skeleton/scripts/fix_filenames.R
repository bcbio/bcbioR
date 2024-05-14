library(stringr)
fsq=list.files("data/Athe_TEA_Seq_3rdbatch/240417_GEX_AT12013_fastq/", pattern = "_S[1-9]")


dir.create("data/Athe_TEA_Seq_3rdbatch/240417_GEX_AT12013_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")  
  fnnew=paste0("20240417_GEX",ending)
  fs::link_create(path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_GEX_AT12013_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_GEX_AT12013_fastq_fixed",fnnew)))  
})


fsq=list.files("data/Athe_TEA_Seq_3rdbatch/240417_ATAC_AT12013_fastq", pattern = "_S[1-9]")
dir.create("data/Athe_TEA_Seq_3rdbatch/240417_ATAC_AT12013_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")  
  fnnew=paste0("20240417_ATAC",ending)
  fs::link_create(path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_ATAC_AT12013_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_ATAC_AT12013_fastq_fixed",fnnew)))  
})


fsq=list.files("data/Athe_TEA_Seq_3rdbatch/240417_ADT_AT12013_fastq/", pattern = "_S[1-9]")
dir.create("data/Athe_TEA_Seq_3rdbatch/240417_ADT_AT12013_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")  
  fnnew=paste0("20240417_ADT",ending)
  fs::link_create(path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_ADT_AT12013_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_ADT_AT12013_fastq_fixed",fnnew)))  
})

fsq=list.files("data/Athe_TEA_Seq_3rdbatch/240417_HTO_AT12013_fastq/", pattern = "_S[1-9]")
dir.create("data/Athe_TEA_Seq_3rdbatch/240417_HTO_AT12013_fastq_fixed")
sapply(fsq, function(fn){
  ending=str_extract(fn, "_S.*.gz")  
  fnnew=paste0("20240417_HTO",ending)
  fs::link_create(path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_HTO_AT12013_fastq",fn)),
                  new_path = fs::path_abs(file.path("data/Athe_TEA_Seq_3rdbatch/240417_HTO_AT12013_fastq_fixed",fnnew)))  
})
