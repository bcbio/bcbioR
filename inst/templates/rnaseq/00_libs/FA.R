library(msigdbr)
library(clusterProfiler)
source <- "https://github.com/bcbio/resources/raw/refs/heads/main/gene_sets/gene_sets/20240904"
get_databases_v2=function(sps="human"){
  gmt.files=list(human=c("h.all.v2024.1.Hs.entrez.gmt",
                         #"c5.go.v2024.1.Hs.entrez.gmt",
                         "c5.go.mf.v2024.1.Hs.entrez.gmt",
                         "c5.go.cc.v2024.1.Hs.entrez.gmt",
                         "c5.go.bp.v2024.1.Hs.entrez.gmt",
                         "c2.cp.reactome.v2024.1.Hs.entrez.gmt",
                         "c2.cp.kegg_legacy.v2024.1.Hs.entrez.gmt"),
                 mouse=c("mh.all.v2024.1.Mm.entrez.gmt",
                         #"m5.go.v2024.1.Mm.entrez.gmt",
                         "m5.go.mf.v2024.1.Mm.entrez.gmt",
                         "m5.go.cc.v2024.1.Mm.entrez.gmt",
                         "m5.go.bp.v2024.1.Mm.entrez.gmt",
                         "m2.cp.reactome.v2024.1.Mm.entrez.gmt",
                         "m2.cp.kegg_legacy.v2024.1.Mm.entrez.gmt"))
  all_in_life=lapply(gmt.files[[sps]], function(gmt){
    df=read.gmt(file.path(source,sps,gmt))
    names(df)=c("gs_name", "entrez_gene")
    df
  })
  names(all_in_life) = str_remove(gmt.files[[sps]], ".v2024.*$")
  all_in_life
}

get_databases=function(sps="human"){
  all_in_life=list(
    msigdbr(species = sps, category = "H") %>% mutate(gs_subcat="Hallmark"),
    #  msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME"),
    msigdbr(species = sps, category = "C2", subcategory = "CP:KEGG"),
    #  msigdbr(species = "human", category = "C2", subcategory = "CP:PID"),
    msigdbr(species = sps, category = "C5", subcategory = "GO:BP"),
    msigdbr(species = sps, category = "C5", subcategory = "GO:MF")
    #  msigdbr(species = "human", category = "C5", subcategory = "HPO"),
    #  msigdbr(species = "human", category = "C3", subcategory = "TFT:GTRD"),
    #  msigdbr(species = "human", category = "C6") %>% mutate(gs_subcat="Oncogenic")
  )
  all_in_life
}

run_fora_v2=function(input, uni, all_in_life){
  total_deg=length(unique(input$ENTREZID))/length(unique(uni$ENTREZID))
  pathways_ora_all = lapply(names(all_in_life), function(database){
    p = all_in_life[[database]]
    #browser()
    pathway = split(x = p$entrez_gene, f = p$gs_name)
    db_name = database
    respath <- fora(pathways = pathway,
                    genes = unique(input$ENTREZID),
                    universe = unique(uni$ENTREZID),
                    minSize  = 15,
                    maxSize  = 500)
    respath  %>%
      mutate(database=db_name, NES=(overlap/size)/(total_deg))
  }) %>% bind_rows() %>%
    mutate(analysis="ORA")
  ora_tb = pathways_ora_all %>% unnest(overlapGenes) %>%
    group_by(pathway) %>%
    left_join(uni, by =c("overlapGenes"="ENTREZID")) %>%
    dplyr::select(pathway, padj, NES, SYMBOL, analysis,
                  database) %>%
    group_by(pathway,padj,NES,database,analysis) %>%
    summarise(genes=paste(SYMBOL,collapse = ","))
  ora_tb

}

run_fgsea_v2=function(input, all_in_life){
  # browser()
  input_gsea <- input$lfc
  names(input_gsea) <- input$ENTREZID
  pathways_all = lapply(names(all_in_life), function(database){
    p = all_in_life[[database]]
    pathway = split(x = p$entrez_gene, f = p$gs_name)
    db_name = database
    respath <- fgsea(pathways = pathway,
                     stats = input_gsea,
                     minSize  = 15,
                     maxSize  = 500)

    as_tibble(respath)  %>%
      mutate(database=db_name)
  }) %>% bind_rows() %>%
    mutate(analysis="GSEA")
  tb = pathways_all %>% unnest(leadingEdge) %>%
    group_by(pathway) %>%
    left_join(input, by =c("leadingEdge"="ENTREZID")) %>%
    dplyr::select(pathway, padj, size, NES, SYMBOL, analysis,
                  database) %>%
    group_by(pathway, padj, size, NES, database, analysis) %>%
    summarise(genes=paste(SYMBOL,collapse = ","))
  tb

}

run_fora=function(input, uni,all_in_life){
  # browser()
  total_deg=length(unique(input))/length(unique(uni$ENTREZID))
  pathways_ora_all = lapply(all_in_life, function(p){
    pathway = split(x = p$entrez_gene, f = p$gs_name)
    db_name = paste(p$gs_cat[1], p$gs_subcat[1],sep=":")
    respath <- fora(pathways = pathway,
                    genes = unique(input$ENTREZID),
                    universe = unique(uni$ENTREZID),
                    minSize  = 15,
                    maxSize  = 500)
    # coll_respath = collapsePathwaysORA(respath[order(pval)][padj < 0.1],
    #                                    pathway, unique(input$ENTREZID), unique(uni$ENTREZID))
    as_tibble(respath)  %>%
      mutate(database=db_name, NES=(overlap/size)/(total_deg))
  }) %>% bind_rows() %>%
    mutate(analysis="ORA")
  ora_tb = pathways_ora_all %>% unnest(overlapGenes) %>%
    group_by(pathway) %>%
    left_join(uni, by =c("overlapGenes"="ENTREZID")) %>%
    dplyr::select(pathway, padj, NES, SYMBOL, analysis,
                  database) %>%
    group_by(pathway,padj,NES,database,analysis) %>%
    summarise(genes=paste(SYMBOL,collapse = ","))
  ora_tb

}

run_fgsea=function(input, all_in_life){
  # browser()
  input_gsea <- input$lfc
  names(input_gsea) <- input$ENTREZID
  pathways_all = lapply(all_in_life, function(p){
    pathway = split(x = p$entrez_gene, f = p$gs_name)
    db_name = paste(p$gs_cat[1], p$gs_subcat[1],sep=":")
    respath <- fgsea(pathways = pathway,
                     stats = input_gsea,
                     minSize  = 15,
                     maxSize  = 500)

    as_tibble(respath)  %>%
      mutate(database=db_name)
  }) %>% bind_rows() %>%
    mutate(analysis="GSEA")
  tb = pathways_all %>% unnest(leadingEdge) %>%
    group_by(pathway) %>%
    left_join(input, by =c("leadingEdge"="ENTREZID")) %>%
    dplyr::select(pathway, padj, size, NES, SYMBOL, analysis,
                  database) %>%
    group_by(pathway, padj, size, NES, database, analysis) %>%
    summarise(genes=paste(SYMBOL,collapse = ","))
  tb

}
