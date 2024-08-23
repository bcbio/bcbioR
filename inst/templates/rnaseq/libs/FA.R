# library(msigdb)
# msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.5')
#
# get_databases_v2=function(){
#   all_in_life=list(
#     GOBP=subsetCollection(msigdb.hs, 'c5', 'GO:BP'),
#     GOMF=subsetCollection(msigdb.hs, 'c5', 'GO:MF'),
#     HALLMARK=subsetCollection(msigdb.hs, 'h'),
#     KEGG=subsetCollection(msigdb.hs, 'c2', 'CP:KEGG')
#   ) %>% lapply(., function(geneset){
#     gs=lapply(geneset, function(x){
#       geneIds(x)
#     })
#     names(gs)=sapply(geneset, setName)
#     gs
#   })
#   all_in_life
# }

get_databases=function(){
  all_in_life=list(
    msigdbr(species = "human", category = "H") %>% mutate(gs_subcat="Hallmark"),
    #  msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME"),
    msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG"),
    #  msigdbr(species = "human", category = "C2", subcategory = "CP:PID"),
    msigdbr(species = "human", category = "C5", subcategory = "GO:BP"),
    msigdbr(species = "human", category = "C5", subcategory = "GO:MF")
    #  msigdbr(species = "human", category = "C5", subcategory = "HPO"),
    #  msigdbr(species = "human", category = "C3", subcategory = "TFT:GTRD"),
    #  msigdbr(species = "human", category = "C6") %>% mutate(gs_subcat="Oncogenic")
  )
all_in_life
}

run_fora_v2=function(input, uni, all_in_life){
  # browser()
  total_deg=length(unique(input$SYMBOL))/length(unique(uni$SYMBOL))
  pathways_ora_all = lapply(names(all_in_life), function(database){
    pathway = all_in_life[[database]]
    #pathway = split(x = p$entrez_gene, f = p$gs_name)
    #db_name = paste(p$gs_cat[1], p$gs_subcat[1],sep=":")
    respath <- fora(pathways = pathway,
                    genes = unique(input$SYMBOL),
                    universe = unique(uni$SYMBOL),
                    minSize  = 15,
                    maxSize  = 500)
    coll_respath = collapsePathwaysORA(respath[order(pval)][padj < 0.1],
                                       pathway, unique(input$SYMBOL), unique(uni$SYMBOL))
    as_tibble(respath[pathway %in% coll_respath$mainPathways])  %>%
      mutate(database=db_name, NES=(overlap/size)/(total_deg))
  }) %>% bind_rows() %>%
    mutate(analysis="ORA")
  ora_tb = pathways_ora_all %>% unnest(overlapGenes) %>%
    group_by(pathway) %>%
    left_join(uni, by =c("overlapGenes"="SYMBOL")) %>%
    dplyr::select(pathway, padj, NES, SYMBOL, analysis,
                  database) %>%
    group_by(pathway,padj,NES,database,analysis) %>%
    summarise(genes=paste(SYMBOL,collapse = ","))
  ora_tb

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
    coll_respath = collapsePathwaysORA(respath[order(pval)][padj < 0.1],
                                       pathway, unique(input$ENTREZID), unique(uni$ENTREZID))
    as_tibble(respath[pathway %in% coll_respath$mainPathways])  %>%
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
