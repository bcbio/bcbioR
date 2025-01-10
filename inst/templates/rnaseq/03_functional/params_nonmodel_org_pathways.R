# set species
# example: Candida
# species_name       = "Candida albicans"
# species_id         = "cal"
# example: Salmonella
species_name     = "Salmonella enterica subsp. enterica serovar Typhimurium 14028S"
species_id       = "seo"

# set input data: DE output (1 row per gene or protein with LFC and adjusted p-value)
# example: Candida
# input_file = "Candida_DEproteins.csv"
# example: Salmonella
input_file = "Salmonella_DEgenes.csv"

# set columns in the input_file for:
#   log2 fold change (lfc)
#   adjusted p-value
#   Uniprot ID (or other mappable ID)
#   gene name
# example: Candida
# colname_lfc = "logFC"
# colname_padj = "adj.P.Val"
# colname_uniprot = "UniprotID"
# colname_gene = "GeneSymbol"
# example: Salmonella
colname_lfc = "Log2FC_Reads"
colname_padj = "P_Value"
colname_uniprot = "Locus_Tag"
colname_gene = "Gene_ID"

# set Uniprot reference data
# get data from UniProt
# - Go to UniProt: https://www.uniprot.org/
# - Search for genes with format "STM14_####"
# - Determine that this is likely Salmonella typhimurium (strain 14028s / SGSC 2262)
# - Search UniProt for taxonomy ID 588858
# - Click "Customize columns"
# - Make sure at least the following are selected:
#     - Names & Taxonomy > Entry Name, Gene Names, Organism, Protein names
#     - External Resources > Sequence > RefSeq
#     - External Resources > Genome annotation > GeneID
#     - Gene Ontology > all
# - Click "Download"
# - Select "Download all", Format = "Excel"
#     - DO NOT CHOOSE "GFF" as this will not give you all the columns you just selected
#     - DO NOT CHOOSE "TSV" as a lot of the columns get messed up for unknown reasons
#     - "Excel" is an odd choice but at least it works
# - Save this file in the "reference" folder
# example: Candida
# uniprot_ref_file = "Candida_reference_uniprotkb_taxonomy_id_237561_2024_07_31.xlsx"
# example: Salmonella
uniprot_ref_file = "Salmonella_reference_uniprotkb_taxonomy_id_588858_2024_09_11.xlsx"

# set column for Uniprot ID (or other mappable ID) in reference data
# example: Candida
# colname_ref = "Entry"
# example: Salmonella
colname_ref = "KEGG"