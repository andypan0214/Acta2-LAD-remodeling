# The purpose of the Rsciprt is to store useful function for cut&Tag project 

library(GenomicRanges) # library for handling genomic ranges
library(ChIPseeker) # library for ChIP-seq data analysis, including peak annotation and visualization
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # library for mouse genome annotation, providing gene models and annotations for the mm10 genome assembly
library(org.Mm.eg.db) # library for mouse gene annotation, providing gene information and annotations for mouse genes
library(clusterProfiler) # library for functional enrichment analysis, including gene ontology and pathway analysis
library(ggplot2) # library for data visualization, providing a powerful and flexible system for creating complex and informative graphics
library(enrichplot) # library for visualizing enrichment results, providing functions for creating various types of plots to visualize the results of functional enrichment analysis



# The function is designed to filter the masterPeak, the data with all information, by selected unique/strick unique peak from gopeak: 
#   1. unique peak indicate those peak that select amount one type of replicate only, like 3 replicate 
#   2. strick peak indicate those peak that select amount one type of replicate only, like 3 replicate and the peak is not found in other condition

# @param res_table: the data frame with all information, including the peak information and the count information
# @param bed_path: the path of the bed file that contains the unique/strick unique peak information
# @param master_gr: the GRanges object that contains the information of all peaks, including the unique/strick unique peak and the non-unique/strick unique peak
# @param select_index: the index of the column in the res_table that contains the count information, which will be used to filter the master_gr removing the low signal range,
#  usually when signal < 5

extract_target_peak<- function(res_table, bed_path, master_gr, select_index){
  
  # 1. Load the Unique Bed 
  data<- read.table(bed_path, header = FALSE)
  grInfo<- GRanges(seqnames = data$V1,
                   ranges = IRanges(start = data$V2, end = data$V3))
  
  # 2. Align coordinates: get the Granges used in the results_table this use the selectR index to match the peak
  analysis_gr<- master_gr[select_index]
  
  # 3. perfom genomic Overlap 
  # find which of the results match the unique/strict peak
  overlap<- findOverlaps(analysis_gr, grInfo)
  
  # 4. Extract the quantitative data
  results_overlap<- res_table[queryHits(overlap), ]
  
  message("Sucessfully matched ", nrow(results_overlap), " out of ", nrow(data))
  return(results_overlap)
  
}

## inverse Match Execution Function 
# @param peakInfo: the data frame with all information, including the peak information and the count information
# @param tssRegion: the region around the transcription start site (TSS) to consider for annotation, default is c(-5000, 1000), which means 5kb upstream and 1kb downstream of the TSS
# @param TxDb: the TxDb object for gene annotation, default is TxDb.Mmusculus.UCSC.mm10.knownGene, which provides gene models and annotations for the mm10 genome assembly
# @param annoDb: the annotation database to use for gene annotation, default is "org.Mm.eg.db", which provides gene information and annotations for mouse genes
chipSeekFunc<- function(peakInfo, tssRegion = c(-5000, 1000), 
                     TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                     annoDb = "org.Mm.eg.db"){
  peak_anno<- annotatePeak(makeGRangesFromDataFrame(peakInfo), 
                          tssRegion = tssRegion, 
                          TxDb = TxDb, 
                          annoDb = annoDb)
  annotated_df<- as.data.frame(peak_anno)
  
  # the following function is design for mouse gene annotation database, if the user want to use other annotation database, like human gene annotation database, then the function will not perform the ID conversion, but return the annotated data frame with ENSEMBL ID, which can be used for downstream analysis
  if(annoDb == "org.Mm.eg.db"){
    id_map <- bitr(annotated_df$ENSEMBL,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
    
    annotated_df <- merge(annotated_df, id_map,
                               by = "ENSEMBL",
                               all.x = TRUE)
    colnames(annotated_df)
    
    return(annotated_df)    
  }
  
  # the following function is design if the user want to use other annotation database, like human gene annotation database, then the function will not perform the ID conversion, but return the annotated data frame with ENSEMBL ID, which can be used for downstream analysis
}

## GO pathway analysis function 
# @param annotated_df: the data frame with annotated peak information, which should contain the ENTREZID column for gene ID
# @param pvalueCutoff: the p-value cutoff for significant GO terms, default is 0.05
# @param method: the ontology category to analyze, default is "BP" for Biological Process, other options include "MF" for Molecular Function and "CC" for Cellular Component
# @param DB: the annotation database to use for GO analysis, default is org.Mm.eg.db for mouse gene annotation, which provides gene information and annotations for mouse genes. If the user want to use other annotation database, like human gene annotation database, then the function will not perform the ID conversion, but return the annotated data frame with ENSEMBL ID, which can be used for downstream analysis

GOPathFunc <- function(annotated_df, pvalueCutoff = 0.05, method = "BP", DB = org.Mm.eg.db){
  
  # 1. Run Enrichment - Keep this as an enrichResult object for plotting
  go_enrich_obj <- clusterProfiler::enrichGO(
    gene          = annotated_df$ENTREZID,
    OrgDb         = DB,
    keyType       = "ENTREZID",
    ont           = method,
    pAdjustMethod = "BH",
    pvalueCutoff  = pvalueCutoff,
    readable      = TRUE)
  
  # 2. Check if any pathways were found to avoid errors in plotting
  if (is.null(go_enrich_obj) || nrow(as.data.frame(go_enrich_obj)) == 0) {
    message("No significant GO terms found.")
    return(NULL)
  }
  
  # 3. Generate the Dotplot
  # Crucial: Use the 'go_enrich_obj', NOT the data frame
  print(
    enrichplot::dotplot(go_enrich_obj, showCategory = 10) + 
      ggtitle(paste("GO", method, "Enrichment")) +
      theme_minimal()
  )
  
  # 4. Return the results as a data frame for the user's records
  return(as.data.frame(go_enrich_obj))
}



