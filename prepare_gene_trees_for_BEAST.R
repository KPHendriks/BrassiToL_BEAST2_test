

#Read in species tree

library(ape)
species_tree<-ape::read.tree(file = "iqtree_ML_tree.tree")

#Get tip names. These should all be present in gene trees also.
tips<-species_tree$tip.label

#Read list of fasta files.
alignments<-list.files(path= "data_gene_trees_clocklike/",
                       pattern = "_taper_final_inclusive.fasta")

#Loop through the alignments to check for sample presence,
#add empty strings for missing samples,
#then export as nexus file that can be read into BEAST2.
library("seqinr")
library("chopper")
for (a in 1:length(alignments)){
  gene<-strsplit(alignments[a], "_")[[1]][1]
  alignment<-seqinr::read.fasta(paste0("data_gene_trees_clocklike/", alignments[a]), 
                                seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
  tips_from_alignment<-names(alignment)
  #Get tips not yet in alignment but in species tree.
  tips_to_add<-tips[!tips %in% tips_from_alignment]
  
  #Quick check if any data in alignment that's not in species tree.
  if(sum(!tips_from_alignment %in% tips)>0){
    print("CAUTION! ALIGNMENTS CONTAINS SAMPLES NOT IN SPECIES TREE!")
  } else {
    print("All samples in alignment found from species tree.")
  }
  
  #Second loop to add strings of "N" to alignment for missing samples.
  for (t in 1:length(tips_to_add)){
    alignment<-c(alignment,
                 strrep("N", nchar(alignment[1])))
    names(alignment)[length(names(alignment))]<-tips_to_add[t]
  }
  #Export alignment as nexus file for BEAST2 to use.
  seqinr::write.fasta(sequences = alignment, 
                      names = names(alignment),
                      open = "w",
                      file.out = paste0("data_gene_trees_clocklike/", gene, "_ready_for_BEAST2.fasta"),
                      nbchar = 10000)
  
  chopper::alg2nex(file = paste0("data_gene_trees_clocklike/", gene, "_ready_for_BEAST2.fasta"),
                   format = "fasta")
  
}
