GOannotate<-function(vecOfGeneNames, whichOnto){

    if(whichOnto == "biological_process"){
        goIDs<-"go_biological_process_id"
    }
    
    if(whichOnto == "molecular_function"){
        goIDs<-"go_molecular_function_id"
    }
    
    if(whichOnto == "cellular_component"){
        goIDs<-"go_cellular_component_id"
    }

    #ensembl<-useMart("ensembl")
    #ensembl<-useDataset("hsapiens_gene_ensembl", mart=ensembl)
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="useast.ensembl.org")


#    resTable<-getBM(attributes=c("hgnc_symbol", goIDs), filters="hgnc_symbol", 
    resTable<-try(getBM(attributes=c("hgnc_symbol", "go_id"), filters="hgnc_symbol",values=vecOfGeneNames, mart=ensembl))

    GOlist<-list("platzhalter")
    newVecOfGeneNames<-list("platzhalter")
    counter<-0
   
    if (class(resTable)=="try-error"){
      invisible(resTable)
    }else{
      for (gene in vecOfGeneNames){    
        if (nrow(resTable[resTable[, 1] == gene,]) > 0){        
            subSet<-resTable[resTable[, 1] == gene, ]
            GoVec<-subSet[, 2]
            GOterms<-GoVec[GoVec!=""]
            if (!is.na(GOterms[1])){
                counter<-counter+1
	        GOlist[[counter]]<-GOterms
                newVecOfGeneNames[[counter]]<-gene
            }
        }
      }
      names(GOlist)<-unlist(newVecOfGeneNames)
      invisible(GOlist)
    }
}



geneNameAnno<-function(whichOnto, feasibleGenes=NULL, gene2GO){

    allGO<-unlist(gene2GO, use.names=FALSE)
    allGO<-unique(allGO)
    
    returnList<-list("platzhalter")
    geneCounter<-0
    
    for (i in 1:length(allGO)){
        geneCounter<-0
        
        for (j in 1:length(gene2GO)){
        
    
            if (allGO[i] %in% gene2GO[[j]]){
                geneCounter<-geneCounter+1
                if(geneCounter == 1){
                    returnList[[i]]<-names(gene2GO)[j]
                }
                if(geneCounter>1){
                    returnList[[i]]<-c(returnList[[i]], names(gene2GO)[j])
                }
            }
        }
    
    }
    names(returnList)<-allGO
    
    ontoGO<-get(paste("GO", whichOnto, "Term", sep = ""))
    goodGO<-intersect(ls(ontoGO), allGO)
    return(returnList[goodGO])
}
