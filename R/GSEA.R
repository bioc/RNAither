gseaAnalysis<-function(hitVector, whichOnto){

    hitVector<-factor(hitVector)
    
    GOlist<-GOannotate(names(hitVector), whichOnto)
    
    if(whichOnto == "biological_process"){
        whichOnto<-"BP"
    }
    if(whichOnto == "molecular_function"){
        whichOnto<-"MF"    
    }
    if(whichOnto == "cellular_component"){
        whichOnto<-"CC"    
    }
    
    GOdata<-new("topGOdata", ontology=whichOnto, allGenes=hitVector, 
    annot=geneNameAnno, gene2GO=GOlist)

    testStat<-new("weightCount", testStatistic=GOFisherTest, name="Fisher test", 
    sigRatio="ratio")
    resultWeight<-getSigGroups(GOdata, testStat)

    l <- GenTable(GOdata, resultWeight, topNodes=length(usedGO(GOdata)))
}



gseaAnalysisPt1<-function(hitVector, whichOnto){
    GOlist<-GOannotate(names(hitVector), whichOnto)
}


gseaAnalysisPt2<-function(hitVector, GOlist, whichOnto){

    hitVector<-factor(hitVector)

    if(whichOnto == "biological_process"){
        whichOnto<-"BP"
    }
    if(whichOnto == "molecular_function"){
        whichOnto<-"MF"    
    }
    if(whichOnto == "cellular_component"){
        whichOnto<-"CC"    
    }
    
    GOdata<-new("topGOdata", ontology=whichOnto, allGenes=hitVector, 
    annot=geneNameAnno, gene2GO=GOlist)

    testStat<-new("weightCount", testStatistic=GOFisherTest, name="Fisher test", 
    sigRatio="ratio")
    resultWeight<-getSigGroups(GOdata, testStat)

    l <- GenTable(GOdata, resultWeight, topNodes=length(usedGO(GOdata)))
}


