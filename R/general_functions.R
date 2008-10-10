createSubset<-function(dataset, listIDs, equalTo){
    subset<-dataset[which(listIDs == equalTo), ]
    invisible(subset)
}

indexSubset<-function(listIDs, equalTo){
    ixx<-which(listIDs == equalTo)
    invisible(ixx)
}


findReplicates<-function(dataset, whichCol, replicateID){
    invisible(which(dataset[[get("whichCol")]] == replicateID))
}


orderGeneIDs<-function(dataset, ID1){
    Indiz<-order(dataset[[get("ID1")]])
    geordneteGenListe<-dataset
    for (i in 1:length(dataset[[get("ID1")]])){
        geordneteGenListe[i, ]<-dataset[Indiz[i], ]
    }
    invisible(geordneteGenListe)
}

rms<-function(Ivec, na.rm=T){
    Ivec <- Ivec[!is.na(Ivec)]
    invisible(sqrt(sum(Ivec^2)/length(Ivec)))
}

trim<-function(Ivec, na.rm=T){
    invisible(mean(Ivec, trim=0.05, na.rm=T))
}

closestToZero<-function(Ivec, na.rm=T){
    wh<-which(Ivec == min(abs(Ivec), na.rm=T))
    invisible(Ivec[wh])
}

furthestFromZero<-function(Ivec, na.rm=T){
    wh<-which(Ivec == max(abs(Ivec), na.rm=T))
    invisible(Ivec[wh])
}

divideChannels<-function(ch1, ch2){
    invisible(ch1/ch2)
}

eraseDataSetColumn<-function(dataset, colname){

    if (is.na(colname)){
        invisible(dataset)
    }else{
        ##erase old column ch2:
        numToErase<-which(colnames(dataset) == get("colname"))
        counter<-0
        merken<-0
        for (i in 1:length(colnames(dataset))){
            if (i!=numToErase){
                counter<-counter+1
                if (counter == 1){
                    newDataset<-dataset[[i]]
                    merken<-i
                }
                if (counter>1){
                    newDataset<-data.frame(newDataset, dataset[[i]])
                    colnames(newDataset)[counter]<-colnames(dataset[i])
                }
            }
        }
        colnames(newDataset)[1]<-colnames(dataset)[merken]
        invisible(newDataset)
    }
}

generateReplicateMat<-function(data, minNbReps, IndexOrInt, col4val, col4anno){

    workdata<-data
    workdata[[get("col4val")]][which(workdata$SpotType == -1)]<-NA_character_

    ##checks the maximum number of replicates:
    maxNumRep<-0
    IDs<-unique(workdata[[get("col4anno")]])

    ##gene IDs with more than one replicate:
    newIDs<-IDs
    for (i in 1:length(IDs)){
        repIndex<-findReplicates(workdata, col4anno, IDs[i])
        if(length(repIndex)>maxNumRep){
            maxNumRep<-length(repIndex)
        }
        if (minNbReps == 2){
            if(length(repIndex)<2){
                newIDs[i]<-NA_character_
            }
        }
    }

    newIDs<-newIDs[!is.na(newIDs)]

    CVmatrix<-matrix(as.double(NA_character_), length(newIDs), maxNumRep)
    indexPosControls<-rep(0, nrow(CVmatrix))
    indexNegControls<-rep(0, nrow(CVmatrix))

    posCounter<-1
    negCounter<-1
    for (i in 1:length(newIDs)){
        repIndex<-findReplicates(workdata, col4anno, newIDs[i])

        if (IndexOrInt == "Intensities"){
            colCounter<-0
            for (j in 1:length(repIndex)){
                if (!is.na(workdata[[get("col4val")]][repIndex[j]])){
                    colCounter<-colCounter+1
                    CVmatrix[i, colCounter]<-as.double(workdata[[get("col4val")]][repIndex[j]])
                }
            }
        }

        if (IndexOrInt == "Index"){

            zwischenVek<-rep(as.integer(NA_character_), length(repIndex))
            for (j in 1:length(repIndex)){
                if (!is.na(workdata[[get("col4val")]][repIndex[j]])){
                    zwischenVek[j]<-as.integer(repIndex[j])
                }
            }
            zwischenVek <- zwischenVek[!is.na(zwischenVek)]

            if (length(zwischenVek)>0){
                CVmatrix[i, 1:length(zwischenVek)]<-as.integer(zwischenVek)
            }
        }


        if (workdata$SpotType[repIndex[1]] == 1){
            indexPosControls[posCounter]<-i
            posCounter<-posCounter+1
        }
        if (workdata$SpotType[repIndex[1]] == 0){
            indexNegControls[negCounter]<-i
            negCounter<-negCounter+1
        }
    }


    ##filter out rows containing only NAs:
    rowCount<-0
    for (i in 1:nrow(CVmatrix)){
        if(sum(!is.na(CVmatrix[i, ]))>0){
            rowCount<-rowCount+1
            if (rowCount == 1){
                goodRows<-i
            }
            if (rowCount>1){
                goodRows<-c(goodRows, i)
            }
        }else{
            if (length(which(indexPosControls == i))>0){
                indexPosControls[which(indexPosControls == i)]<-0
            }
            if (length(which(indexNegControls == i))>0){
                indexNegControls[which(indexNegControls == i)]<-0
            }
            ##adapt rownumbers in indexNegControls and indexPosControls:
            indexPosControls[indexPosControls>i]<-indexPosControls[indexPosControls>i]-1
            indexNegControls[indexNegControls>i]<-indexNegControls[indexNegControls>i]-1

        }
    }

    indexPosControls<-indexPosControls[which(indexPosControls!=0)]
    indexNegControls<-indexNegControls[which(indexNegControls!=0)]

    rownames(CVmatrix)<-newIDs

    if(rowCount != 0){
        CVmatrix<-CVmatrix[goodRows, ]
    }else{
        print("#############################################")
        print("#Warning. Replicate matrix contains only NA!#")
        print("#############################################")
    }

    invisible(list(CVmatrix, indexPosControls, indexNegControls))
}


generateRepMatNoFilter<-function(data, minNbReps, IndexOrInt, col4val, col4anno){

    maxNumRep<-0
    IDs<-unique(data[[get("col4anno")]])
    newIDs<-IDs
    for (i in 1:length(IDs)){
        repIndex<-findReplicates(data, col4anno, IDs[i])
        if(length(repIndex)>maxNumRep){
            maxNumRep<-length(repIndex)
        }
        if (minNbReps == 2){
            if(length(repIndex)<2){
                newIDs[i]<-NA_character_
            }
        }
    }

    newIDs<-newIDs[!is.na(newIDs)]

    CVmatrix<-matrix(as.double(NA_character_), length(newIDs), maxNumRep)
    indexPosControls<-rep(0, nrow(CVmatrix))
    indexNegControls<-rep(0, nrow(CVmatrix))

    posCounter<-1
    negCounter<-1
    for (i in 1:length(newIDs)){
        repIndex<-findReplicates(data, col4anno, newIDs[i])

        if (IndexOrInt == "Intensities"){
            for (j in 1:length(repIndex)){
                CVmatrix[i, j]<-as.double(data[[get("col4val")]][repIndex[j]])

            }
        }

        if (IndexOrInt == "Index"){
            for (j in 1:length(repIndex)){
                CVmatrix[i, j]<-as.integer(repIndex[j])
            }
        }


        if (data$SpotType[repIndex[1]] == 1){
            indexPosControls[posCounter]<-i
            posCounter<-posCounter+1
        }
        if (data$SpotType[repIndex[1]] == 0){
            indexNegControls[negCounter]<-i
            negCounter<-negCounter+1
        }


    }

    indexPosControls<-indexPosControls[which(indexPosControls!=0)]
    indexNegControls<-indexNegControls[which(indexNegControls!=0)]

    rownames(CVmatrix)<-newIDs

    invisible(list(CVmatrix, indexPosControls, indexNegControls))

}


summarizeReps<-function(data, funSum, col4val, col4anno, cols2del){

    for (a in 1:length(cols2del)){
        data<-eraseDataSetColumn(data, cols2del[a])
    }

    data<-data[which(data$SpotType!=-1), ]

    tp<-generateRepMatNoFilter(data, 1, "Index", col4val[1], col4anno)
    replicateMatrix<-tp[[1]]

    Spotnumber<-rep(0, nrow(replicateMatrix))
    SpotType<-rep(0, nrow(replicateMatrix))
    Internal_GeneID<-rep(NA_character_, nrow(replicateMatrix))
    GeneName<-rep(NA_character_, nrow(replicateMatrix))
    LabtekNb<-rep(0, nrow(replicateMatrix))
    RowNb<-rep(0, nrow(replicateMatrix))
    ColNb<-rep(0, nrow(replicateMatrix))
    ScreenNb<-rep(0, nrow(replicateMatrix))

    listForFlexDataCols<-list("platzhalter")
    for (b in 1:length(col4val)){
        listForFlexDataCols[[b]]<-rep(0, nrow(replicateMatrix))
    }

    for (i in 1:nrow(replicateMatrix)){
        tp<-replicateMatrix[i, !is.na(replicateMatrix[i, ])]

        Spotnumber<-""
        LabtekNb<-""
        RowNb<-""
        ColNb<-""
        ScreenNb<-""

        for (j in 1:length(tp)){

            Spotnumber<-paste(Spotnumber, data$Spotnumber[tp[j]], sep=",")
            LabtekNb<-paste(LabtekNb, data$LabtekNb[tp[j]], sep=",")
            RowNb<-paste(RowNb, data$RowNb[tp[j]], sep=",")
            ColNb<-paste(ColNb, data$ColNb[tp[j]], sep=",")
            ScreenNb<-paste(ScreenNb, data$ScreenNb[tp[j]], sep=",")
        }
        
        Spotnumber<-substr(Spotnumber, 2, nchar(Spotnumber))
        LabtekNb<-substr(LabtekNb, 2, nchar(LabtekNb))
        RowNb<-substr(RowNb, 2, nchar(RowNb))
        ColNb<-substr(ColNb, 2, nchar(ColNb))
        ScreenNb<-substr(ScreenNb, 2, nchar(ScreenNb))


        Spotnumber[i]<-Spotnumber
        LabtekNb[i]<-LabtekNb
        RowNb[i]<-RowNb
        ColNb[i]<-ColNb
        ScreenNb[i]<-ScreenNb


        SpotType[i]<-data$SpotType[tp[1]]
        Internal_GeneID[i]<-as.character(data$Internal_GeneID[tp[1]])
        GeneName[i]<-as.character(data$GeneName[tp[1]])

        if (length(unique(data$SpotType[tp]))>1){
            s1<-"Replicates for siRNA"
            s2<-"have different spot types:"
            s<-paste(s1, data$Internal_GeneID[tp[1]], s2, data$SpotType[tp], sep=" ")
            print(s)
        }
        
        if (length(unique(data$Internal_GeneID[tp]))>1){
        
            s1<-"Replicates for gene"
            s2<-"have different Internal IDs:"
            s<-paste(s1, data$GeneName[tp[1]], s2, data$Internal_GeneID[tp], sep=" ")
            print(s)

            Internal_GeneID[i]<-""
            idss<-unique(data$Internal_GeneID[tp])
            for (g in 1:length(idss)){
                Internal_GeneID[i]<-paste(Internal_GeneID[i], idss[g], sep=",")
            }
            
            blop<-substr(Internal_GeneID[i], 2, nchar(Internal_GeneID[i]))
            Internal_GeneID[i]<-blop
        }
        
        if (length(unique(data$GeneName[tp]))>1){
            s1<-"Replicates for siRNA"
            s2<-"have different NCBI IDs:"
            s<-paste(s1, data$Internal_GeneID[tp[1]], s2, data$GeneName[tp], sep=" ")
            print(s)
        }

        if (length(tp)>0){

            for (h in 1:length(col4val)){
                currColName<-col4val[h]
                if (sum(!is.na(data[[get("currColName")]][tp])>0)){
                    listForFlexDataCols[[h]][i]<-funSum(data[[get("currColName")]][tp], na.rm=T)
                }else{
                    listForFlexDataCols[[h]][i]<-as.double(NA_integer_)
                }
            }

        }

        if (length(tp) == 0){

            for (h in 1:length(col4val)){
                currColName<-col4val[h]
                listForFlexDataCols[[h]][i]<-as.double(NA_character_)
            }
        }
    }




    summarizedDataset<-data.frame(Spotnumber, Internal_GeneID, 
    GeneName, SpotType, LabtekNb, RowNb, ColNb, ScreenNb, 
    stringsAsFactors=F)
    
    for (c in 1:length(col4val)){
        summarizedDataset<-data.frame(summarizedDataset, listForFlexDataCols[[c]])
    }
    colNums<-length(colnames(summarizedDataset))
    colnames(summarizedDataset)[(colNums-length(col4val)+1):colNums]<-col4val

    invisible(summarizedDataset)
}




summarizeRepsNoFiltering<-function(data, funSum, col4val, col4anno, cols2del){

    for (a in 1:length(cols2del)){
        data<-eraseDataSetColumn(data, cols2del[a])
    }

    for (i in 1:length(col4val)){
        currCol<-col4val[i]
        data[[get("currCol")]][which(data$SpotType == -1)]<-as.double(NA_character_)
    }

    tp<-generateRepMatNoFilter(data, 1, "Index", col4val[1], col4anno)
    replicateMatrix<-tp[[1]]

    Spotnumber<-rep(0, nrow(replicateMatrix))
    SpotType<-rep(0, nrow(replicateMatrix))
    Internal_GeneID<-rep(NA_character_, nrow(replicateMatrix))
    GeneName<-rep(NA_character_, nrow(replicateMatrix))
    LabtekNb<-rep(0, nrow(replicateMatrix))
    RowNb<-rep(0, nrow(replicateMatrix))
    ColNb<-rep(0, nrow(replicateMatrix))
    ScreenNb<-rep(0, nrow(replicateMatrix))
    listForFlexDataCols<-list("platzhalter")

    for (b in 1:length(col4val)){
        listForFlexDataCols[[b]]<-rep(0, nrow(replicateMatrix))
    }

    for (i in 1:nrow(replicateMatrix)){
        tp<-replicateMatrix[i, !is.na(replicateMatrix[i, ])]

        Spotnumber<-""
        LabtekNb<-""
        RowNb<-""
        ColNb<-""
        ScreenNb<-""

        for (j in 1:length(tp)){

            Spotnumber<-paste(Spotnumber, data$Spotnumber[tp[j]], sep=",")
            LabtekNb<-paste(LabtekNb, data$LabtekNb[tp[j]], sep=",")
            RowNb<-paste(RowNb, data$RowNb[tp[j]], sep=",")
            ColNb<-paste(ColNb, data$ColNb[tp[j]], sep=",")
            ScreenNb<-paste(ScreenNb, data$ScreenNb[tp[j]], sep=",")
        }
        
        Spotnumber<-substr(Spotnumber, 2, nchar(Spotnumber))
        LabtekNb<-substr(LabtekNb, 2, nchar(LabtekNb))
        RowNb<-substr(RowNb, 2, nchar(RowNb))
        ColNb<-substr(ColNb, 2, nchar(ColNb))
        ScreenNb<-substr(ScreenNb, 2, nchar(ScreenNb))

        Spotnumber[i]<-Spotnumber
        LabtekNb[i]<-LabtekNb
        RowNb[i]<-RowNb
        ColNb[i]<-ColNb
        ScreenNb[i]<-ScreenNb


        SpotType[i]<-data$SpotType[tp[1]]
        Internal_GeneID[i]<-as.character(data$Internal_GeneID[tp[1]])
        GeneName[i]<-as.character(data$GeneName[tp[1]])


        if (length(unique(data$Internal_GeneID[tp]))>1){
            s1<-"Replicates for gene"
            s2<-"have different Internal IDs:"
            s<-paste(s1, data$GeneName[tp[1]], s2, data$Internal_GeneID[tp], sep=" ")
            print(s)

            Internal_GeneID[i]<-""
            idss<-unique(data$Internal_GeneID[tp])
            for (g in 1:length(idss)){
                Internal_GeneID[i]<-paste(Internal_GeneID[i], idss[g], sep=",")
            }
            
            blop<-substr(Internal_GeneID[i], 2, nchar(Internal_GeneID[i]))
            Internal_GeneID[i]<-blop
        }
        
        if (length(unique(data$GeneName[tp]))>1){
            s1<-"Replicates for siRNA"
            s2<-"have different NCBI IDs:"
            s<-paste(s1, data$Internal_GeneID[tp[1]], s2, data$GeneName[tp], sep=" ")
            print(s)
        }

        if (length(tp)>0){

            for (h in 1:length(col4val)){
                currColName<-col4val[h]
                if (sum(!is.na(data[[get("currColName")]][tp])>0)){
                    listForFlexDataCols[[h]][i]<-funSum(data[[get("currColName")]][tp], na.rm=T)
                }else{
                    listForFlexDataCols[[h]][i]<-as.double(NA_integer_)
                }
            }
        }

        if (length(tp) == 0){

            for (h in 1:length(col4val)){
                currColName<-col4val[h]
                listForFlexDataCols[[h]][i]<-as.double(NA_character_)
            }
        }
    }


    summarizedData<-data.frame(Spotnumber, Internal_GeneID, GeneName, 
    SpotType, LabtekNb, RowNb, ColNb, ScreenNb, 
    stringsAsFactors=F)
    
    for (c in 1:length(col4val)){
        summarizedData<-data.frame(summarizedData, listForFlexDataCols[[c]])
    }
    colNums<-length(colnames(summarizedData))
    colnames(summarizedData)[(colNums-length(col4val)+1):colNums]<-col4val

    invisible(summarizedData)
}


sumChannels<-function(header, dataset, funName, colname4ch1, colname4ch2){

    dataset<-saveOldIntensityColumns(dataset, colname4ch1)
    dataset<-saveOldIntensityColumns(dataset, colname4ch2)

    blop<-funName(dataset[[get("colname4ch1")]], dataset[[get("colname4ch2")]])
    dataset[[get("colname4ch1")]]<-blop
    newDataset<-eraseDataSetColumn(dataset, colname4ch2)

    no<-which(colnames(newDataset) == colname4ch1)
    colnames(newDataset)[[no]]<-paste(colname4ch1, colname4ch2, sep="_")

    header[3]<-gsub("NA", "", header[3])
    s1<-"Summarization of channel "
    s2<-" and "
    header[3]<-paste(header[3], s1, colname4ch1, s2, colname4ch2, ",", sep="")

    invisible(list(header, newDataset))
}

