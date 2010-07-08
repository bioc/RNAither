LiWongRank<-function(header, dataset, listOfArgs){

    col4val<-listOfArgs[[1]]
    col4anno<-listOfArgs[[2]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    numOfPlates=max(dataset$LabtekNb)

    if (length(unique(dataset$ScreenNb)) == 1){
        stop("Comparison not possible - Only one screen in dataset")
    }

    ##list of genes occuring only once on the plate:
    geneList<-list("platzhalter")

    ##test if requirements for normalization are met:
    for (i in 1:numOfPlates){
        if (length(which(dataset$LabtekNb == i))>0){
        
            tempSub<-createSubset(dataset, dataset$LabtekNb, i)

            ##test if the plate layout is the same:
            for (j in minOfScreens:(numOfScreens-1)){
            
                if (length(which(tempSub$ScreenNb == j))>0 
                & length(which(tempSub$ScreenNb == (j+1)))>0){
                
                    Screen1<-createSubset(tempSub, tempSub$ScreenNb, j)
                    Screen2<-createSubset(tempSub, tempSub$ScreenNb, (j+1))
                    
                    if (sum(Screen1[[get("col4anno")]]!=Screen2[[get("col4anno")]])>0 
                    | length(Screen1[[get("col4anno")]])!=length(Screen2[[get("col4anno")]])){
                    
                        s1<-"Plate layout of Screen"
                        s2<-"Plate"
                        s3<-"is not the same as plate layout of Screen"
                        stop(paste(s1, j, s2, i, s3, (j+1), s2, i, sep=" "))
                    }            
                }
            }
        
            ##test if there is more than one replicate of each siRNA per plate type:
            testPlate<-createSubset(tempSub, tempSub$ScreenNb, 1)
            temp<-generateRepMatNoFilter(testPlate, 1, "Index", col4val, col4anno)
            testReplicaMatrix<-temp[[1]]
            
            if (sum(!is.na(testReplicaMatrix[, 2]))/nrow(testReplicaMatrix)>0.2){
                stop("Too many replicates on same plate - cannot normalize Screen")
            }else{
                Indexes<-which(is.na(testReplicaMatrix[, 2]) == T)
                namesVec<-rownames(testReplicaMatrix)
                geneList[[i]]<-namesVec[Indexes]
            }
        }            
    }


    datasetForComputing<-dataset[which(dataset$SpotType!=-1), ]
    
    for (i in 1:numOfPlates){
        if (length(which(datasetForComputing$LabtekNb == i))>0){
        
            matRow<-length(geneList[[i]])
            matCol<-length(unique(datasetForComputing$LabtekNb))
            
            matOfRanksForPlateTypei<-matrix(0, matRow, matCol)
            rownames(matOfRanksForPlateTypei)<-geneList[[i]]
            ##contains for each plate of plate layout/type i (columns) the rank of each 
            ##gene (rows) on the plate
        
            tempSub<-createSubset(datasetForComputing, datasetForComputing$LabtekNb, i)
            tempSubset<-createSubset(dataset, dataset$LabtekNb, i) 
            IXXi<-indexSubset(dataset$LabtekNb, i)

            screenCount<-0
            for (j in minOfScreens:numOfScreens){
                if (length(which(tempSub$ScreenNb == j))>0){
                    screenCount<-screenCount+1
                    currPlate<-createSubset(tempSub, tempSub$ScreenNb, j)

                    currPlate<-currPlate[which(currPlate[[get("col4anno")]] %in% geneList[[i]]), ]
                    tempp<-orderGeneIDs(currPlate, col4val)
                    orderedListCurr<-tempp[[get("col4anno")]]

                    count1<-0
                    for (gene in geneList[[i]]){
                        count1<-count1+1
                        if (length(which(orderedListCurr == gene))!=0){
                            matOfRanksForPlateTypei[count1, screenCount]<-which(orderedListCurr == gene)
                        }else{
                            matOfRanksForPlateTypei[count1, screenCount]<-as.double(NA_character_)
                        }
                    }
                }
            }
        }

        rankSDpersiRNA<-rep(0, nrow(matOfRanksForPlateTypei))
        names(rankSDpersiRNA)<-rownames(matOfRanksForPlateTypei)
        for (k in 1:length(rankSDpersiRNA)){
            rankSDpersiRNA[k]<-sd(matOfRanksForPlateTypei[k, ], na.rm=T)
        }

        hist(rankSDpersiRNA, breaks<-20)
        
        if(interactive()){
        	s1<-"Plate type"
	        s2<-"- Set a threshold to print out invariant probe set genes\r\n"
	        userInput<-readline(paste(s1, i, s2, sep=" "))
	}else{
	    if(i==1){
	        userInput<-0.5
            }
            if(i==2){
                userInput<-2
            }
	}
        
        tempVecc<-geneList[[i]]

        genesUnderUserThresh<-tempVecc[rankSDpersiRNA<as.double(userInput)]
        
        if (length(genesUnderUserThresh)!=0){
            c1<-as.character(genesUnderUserThresh)
            c2<-as.double(rankSDpersiRNA[rankSDpersiRNA<as.double(userInput)])
            print(cbind(c1, c2))
            
            if(interactive()){
                s1<-"Choose one gene for normalization - Input row number above.\r\n"
                userInput2<-readline(s1)
            }else{
                userInput2<-3
            }
        }else{
            stop("No genes under threshold. Normalization aborted.")
        }
        

        for (j in minOfScreens:numOfScreens){
            if (length(which(tempSubset$ScreenNb == j))>0){
            
                currPlate<-createSubset(tempSubset, tempSubset$ScreenNb, j)
                IXXj<-indexSubset(tempSubset$ScreenNb, j)

                c3<-genesUnderUserThresh[as.integer(userInput2)]
                c4<-currPlate[[get("col4anno")]]
                lineOfGene<-which(c3 == c4)
                
                c5<-currPlate[[get("col4val")]]
                c6<-currPlate[[get("col4val")]][lineOfGene]
                currPlate[[get("col4val")]]<-c5/c6

                tempSubset[[get("col4val")]][IXXj]<-currPlate[[get("col4val")]]
            }
        }
        dataset[[get("col4val")]][IXXi]<-tempSubset[[get("col4val")]]
    }
    header[3]<-gsub("NA", "", header[3])
    header[3]<-paste(header[3], "LiWongRank,")
    invisible(list(header, dataset))
}



varAdjust<-function(header, dataset, listOfArgs){

    flag<-listOfArgs[[1]]
    col4val<-listOfArgs[[2]]
    excludeControlsFlag<-listOfArgs[[3]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    if (flag == 1){
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
            
                subset<-createSubset(dataset, dataset$ScreenNb, i)
                IXX<-indexSubset(dataset$ScreenNb, i)

                specialSubset<-subset[which(subset$SpotType!=-1), ]
                if (excludeControlsFlag == 1){
                    specialSubset<-specialSubset[which(specialSubset$SpotType == 2), ]
                }
                
                c1<-subset[[get("col4val")]]
                c2<-mad(specialSubset[[get("col4val")]], na.rm=T)
                subset[[get("col4val")]]<-c1/c2
                dataset[[get("col4val")]][IXX]<-subset[[get("col4val")]]
            }
        }
    }
    if (flag == 2){
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
            
                subset<-createSubset(dataset, dataset$ScreenNb, i)
                IXXi<-indexSubset(dataset$ScreenNb, i)
                numOfPlates<-max(subset$LabtekNb)
                minOfPlates<-min(subset$LabtekNb)
                
                for (j in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == j))>0){

                        subsubset<-createSubset(subset, subset$LabtekNb, j)
                        IXXj<-indexSubset(subset$LabtekNb, j)

                        specialSubset<-subsubset[which(subsubset$SpotType!=-1), ]

                        if (excludeControlsFlag == 1){
                            specialSubset<-specialSubset[which(specialSubset$SpotType == 2), ]
                        }
                        
                        c1<-subsubset[[get("col4val")]]
                        c2<-mad(specialSubset[[get("col4val")]], na.rm=T)
                        subsubset[[get("col4val")]]<-c1/c2
                        
                        subset[[get("col4val")]][IXXj]<-subsubset[[get("col4val")]]
                    }
                }
                dataset[[get("col4val")]][IXXi]<-subset[[get("col4val")]]
            }
        }
    }


    header[3]<-gsub("NA", "", header[3])
    if (flag == 1){
        header[3]<-paste(header[3], "VarAdjust", col4val, "(by screen/exp),")
    }
    if (flag == 2){
        header[3]<-paste(header[3], "VarAdjust", col4val, "(by plate),")
    }
    invisible(list(header, dataset))
}




divNorm <- function(header, dataset, listOfArgs){
    
    
    funName<-listOfArgs[[1]]
    flag<-listOfArgs[[2]]
    flagForNorm<-listOfArgs[[3]]
    col4val<-listOfArgs[[4]]
    excludeControlsFlag<-listOfArgs[[5]]
    
    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    if (flag == 1){
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-createSubset(dataset, dataset$ScreenNb, i)
                IXX<-indexSubset(dataset$ScreenNb, i)
                
                specialSubset<-subset[which(subset$SpotType!=-1), ]
                
                if (excludeControlsFlag == 1){
                    specialSubset<-specialSubset[which(specialSubset$SpotType == 2), ]
                }

                if (flagForNorm == 1){
                    c1<-subset[[get("col4val")]]
                    c2<-funName(specialSubset[[get("col4val")]], na.rm=T)
                    subset[[get("col4val")]]<-c1/c2
                }
                if (flagForNorm == 2){
                    c1<-subset[[get("col4val")]]
                    c2<-funName(specialSubset[[get("col4val")]], na.rm=T)
                    subset[[get("col4val")]]<-c1-c2
                }
                dataset[[get("col4val")]][IXX]<-subset[[get("col4val")]]
            }
        }
    }
    if (flag == 2){
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
                subset<-createSubset(dataset, dataset$ScreenNb, i)
                IXXi<-indexSubset(dataset$ScreenNb, i)
                numOfPlates<-max(subset$LabtekNb)
                minOfPlates<-min(subset$LabtekNb)
                
                for (j in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == j))>0){
                
                        subsubset<-createSubset(subset, subset$LabtekNb, j)
                        IXXj<-indexSubset(subset$LabtekNb, j)

                        specialSubset<-subsubset[which(subsubset$SpotType!=-1), ]

                        if (excludeControlsFlag == 1){
                            specialSubset<-specialSubset[which(specialSubset$SpotType == 2), ]
                        }

                        if (flagForNorm == 1){
                            c1<-subsubset[[get("col4val")]]
                            c2<-funName(specialSubset[[get("col4val")]], na.rm=T)
                            subsubset[[get("col4val")]]<-c1/c2
                        }
                        if (flagForNorm == 2){
                            c1<-subsubset[[get("col4val")]]
                            c2<-funName(specialSubset[[get("col4val")]], na.rm=T)
                            subsubset[[get("col4val")]]<-c1-c2
                        }
                        subset[[get("col4val")]][IXXj]<-subsubset[[get("col4val")]]
                    }
                }
                dataset[[get("col4val")]][IXXi]<-subset[[get("col4val")]]
            }
        }
    }

    header[3]<-gsub("NA", "", header[3])
    if (flag == 1){
        if (round(funName(c(2, 3, 5)), digits=2) == 3.33){
            if (flagForNorm == 1){
                s1<-"normalized on screen/exp mean (division),"
                
            }
            if (flagForNorm == 2){
                s1<-"normalized on screen/exp mean (substraction),"
            }
        }
        if (funName(c(2, 3, 5)) == 3){
            if (flagForNorm == 1){
                s1<-"normalized on screen/exp median (division),"
            }
            if (flagForNorm == 2){
                s1<-"normalized on screen/exp median (substraction),"
            }
        }
    }
    if (flag == 2){
        if (round(funName(c(2, 3, 5)), digits=2) == 3.33){
            if (flagForNorm == 1){
                s1<-"normalized on plate mean (division),"
            }
            if (flagForNorm == 2){
                s1<-"normalized on plate mean (substraction),"
            }
        }
        if (funName(c(2, 3, 5)) == 3){
            if (flagForNorm == 1){
                s1<-"normalized on plate median (division),"
            }
            if (flagForNorm == 2){
                s1<-"normalized on plate median (substraction),"
            }
        }
    }
    header[3]<-paste(header[3], col4val, s1)
    invisible(list(header, dataset))
}

quantileNormalization<-function(header, dataset, listOfArgs){
    
    flag<-listOfArgs[[1]]
    col4val<-listOfArgs[[2]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    dataset[[get("col4val")]][which(dataset$SpotType == -1)]<-as.double(NA)

    if (flag == 1){
        numOfScreens<-max(dataset$ScreenNb)
        minOfScreens<-min(dataset$ScreenNb)
        maxIntensVals<-0
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
            
                subset<-createSubset(dataset, dataset$ScreenNb, i)
            
                if (length(subset[[get("col4val")]])>maxIntensVals){
                    maxIntensVals<-length(subset[[get("col4val")]])
                }
            }
        }

        matrixForNorm<-matrix(0, maxIntensVals, length(unique(dataset$ScreenNb)))

        NAix<-which(is.na(dataset[[get("col4val")]]))

        screenCounter<-0
        for (j in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == j))>0){
                screenCounter<-screenCounter+1
        
                subset<-createSubset(dataset, dataset$ScreenNb, j)

                if (length(subset[[get("col4val")]])<maxIntensVals){
                    c1<-subset[[get("col4val")]]
                    c2<-rep(as.double(NA_character_), (maxIntensVals-length(subset[[get("col4val")]])))
                    tempSigIntensity<-c(c1, c2)
                    
                }else{
                    tempSigIntensity<-subset[[get("col4val")]]
                }
                matrixForNorm[, screenCounter]<-tempSigIntensity
            }
        }

        normalizedMatrix<-normalizeBetweenArrays(matrixForNorm, method="quantile")
        
        zurueckFuehrung<-as.vector(normalizedMatrix)
        zurueckFuehrung<-zurueckFuehrung[!is.na(zurueckFuehrung)]
        naCounter<-0

        for (i in 1:length(dataset[[get("col4val")]])){
            if (i %in% NAix){
                dataset[[get("col4val")]][i]<-as.double(NA_character_)
            }else{
                naCounter<-naCounter+1
                dataset[[get("col4val")]][i]<-zurueckFuehrung[naCounter]
            }
        }
    }
    
    if (flag == 2){
        numOfScreens<-max(dataset$ScreenNb)
        minOfScreens<-min(dataset$ScreenNb)
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-createSubset(dataset, dataset$ScreenNb, i)
                IXXi<-indexSubset(dataset$ScreenNb, i)
                NAix<-which(is.na(subset[[get("col4val")]]))
                
                numOfPlates<-max(subset$LabtekNb)
                minOfPlates<-min(subset$LabtekNb)
                
                for (j in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == j))>0){
                
                        subsubset<-createSubset(subset, subset$LabtekNb, j)
                        maxIntensVals<-0

                        if (length(subsubset[[get("col4val")]])>maxIntensVals){
                            maxIntensVals<-length(subsubset[[get("col4val")]])
                        }
                    }
                }

                matrixForNorm<-matrix(0, maxIntensVals, length(unique(subset$LabtekNb)))

                for (k in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == k))>0){
                
                        subsubset<-createSubset(subset, subset$LabtekNb, k)

                        if (length(subsubset[[get("col4val")]])<maxIntensVals){
                            c1<-subsubset[[get("col4val")]]
                            c2<-rep(as.double(NA_character_), (maxIntensVals-length(subsubset[[get("col4val")]])))
                            tempSigIntensity<-c(c1, c2)
                            
                        }else{
                            tempSigIntensity<-subsubset[[get("col4val")]]
                        }
                        matrixForNorm[, k]<-tempSigIntensity
                    }
                }

                normalizedMatrix<-normalizeBetweenArrays(matrixForNorm, method="quantile")
        
                zurueckFuehrung<-as.vector(normalizedMatrix)
                zurueckFuehrung<-zurueckFuehrung[!is.na(zurueckFuehrung)]
                naCounter<-0

                for (m in 1:length(subset[[get("col4val")]])){
                    if (m %in% NAix){
                        subset[[get("col4val")]][m]<-as.double(NA_character_)
                    }else{
                        naCounter<-naCounter+1
                        subset[[get("col4val")]][m]<-zurueckFuehrung[naCounter]
                    }
                }
                dataset[[get("col4val")]][IXXi]<-subset[[get("col4val")]]
            }
        }
    }

    header[3]<-gsub("NA", "", header[3])
    if (flag == 1){
        header[3]<-paste(header[3], col4val, "quantile normalization (by screen/exp),")
    }
    if (flag == 2){
        header[3]<-paste(header[3], col4val, "quantile normalization (by plate),")
    }
    invisible(list(header, dataset))
}




BScore <- function(header, dataset, listOfArgs){
    
    col4val<-listOfArgs[[1]]
    excludeControlsFlag<-listOfArgs[[2]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    dataset[[get("col4val")]][which(dataset$SpotType == -1)]<-as.double(NA)
    
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
    
            subset<-createSubset(dataset, dataset$ScreenNb, i)
            IXXi<-indexSubset(dataset$ScreenNb, i)
            numOfPlates<-max(subset$LabtekNb)
            minOfPlates<-min(subset$LabtekNb)

            for (j in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == j))>0){
            
                    subsubset<-createSubset(subset, subset$LabtekNb, j)
                    IXXj<-indexSubset(subset$LabtekNb, j)

                    if (excludeControlsFlag == 1){
                        posControlRows<-subsubset[which(subsubset$SpotType == 1), ]
                        negControlRows<-subsubset[which(subsubset$SpotType == 0), ]
                        subsubset[[get("col4val")]][which(subsubset$SpotType!=2)]<-as.double(NA)
                    }

                    c1<-subsubset[[get("col4val")]]
                    c2<-max(subsubset$RowNb)
                    c3<-max(subsubset$ColNb)
                    
                    
                    ####added on 15.02.2010:
                    medPolishMatrix<-matrix(NA, c2, c3, byrow=T)
                    for (p in 1:length(c1)){
                        medPolishMatrix[subsubset$RowNb[p],subsubset$ColNb[p]]<-c1[p]
                    }
                    ##############

                    fitResults<-medpolish(medPolishMatrix, maxiter=100, trace.iter=FALSE, na.rm=TRUE)
                    BscoreMat<-fitResults$residuals/mad(fitResults$residuals, na.rm=TRUE)
                    
                    if (excludeControlsFlag == 1){
                        
                        for (u in 1:nrow(posControlRows)){
                            c1<-fitResults$overall
                            c2<-fitResults$row[(posControlRows$RowNb[u])]
                            c3<-fitResults$col[(posControlRows$ColNb[u])]
                            normedValsPosContr<-posControlRows[[get("col4val")]][u]-c1-c2-c3
                            normedValsPosContr<-(normedValsPosContr)/mad(fitResults$residuals, na.rm=T)
                            
                            BscoreMat[posControlRows$RowNb[u], posControlRows$ColNb[u]]<-normedValsPosContr
                        }
                        for (u in 1:nrow(negControlRows)){
                            c1<-fitResults$overall
                            c2<-fitResults$row[(negControlRows$RowNb[u])]
                            c3<-fitResults$col[(negControlRows$ColNb[u])]
                            normedValsNegContr<-negControlRows[[get("col4val")]][u]-c1-c2-c3
                            normedValsNegContr<-(normedValsNegContr)/mad(fitResults$residuals, na.rm=T)
                            
                            BscoreMat[negControlRows$RowNb[u], negControlRows$ColNb[u]]<-normedValsNegContr
                        }
                    }
                   
                    #cou<-0
                    #for (k in seq(min(IXXj), max(IXXj), max(subsubset$ColNb))){
                    #    cou<-cou+1
                    #    subset[[get("col4val")]][k:(k+max(subsubset$ColNb)-1)]<-BscoreMat[cou, ]
                    #}
                    
                    ####added on 15.02.2010:
#                    for (k in 1:length(c1)){	#1 bis 384 (1 plate)
                    	#print("BscoreMat[subsubset$RowNb,subsubset$ColNb]")
                    	#print(BscoreMat[subsubset$RowNb,subsubset$ColNb])
                    	
#                        subset[[get("col4val")]][k]<-BscoreMat[subsubset$RowNb,subsubset$ColNb]
#                    }

#Change NR of what was added 15.02.2010; (08.07.2010):
                    for (k in IXXj){                    	
                        subset[[get("col4val")]][k]<-BscoreMat[subset$RowNb[k],subset$ColNb[k]]
                    }                    
#################
                }
            }
            dataset[[get("col4val")]][IXXi]<-subset[[get("col4val")]]
        }
    }

    header[3]<-gsub("NA", "", header[3])
    header[3]<-paste(header[3], col4val, "BScore normalization,")
    invisible(list(header, dataset))
}


ZScore<-function(header, dataset, listOfArgs){
    
    col4val<-listOfArgs[[1]]
    excludeControlsFlag<-listOfArgs[[2]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
    
            subset<-createSubset(dataset, dataset$ScreenNb, i)
            IXXi<-indexSubset(dataset$ScreenNb, i)
            numOfPlates<-max(subset$LabtekNb)
            minOfPlates<-min(subset$LabtekNb)

            for (j in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == j))>0){
            
                    subsubset<-createSubset(subset, subset$LabtekNb, j)
                    IXXj<-indexSubset(subset$LabtekNb, j)
    
                    specialSubset<-subsubset[which(subsubset$SpotType!=-1), ]
                    if (excludeControlsFlag == 1){
                        specialSubset<-specialSubset[which(specialSubset$SpotType == 2), ]
                    }

                    if (sum(!is.na(specialSubset[[get("col4val")]]))!=0){
                        if (mad(specialSubset[[get("col4val")]], na.rm=T)!=0){
                            c1<-subsubset[[get("col4val")]]
                            c2<-median(specialSubset[[get("col4val")]], na.rm=T)
                            c3<-mad(specialSubset[[get("col4val")]], na.rm=T)
                            temp<-(c1-c2)/c3
                        }else{
                            temp<-rep(NA_integer_, nrow(subsubset))
                        }
                    }else{
                        temp<-rep(NA_integer_, nrow(subsubset))
                    }
                    subset[[get("col4val")]][IXXj]<-temp
                }
            }
            dataset[[get("col4val")]][IXXi]<-subset[[get("col4val")]]
        }
    }

    header[3]<-gsub("NA", "", header[3])
    header[3]<-paste(header[3], col4val, "ZScore normalization,")
    invisible(list(header, dataset))
}



ZScorePerScreen<-function(header, dataset, listOfArgs){

    col4val<-listOfArgs[[1]]
    excludeControlsFlag<-listOfArgs[[2]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
    
            subset<-createSubset(dataset, dataset$ScreenNb, i)
            IXXi<-indexSubset(dataset$ScreenNb, i)
        
        
            specialSubset<-subset[which(subset$SpotType!=-1), ]

            if (excludeControlsFlag == 1){
                specialSubset<-specialSubset[which(specialSubset$SpotType == 2), ]
            }

            if (sum(!is.na(specialSubset[[get("col4val")]]))!=0){
                if (mad(specialSubset[[get("col4val")]], na.rm=T)!=0){
                    c1<-subset[[get("col4val")]]
                    c2<-median(specialSubset[[get("col4val")]], na.rm=T)
                    c3<-mad(specialSubset[[get("col4val")]], na.rm=T)
                    temp<-(c1-c2)/c3
                }else{
                    temp<-rep(NA_integer_, nrow(subset))
                }
            }else{
                temp<-rep(NA_integer_, nrow(subset))
            }
        
            dataset[[get("col4val")]][IXXi]<-temp
        }
    }

    header[3]<-gsub("NA", "", header[3])
    header[3]<-paste(header[3], col4val, "ZScore per screen normalization,")

    invisible(list(header, dataset))
}




subtractBackground<-function(header, dataset, listOfArgs){
    
    col4val<-listOfArgs[[1]]
    col4BG<-listOfArgs[[2]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    dataset[[get("col4val")]]<-dataset[[get("col4val")]]-dataset[[get("col4BG")]]

    header[3]<-gsub("NA", "", header[3])
    s1<-" Background substracted from "
    s2<-","
    header[3]<-paste(header[3], s1, col4val, s2, sep="")

    invisible(list(header, dataset))
}

lowessNorm<-function(header, dataset, listOfArgs){

    col4ch1<-listOfArgs[[1]]
    col4ch2<-listOfArgs[[2]]
    if (length(listOfArgs)>2){
    	smSpan<-listOfArgs[[3]]
    }else{
        smSpan<-(2/3)
    }

    dataset<-saveOldIntensityColumns(dataset, col4ch1)
    dataset<-saveOldIntensityColumns(dataset, col4ch2)  

    dataset[[get("col4ch1")]][which(dataset$SpotType == -1)]<-as.double(NA)
    dataset[[get("col4ch2")]][which(dataset$SpotType == -1)]<-as.double(NA)  

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
    
            subset<-createSubset(dataset, dataset$ScreenNb, i)
            IXXi<-indexSubset(dataset$ScreenNb, i)
            numOfPlates<-max(subset$LabtekNb)
            minOfPlates<-min(subset$LabtekNb)

            for (j in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == j))>0){
            
                    subsubset<-createSubset(subset, subset$LabtekNb, j)
                    IXXj<-indexSubset(subset$LabtekNb, j)
    

                    ch1<-subsubset[[get("col4ch1")]]
                    ch2<-subsubset[[get("col4ch2")]]
                    retX<-ch1
                    retY<-ch2
                    use<-(!is.na(ch1))&(!is.na(ch2))
                    if (sum(use)!=0){
                        s<-sort.int(ch1[use], index.return=T)
                        s2<-sort.int(s$ix, index.return=T)
                        z<-lowess(ch1[use], ch2[use], f=smSpan)
                        ch2Lowess<-ch2[use]-z$y[s2$ix]
                        retY[use]<-ch2Lowess
                    }
                    retY[!use]<-NA
                    subsubset[[get("col4ch1")]]<-retX
                    subsubset[[get("col4ch2")]]<-retY

                    subset[[get("col4ch1")]][IXXj]<-subsubset[[get("col4ch1")]]
                    subset[[get("col4ch2")]][IXXj]<-subsubset[[get("col4ch2")]]
                }
            }
            dataset[[get("col4ch1")]][IXXi]<-subset[[get("col4ch1")]]
            dataset[[get("col4ch2")]][IXXi]<-subset[[get("col4ch2")]]
        }
    }

    header[3]<-gsub("NA", "", header[3])
    s1<-"Lowess normalization, ch1 "
    s2<-", ch2 "
    header[3]<-paste(header[3], s1, col4ch1, s2, col4ch2, ",", sep="")

    invisible(list(header, dataset))
}



controlNorm<-function(header,dataset,listOfArgs){

    flag<-listOfArgs[[1]]
    flag2<-listOfArgs[[2]]
    col4val<-listOfArgs[[3]]
    flag3<-listOfArgs[[4]]

    dataset<-saveOldIntensityColumns(dataset, col4val)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)


    if (flag==1){
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb==i))>0){

                subset<-createSubset(dataset, dataset$ScreenNb,i)
                ixx<-indexSubset(dataset$ScreenNb, i)

                if (flag2==0){
                    normSubset<-createSubset(subset, subset$SpotType, 0)
                }
                if (flag2==1){
                    normSubset<-createSubset(subset, subset$SpotType, 1)
                }
                if (typeof(flag2)=="character"){
                    normSubset<-createSubset(subset, subset$GeneName, flag2)
                }
                c1 <- subset[[get("col4val")]]
                if (flag3==1){
                    subset[[get("col4val")]]<-c1/median(normSubset[[get("col4val")]],na.rm=T)
                }else{
                    subset[[get("col4val")]]<-c1-median(normSubset[[get("col4val")]],na.rm=T)
                }
                    dataset[ixx,]<-subset
            }
        }
    }

    if (flag==2){
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb==i))>0){

                subset<-createSubset(dataset, dataset$ScreenNb,i)
                ixx<-indexSubset(dataset$ScreenNb, i)
				
                minOfPlates=min(subset$LabtekNb)
                numOfPlates=max(subset$LabtekNb)

                for (j in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb==j))>0){
				
                        subsubset<-createSubset(subset, subset$LabtekNb, j)
                        iyy<-indexSubset(subset$LabtekNb, j)

                        if (flag2==0){
                            normSubset<-createSubset(subsubset, subsubset$SpotType, 0)
                        }
                        if (flag2==1){
                            normSubset<-createSubset(subsubset, subsubset$SpotType, 1)
                        }
                        if (typeof(flag2)=="character"){
                            normSubset<-createSubset(subset, subset$GeneName, flag2)
                        }
                        c1 <- subsubset[[get("col4val")]]
                        if (flag3==1){
                            subsubset[[get("col4val")]]<-c1/median(normSubset[[get("col4val")]],na.rm=T)
                        }else{
                            subsubset[[get("col4val")]]<-c1-median(normSubset[[get("col4val")]],na.rm=T)
                        }
                        subset[iyy,]<-subsubset
                    }
                }
                dataset[ixx,]<-subset
            }
        }
    }

    header[3]<-gsub("NA", "", header[3])
    s1 <- "Normalization of"
    if (flag==1){
        if (flag2==0){
            if (flag3==1){
                s2 <- "on negative controls (by screen/exp)(div),"
            }else{
                s2 <- "on negative controls (by screen/exp)(substr),"
            }
        }
        if (flag2==1){
            if (flag3==1){
                s2 <- "on positive controls (by screen/exp)(div),"
            }else{
                s2 <- "on positive controls (by screen/exp)(substr),"
            }
        }
        if (typeof(flag2)=="character"){
            if (flag3==1){
                s2 <- paste("on control", flag2, "(by screen/exp)(div),")
            }else{
                s2 <- paste("on control", flag2, "(by screen/exp)(substr),")
            }
        }
    }
    if (flag==2){
        if (flag2==0){
            if (flag3==1){
                s2 <- "on negative controls (by plate)(div),"
            }else{
                s2 <- "on negative controls (by plate)(substr),"
            }
        }
        if (flag2==1){
            if (flag3==1){
                s2 <- "on positive controls (by plate)(div),"
            }else{
                s2 <- "on positive controls (by plate)(substr),"
            }
        }
        if (typeof(flag2)=="character"){
            if (flag3==1){
                s2 <- paste("on control", flag2, "(by plate)(div),")
            }else{
                s2 <- paste("on control", flag2, "(by plate)(substr),")
            }
        }
    }
    header[3] <- paste(header[3], s1, col4val, s2)
    invisible(list(header, dataset))
}



saveDataset<-function(header, data, dataSetFile){

    write.table(header, file=dataSetFile, quote=F, col.names=F, row.names=F)
    suppressWarnings({write.table(data, file=dataSetFile, sep="\t", quote=F, col.names=T, append=T)})

}


saveOldIntensityColumns<-function(dataset, col4val){

    if (length(grep(col4val, colnames(dataset)))>0){

        lengthVal<-length(grep(col4val, colnames(dataset)))
        newColName<-paste(col4val, ".old", lengthVal, sep="")
        dataset<-data.frame(dataset, dataset[[get("col4val")]], stringsAsFactors=F)
        lastColNum<-ncol(dataset)
        colnames(dataset)[lastColNum]<-newColName

    }else{
        newColName<-paste(col4val, ".old", sep="")
        dataset<-data.frame(dataset, dataset[[get("col4val")]], stringsAsFactors=F)
        lastColNum<-ncol(dataset)
        colnames(dataset)[lastColNum]<-newColName
    }
    invisible(dataset)
}




