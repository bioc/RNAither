ZScorePlot<-function(header, dataset, flag, col4plot, col4anno, plotTitle, 
showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    if(flag == 1){
        numOfScreens<-max(dataset$ScreenNb)
        minOfScreens<-min(dataset$ScreenNb)
        vekForPlot<-rep(0, length(dataset[[get("col4plot")]]))
        count<-0
        for (j in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == j))>0){
            
                subset<-createSubset(dataset, dataset$ScreenNb, j)
                minOfPlates<-min(subset$LabtekNb)
                numOfPlates<-max(subset$LabtekNb)
        
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
            
                        count<-count+1
                        vekForPlot[count]<-paste(j, i, sep="_")
                    }
                }
            }
        }
    }else{
        vekForPlot<-dataset[[get("col4anno")]]
    }
    
    currData<-dataset[[get("col4plot")]]
    sd2a<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
    sd2b<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
    sd1a<-mean(currData, na.rm=T)-sd(currData, na.rm=T)
    sd1b<-mean(currData, na.rm=T)+sd(currData, na.rm=T)

    if (sd2a<min(currData, na.rm=T)){
        yymin<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
    }else{
        yymin<-min(currData, na.rm=T)
    }
    if (sd2b>max(currData, na.rm=T)){
        yymax<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
    }else{
        yymax<-max(currData, na.rm=T)
    }
    
    if (showPlot == 1){
        if(interactive()){
            plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, main=plotTitle, 
            cex.main=0.8)
        
            if (flag == 1){
                axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", 
                ylab=col4plot, cex.axis=0.6, las=2)
            }else{
                axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="siRNA", 
                ylab=col4plot, cex.axis=0.6, las=2)
            }
        
            abline(mean(currData, na.rm=T), 0)
            abline(sd2b, 0, col="red")
            abline(sd2a, 0, col="red")        
            abline(sd1b, 0, col="green")
            abline(sd1a, 0, col="green")

            identify(dataset[[get("col4plot")]], labels=dataset[[get("col4anno")]])
        }
    }
    
    headerTemp<-strsplit(header[1], ", ")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(plotName, ".pdf", sep=""))
    
        plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, main=plotTitle, 
        cex.main=0.8)
        
        if (flag == 1){
            axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", 
            ylab=col4plot, cex.axis=0.6, las=2)
        }else{
            axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="siRNA", 
            ylab=col4plot, cex.axis=0.6, las=2)
        }
        
        abline(mean(currData, na.rm=T), 0)
        abline(sd2b, 0, col="red")
        abline(sd2a, 0, col="red")        
        abline(sd1b, 0, col="green")
        abline(sd1a, 0, col="green")
        
        identify(dataset[[get("col4plot")]], labels=dataset[[get("col4anno")]])    
    dev.off()
    
    png(paste(plotName, ".png", sep=""))
    
        plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, main=plotTitle, 
        cex.main=0.8)
        
        if (flag == 1){
            axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", 
            ylab=col4plot, cex.axis=0.6, las=2)
        }else{
            axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="siRNA", 
            ylab=col4plot, cex.axis=0.6, las=2)
        }
        
        abline(mean(currData, na.rm=T), 0)
        abline(sd2b, 0, col="red")
        abline(sd2a, 0, col="red")        
        abline(sd1b, 0, col="green")
        abline(sd1a, 0, col="green")
        
        identify(dataset[[get("col4plot")]], labels=dataset[[get("col4anno")]])    
    dev.off()
    invisible(plotName)
}


Ttest<-function(dataset, listofargs){

    testType<-listofargs[[1]]
    reference<-listofargs[[2]]
    col4val<-listofargs[[3]]
    col4anno<-listofargs[[4]]
    
    temp<-generateReplicateMat(dataset, 1, "Intensities", col4val, col4anno)
    replicaMatrix<-temp[[1]]
    pValVec<-rep(0, nrow(replicaMatrix))
    names(pValVec)<-rownames(replicaMatrix)


    for (i in 1:nrow(replicaMatrix)){

        if (var(replicaMatrix[i, ], na.rm=T)!=0 
        & !is.na(var(replicaMatrix[i, ], na.rm=T))){
        
            if(typeof(reference)!="character"){

                stat<-t.test(replicaMatrix[i, ], alt=testType, mu=reference)
            
            }else{
                rowNum<-which(rownames(replicaMatrix) == reference)
                stat<-t.test(replicaMatrix[i, ], replicaMatrix[rowNum, ], alt=testType)
            }
    
            pValVec[i]<-stat$p.value

        }else{
            pValVec[i]<-as.double(NA_character_)
        }
    }
    dataset<-incorporatepValVec(dataset, pValVec, replicaMatrix, col4anno, 
    paste("pValue.ttest", testType, sep="_"))

    invisible(list(pValVec, dataset, paste("pValue.ttest", testType, sep="_"), 
    "t test"))
}



MannWhitney<-function(dataset, listofargs){

    testType<-listofargs[[1]]
    reference<-listofargs[[2]]
    col4val<-listofargs[[3]]
    col4anno<-listofargs[[4]]

    temp<-generateReplicateMat(dataset, 1, "Intensities", col4val, col4anno)
    replicaMatrix<-temp[[1]]
    pValVec<-rep(0, nrow(replicaMatrix))
    names(pValVec)<-rownames(replicaMatrix)


    for (i in 1:nrow(replicaMatrix)){
        if (var(replicaMatrix[i, ], na.rm=T)!=0 
        & !is.na(var(replicaMatrix[i, ], na.rm=T))){
        
            if(typeof(reference)!="character"){

                stat<-wilcox.test(replicaMatrix[i, ], alt=testType, mu=reference)
            
            }else{
                rowNum<-which(rownames(replicaMatrix) == reference)
                stat<-wilcox.test(replicaMatrix[i, ], replicaMatrix[rowNum, ], alt=testType)
            }
    
            pValVec[i]<-stat$p.value

        }else{
            pValVec[i]<-as.double(NA_character_)
        }
    }
    dataset<-incorporatepValVec(dataset, pValVec, replicaMatrix, col4anno, 
    paste("pValue.mannwhitney", testType, sep="_"))
    
    invisible(list(pValVec, dataset, paste("pValue.mannwhitney", testType, sep="_"), 
    "Mann-Whitney test"))
}




RankProduct<-function(dataset, listofargs){


    permutations<-listofargs[[1]]
    flag<-listofargs[[2]]
    col4val<-listofargs[[3]]
    col4anno<-listofargs[[4]]

    if (flag == 1){
        testType<-"l"
    }else{
        testType<-"g"
    }
    
    temp<-generateReplicateMat(dataset, 1, "Intensities", col4val, col4anno)
    replicaMatrix<-temp[[1]]

    res<-RP(replicaMatrix, rep(1, ncol(replicaMatrix)), num.perm=permutations, 
    na.rm=F, gene.names=rownames(replicaMatrix), plot=T)


    dataset<-incorporatepValVec(dataset, res$pfp[, flag], replicaMatrix, 
    col4anno, paste("pValue.rankproduct", testType, sep="_"))

    pValVec<-res$pfp[, flag]
    names(pValVec)<-rownames(replicaMatrix)

    tt<-paste("pValue.rankproduct", testType, sep="_")
    invisible(list(pValVec, dataset, tt, "Rank product test"))
}



multTestAdjust<-function(pValVec, adjustMethod){

    correctedPvals<-p.adjust(pValVec, method=adjustMethod)    
    invisible(pValVec)
}


savepValVec<-function(pValVec, filename){
    pValVec<-sort(pValVec)
    write.table(pValVec, file=filename, quote=F)
}

incorporatepValVec<-function(dataset, pValVec, replicaMatrix, col4anno, 
colname4pval){

    pValue<-rep(NA_integer_, nrow(dataset))

    dataset<-data.frame(dataset, pValue, stringsAsFactors=F)
    lastColNum<-ncol(dataset)
    colnames(dataset)[lastColNum]<-colname4pval

    for (j in 1:nrow(replicaMatrix)){

        indexes<-which(dataset[[get("col4anno")]] == rownames(replicaMatrix)[j])
        dataset[[get("colname4pval")]][indexes]<-pValVec[j]
    }
    invisible(dataset)
}


hitselectionPval<-function(dataset, pValVec, col4val, col4sel, thresh, 
col4anno, file4hits){

    temp<-generateReplicateMat(dataset, 1, "Intensities", col4val, col4anno)
    replicaMatrix<-temp[[1]]
    
    indexes<-which(pValVec<thresh)

    hitVector<-rep(0, nrow(replicaMatrix))
    names(hitVector)<-rownames(replicaMatrix)
    
    if (length(indexes) == 0){
        thresh<-min(pValVec, na.rm=T)+0.00001
        s1<-"No pvalues under threshold. Threshold was increased to"
        s2<-"for hit selection."
        print(paste(s1, thresh, s2))
        indexes<-which(pValVec<thresh)
    }
        
    hitVector[indexes]<-1
        
    dataset<-incorporatepValVec(dataset, hitVector, replicaMatrix, col4anno, 
    col4sel)
    
    dataset[[get("col4sel")]][which(dataset$SpotType == -1)]<-0
    
    
    toWrite<-replicaMatrix[which(hitVector == 1), ]
    
    
    if (length(which(hitVector == 1)) == 1){
        toWrite <- toWrite[!is.na(toWrite)]
        output<-c(pValVec[indexes], median(toWrite, na.rm=T), toWrite)
        names(output)<-rep(col4val, length(output))
        names(output)[1]<-"pvalue"
        names(output)[2]<-"median"
        
        write.table(t(output), file=file4hits, 
        row.names=names(hitVector[which(hitVector == 1)]), sep="\t", quote=F)
    }else{
    
        medianVecToWrite<-rep(0, nrow(toWrite))
        for (b in 1:nrow(toWrite)){
            medianVecToWrite[b]<-median(toWrite[b, ], na.rm=T)
        }
    
        for (m in 1:ncol(toWrite)){
            if (sum(!is.na(toWrite[, m])) == 0){
                break
            }
        }

        output<-cbind(pValVec[indexes], medianVecToWrite, toWrite[, 1:(m-1)])
        colnames(output)<-rep(col4val, ncol(output))
        colnames(output)[1]<-"pvalue"
        colnames(output)[2]<-"median"
        indexesSort<-order(output[, 1])

        write.table(output[indexesSort, ], file=file4hits, sep="\t", quote=F)
    }

    invisible(list(dataset, hitVector, replicaMatrix, thresh))
}


hitselectionZscore<-function(dataset, col4zscore, col4sel, thresh, flag, 
flag2, col4anno, sumFunc, file4hits){
    
    if (flag == 1){
    
        if (flag2 == 1){
    
            dataset<-orderGeneIDs(dataset, col4zscore)
            hitVector<-rep(0, nrow(dataset))
            
            if (2*abs(thresh)>nrow(dataset)){
                s1<-"Please choose smaller threshold for Z-score hit selection. Only"
                s2<-"values available."
                print(paste(s1, nrow(dataset), s2))
                write.table(colnames(dataset), file=file4hits, sep="\t", quote=F)
                return(list(dataset, hitVector))
            }
    
        
            if (thresh<0){
        
                count<-0
                for (i in nrow(dataset):1){
            
                    if (dataset$SpotType[i]!=-1 & !is.na(dataset[[get("col4zscore")]][i])){
                        count<-count+1
                        hitVector[i]<-1
            
                        if (count == abs(thresh)){
                            break
                        }
                    }
                }
            }
    
            if (thresh>0){
    
                count<-0
                for (i in 1:nrow(dataset)){
                    if (dataset$SpotType[i]!=-1 & !is.na(dataset[[get("col4zscore")]][i])){
                        count<-count+1
                        hitVector[i]<-1
                        
                        if (count == abs(thresh)){
                            break
                        }
                    }
                }
            }
        
            if (thresh == 0){
                hitVector<-rep(1, nrow(dataset))
            }

            dataset<-data.frame(dataset, hitVector, stringsAsFactors=F)
            lastColNum<-ncol(dataset)
            colnames(dataset)[lastColNum]<-col4sel
        
            names(hitVector)<-dataset[[get("col4anno")]]
            
            write.table(dataset[which(hitVector == 1), ], file=file4hits, sep="\t", quote=F)
        }
        
        if (flag2 == (-2)){
            
            indexes<-which(dataset[[get("col4zscore")]]<thresh)
            
            hitVector<-rep(0, nrow(dataset))
            names(hitVector)<-dataset[[get("col4anno")]]
            
            if (length(indexes) == 0){
                thresh<-min(dataset[[get("col4zscore")]], na.rm=T)+0.000001
                s1<-"No Z-scores under threshold. Threshold was increased to"
                s2<-"for hit selection."
                print(paste(s1, thresh, s2))
                indexes<-which(dataset[[get("col4zscore")]]<thresh)
            }
            
            hitVector[indexes]<-1
            hitVector[which(dataset$SpotType == -1)]<-0
            
            dataset<-data.frame(dataset, hitVector, stringsAsFactors=F)
            lastColNum<-ncol(dataset)
            colnames(dataset)[lastColNum]<-col4sel
            
            datasetToWrite<-orderGeneIDs(dataset, col4zscore)
            
            write.table(datasetToWrite[which(datasetToWrite[[get("col4sel")]] == 1), ], 
            file=file4hits, sep="\t", quote=F)
        
        }
        
        if (flag2 == 2){
                    
            indexes<-which(dataset[[get("col4zscore")]]>thresh)
                    
            hitVector<-rep(0, nrow(dataset))
            names(hitVector)<-dataset[[get("col4anno")]]
            
            if (length(indexes) == 0){
                thresh<-max(dataset[[get("col4zscore")]], na.rm=T)-0.000001
                s1<-"No Z-scores over threshold. Threshold was decreased to"
                s2<-"for hit selection."
                print(paste(s1, thresh, s2))
                indexes<-which(dataset[[get("col4zscore")]]>thresh)
            }
            
            hitVector[indexes]<-1
            hitVector[which(dataset$SpotType == 1)]<-0
                    
            dataset<-data.frame(dataset, hitVector, stringsAsFactors=F)
            lastColNum<-ncol(dataset)
            colnames(dataset)[lastColNum]<-col4sel
                    
            datasetToWrite<-orderGeneIDs(dataset, col4zscore)
            
            write.table(datasetToWrite[which(datasetToWrite[[get("col4sel")]] == 1), ], 
            file=file4hits, sep="\t", quote=F)
        }
    }
    
    if (flag == 2){
    
        if (flag2 == 1){
    
            temp<-generateReplicateMat(dataset, 1, "Intensities", col4zscore, col4anno)
            replicaMatrix<-temp[[1]]    
        
            summedZscores<-rep(0, nrow(replicaMatrix))
            for (i in 1:nrow(replicaMatrix)){
            
                summedZscores[i]<-sumFunc(replicaMatrix[i, ], na.rm=T)
            
            }
            setToOrder<-data.frame(rownames(replicaMatrix), summedZscores, 
            stringsAsFactors=F)
            
            colnames(setToOrder)[1]<-"siRNAs"
            colnames(setToOrder)[2]<-"summedZscores"
        
            indexes<-order(setToOrder[, 2])

            hitVector<-rep(0, nrow(replicaMatrix))
            names(hitVector)<-rownames(replicaMatrix)
        
        
            if (2*abs(thresh)>nrow(replicaMatrix)){
                s1<-"Please choose smaller threshold for Z-score hit selection. Only"
                s2<-"values available."
                print(paste(s1, nrow(replicaMatrix), s2))
                write.table(c("summarizedZScore", col4zscore), file=file4hits, sep="\t", 
                quote=F)
                return(list(dataset, hitVector))
            }
        
        
            if (thresh<0){
                indexSubset<-indexes[(length(indexes)-abs(thresh)+1):length(indexes)]
            }
    
            if (thresh>0){    
                indexSubset<-indexes[1:abs(thresh)]            
            }
        
            if (thresh == 0){
                indexSubset<-indexes
            }
        
            hitVector[indexSubset]<-1
            
            dataset<-incorporatepValVec(dataset, hitVector, replicaMatrix, col4anno, 
            col4sel)
            
            dataset[[get("col4sel")]][which(dataset$SpotType == -1)]<-0
        
            toWrite<-replicaMatrix[indexSubset, ]
    
            if (length(which(hitVector == 1)) == 1){
                toWrite <- toWrite[!is.na(toWrite)]
                output<-c(summedZscores[indexSubset], toWrite)
                names(output)<-rep(col4zscore, length(output))
                names(output)[1]<-"summarizedZScore"
                write.table(t(output), file=file4hits, 
                row.names=names(hitVector[which(hitVector == 1)]), sep="\t", quote=F)

            }else{
            
                for (m in 1:ncol(toWrite)){
                    if (sum(!is.na(toWrite[, m])) == 0){
                        break
                    }
                }
                output<-cbind(summedZscores[indexSubset], toWrite[, 1:(m-1)])
                colnames(output)<-rep(col4zscore, ncol(output))
                colnames(output)[1]<-"summarizedZScore"
                write.table(output, file=file4hits, sep="\t", quote=F)
            }
        }
        
        if (flag2 == (-2)){
        
            temp<-generateReplicateMat(dataset, 1, "Intensities", col4zscore, col4anno)
            replicaMatrix<-temp[[1]]
            
            summedZscores<-rep(0, nrow(replicaMatrix))
            for (i in 1:nrow(replicaMatrix)){
                summedZscores[i]<-sumFunc(replicaMatrix[i, ], na.rm=T)            
            }
            
            indexes<-which(summedZscores<thresh)
                
            hitVector<-rep(0, nrow(replicaMatrix))
            names(hitVector)<-rownames(replicaMatrix)
        
        
            if (length(indexes) == 0){
                thresh<-min(summedZscores, na.rm=T)+0.000001
                s1<-"No Z-scores under threshold. Threshold was increased to"
                s2<-"for hit selection."
                print(paste(s1, thresh, s2))
                indexes<-which(summedZscores<thresh)
            }
        
            hitVector[indexes]<-1
                
            dataset<-incorporatepValVec(dataset, hitVector, replicaMatrix, col4anno, 
            col4sel)
            dataset[[get("col4sel")]][which(dataset$SpotType == -1)]<-0
            
            
            toWrite<-replicaMatrix[which(hitVector == 1), ]
            
            if (length(which(hitVector == 1)) == 1){
                toWrite <- toWrite[!is.na(toWrite)]
                output<-c(summedZscores[indexes], toWrite)
                names(output)<-rep(col4zscore, length(output))
                names(output)[1]<-"summarizedZScore"
                
                write.table(t(output), file=file4hits, 
                row.names=names(hitVector[which(hitVector == 1)]), sep="\t", quote=F)
                
            }else{

                for (m in 1:ncol(toWrite)){
                    if (sum(!is.na(toWrite[, m])) == 0){
                        break
                    }
                }
                output<-cbind(summedZscores[indexes], toWrite[, 1:(m-1)])
                colnames(output)<-rep(col4zscore, ncol(output))
                colnames(output)[1]<-"summarizedZScore"
                indexesSort<-order(output[, 1])

                write.table(output[indexesSort, ], file=file4hits, sep="\t", quote=F)
            }

        }

        if (flag2 == 2){

            temp<-generateReplicateMat(dataset, 1, "Intensities", col4zscore, col4anno)
            replicaMatrix<-temp[[1]]
            
            summedZscores<-rep(0, nrow(replicaMatrix))
            for (i in 1:nrow(replicaMatrix)){
                summedZscores[i]<-sumFunc(replicaMatrix[i, ], na.rm=T)            
            }
            
            indexes<-which(summedZscores>thresh)

            hitVector<-rep(0, nrow(replicaMatrix))
            names(hitVector)<-rownames(replicaMatrix)
            
            if (length(indexes) == 0){
                thresh<-max(summedZscores, na.rm=T)-0.00001
                s1<-"No Z-scores over threshold. Threshold was decreased to"
                s2<-"for hit selection."
                print(paste(s1, thresh, s2))
                indexes<-which(summedZscores>thresh)
            }

            hitVector[indexes]<-1
                
            dataset<-incorporatepValVec(dataset, hitVector, replicaMatrix, col4anno, 
            col4sel)
            
            dataset[[get("col4sel")]][which(dataset$SpotType == -1)]<-0
            
            
            toWrite<-replicaMatrix[which(hitVector == 1), ]
            
            if (length(which(hitVector == 1)) == 1){
                toWrite <- toWrite[!is.na(toWrite)]
                output<-c(summedZscores[indexes], toWrite)
                names(output)<-rep(col4zscore, length(output))
                names(output)[1]<-"summarizedZScore"
                
                write.table(t(output), file=file4hits, 
                row.names=names(hitVector[which(hitVector == 1)]), sep="\t", quote=F)
                
            }else{
            
                for (m in 1:ncol(toWrite)){
                    if (sum(!is.na(toWrite[, m])) == 0){
                        break
                    }
                }
                output<-cbind(summedZscores[indexes], toWrite[, 1:(m-1)])
                colnames(output)<-rep(col4zscore, ncol(output))
                colnames(output)[1]<-"summarizedZScore"
                indexesSort<-order(output[, 1])

                write.table(output[indexesSort, ], file=file4hits, sep="\t", quote=F)
            }
        }
    }
    
    invisible(list(dataset, hitVector, thresh))
}



hitselectionZscorePval<-function(dataset, pValVec, col4zscore, col4sel, 
thresh, thresh2, flag2, col4anno, sumFunc, file4hits){

    if (flag2 == (-2)){

        temp<-generateReplicateMat(dataset, 1, "Intensities", col4zscore, col4anno)
        replicaMatrix<-temp[[1]]

        summedZscores<-rep(0, nrow(replicaMatrix))
        for (i in 1:nrow(replicaMatrix)){
            summedZscores[i]<-sumFunc(replicaMatrix[i, ], na.rm=T)            
        }

        indexesZscore<-which(summedZscores<thresh)
        indexesPval<-which(pValVec<thresh2)

        hitVector<-rep(0, nrow(replicaMatrix))
        names(hitVector)<-rownames(replicaMatrix)
        hitVectorZscore<-rep(0, nrow(replicaMatrix))
        hitVectorPval<-rep(0, nrow(replicaMatrix))

        if (length(indexesZscore) == 0){
            thresh<-min(summedZscores, na.rm=T)+0.000001
            s1<-"No Z-scores under threshold. Threshold was increased to"
            s2<-"for combined Z-score/p-value hit selection."
            print(paste(s1, thresh, s2))
            indexesZscore<-which(summedZscores<thresh)

        }
        if (length(indexesPval) == 0){
            thresh2<-min(pValVec, na.rm=T)+0.000001
            s1<-"No p-values under threshold. Threshold was increased to"
            s2<-"for combined Z-score/p-value hit selection."
            print(paste(s1, thresh2, s2))
            indexesPval<-which(pValVec<thresh2)
        }

        hitVectorZscore[indexesZscore]<-1
        hitVectorPval[indexesPval]<-1
        
        hitVector<-hitVectorZscore&hitVectorPval
        hitVector<-as.numeric(hitVector)
        names(hitVector)<-rownames(replicaMatrix)
        indexes<-which(hitVector == 1)
        
        if (sum(hitVector) == 0){
            s1<-"No scores under combined threshold of p-value <"
            s2<-"and ZScore <"
            s3<-". Z-score threshold was increased to"
            s4<-"for combined Z-score/p-value hit selection."
            c1<-min(summedZscores[indexesPval], na.rm=T)+0.00001
            print(paste(s1, thresh2, s2, thresh, s3, c1, s4))
            
            thresh<-min(summedZscores[indexesPval], na.rm=T)+0.00001
            indexesZscore<-which(summedZscores<thresh)
            hitVectorZscore[indexesZscore]<-1
            hitVector<-hitVectorZscore&hitVectorPval
            hitVector<-as.numeric(hitVector)
            names(hitVector)<-rownames(replicaMatrix)
            indexes<-which(hitVector == 1)
        }

        dataset<-incorporatepValVec(dataset, hitVector, replicaMatrix, col4anno, 
        col4sel)
        dataset[[get("col4sel")]][which(dataset$SpotType == -1)]<-0


        toWrite<-replicaMatrix[which(hitVector == 1), ]
        
        if (length(which(hitVector == 1)) == 1){
            toWrite <- toWrite[!is.na(toWrite)]
            output<-c(summedZscores[indexes], pValVec[indexes], toWrite)
            names(output)<-rep(col4zscore, length(output))
            names(output)[1]<-"summarizedZScore"
            names(output)[2]<-"p-value"
            
            write.table(t(output), file=file4hits, 
            row.names=names(hitVector[which(hitVector == 1)]), sep="\t", quote=F)

        }else{
        
            for (m in 1:ncol(toWrite)){
                if (sum(!is.na(toWrite[, m])) == 0){
                    break
                }
            }
            output<-cbind(summedZscores[indexes], pValVec[indexes], toWrite[, 1:(m-1)])
            colnames(output)<-rep(col4zscore, ncol(output))
            colnames(output)[1]<-"summarizedZScore"
            colnames(output)[2]<-"p-value"
            indexesSort<-order(output[, 1])

            write.table(output[indexesSort, ], file=file4hits, sep="\t", quote=F)
        }

    }

    if (flag2 == 2){

        temp<-generateReplicateMat(dataset, 1, "Intensities", col4zscore, col4anno)
        replicaMatrix<-temp[[1]]

        summedZscores<-rep(0, nrow(replicaMatrix))
        for (i in 1:nrow(replicaMatrix)){
            summedZscores[i]<-sumFunc(replicaMatrix[i, ], na.rm=T)            
        }

        indexesZscore<-which(summedZscores>thresh)
        indexesPval<-which(pValVec<thresh2)

        hitVector<-rep(0, nrow(replicaMatrix))
        names(hitVector)<-rownames(replicaMatrix)
        hitVectorZscore<-rep(0, nrow(replicaMatrix))
        hitVectorPval<-rep(0, nrow(replicaMatrix))


        if (length(indexesZscore) == 0){
            thresh<-max(summedZscores, na.rm=T)-0.00001
            s1<-"No Z-scores over threshold. Threshold was decreased to"
            s2<-"for combined Z-score/p-value hit selection."
            print(paste(s1, thresh, s2))
            indexesZscore<-which(summedZscores>thresh)
        }
        if (length(indexesPval) == 0){
            thresh2<-min(pValVec, na.rm=T)-0.00001
            s1<-"No p-values under threshold. Threshold was increased to"
            s2<-"for combined Z-score/p-value hit selection."
            print(paste(s1, thresh2, s2))
            indexesPval<-which(pValVec<thresh2)
        }

        hitVectorZscore[indexesZscore]<-1
        hitVectorPval[indexesPval]<-1

        hitVector<-hitVectorZscore&hitVectorPval
        hitVector<-as.numeric(hitVector)
        names(hitVector)<-rownames(replicaMatrix)
        indexes<-which(hitVector == 1)

        if (sum(hitVector) == 0){
            s1<-"No scores over combined threshold of p-value <"
            s2<-"and ZScore >"
            s3<-". Z-score threshold was decreased to"
            s4<-"for combined Z-score/p-value hit selection."
            c1<-max(summedZscores[indexesPval], na.rm=T)-0.00001
            print(paste(s1, thresh2, s2, thresh, s3, c1, s4))
            
            thresh<-max(summedZscores[indexesPval], na.rm=T)-0.00001
            indexesZscore<-which(summedZscores>thresh)
            hitVectorZscore[indexesZscore]<-1
            hitVector<-hitVectorZscore&hitVectorPval
            hitVector<-as.numeric(hitVector)
            names(hitVector)<-rownames(replicaMatrix)
            indexes<-which(hitVector == 1)
        }

        dataset<-incorporatepValVec(dataset, hitVector, replicaMatrix, col4anno, 
        col4sel)
        
        dataset[[get("col4sel")]][which(dataset$SpotType == -1)]<-0


        toWrite<-replicaMatrix[which(hitVector == 1), ]
        
        if (length(which(hitVector == 1)) == 1){
            toWrite <- toWrite[!is.na(toWrite)]
            output<-c(summedZscores[indexes], pValVec[indexes], toWrite)
            names(output)<-rep(col4zscore, length(output))
            names(output)[1]<-"summarizedZScore"
            names(output)[2]<-"p-value"
            
            write.table(t(output), file=file4hits, 
            row.names=names(hitVector[which(hitVector == 1)]), sep="\t", quote=F)

        }else{

            for (m in 1:ncol(toWrite)){
                if (sum(!is.na(toWrite[, m])) == 0){
                    break
                }
            }
            output<-cbind(summedZscores[indexes], pValVec[indexes], toWrite[, 1:(m-1)])
            colnames(output)<-rep(col4zscore, ncol(output))
            colnames(output)[1]<-"summarizedZScore"
            colnames(output)[2]<-"p-value"
            indexesSort<-order(output[, 1])

            write.table(output[indexesSort, ], file=file4hits, sep="\t", quote=F)
        }
    }    
    invisible(list(dataset, hitVector, thresh, thresh2))
}





spatialDistribHits<-function(header, dataset, plotTitle, col4hits, col4anno, 
showPlot){
    spatialDistrib(header, dataset, plotTitle, col4hits, col4anno, showPlot)
}



volcanoPlot<-function(header, dataset, col4plotx, col4ploty, col4anno, plotTitle, 
sigLevel, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    currDataX<-dataset[[get("col4plotx")]]
    currDataY<-dataset[[get("col4ploty")]]
    
    if (showPlot == 1){
        if(interactive()){
            plot(currDataX, -log10(currDataY), main=plotTitle, xlab=col4plotx, 
            ylab=paste("-log(", col4ploty, ")", sep=""))
        
            if (length(sigLevel) == 1){
                abline(h=-log10(sigLevel), col="green")
            }
            if (length(sigLevel) == 2){
                abline(h=-log10(sigLevel), col="green")
                abline(v=sigLevel[2], col="red")
            }
            if (length(sigLevel) == 3){
                abline(h=-log10(sigLevel), col="green")
                abline(v=sigLevel[2], col="red")
                abline(v=sigLevel[3], col="red")
            }
            identify(currDataX, -log10(currDataY), labels=dataset[[get("col4anno")]])
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    pdf(paste(plotName, ".pdf", sep=""))
        plot(currDataX, -log10(currDataY), main=plotTitle, xlab=col4plotx, 
        ylab=paste("-log(", col4ploty, ")", sep=""))
        
        if (length(sigLevel) == 1){
            abline(h=-log10(sigLevel), col="green")
        }
        if (length(sigLevel) == 2){
            abline(h=-log10(sigLevel), col="green")
            abline(v=sigLevel[2], col="red")
        }
        if (length(sigLevel) == 3){
            abline(h=-log10(sigLevel), col="green")
            abline(v=sigLevel[2], col="red")
            abline(v=sigLevel[3], col="red")
        }
    dev.off()
    png(paste(plotName, ".png", sep=""))
        plot(currDataX, -log10(currDataY), main=plotTitle, xlab=col4plotx, 
        ylab=paste("-log(", col4ploty, ")", sep=""))
        
        if (length(sigLevel) == 1){
            abline(h=-log10(sigLevel), col="green")
        }
        if (length(sigLevel) == 2){
            abline(h=-log10(sigLevel), col="green")
            abline(v=sigLevel[2], col="red")
        }
        if (length(sigLevel) == 3){
            abline(h=-log10(sigLevel), col="green")
            abline(v=sigLevel[2], col="red")
            abline(v=sigLevel[3], col="red")
        }
    dev.off()
    invisible(plotName)
}



vennDiag<-function(header, listOfCols, listOfNames, plotTitle, showPlot){

    if (length(listOfCols) == 2){
    
        mat<-cbind(listOfCols[[1]], listOfCols[[2]])

        obj<-vennCounts(mat)
        
        if (showPlot == 1){
            if(interactive()){           
                vennDiagram(obj, names=c(listOfNames[[1]], listOfNames[[2]]), 
                circle.col=c("blue", "red"), cex=0.8, main=plotTitle, mar=rep(1, 4), 
                cex.main=0.8)
            }
        }

        headerTemp<-strsplit(header[1], ",")
        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
        pdf(paste(plotName, ".pdf", sep=""))
            vennDiagram(obj, names=c(listOfNames[[1]], listOfNames[[2]]), 
            circle.col=c("blue", "red"), cex=0.8, main=plotTitle, 
            mar=rep(1, 4), cex.main=0.8)
        dev.off()
        
        png(paste(plotName, ".png", sep=""))
            vennDiagram(obj, names=c(listOfNames[[1]], listOfNames[[2]]), 
            circle.col=c("blue", "red"), cex=0.8, main=plotTitle, mar=rep(1, 4), 
            cex.main=0.8)
        dev.off()
    
    }else{
    
        mat<-cbind(listOfCols[[1]], listOfCols[[2]], listOfCols[[3]])

        obj<-vennCounts(mat)
        
        if (showPlot == 1){
            if(interactive()){
                vennDiagram(obj, names=c(listOfNames[[1]], listOfNames[[2]], listOfNames[[3]]), 
                circle.col=c("blue", "green", "red"), cex=0.8, main=plotTitle, mar=rep(1, 4), 
                cex.main=0.8)
            }
        }

        headerTemp<-strsplit(header[1], ",")
        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
        pdf(paste(plotName, ".pdf", sep=""))
            vennDiagram(obj, names=c(listOfNames[[1]], listOfNames[[2]], listOfNames[[3]]), 
            circle.col=c("blue", "green", "red"), cex=0.8, main=plotTitle, mar=rep(1, 4), 
            cex.main=0.8)
        dev.off()
        png(paste(plotName, ".png", sep=""))
            vennDiagram(obj, names=c(listOfNames[[1]], listOfNames[[2]], listOfNames[[3]]), 
            circle.col=c("blue", "green", "red"), cex=0.8, main=plotTitle, mar=rep(1, 4), 
            cex.main=0.8)
        dev.off()
    }
    invisible(plotName)
}


compareHits<-function(hitVec1, hitVec2, namesHitVec1, namesHitVec2){

    if (length(hitVec1)!=length(hitVec2)){
        s1<-"Hit vectors have different lengths"
        s2<-"and"
        stop(paste(s1, length(hitVec1), s2, length(hitVec2)))
    }
    if (length(hitVec1)!=length(namesHitVec1)){
        s1<-"Hit vector and annotation vector have different lengths"
        s2<-"and"
        stop(paste(s1, length(hitVec1), s2, length(namesHitVec1)))
    }

    setToOrder1<-data.frame(namesHitVec1, hitVec1, stringsAsFactors=F)
    indexes1<-order(setToOrder1[, 1])
    hitVec1<-hitVec1[indexes1]
    namesHitVec1<-namesHitVec1[indexes1]
    
    setToOrder2<-data.frame(namesHitVec2, hitVec2, stringsAsFactors=F)
    indexes2<-order(setToOrder2[, 1])
    hitVec2<-hitVec2[indexes2]
    namesHitVec2<-namesHitVec2[indexes2]
    
    if (all(namesHitVec2!=namesHitVec1)){
        print(cbind(namesHitVec1, namesHitVec2))
        stop("Annotation vectors are different")
    }
    
    newHitVec<-hitVec1[which(hitVec1 == hitVec2)]
    newNameVec<-namesHitVec1[which(hitVec1 == hitVec2)]

    finalNameVec<-newNameVec[which(newHitVec == 1)]
    finalNameVec<-unique(finalNameVec)
}