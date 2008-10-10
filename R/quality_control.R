discardLabtek<-function(data, screenNr, labtekNr){
    subset<-createSubset(data, data$ScreenNb, screenNr)
    IX<-indexSubset(data$ScreenNb, screenNr)
    subset$SpotType[which(subset$LabtekNb == labtekNr)]<-(-1)
    data[IX, ]<-subset
    invisible(data)
}


discardWells<-function(data, screenNr, labtekNr, vecPositions){
    subset<-createSubset(data, data$ScreenNb, screenNr)
    IX<-indexSubset(data$ScreenNb, screenNr)

    subsubset<-createSubset(subset, subset$LabtekNb, labtekNr)
    IXX<-indexSubset(subset$LabtekNb, labtekNr)

    subsubset$SpotType[vecPositions]<-(-1)

    subset[IXX, ]<-subsubset
    data[IX, ]<-subset
    
    invisible(data)
}


channelPlot <- function(header, dataset, vecOfChannels, flag, plotTitle,
showPlot){

  dataset<-dataset[which(dataset$SpotType!=-1), ]

  if (flag == 0){

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

    combiCounter<-0
    for (j in 1:(length(vecOfChannels)-1)){
        for (k in (j+1):length(vecOfChannels)){
            combiCounter<-combiCounter+1
            ch1<-vecOfChannels[j]
            ch2<-vecOfChannels[k]

            dataset<-dataset[which(!is.na(dataset[[get("ch1")]])), ]
            dataset<-dataset[which(!is.na(dataset[[get("ch2")]])), ]

            tt<-paste(plotTitle, " - Channel ", j, " vs. Channel ", k, sep="")
            x<-dataset[[get("ch1")]]
            y<-dataset[[get("ch2")]]

            if (showPlot == 1){
                if(interactive()){
                    x11()
                    plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                    lines(lowess(x, y), col="red")
                }
            }
            pdf(paste(plotName, "(", combiCounter, ").pdf", sep=""))
                plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                lines(lowess(x, y), col="red")
            dev.off()
            png(paste(plotName, "(", combiCounter, ").png", sep=""), width=300, height=300)
                plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                lines(lowess(x, y), col="red")
            dev.off()
        }
    }
    return(plotName)
  }

  if (flag == 1){

    minOfScreens<-min(dataset$ScreenNb)
    numOfScreens<-max(dataset$ScreenNb)

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
        
            subset<-dataset[which(dataset$ScreenNb == i), ]

            combiCounter<-0
            for (j in 1:(length(vecOfChannels)-1)){
                for (k in (j+1):length(vecOfChannels)){
                    combiCounter<-combiCounter+1
                    ch1<-vecOfChannels[j]
                    ch2<-vecOfChannels[k]

                    subset<-subset[which(!is.na(subset[[get("ch1")]])), ]
                    subset<-subset[which(!is.na(subset[[get("ch2")]])), ]
                    
                    tt<-paste(plotTitle, " - Channel ", j, " vs. Channel ", k, sep="")
		    x<-subset[[get("ch1")]]
                    y<-subset[[get("ch2")]]

                    if (showPlot == 1){
                        if(interactive()){
                            x11()
                            plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                            lines(lowess(x, y), col="red")
                        }
                    }
                    pdf(paste(plotName, "(", combiCounter, ") (Exp. ", i, ").pdf", sep=""))
                        plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                        lines(lowess(x, y), col="red")
                    dev.off()
                    png(paste(plotName, "(", combiCounter, ") (Exp. ", i, ").png", sep=""), 
                    width=300, height=300)
                        plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                        lines(lowess(x, y), col="red")
                    dev.off()
                }
            }
        }
    }
    return(list(plotName, minOfScreens, numOfScreens))
  }

  if (flag == 2){

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
        
            subset<-dataset[which(dataset$ScreenNb == j), ]
            
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)

            headerTemp<-strsplit(header[1], ",")
            plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

            for (i in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == i))>0){
                    
                    subsubset<-subset[which(subset$LabtekNb == i), ]
                    
                    combiCounter<-0
                    for (o in 1:(length(vecOfChannels)-1)){
                        for (p in (o+1):length(vecOfChannels)){
                            combiCounter<-combiCounter+1
                            ch1<-vecOfChannels[o]
                            ch2<-vecOfChannels[p]

                            subsubset<-subsubset[which(!is.na(subsubset[[get("ch1")]])), ]
                            subsubset<-subsubset[which(!is.na(subsubset[[get("ch2")]])), ]
                            
                            tt<-paste(plotTitle, " - Channel ", o, " vs. Channel ", p, sep="")
   		            x<-subsubset[[get("ch1")]]
                            y<-subsubset[[get("ch2")]]

                            if (showPlot == 1){
                                if(interactive()){
                                    x11()
                                    plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                                    lines(lowess(x, y), col="red")
                                }
                            }
                            pdf(paste(plotName, "(", combiCounter, ")_Exp_", j, "_Plate_", i, "_.pdf", 
                            sep=""))
                                plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                                lines(lowess(x, y), col="red")
                            dev.off()
                            png(paste(plotName, "(", combiCounter, ")_Exp_", j, "_Plate_", i, "_.png", 
                            sep=""), width=300, height=300)
                                plot(x, y, main=tt, xlab="Ch1", ylab="Ch2", cex.main=0.7)
                                lines(lowess(x, y), col="red")
                            dev.off()
                        }
                    }
                }
            }
        }
    }
    return(list(plotName, c(minOfScreens, numOfScreens), 
    c(minOfPlates, numOfPlates)))
  }
}


plotBarGeruest <- function(dataset, col4val){

  wells <- nrow(dataset)
  pos<-which(dataset$SpotType == 1)
  neg<-which(dataset$SpotType == 0)
  colors <- rep("black", wells)
  colors[pos] <- "green"
  colors[neg] <- "red"
  data<-dataset[[get("col4val")]]
  med <- median(data, na.rm=T)
  dev <- mad(data, na.rm=T)
  if (is.na(med))
    med <- 0
  if (is.na(dev))
    dev <- 1
  minim <- max(min(as.numeric(data), na.rm=TRUE), med-8*dev)
  maxim <- min(max(as.numeric(data), na.rm=TRUE), med+8*dev)
  if (minim  ==  Inf)
    minim <- -5
  if (maxim  ==  -Inf)
    maxim <- -5

  invisible(list(minim, maxim, data, colors, med, dev))
}

plotBar <- function(header, dataset, col4val, flag, plotTitle, showPlot){

  dataset<-dataset[which(dataset$SpotType!=-1), ]


  if (flag == 0){
      plotData<-plotBarGeruest(dataset, col4val)

    minim<-plotData[[1]]
    maxim<-plotData[[2]]
    data<-plotData[[3]]
    colors<-plotData[[4]]
    med<-plotData[[5]]
    dev<-plotData[[6]]
    
    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    if (showPlot == 1){
        if(interactive()){
            boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
            abline(med, 0, col="blue")
            abline(med+dev, 0, col="green")
            abline(med-dev, 0, col="green")
            abline(med+2*dev, 0, col="red")
            abline(med-2*dev, 0, col="red")
        }
    }
    pdf(paste(plotName, ".pdf", sep=""))
        boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
        abline(med, 0, col="blue")
        abline(med+dev, 0, col="green")
        abline(med-dev, 0, col="green")
        abline(med+2*dev, 0, col="red")
        abline(med-2*dev, 0, col="red")
    dev.off()
    png(paste(plotName, ".png", sep=""), width=300, height=300)
        boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
        abline(med, 0, col="blue")
        abline(med+dev, 0, col="green")
        abline(med-dev, 0, col="green")
        abline(med+2*dev, 0, col="red")
        abline(med-2*dev, 0, col="red")
    dev.off()

    return(plotName)
  }

  if (flag == 1){

    minOfScreens<-min(dataset$ScreenNb)
    numOfScreens<-max(dataset$ScreenNb)

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
        
            subset<-dataset[which(dataset$ScreenNb == i), ]

            plotData<-plotBarGeruest(subset, col4val)
            minim<-plotData[[1]]
            maxim<-plotData[[2]]
            data<-plotData[[3]]
            colors<-plotData[[4]]
            med<-plotData[[5]]
            dev<-plotData[[6]]

            if (showPlot == 1){
                if(interactive()){
                    x11()
                    boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
                    abline(med, 0, col="blue")
                    abline(med+dev, 0, col="green")
                    abline(med-dev, 0, col="green")
                    abline(med+2*dev, 0, col="red")
                    abline(med-2*dev, 0, col="red")
                }
            }

            pdf(paste(plotName, "(Exp. ", i, ").pdf", sep=""))
                boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
                abline(med, 0, col="blue")
                abline(med+dev, 0, col="green")
                abline(med-dev, 0, col="green")
                abline(med+2*dev, 0, col="red")
                abline(med-2*dev, 0, col="red")
            dev.off()
            png(paste(plotName, "(Exp. ", i, ").png", sep=""), width=300, height=300)
                boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
                abline(med, 0, col="blue")
                abline(med+dev, 0, col="green")
                abline(med-dev, 0, col="green")
                abline(med+2*dev, 0, col="red")
                abline(med-2*dev, 0, col="red")
            dev.off()
        }
    }
    return(list(plotName, minOfScreens, numOfScreens))
  }

  if (flag == 2){

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
        
            subset<-dataset[which(dataset$ScreenNb == j), ]
            
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)

            headerTemp<-strsplit(header[1], ",")
            plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

            for (i in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == i))>0){
                    
                    subsubset<-subset[which(subset$LabtekNb == i), ]
                    plotData<-plotBarGeruest(subsubset, col4val)
                    minim<-plotData[[1]]
                    maxim<-plotData[[2]]
                    data<-plotData[[3]]
                    colors<-plotData[[4]]
                    med<-plotData[[5]]
                    dev<-plotData[[6]]

                    if (showPlot == 1){
                        if(interactive()){
                            x11()
                            boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
                            abline(med, 0, col="blue")
                            abline(med+dev, 0, col="green")
                            abline(med-dev, 0, col="green")
                            abline(med+2*dev, 0, col="red")
                            abline(med-2*dev, 0, col="red")
                        }
                    }

                    pdf(paste(plotName, "Exp", j, "Plate", i, ".pdf", sep="_"))
                        boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
                        abline(med, 0, col="blue")
                        abline(med+dev, 0, col="green")
                        abline(med-dev, 0, col="green")
                        abline(med+2*dev, 0, col="red")
                        abline(med-2*dev, 0, col="red")
                    dev.off()
                    png(paste(plotName, "Exp", j, "Plate", i, ".png", sep="_"), width=300, height=300)
                        boxplot(data.frame(t(data)), border=colors, col=colors, ylim=c(minim, maxim))
                        abline(med, 0, col="blue")
                        abline(med+dev, 0, col="green")
                        abline(med-dev, 0, col="green")
                        abline(med+2*dev, 0, col="red")
                        abline(med-2*dev, 0, col="red")
                    dev.off()
                }
            }
        }
    }
    return(list(plotName, c(minOfScreens, numOfScreens), c(minOfPlates, 
    numOfPlates)))
  }
}




ZScorePlotTwo<-function(header, dataset, flag, flag2, col4plot, col4anno, 
plotTitle, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    if (flag == 0){

        yy<-ylabCompute(dataset, col4plot, flag2)
        yymin<-yy[[1]]
        yymax<-yy[[2]]
        currData<-dataset[[get("col4plot")]]
        
        if (flag2 == 0){
            
            sd2a<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
            sd2b<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
            sd1a<-mean(currData, na.rm=T)+sd(currData, na.rm=T)
            sd1b<-mean(currData, na.rm=T)-sd(currData, na.rm=T)
        }
        if (flag2 == 1){
            
            sd2a<-median(currData, na.rm=T)+2*mad(currData, na.rm=T)
            sd2b<-median(currData, na.rm=T)-2*mad(currData, na.rm=T)
            sd1a<-median(currData, na.rm=T)+mad(currData, na.rm=T)
            sd1b<-median(currData, na.rm=T)-mad(currData, na.rm=T)
        }
        
    
        if (showPlot == 1){
            if(interactive()){
                plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                main=plotTitle, cex.main=0.8)
            
                if (flag2 == 0){
                    abline(mean(currData, na.rm=T),0,col="blue")
                }

                if (flag2 == 1){
                    abline(median(currData, na.rm=T),0,col="blue")
                }
            
                abline(sd2a,0,col="red")
                abline(sd2b,0,col="red")
                abline(sd1a,0,col="green")
                abline(sd1b,0,col="green")
                identify(dataset[[get("col4plot")]], labels=dataset[[get("col4anno")]])
            }
        }
        
        headerTemp<-strsplit(header[1], ",")
        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
        pdf(paste(plotName, ".pdf", sep=""))
            plot(dataset[[get("col4plot")]], xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
            main=plotTitle, cex.main=0.8)
            
            if (flag2 == 0){
                abline(mean(currData, na.rm=T),0,col="blue")
            }
            if (flag2 == 1){
                abline(median(currData, na.rm=T),0,col="blue")
            }
            
            abline(sd2a,0,col="red")
            abline(sd2b,0,col="red")
            abline(sd1a,0,col="green")
            abline(sd1b,0,col="green")
        dev.off()
        
        png(paste(plotName, ".png", sep=""), width=300, height=300)
            plot(dataset[[get("col4plot")]], xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
            main=plotTitle, cex.main=0.8)
            
            if (flag2 == 0){
                abline(mean(currData, na.rm=T),0,col="blue")
            }
            if (flag2 == 1){
                abline(median(currData, na.rm=T),0,col="blue")
            }
            
            abline(sd2a,0,col="red")
            abline(sd2b,0,col="red")
            abline(sd1a,0,col="green")
            abline(sd1b,0,col="green")
        dev.off()
        return(plotName)
    }


    if (flag == 1){
    
        minOfScreens<-min(dataset$ScreenNb)
        numOfScreens<-max(dataset$ScreenNb)
    
        headerTemp<-strsplit(header[1], ",")
        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
            
                subset<-dataset[which(dataset$ScreenNb == i), ]
                currData<-subset[[get("col4plot")]]

                yy<-ylabCompute(subset, col4plot, flag2)
                yymin<-yy[[1]]
                yymax<-yy[[2]]
                
                if (flag2 == 0){
                    
                    sd2a<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
		    sd2b<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
		    sd1a<-mean(currData, na.rm=T)+sd(currData, na.rm=T)
                    sd1b<-mean(currData, na.rm=T)-sd(currData, na.rm=T)
                }
                if (flag2 == 1){
                    
                    sd2a<-median(currData, na.rm=T)+2*mad(currData, na.rm=T)
		    sd2b<-median(currData, na.rm=T)-2*mad(currData, na.rm=T)
		    sd1a<-median(currData, na.rm=T)+mad(currData, na.rm=T)
		    sd1b<-median(currData, na.rm=T)-mad(currData, na.rm=T)
		}
                
            
                if (showPlot == 1){
                    if(interactive()){
                        plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                        main=plotTitle, cex.main=0.8)
                    
                        if (flag2 == 0){                        
                            abline(mean(currData, na.rm=T),0,col="blue")
                        }
        
                        if (flag2 == 1){
                            abline(median(currData, na.rm=T),0,col="blue")
                        }
                    
                        abline(sd2a,0,col="red")
                        abline(sd2b,0,col="red")
                        abline(sd1a,0,col="green")
                        abline(sd1b,0,col="green")
                        identify(subset[[get("col4plot")]], labels=subset[[get("col4anno")]])
                    }
                }
                
                headerTemp<-strsplit(header[1], ", ")
                plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
                
                pdf(paste(plotName, "(Exp. ", i, ").pdf", sep=""))
                    plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                    main=plotTitle, cex.main=0.8)
                    
                    if (flag2 == 0){
                        abline(mean(currData, na.rm=T),0,col="blue")
                    }
                    if (flag2 == 1){
                        abline(median(currData, na.rm=T),0,col="blue")
                    }
                    abline(sd2a,0,col="red")
		    abline(sd2b,0,col="red")
		    abline(sd1a,0,col="green")
                    abline(sd1b,0,col="green")
                dev.off()
                
                png(paste(plotName, "(Exp. ", i, ").png", sep=""), width=300, height=300)
                    plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                    main=plotTitle, cex.main=0.8)
                    
                    if (flag2 == 0){
                        abline(mean(currData, na.rm=T),0,col="blue")
                    }
                    if (flag2 == 1){
                        abline(median(currData, na.rm=T),0,col="blue")
                    }
                    abline(sd2a,0,col="red")
		    abline(sd2b,0,col="red")
		    abline(sd1a,0,col="green")
                    abline(sd1b,0,col="green")
                dev.off()
            }
        }
        return(list(plotName, minOfScreens, numOfScreens))
    }
    
    if (flag == 2){
    
        numOfScreens<-max(dataset$ScreenNb)
        minOfScreens<-min(dataset$ScreenNb)
    
        for (j in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == j))>0){
            
                subset<-dataset[which(dataset$ScreenNb == j), ]
                
                minOfPlates<-min(subset$LabtekNb)
                numOfPlates<-max(subset$LabtekNb)
    
                headerTemp<-strsplit(header[1], ",")
                plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                        
                        subsubset<-subset[which(subset$LabtekNb == i), ]
                        currData<-subsubset[[get("col4plot")]]

                        yy<-ylabCompute(subsubset, col4plot, flag2)
                        yymin<-yy[[1]]
                        yymax<-yy[[2]]
                        
                        if (flag2 == 0){
                            sd2a<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
		            sd2b<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
		            sd1a<-mean(currData, na.rm=T)+sd(currData, na.rm=T)
                            sd1b<-mean(currData, na.rm=T)-sd(currData, na.rm=T)
                        }
                        if (flag2 == 1){
                            sd2a<-median(currData, na.rm=T)+2*mad(currData, na.rm=T)
 		            sd2b<-median(currData, na.rm=T)-2*mad(currData, na.rm=T)
 		            sd1a<-median(currData, na.rm=T)+mad(currData, na.rm=T)
 		            sd1b<-median(currData, na.rm=T)-mad(currData, na.rm=T)
 		        }
                    
                        if (showPlot == 1){
                            if(interactive()){
                                plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                                main=plotTitle, cex.main=0.8)
                            
                                if (flag2 == 0){
                                    abline(mean(currData, na.rm=T),0,col="blue")
                                }
                
                                if (flag2 == 1){
                                    abline(median(currData, na.rm=T),0,col="blue")
                                }
                                abline(sd2a,0,col="red")
	        	        abline(sd2b,0,col="red")
        		        abline(sd1a,0,col="green")
                                abline(sd1b,0,col="green")
                                identify(subsubset[[get("col4plot")]], labels=subsubset[[get("col4anno")]])
                            }
                        }
                        
                        headerTemp<-strsplit(header[1], ", ")
                        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
                        
                        pdf(paste(plotName, "Exp", j, "Plate", i, ".pdf", sep="_"))

                            plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                            main=plotTitle, cex.main=0.8)
                            
                            if (flag2 == 0){
                                abline(mean(currData, na.rm=T),0,col="blue")
                            }
                            if (flag2 == 1){
                                abline(median(currData, na.rm=T),0,col="blue")
                            }
                            abline(sd2a,0,col="red")
		            abline(sd2b,0,col="red")
    		            abline(sd1a,0,col="green")
                            abline(sd1b,0,col="green")
                        dev.off()
                        
                        png(paste(plotName, "Exp", j, "Plate", i, ".png", sep="_"), 
                        width=300, height=300)
                        
                            plot(currData, xaxt='n', ylim=c(yymin, yymax), ylab=col4plot, 
                            main=plotTitle, cex.main=0.8)
                            
                            if (flag2 == 0){
                                abline(mean(currData, na.rm=T),0,col="blue")
                            }
                            if (flag2 == 1){
                                abline(median(currData, na.rm=T),0,col="blue")
                            }
                            abline(sd2a,0,col="red")
		            abline(sd2b,0,col="red")
    		            abline(sd1a,0,col="green")
                            abline(sd1b,0,col="green")
                        dev.off()
                    }
                }
            }
        }
        return(list(plotName, c(minOfScreens, numOfScreens), 
        c(minOfPlates, numOfPlates)))
    }
}


ylabCompute<-function(dataset, col4plot, flag2){

    currData<-dataset[[get("col4plot")]]

    if (flag2 == 0){
            
	dev1<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
	dev2<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
	
        if (dev1<min(currData, na.rm=T)){
            yymin<-mean(currData, na.rm=T)-2*sd(currData, na.rm=T)
        }else{
            yymin<-min(currData, na.rm=T)
        }
        
        if (dev2>max(currData, na.rm=T)){
            yymax<-mean(currData, na.rm=T)+2*sd(currData, na.rm=T)
        }else{
            yymax<-max(currData, na.rm=T)
        }
    }

    if (flag2 == 1){
    
        dev1<-median(currData, na.rm=T)-2*mad(currData, na.rm=T)
        dev2<-median(currData, na.rm=T)+2*mad(currData, na.rm=T)

        if (dev1<min(currData, na.rm=T)){
            yymin<-median(currData, na.rm=T)-2*mad(currData, na.rm=T)
        }else{
            yymin<-min(currData, na.rm=T)
        }
        if (dev2>max(currData, na.rm=T)){
            yymax<-median(currData, na.rm=T)+2*mad(currData, na.rm=T)
        }else{
            yymax<-max(currData, na.rm=T)
        }
    }

    invisible(list(yymin, yymax))
}




ZPrime<-function(posControls, negControls){

    ZPrimeFactor<-1-3*((mad(posControls, na.rm=TRUE)+mad(negControls, na.rm=TRUE))/
    (abs(median(posControls, na.rm=TRUE)-median(negControls, na.rm=TRUE))))

}

dynamicRange<-function(dataset, channel){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    totNbPlates4allScreens<-0
    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
            subset<-createSubset(dataset, dataset$ScreenNb, j)
            totNbPlates4allScreens<-totNbPlates4allScreens+length(unique(subset$LabtekNb))
        }
    }

    DR<-rep(0, totNbPlates4allScreens)

    count<-0
    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
            partSet<-dataset[which(dataset$ScreenNb == j), ]
            minOfPlates<-min(partSet$LabtekNb)
            numOfPlates<-max(partSet$LabtekNb)
            
            for (i in minOfPlates:numOfPlates){
                if (length(which(partSet$LabtekNb == i))>0){
                    count<-count+1    
                    subset<-partSet[which(partSet$LabtekNb == i), ]
                    DR[count]<-mean(subset[[get("channel")]][which(subset$SpotType == 0)])/
                    mean(subset[[get("channel")]][which(subset$SpotType == 1)])
                }
            }
        }
    }
    invisible(DR)
}


numCellQualControl<-function(DataSetFile, nbLinesHeader, plotTitle){

    header<-readLines(DataSetFile, nbLinesHeader)
    data<-read.table(DataSetFile, skip=nbLinesHeader, colClasses=c(NA, NA, NA, NA, 
    "factor", NA, NA, NA, NA, NA, NA, NA, NA, NA), stringsAsFactors=FALSE)
    
    threshNumCells<-mean(data$NbCells, na.rm=T)-3*sd(data$NbCells, na.rm=T)
    upperThreshNumCells<-mean(data$NbCells, na.rm=T)+3*sd(data$NbCells, na.rm=T)

    linesUnderThreshNumCells<-data[which(data$NbCells<threshNumCells), ]
    lineNbsUnderThreshNumCells<-which(data$NbCells<threshNumCells)

    linesOverUpperThreshNumCells<-data[which(data$NbCells>upperThreshNumCells), ]
    lineNbsOverUpperThreshNumCells<-which(data$NbCells>upperThreshNumCells)

    
    s1<-"Number of cells per well; Data for wells under lower threshold:"
    print(paste(s1, threshNumCells, sep=" "), quote=F)
    print("", quote=F)
    
    xlabHist<-"Number of cells per well"
    plotTitleHist<-"Distribution of the number of cells"
    
    if(interactive()){
        showThreshHisto<-hist(data$NbCells, breaks=20, xlab=xlabHist, 
        main=plotTitleHist, cex.main=0.8)
        lines(c(threshNumCells, threshNumCells), c(0, max(showThreshHisto$counts)), 
        col="red")
        lines(c(upperThreshNumCells, upperThreshNumCells), 
        c(0, max(showThreshHisto$counts)), col="red")
    }

    if (nrow(linesUnderThreshNumCells) == 0){
        print("No wells under lower threshold.")
    }else{
        print(cbind(lineNbsUnderThreshNumCells, linesUnderThreshNumCells))
    }

    s3<-"\r\nTo set new lower threshold enter value, otherwise \"n\".\r\n\r\n"
    userInput1<-readline(s3)
    
    if (!is.na(as.integer(userInput1))){
    ##set new user-defined threshold:
        linesUnderThreshNumCells<-data[which(data$NbCells<as.double(userInput1)), ]
        lineNbsUnderThreshNumCells<-which(data$NbCells<as.double(userInput1))
        threshNumCells<-as.double(userInput1)

        s4<-"Number of cells per well; Data for wells under new lower threshold:"
        print(paste(s4, threshNumCells, sep=" "), quote=F)
        if (nrow(linesUnderThreshNumCells) == 0){
            print("No wells under lower threshold.")
        }else{
            print(cbind(lineNbsUnderThreshNumCells, linesUnderThreshNumCells))
        }

    }else{
        ##generates warning "NAs introduced by coercion"
        if(is.na(match(userInput1, "n")) | match(userInput1, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    s5<-"\r\nDiscard all wells under lower threshold? [y-n]\r\n\r\n"
    userInput2<-readline(s5)
    if (!is.na(match(userInput2, "y")) & match(userInput2, "y") == 1){
        
        data$SpotType[lineNbsUnderThreshNumCells]<-(-1)

        ##save list of discarded siRNA values:
        headerTemp<-strsplit(header[1], ",")
        s6<-"numCellQualControl_discarded_lower.txt"
        filnam<-paste(headerTemp[[1]][2], s6, sep="_")
        makeHeader<-""
        for (h in 1:length(colnames(data))){
            if (h == 1){
                makeHeader<-colnames(data)[h]
            }else{
                makeHeader<-paste(makeHeader, colnames(data)[h], sep="\t")
            }
        }
        write.table(makeHeader, file=filnam, quote=F, col.names=F, row.names=F)
        write.table(cbind(lineNbsUnderThreshNumCells, linesUnderThreshNumCells), 
        file=filnam, quote=F, sep="\t", col.names=F, row.names=F, append=T)

    }else{
        if(is.na(match(userInput2, "n")) | match(userInput2, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    s7<-"Number of cells per well; Data for wells over upper threshold:"
    print(paste(s7, upperThreshNumCells, sep=" "), quote=F)
    print("", quote=F)
    
    if (nrow(linesOverUpperThreshNumCells) == 0){
        print("No wells over upper threshold.")
    }else{
        print(cbind(lineNbsOverUpperThreshNumCells, linesOverUpperThreshNumCells))
    }

    s8<-"\r\nTo set new upper threshold enter value, otherwise \"n\".\r\n\r\n"
    userInput3<-readline(s8)
    
    if (!is.na(as.integer(userInput3))){
        linesOverUpperThreshNumCells<-data[which(data$NbCells>as.double(userInput3)), ]
        lineNbsOverUpperThreshNumCells<-which(data$NbCells>as.double(userInput3))
        upperThreshNumCells<-as.double(userInput3)


        s9<-"Number of cells per well; Data for wells over new upper threshold:"
        print(paste(s9, upperThreshNumCells, sep=" "), quote=F)
        if (nrow(linesOverUpperThreshNumCells) == 0){
            print("No wells over upper threshold.")
        }else{
            print(cbind(lineNbsOverUpperThreshNumCells, linesOverUpperThreshNumCells))
        }

    }else{
        ##here warning "NAs introduced by coercion"?
        if(is.na(match(userInput3, "n")) | match(userInput3, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    s1<-"\r\nDiscard all wells over upper threshold? [y-n]\r\n\r\n"
    userInput4<-readline(s1)
    if (!is.na(match(userInput4, "y")) & match(userInput4, "y") == 1){

        data$SpotType[lineNbsOverUpperThreshNumCells]<-(-1)

        headerTemp<-strsplit(header[1], ",")
        histoName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")
        
        pdf(histoName)
            if (!is.na(match(userInput2, "y")) & match(userInput2, "y") == 1){
                showThreshHisto<-hist(data$NbCells, breaks=20, xlab="Number of cells per well", 
                main=plotTitle, cex.main=0.8)
                
                lines(c(upperThreshNumCells, upperThreshNumCells), 
                c(0, max(showThreshHisto$counts)), col="red")
                lines(c(threshNumCells, threshNumCells), c(0, max(showThreshHisto$counts)), 
                col="red")
            }
            if (!is.na(match(userInput2, "n")) & match(userInput2, "n") == 1){
                showThreshHisto<-hist(data$NbCells, breaks=20, xlab="Number of cells per well", 
                main=plotTitle, cex.main=0.8)
                lines(c(upperThreshNumCells, upperThreshNumCells), 
                c(0, max(showThreshHisto$counts)), col="red")
            }
        dev.off()
        ##save list of discarded siRNA values:
        filnam<-paste(headerTemp[[1]][2], "numCellQualControl_discarded_higher.txt", sep="_")
        makeHeader<-""
        for (h in 1:length(colnames(data))){
            if (h == 1){
                makeHeader<-colnames(data)[h]
            }else{
                makeHeader<-paste(makeHeader, colnames(data)[h], sep="\t")
            }
        }
        write.table(makeHeader, file=filnam, quote=F, col.names=F, row.names=F)
        write.table(cbind(lineNbsOverUpperThreshNumCells, linesOverUpperThreshNumCells), 
        file=filnam, quote=F, sep="\t", col.names=F, row.names=F, append=T)

    }else{

        headerTemp<-strsplit(header[1], ",")
        histoName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")
        pdf(histoName)
            if (!is.na(match(userInput2, "y")) & match(userInput2, "y") == 1){
                showThreshHisto<-hist(data$NbCells, breaks=20, xlab="Number of cells per well", 
                main=plotTitle, cex.main=0.8)
                lines(c(threshNumCells, threshNumCells), c(0, max(showThreshHisto$counts)), 
                col="red")
            }
            if (!is.na(match(userInput2, "n")) & match(userInput2, "n") == 1){
                showThreshHisto<-hist(data$NbCells, breaks=20, xlab="Number of cells per well", 
                main=plotTitle, cex.main=0.8)
            }
        dev.off()

        if(is.na(match(userInput4, "n")) | match(userInput4, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    write.table(header, file=DataSetFile, quote=F, col.names=F, row.names=F)
    write.table(data, file=DataSetFile, sep="\t", quote=F, col.names=T, append=T)
    invisible(histoName)
}


percCellQualControl<-function(DataSetFile, nbLinesHeader, plotTitle){

    header<-readLines(DataSetFile,nbLinesHeader)
    data<-read.table(DataSetFile, skip=nbLinesHeader, colClasses=c(NA, NA, NA, NA, 
    "factor", NA, NA, NA, NA, NA, NA, NA, NA, NA), stringsAsFactors=FALSE)
	
    threshPcCells<-mean(data$PercCells, na.rm=T)-3*sd(data$PercCells, na.rm=T)
    upperThreshPcCells<-mean(data$PercCells, na.rm=T)+3*sd(data$PercCells, na.rm=T)
    
    linesUnderthreshPcCells<-data[which(data$PercCells<threshPcCells), ]
    lineNbsUnderthreshPcCells<-which(data$PercCells<threshPcCells)
    
    linesOverUpperthreshPcCells<-data[which(data$PercCells>upperThreshPcCells), ]
    lineNbsOverUpperthreshPcCells<-which(data$PercCells>upperThreshPcCells)

    print("", quote=F)
    s1<-"Percentage of cells per well; Data for wells under lower threshold:"
    print(paste(s1, threshPcCells, sep=" "), quote=F)
    print("", quote=F)
    
    xlabs<-"Percentage of cells per well"
    tit<-"Distribution of the percentage of cells"
    if(interactive()){
        showThreshHisto<-hist(data$PercCells, breaks=20, xlab=xlabs, main=tit, 
        cex.main=0.8)
        lines(c(threshPcCells, threshPcCells), c(0, max(showThreshHisto$counts)), 
        col="red")
    }

    if (nrow(linesUnderthreshPcCells) == 0){
        print("No wells under lower threshold.")
    }else{
        print(cbind(lineNbsUnderthreshPcCells, linesUnderthreshPcCells))
    }

    s2<-"\r\n\r\nTo set new lower threshold enter value, otherwise \"n\"\r\n\r\n"
    userInput1<-readline(s2)
    
    if (!is.na(as.integer(userInput1))){
        linesUnderthreshPcCells<-data[which(data$PercCells<as.double(userInput1)), ]
        lineNbsUnderthreshPcCells<-which(data$PercCells<as.double(userInput1))
        threshPcCells<-as.double(userInput1)

        s3<-"Percentage of cells per well; Data for wells under new lower threshold:"
        print(paste(s3, threshPcCells, sep=" "), quote=F)
        if (nrow(linesUnderthreshPcCells) == 0){
            print("No wells under lower threshold.")
        }else{
            print(cbind(lineNbsUnderthreshPcCells, linesUnderthreshPcCells))
        }


    }else{
        if(is.na(match(userInput1, "n")) | match(userInput1, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    s4<-"\r\nDiscard all wells under lower threshold? [y-n]\r\n\r\n"
    userInput2<-readline(s4)
    if (!is.na(match(userInput2, "y")) & match(userInput2, "y") == 1){

        data$SpotType[lineNbsUnderthreshPcCells]<-(-1)

        ##save list of discarded siRNA values:
        headerTemp<-strsplit(header[1], ",")
        s5<-"percCellQualControl_discarded_lower.txt"
        filnam<-paste(headerTemp[[1]][2], s5, sep="_")
        makeHeader<-""
        for (h in 1:length(colnames(data))){
            if (h == 1){
                makeHeader<-colnames(data)[h]
            }else{
                makeHeader<-paste(makeHeader, colnames(data)[h], sep="\t")
            }
        }
        write.table(makeHeader, file=filnam, quote=F, col.names=F, row.names=F)
        write.table(cbind(lineNbsUnderthreshPcCells, linesUnderthreshPcCells), 
        file=filnam, quote=F, sep="\t", col.names=F, row.names=F, append=T)

    }else{
        if(is.na(match(userInput2, "n")) | match(userInput2, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    s6<-"Percentage of cells per well; Data for wells over upper threshold:"
    print(paste(s6, upperThreshPcCells, sep=" "), quote=F)
    print("", quote=F)
    
    if (nrow(linesOverUpperthreshPcCells) == 0){
        print("No wells over upper threshold.")
    }else{
        print(cbind(lineNbsOverUpperthreshPcCells, linesOverUpperthreshPcCells))
    }

    s7<-"\r\nTo set new upper threshold enter value, otherwise \"n\".\r\n\r\n"
    userInput3<-readline(s7)
    
    if (!is.na(as.integer(userInput3))){
        linesOverUpperthreshPcCells<-data[which(data$PercCells>as.double(userInput3)), ]
        lineNbsOverUpperthreshPcCells<-which(data$PercCells>as.double(userInput3))
        upperThreshPcCells<-as.double(userInput3)

        s8<-"Percentage of cells per well; Data for wells over new upper threshold:"
        print(paste(s8, upperThreshPcCells, sep=" "), quote=F)
        if (nrow(linesOverUpperthreshPcCells) == 0){
            print("No wells over upper threshold.")
        }else{
            print(cbind(lineNbsOverUpperthreshPcCells, linesOverUpperthreshPcCells))
        }

    }else{
        ##generates warning "NAs introduced by coercion"
        if(is.na(match(userInput3, "n")) | match(userInput3, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    s9<-"\r\nDiscard all wells over upper threshold? [y-n]\r\n\r\n"
    userInput4<-readline(s9)
    if (!is.na(match(userInput4, "y")) & match(userInput4, "y") == 1){
        data$SpotType[lineNbsOverUpperthreshPcCells]<-(-1)

        headerTemp<-strsplit(header[1], ",")
        histoName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")
        pdf(histoName)
            if (!is.na(match(userInput2, "y")) & match(userInput2, "y") == 1){
            
                showThreshHisto<-hist(data$PercCells, breaks=20, 
                xlab="Percentage of cells per well", main=plotTitle, cex.main=0.8)
                
                lines(c(upperThreshPcCells, upperThreshPcCells), 
                c(0, max(showThreshHisto$counts)), col="red")
                lines(c(threshPcCells, threshPcCells), c(0, max(showThreshHisto$counts)), 
                col="red")
            }
            if (!is.na(match(userInput2, "n")) & match(userInput2, "n") == 1){
                showThreshHisto<-hist(data$PercCells, breaks=20, 
                xlab="Percentage of cells per well", main=plotTitle, cex.main=0.8)
                
                lines(c(upperThreshPcCells, upperThreshPcCells), 
                c(0, max(showThreshHisto$counts)), col="red")
            }
        dev.off()
        ##save list of discarded siRNA values:
        s1<-"percCellQualControl_discarded_higher.txt"
        filnam<-paste(headerTemp[[1]][2], s1, sep="_")
        makeHeader<-""
        for (h in 1:length(colnames(data))){
            if (h == 1){
                makeHeader<-colnames(data)[h]
            }else{
                makeHeader<-paste(makeHeader, colnames(data)[h], sep="\t")
            }
        }
        write.table(makeHeader, file=filnam, quote=F, col.names=F, row.names=F)
        write.table(cbind(lineNbsOverUpperthreshPcCells, linesOverUpperthreshPcCells), 
        file=filnam, quote=F, sep="\t", col.names=F, row.names=F, append=T)

    }else{

        headerTemp<-strsplit(header[1], ",")
        histoName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")
        pdf(histoName)
            if (!is.na(match(userInput2, "y")) & match(userInput2, "y") == 1){
                showThreshHisto<-hist(data$PercCells, breaks=20, 
                xlab="Percentage of cells per well", main=plotTitle, cex.main=0.8)
                
                lines(c(threshPcCells, threshPcCells), c(0, max(showThreshHisto$counts)), 
                col="red")
            }
            if (!is.na(match(userInput2, "n")) & match(userInput2, "n") == 1){
                showThreshHisto<-hist(data$PercCells, breaks=20, 
                xlab="Percentage of cells per well", main=plotTitle, cex.main=0.8)
            }
        dev.off()

        if(is.na(match(userInput4, "n")) | match(userInput4, "n")!=1){
            print("Error. Wrong input.")
        }
    }

    write.table(header, file=DataSetFile, quote=F, col.names=F, row.names=F)
    write.table(data, file=DataSetFile, sep="\t", quote=F, col.names=T, append=T)
    invisible(histoName)
}


ZPRIMEQualControl<-function(header, data, channel, plotTitle, showPlot){

    data<-data[which(data$SpotType!=-1), ]

    numOfScreens<-max(data$ScreenNb)
    minOfScreens<-min(data$ScreenNb)

    totalNbOfPlates4allScreens<-0
    for (j in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == j))>0){
            subset<-createSubset(data, data$ScreenNb, j)
            lenData<-length(unique(subset$LabtekNb))
            totalNbOfPlates4allScreens<-totalNbOfPlates4allScreens+lenData
        }
    }

    ZPRIME<-rep(0, totalNbOfPlates4allScreens)

    count<-0
    for (j in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == j))>0){
        
            smallSet<-data[which(data$ScreenNb == j), ]
            
            minOfPlates<-min(smallSet$LabtekNb)
            numOfPlates<-max(smallSet$LabtekNb)
        
            for (i in minOfPlates:numOfPlates){
                if (length(which(smallSet$LabtekNb == i))>0){
                    count<-count+1
                
                    subset<-smallSet[which(smallSet$LabtekNb == i), ]
                    posContr<-subset[[get("channel")]][which(subset$SpotType == 1)]
                    negContr<-subset[[get("channel")]][which(subset$SpotType == 0)]
                    ZPRIME[count]<-ZPrime(posContr, negContr)
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    filnam<-paste(headerTemp[[1]][2], "Z'Scores.txt", sep="_")

    if (showPlot == 1){
        print("", quote=F)
    }
    vekForPlot<-rep(0, length(ZPRIME))
    count<-0
    Exp<-rep(0, length(ZPRIME))
    Pla<-rep(0, length(ZPRIME))
    for (j in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == j))>0){
    
            smallSet<-data[which(data$ScreenNb == j), ]
                
            minOfPlates<-min(smallSet$LabtekNb)
            numOfPlates<-max(smallSet$LabtekNb)
    
            for (i in minOfPlates:numOfPlates){
                if (length(which(smallSet$LabtekNb == i))>0){
            
                    count<-count+1
                    s2<-"The Z\'Score for Experiment"
                    ausgabe1<-paste(s2, j, "Plate", i, "is", ZPRIME[count], sep=" ")
                    Exp[count]<-j
                    Pla[count]<-i
                    if (showPlot == 1){
                        print(ausgabe1, quote=F)
                    }
                    vekForPlot[count]<-paste(j, i, sep="_")
                }
            }
        }
    }
    if (showPlot == 1){
        print("", quote=F)
    }
    
    ausgabe2<-data.frame(Exp, Pla, ZPRIME)
    colnames(ausgabe2)<-c("Experiment", "Plate", "Z_Prime_Score")
    write.table(ausgabe2, file=filnam, quote=F, row.names=F)
    
    if (showPlot == 1){
        if(interactive()){
            plot(ZPRIME, xaxt='n', main=plotTitle, cex.main=0.8)
            axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", cex.axis=0.6, 
            las=2)
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(plotName, ".pdf", sep=""))
        plot(ZPRIME, xaxt='n', main=plotTitle, cex.main=0.8)
        axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", cex.axis=0.6, 
        las=2)
    dev.off()
    png(paste(plotName, ".png", sep=""), width=300, height=300)
        plot(ZPRIME, xaxt='n', main=plotTitle, cex.main=0.8, cex.lab=0.8, cex.axis=0.8)
        axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", cex.axis=0.6, 
        las=2)
    dev.off()

    ZPrimeTabelle<-read.table(filnam)
    invisible(list(plotName, ZPrimeTabelle))
}


SNRQualControl<-function(dataSetFile, nbLinesHeader, channel, noise, plotTitle){

    header<-readLines(dataSetFile, nbLinesHeader)
    data<-read.table(dataSetFile, skip=nbLinesHeader, colClasses=c(NA, NA, NA, NA, 
    "factor", NA, NA, NA, NA, NA, NA, NA, NA, NA), stringsAsFactors=FALSE)

    data<-data[which(data$SpotType!=-1), ]

    SNR<-data[[get("channel")]]/data[[get("noise")]]

    if(interactive()){
        hist(SNR, breaks=20, xlab="SNR per well")
    }

    numOfScreens<-max(data$ScreenNb)
    minOfScreens<-min(data$ScreenNb)

    c1<-ceiling(sqrt(length(unique(data$ScreenNb))))
    c2<-ceiling(length(unique(data$ScreenNb))/c1)
    
    if(interactive()){
        x11()
        par(mfrow=c(c1, c2))
    }
    for (i in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == i))>0){
            subsetSnr<-SNR[which(data$ScreenNb == i)]
            
            if(interactive()){
                hist(subsetSnr, breaks=20, xlab=paste("SNR per well", sep=" "), 
                main=paste(plotTitle, "for Exp.", i, sep=" "))
            }
        }
    }

    for (j in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == j))>0){
        
            subset<-data[which(data$ScreenNb == j), ]
            
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
            
            if(interactive()){
                x11()
                c3<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                c4<-ceiling(length(unique(subset$LabtekNb))/c3)
                par(mfrow=c(c3, c4))
            }
            
            for (i in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == i))>0){
                
                    subsetSnr<-SNR[which(subset$LabtekNb == i)]
                    
                    if(interactive()){
                        hist(subsetSnr, breaks=20, xlab=paste("SNR per well", sep=" "), 
                        main=paste(plotTitle, "for Exp.", j, "Plate", i, sep=" "))
                    }
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    histoName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")
    pdf(histoName)
        hist(SNR, breaks=20, xlab="SNR per well")
    dev.off()

    histoName<-paste(headerTemp[[1]][2], plotTitle, "PerExp.pdf", sep="_")
    pdf(histoName)
        par(mfrow=c(c1, c2))
        for (i in minOfScreens:numOfScreens){
            if (length(which(data$ScreenNb == i))>0){
            
                subsetSnr<-SNR[which(data$ScreenNb == i)]
                
                hist(subsetSnr, breaks=20, xlab=paste("SNR per well", sep=" "), 
                main=paste(plotTitle, "for Exp.", i, sep=" "), cex.main=0.8)
            }
        }
    dev.off()

    for (j in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == j))>0){
            subset<-data[which(data$ScreenNb == j), ]
            
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
            
            histoName<-paste(headerTemp[[1]][2], plotTitle, "PerLabtek_Exp", j, ".pdf", sep="_")
            
            pdf(histoName)
                c3<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                c4<-ceiling(length(unique(subset$LabtekNb))/c3)
                par(mfrow=c(c3, c4))
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                
                        subsetSnr<-SNR[which(subset$LabtekNb == i)]
                        
                        hist(subsetSnr, breaks=20, xlab=paste("SNR per well", sep=" "), 
                        main=paste(plotTitle, "for Exp.", j, "Plate", i, sep=" "), cex.main=0.8)
                    }
                }
            dev.off()
        }
    }
}

DRQualControl<-function(header, data, nbLinesHeader, channel, plotTitle){

    data<-data[which(data$SpotType!=-1), ]

    numOfScreens<-max(data$ScreenNb)
    minOfScreens<-min(data$ScreenNb)

    DR<-dynamicRange(data, channel)

    headerTemp<-strsplit(header[1], ",")
    filnam<-paste(headerTemp[[1]][2], "DR.txt", sep="_")
    write.table(paste("Experiment", "Plate", "DR", sep="\t"), file=filnam, quote=F, 
    col.names=F, row.names=F)

    print("", quote=F)
    vekForPlot<-rep(0, length(DR))
    count<-0
    for (j in minOfScreens:numOfScreens){
        if (length(which(data$ScreenNb == j))>0){

            subset<-createSubset(data, data$ScreenNb, j)
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
        
            for (i in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == i))>0){
            
                    count<-count+1
                    s1<-"The dynamic range for Experiment"
                    ausgabe1<-paste(s1, j, "Plate", i, "is", DR[count], sep=" ")
                    ausgabe2<-paste(j, i, DR[count], sep="\t")
                    print(ausgabe1, quote=F)
                    vekForPlot[count]<-paste(j, i, sep="_")
                    write.table(ausgabe2, file=filnam, quote=F, col.names=F, row.names=F, append=T)
                }
            }
        }
    }
    print("", quote=F)
    
    if(interactive()){
        plot(DR, xaxt='n', main=plotTitle, cex.main=0.8)
        axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", 
        cex.axis=0.6, las=2)
    }
    
    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")
    pdf(plotName)
        plot(DR, xaxt='n', main=plotTitle, cex.main=0.8)
        axis(1, 1:length(vekForPlot), labels=vekForPlot, xlab="Exp_Plate", 
        cex.axis=0.6, las=2)
    dev.off()
}




plotControlHisto<-function(header, dataset, channel, plotTitle, showPlot){

    All<-dataset[[get("channel")]][which(dataset$SpotType!=-1)]

    Neg<-dataset[[get("channel")]][which(dataset$SpotType == 0)]
    Pos<-dataset[[get("channel")]][which(dataset$SpotType == 1)]
    Other<-dataset[[get("channel")]][which(dataset$SpotType == 2)]

    if (sum(!is.na(All)) == 0){
        stop("Cannot plot histogram (only NAs in dataset)")
    }

    if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
        a<-min(All, na.rm=T)
        b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
        d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
        computeHistoBreaks=seq(a, b, d)
    }else{
        a<-min(All, na.rm=T)
        b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
        d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
        computeHistoBreaks=seq(a, b, d)
    }

    all<-hist(Other, plot=F, breaks=computeHistoBreaks)
    pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
    neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

    if (showPlot == 1){
        if(interactive()){
            plot(all, main=plotTitle, xlab=channel, cex.main=0.7)
            par(fg="green")
            lines(pos, angle=45, density=20)
            par(fg="red")
            lines(neg, angle=-45, density=20)
            par(fg="black")
            legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
            angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.8)
        }
    }

    headerTemp<-strsplit(header[1], ",")
    histoName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(histoName, ".pdf", sep=""))
        plot(all, main=plotTitle, xlab=channel, cex.main=0.7)
        par(fg="green")
        lines(pos, angle=45, density=20)
        par(fg="red")
        lines(neg, angle=-45, density=20)
        par(fg="black")
        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.8)
    dev.off()
    png(paste(histoName, ".png", sep=""), width=300, height=300)
        plot(all, main=plotTitle, xlab=channel, cex.main=0.8)
        par(fg="green")
        lines(pos, angle=45, density=20)
        par(fg="red")
        lines(neg, angle=-45, density=20)
        par(fg="black")
        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.7)
    dev.off()
    invisible(histoName)
}


plotControlHistoPerplate<-function(header, dataset, channel, plotTitle, 
plotDesign, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    if (showPlot == 1){
        if(interactive()){
            for (j in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == j))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == j), ]
            
                    minOfPlates<-min(subset$LabtekNb)
                    numOfPlates<-max(subset$LabtekNb)
            
                    x11()
                    if (plotDesign == 1){
                        c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                        c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                        par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
                    }

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                
                            subsubset<-subset[which(subset$LabtekNb == i), ]
                            All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]
                            Neg<-subsubset[[get("channel")]][which(subsubset$SpotType == 0)]
                            Pos<-subsubset[[get("channel")]][which(subsubset$SpotType == 1)]
                            Other<-subsubset[[get("channel")]][which(subsubset$SpotType == 2)]

                            if (sum(!is.na(All)) == 0){
                                s1<-"Cannot plot histogram for Exp."
                                print(paste(s1, j, "plate", i, "(only NAs in dataset)", sep=" "))
                            }else{

                                if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                    d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                    computeHistoBreaks=seq(a, b, d)
                                }else{
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+abs(a+0.5)+((max(All, na.rm=T)-a+1)/20)
                                    d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                    computeHistoBreaks=seq(a, b, d)
                                }

                                all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                                pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                                neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                                plot(all, main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, 
                                cex.main=0.6)
                                par(fg="green")
                                lines(pos, angle=45, density=20)
                                par(fg="red")
                                lines(neg, angle=-45, density=20)
                                par(fg="black")
                                legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                                angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.8)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T, cex=0.8)
                            }
                        }
                    }
                }
            }
        }
    }

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){

            subset<-dataset[which(dataset$ScreenNb == j), ]
            
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)

            headerTemp<-strsplit(header[1], ",")
            histoName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

            if (plotDesign == 1){
            
                pdf(paste(histoName, "Exp", j, "PerPlate", ".pdf", sep="_"))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))

    
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
         
                            All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]    
                            Neg<-subsubset[[get("channel")]][which(subsubset$SpotType == 0)]
                            Pos<-subsubset[[get("channel")]][which(subsubset$SpotType == 1)]
                            Other<-subsubset[[get("channel")]][which(subsubset$SpotType == 2)]

                            if (sum(!is.na(All))!=0){    
                                if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                    d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                    computeHistoBreaks=seq(a, b, d)
                                }else{
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+abs(a+0.5)+((max(All, na.rm=T)-a+1)/20)
                                    d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                    computeHistoBreaks=seq(a, b, d)
                                }

                                all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                                pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                                neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                                plot(all, main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, 
                                cex.main=0.7)
                                par(fg="green")
                                lines(pos, angle=45, density=20)
                                par(fg="red")
                                lines(neg, angle=-45, density=20)
                                par(fg="black")
                                legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                                angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            }
                        }
                    }
                dev.off()
                
                png(paste(histoName, "Exp", j, "PerPlate", ".png", sep="_"))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))

    
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
                            All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]    
                            Neg<-subsubset[[get("channel")]][which(subsubset$SpotType == 0)]
                            Pos<-subsubset[[get("channel")]][which(subsubset$SpotType == 1)]
                            Other<-subsubset[[get("channel")]][which(subsubset$SpotType == 2)]

                            if (sum(!is.na(All))!=0){    
                                if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                    d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                    computeHistoBreaks=seq(a, b, d)
                                }else{
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+abs(a+0.5)+((max(All, na.rm=T)-a+1)/20)
                                    d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                    computeHistoBreaks=seq(a, b, c)
                                }

                                all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                                pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                                neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                                plot(all, main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, 
                                cex.main=0.7)
                                par(fg="green")
                                lines(pos, angle=45, density=20)
                                par(fg="red")
                                lines(neg, angle=-45, density=20)
                                par(fg="black")
                                legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                                angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            }
                        }
                    }
                dev.off()
                
            }else{
                    
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                    
                        subsubset<-subset[which(subset$LabtekNb == i), ]         
                        All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]
                        Neg<-subsubset[[get("channel")]][which(subsubset$SpotType == 0)]
                        Pos<-subsubset[[get("channel")]][which(subsubset$SpotType == 1)]
                        Other<-subsubset[[get("channel")]][which(subsubset$SpotType == 2)]
                        
                        if (sum(!is.na(All))!=0){
    
                            if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                a<-min(All, na.rm=T)
                                b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                computeHistoBreaks=seq(a, b, d)
                            }else{
                                a<-min(All, na.rm=T)
                                b<-max(All, na.rm=T)+abs(a+0.5)+((max(All, na.rm=T)-a+1)/20)
                                d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                computeHistoBreaks=seq(a, b, d)
                            }

                            all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                            pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                            neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                            pdf(paste(histoName, "Exp", j, "PerPlate", i, ".pdf", sep="_"))

                                plot(all, main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, 
                                cex.main=0.7)
                                par(fg="green")
                                lines(pos, angle=45, density=20)
                                par(fg="red")
                                lines(neg, angle=-45, density=20)
                                par(fg="black")
                                legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                                angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            dev.off()
                            
                            png(paste(histoName, "Exp", j, "PerPlate", i, ".png", sep="_"), 
                            width=300, height=300)

                                plot(all, main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, 
                                cex.main=0.5)
                                par(fg="green")
                                lines(pos, angle=45, density=20)
                                par(fg="red")
                                lines(neg, angle=-45, density=20)
                                par(fg="black")
                                legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                                angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            dev.off()
                        }
                    }
                }
            }
        }
    }
    invisible(list(histoName, c(minOfScreens, numOfScreens), 
    c(minOfPlates, numOfPlates)))
}

plotControlHistoPerscreen<-function(header, dataset, channel, plotTitle, 
plotDesign, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    minOfScreens<-min(dataset$ScreenNb)
    numOfScreens<-max(dataset$ScreenNb)
    
    if (showPlot == 1){
        if(interactive()){
    
            if (plotDesign == 1){
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
                c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
                par(mfrow=c(c1, c2))
            }

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
    
                    subset<-dataset[which(dataset$ScreenNb == i), ]     
                    All<-subset[[get("channel")]][which(subset$SpotType!=-1)]
                    Neg<-subset[[get("channel")]][which(subset$SpotType == 0)]
                    Pos<-subset[[get("channel")]][which(subset$SpotType == 1)]
                    Other<-subset[[get("channel")]][which(subset$SpotType == 2)]

                    if (sum(!is.na(All)) == 0){
                        s1<-"Cannot plot histogram for Exp."
                        print(paste(s1, i, "(only NAs in dataset)", sep=" "))
                    }else{

                        if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                            d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                            computeHistoBreaks=seq(a, b, d)
                        }else{
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                            d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                            computeHistoBreaks=seq(a, b, d)
                        }
    
                        all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                        pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                        neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                        plot(all, main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, 
                        cex.main=0.7)
                        par(fg="green")
                        lines(pos, angle=45, density=20)
                        par(fg="red")
                        lines(neg, angle=-45, density=20)
                        par(fg="black")
                        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.8)
                    }
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    histoName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    if (plotDesign == 1){
        pdf(paste(histoName, ".pdf", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]
 
                    All<-subset[[get("channel")]][which(subset$SpotType!=-1)]
                    Neg<-subset[[get("channel")]][which(subset$SpotType == 0)]
                    Pos<-subset[[get("channel")]][which(subset$SpotType == 1)]
                    Other<-subset[[get("channel")]][which(subset$SpotType == 2)]

                    if (sum(!is.na(All))!=0){

                        if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                            d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                            computeHistoBreaks=seq(a, b, d)
                        }else{
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                            d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                            computeHistoBreaks=seq(a, b, d)
                        }
    
                        all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                        pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                        neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                        plot(all, main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, 
                        cex.main=0.5)
                        par(fg="green")
                        lines(pos, angle=45, density=20)
                        par(fg="red")
                        lines(neg, angle=-45, density=20)
                        par(fg="black")
                        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                    }
                }
            }
        dev.off()
        
        png(paste(histoName, ".png", sep=""))

            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
	    c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
            
            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ] 
                    All<-subset[[get("channel")]][which(subset$SpotType!=-1)]
                    Neg<-subset[[get("channel")]][which(subset$SpotType == 0)]
                    Pos<-subset[[get("channel")]][which(subset$SpotType == 1)]
                    Other<-subset[[get("channel")]][which(subset$SpotType == 2)]

                    if (sum(!is.na(All))!=0){

                        if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                            d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                            computeHistoBreaks=seq(a, b, d)
                        }else{
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                            d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                            computeHistoBreaks=seq(a, b, d)
                        }
    
                        all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                        pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                        neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                        plot(all, main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, 
                        cex.main=0.5)
                        par(fg="green")
                        lines(pos, angle=45, density=20)
                        par(fg="red")
                        lines(neg, angle=-45, density=20)
                        par(fg="black")
                        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                    }
                }
            }
        dev.off()

    }else{

        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-dataset[which(dataset$ScreenNb == i), ]

                All<-subset[[get("channel")]][which(subset$SpotType!=-1)]
                Neg<-subset[[get("channel")]][which(subset$SpotType == 0)]
                Pos<-subset[[get("channel")]][which(subset$SpotType == 1)]
                Other<-subset[[get("channel")]][which(subset$SpotType == 2)]

                if (sum(!is.na(All))!=0){

                    if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                        a<-min(All, na.rm=T)
                        b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                        d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                        computeHistoBreaks=seq(a, b, d)
                    }else{
                        a<-min(All, na.rm=T)
                        b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                        d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                        computeHistoBreaks=seq(a, b, d)
                    }
    
                    all<-hist(Other, plot=F, breaks=computeHistoBreaks)
                    pos<-hist(Pos, plot=F, breaks=computeHistoBreaks)
                    neg<-hist(Neg, plot=F, breaks=computeHistoBreaks)

                    pdf(paste(histoName, "(Exp. ", i, ").pdf", sep=""))
                        plot(all, main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, 
                        cex.main=0.5)
                        par(fg="green")
                        lines(pos, angle=45, density=20)
                        par(fg="red")
                        lines(neg, angle=-45, density=20)
                        par(fg="black")
                        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                    dev.off()
                    
                    png(paste(histoName, "(Exp. ", i, ").png", sep=""), width=300, height=300)
                        plot(all, main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, 
                        cex.main=0.5)
                        par(fg="green")
                        lines(pos, angle=45, density=20)
                        par(fg="red")
                        lines(neg, angle=-45, density=20)
                        par(fg="black")
                        legend("topleft", c("Data", "Positive Controls", "Negative Controls"), 
                        angle=c(0, 45, -45), density=c(0, 20, 20), bty="n", cex=0.5)
                    dev.off()
                    
                }
            }
        }
    }
    invisible(list(histoName, minOfScreens, numOfScreens))
}






plotHisto<-function(header, dataset, channel, plotTitle, showPlot){

    All<-dataset[[get("channel")]][which(dataset$SpotType!=-1)]

    if (sum(!is.na(All)) == 0){
        stop("Cannot plot histogram (only NAs in dataset)")
    }

    if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
        a<-min(All, na.rm=T)
        b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
        d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
        computeHistoBreaks=seq(a, b, d)
    }else{
        a<-min(All, na.rm=T)
        b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
        d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
        computeHistoBreaks=seq(a, b, d)
    }

    if (showPlot == 1){
        if(interactive()){
            hist(All, breaks=computeHistoBreaks, main=plotTitle, xlab=channel, cex.main=0.8)
        }
    }
    
    headerTemp<-strsplit(header[1], ",")
    histoName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(histoName, ".pdf", sep=""))
        hist(All, breaks=computeHistoBreaks, main=plotTitle, xlab=channel, cex.main=0.8)
    dev.off()
    
    png(paste(histoName, ".png", sep=""), width=300, height=300)
        hist(All, breaks=computeHistoBreaks, main=plotTitle, xlab=channel, cex.main=0.8)
    dev.off()
    invisible(histoName)
}

plotHistoPerplate<-function(header, dataset, channel, plotTitle, plotDesign, 
showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    if (showPlot == 1){
        if(interactive()){
            for (j in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == j))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == j), ]
                
                    minOfPlates<-min(subset$LabtekNb)
                    numOfPlates<-max(subset$LabtekNb)
            
                    if (plotDesign == 1){
                        x11()
                        c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                        c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                        par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
                    }

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                
                            subsubset<-subset[which(subset$LabtekNb == i), ]     
                            All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]

                            if (sum(!is.na(All)) == 0){
                                s1<-"Cannot plot histogram for Exp."
                                print(paste(s1, j, "plate", i, "(only NAs in dataset)", sep=" "))
                            }else{

                                if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                    d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                    computeHistoBreaks=seq(a, b, d)
                                }else{
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                                    d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                    computeHistoBreaks=seq(a, b, d)
                                }
                        
                                hist(All, breaks=computeHistoBreaks, 
                                main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, cex.main=0.6)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T, cex=0.8)
                            }
                        }
                    }
                }
            }
        }
    }

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){

            subset<-dataset[which(dataset$ScreenNb == j), ]
            
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)

            headerTemp<-strsplit(header[1], ",")
            histoName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

            if (plotDesign == 1){

                pdf(paste(histoName, "Exp", j, "PerPlate", ".pdf", sep="_"))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
    
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
                            All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]
    
                            if (sum(!is.na(All))!=0){
    
                                if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                    d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                    computeHistoBreaks=seq(a, b, d)
                                }else{
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                                    d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                    computeHistoBreaks=seq(a, b, d)
                                }

                                hist(All, breaks=computeHistoBreaks, 
                                main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, cex.main=0.7)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            }
                        }
                    }
                dev.off()
                
                png(paste(histoName, "Exp", j, "PerPlate", ".png", sep="_"))
                    
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
    
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
                            All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]
    
                            if (sum(!is.na(All))!=0){
    
                                if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                    d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                    computeHistoBreaks=seq(a, b, d)
                                }else{
                                    a<-min(All, na.rm=T)
                                    b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                                    d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                    computeHistoBreaks=seq(a, b, d)
                                }

                                hist(All, breaks=computeHistoBreaks, 
                                main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, cex.main=0.7)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            }
                        }
                    }
                dev.off()
                
            }else{
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                    
                        subsubset<-subset[which(subset$LabtekNb == i), ]         
                        All<-subsubset[[get("channel")]][which(subsubset$SpotType!=-1)]
                        
                        if (sum(!is.na(All))!=0){
    
                            if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                                a<-min(All, na.rm=T)
                                b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                                d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                                computeHistoBreaks=seq(a, b, d)
                            }else{
                                a<-min(All, na.rm=T)
                                b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                                d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                                computeHistoBreaks=seq(a, b, d)
                            }

                            pdf(paste(histoName, "Exp", j, "PerPlate", i, ".pdf", sep=""))
                                hist(All, breaks=computeHistoBreaks, 
                                main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, cex.main=0.7)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            dev.off()
                            
                            png(paste(histoName, "Exp", j, "PerPlate", i, ".png", sep=""), 
                            width=300, height=300)
                                hist(All, breaks=computeHistoBreaks, 
                                main=paste(plotTitle, "for plate", i, sep=" "), xlab=channel, cex.main=0.7)
                                mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                            dev.off()
                        }
                    }
                }            
            }
        }
    }
    invisible(list(histoName, c(minOfScreens, numOfScreens), 
    c(minOfPlates, numOfPlates)))
}


plotHistoPerscreen<-function(header, dataset, channel, plotTitle, plotDesign, 
showPlot){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    minOfScreens<-min(dataset$ScreenNb)
    numOfScreens<-max(dataset$ScreenNb)
    
    if (showPlot == 1){
        if(interactive()){
    
            if (plotDesign == 1){
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
                c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
                par(mfrow=c(c1, c2))
            }

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
    
                    subset<-dataset[which(dataset$ScreenNb == i), ]
                    All<-subset[[get("channel")]][which(subset$SpotType!=-1)]

                    if (sum(!is.na(All)) == 0){
                        s1<-"Cannot plot histogram for Exp."
                        print(paste(s1, i, "(only NAs in dataset)", sep=" "))
                    }else{
                        if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                            d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                            computeHistoBreaks=seq(a, b, d)
                        }else{
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                            d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                            computeHistoBreaks=seq(a, b, d)
                        }
    
                        hist(All, breaks=computeHistoBreaks, 
                        main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, cex.main=0.8)
                    }
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    histoName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    if (plotDesign == 1){
        pdf(paste(histoName, ".pdf", sep=""))

            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]
 
                    All<-subset[[get("channel")]][which(subset$SpotType!=-1)]

                    if (sum(!is.na(All))!=0){

                        if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                            d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                            computeHistoBreaks=seq(a, b, d)    
                        }else{
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                            d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                            computeHistoBreaks=seq(a, b, d)
                        }
    
                        hist(All, breaks=computeHistoBreaks, 
                        main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, cex.main=0.6)
                    }
                }
            }
        dev.off()
        
        png(paste(histoName, ".png", sep=""))

            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ] 
                    All<-subset[[get("channel")]][which(subset$SpotType!=-1)]

                    if (sum(!is.na(All))!=0){

                        if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                            d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                            computeHistoBreaks=seq(a, b, d)    
                        }else{
                            a<-min(All, na.rm=T)
                            b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                            d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                            computeHistoBreaks=seq(a, b, d)
                        }
    
                        hist(All, breaks=computeHistoBreaks, 
                        main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, cex.main=0.6)
                    }
                }
            }
        dev.off()
    }else{
    
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-dataset[which(dataset$ScreenNb == i), ]
                All<-subset[[get("channel")]][which(subset$SpotType!=-1)]

                if (sum(!is.na(All))!=0){

                    if (round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)!=0){
                        a<-min(All, na.rm=T)
                        b<-max(All, na.rm=T)+(round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20))
                        d<-round((max(All, na.rm=T)-min(All, na.rm=T)+1)/20)
                        computeHistoBreaks=seq(a, b, d)    
                    }else{
                        a<-min(All, na.rm=T)
                        b<-max(All, na.rm=T)+abs(min(All, na.rm=T)+0.5)
                        d<-(max(All, na.rm=T)-min(All, na.rm=T)+1)/20
                        computeHistoBreaks=seq(a, b, d)
                    }
                    pdf(paste(histoName, "(Exp. ", i, ").pdf", sep=""))
                        hist(All, breaks=computeHistoBreaks, 
                        main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, cex.main=0.6)
                    dev.off()
                    
                    png(paste(histoName, "(Exp. ", i, ").png", sep=""), width=300, height=300)
                        hist(All, breaks=computeHistoBreaks, 
                        main=paste(plotTitle, "for Exp.", i, sep=" "), xlab=channel, cex.main=0.6)
                    dev.off()
                }
            }
        }
    }
    invisible(list(histoName, minOfScreens, numOfScreens))
}




plotQQ<-function(header, dataset, channel, plotTitle, showPlot){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    if (showPlot == 1){
        if(interactive()){
            qqnorm(dataset[[get("channel")]], main=plotTitle)
            qqline(dataset[[get("channel")]])
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(plotName, ".pdf", sep=""))
        qqnorm(dataset[[get("channel")]], main=plotTitle)
        qqline(dataset[[get("channel")]])
    dev.off()
    png(paste(plotName, ".png", sep=""), width=300, height=300)
        qqnorm(dataset[[get("channel")]], main=plotTitle)
        qqline(dataset[[get("channel")]])
    dev.off()
    invisible(plotName)
}


plotQQperplate<-function(header, dataset, channel, plotTitle, plotDesign, 
showPlot){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    if (showPlot == 1){
        if(interactive()){
            for (j in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == j))>0){
    
                    subset<-dataset[which(dataset$ScreenNb == j), ]
            
                    minOfPlates<-min(subset$LabtekNb)
                    numOfPlates<-max(subset$LabtekNb)

                    x11()
                    if (plotDesign == 1){
                        c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                        c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                        par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
                    }
        
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
    
                            qqnorm(subsubset[[get("channel")]], 
                            main=paste(plotTitle, "for plate", i, sep=" "), cex.main=0.8)
                            qqline(subsubset[[get("channel")]])
                            mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                        }
                    }
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
    
            subset<-dataset[which(dataset$ScreenNb == j), ]
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
            
            plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
            
            if (plotDesign == 1){

                pdf(paste(plotName, "Exp", j, ".pdf", sep=""))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
                
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
    
                            qqnorm(subsubset[[get("channel")]], 
                            main=paste(plotTitle, "for plate", i, sep=" "), cex.main=0.8)
                            qqline(subsubset[[get("channel")]])
                            mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                        }
                    }
                dev.off()
                
                png(paste(plotName, "Exp", j, ".png", sep=""))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
                
                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]
    
                            qqnorm(subsubset[[get("channel")]], 
                            main=paste(plotTitle, "for plate", i, sep=" "), cex.main=0.8)
                            qqline(subsubset[[get("channel")]])
                            mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                        }
                    }
                dev.off()
            }else{
                
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                    
                        subsubset<-subset[which(subset$LabtekNb == i), ]
    
                        pdf(paste(plotName, "Exp", j, "_PerPlate", i, ".pdf", sep=""))
                            qqnorm(subsubset[[get("channel")]], 
                            main=paste(plotTitle, "for plate", i, sep=" "), cex.main=0.8)
                            qqline(subsubset[[get("channel")]])
                            mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                        dev.off()
                        
                        png(paste(plotName, "Exp", j, "_PerPlate", i, ".png", sep=""), 
                        width=300, height=300)
                            qqnorm(subsubset[[get("channel")]], 
                            main=paste(plotTitle, "for plate", i, sep=" "), cex.main=0.8)
                            qqline(subsubset[[get("channel")]])
                            mtext(paste(plotTitle, "for Experiment", j, sep=" "), side=3, outer=T)
                        dev.off()
                    }
                }            
            }
        }
    }
    invisible(list(plotName, c(minOfScreens, numOfScreens), 
    c(minOfPlates, numOfPlates)))
}

plotQQperscreen<-function(header, dataset, channel, plotTitle, plotDesign, 
showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    
    if (showPlot == 1){
        if(interactive()){
            if (plotDesign == 1){
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
                c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
                par(mfrow=c(c1, c2))
            }

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
    
                    subset<-dataset[which(dataset$ScreenNb == i), ]
    
                    qqnorm(subset[[get("channel")]], 
                    main=paste(plotTitle, "for screen", i, sep=" "), cex.main=0.8)
                    qqline(subset[[get("channel")]])
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    if (plotDesign == 1){
        pdf(paste(plotName, ".pdf", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
        
            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]
        
                    qqnorm(subset[[get("channel")]], 
                    main=paste(plotTitle, "for screen", i, sep=" "), cex.main=0.8)
                    qqline(subset[[get("channel")]])
                }
            }
        dev.off()
        
        png(paste(plotName, ".png", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
        
            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]
        
                    qqnorm(subset[[get("channel")]], 
                    main=paste(plotTitle, "for screen", i, sep=" "), cex.main=0.8)
                    qqline(subset[[get("channel")]])
                }
            }
        dev.off()
    }else{
        
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-dataset[which(dataset$ScreenNb == i), ]
        
                pdf(paste(plotName, " (Exp.", i, ").pdf", sep=""))
                    qqnorm(subset[[get("channel")]], 
                    main=paste(plotTitle, "for screen", i, sep=" "), cex.main=0.8)
                    qqline(subset[[get("channel")]])
                dev.off()
                
                png(paste(plotName, " (Exp.", i, ").png", sep=""), width=300, height=300)
                    qqnorm(subset[[get("channel")]], 
                    main=paste(plotTitle, "for screen", i, sep=" "), cex.main=0.8)
                    qqline(subset[[get("channel")]])
                dev.off()
            }
        }    
    }
    invisible(list(plotName, minOfScreens, numOfScreens))
}


replicatesSpearmancor<-function(header, dataset, flag, col4val, col4anno, 
fileNameSuffix){

    if (flag == 1){

        numOfScreens<-max(dataset$ScreenNb)
        minOfScreens<-min(dataset$ScreenNb)

        headerTemp<-strsplit(header[1], ",")
        filnam<-paste(headerTemp[[1]][2], fileNameSuffix, "Spearmancor.txt", sep="_")
        colNams<-paste("Exp", "Replicate", "Replicate", "Correlation_coeff", sep="\t")
        write.table(colNams, file=filnam, quote=F, col.names=F, row.names=F)

        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
            
                subset<-dataset[which(dataset$ScreenNb == i), ]

                CVmatrix<-generateReplicateMat(subset, 2, "Intensities", col4val, col4anno)
                CVmatrix<-CVmatrix[[1]]

                print("", quote=F)
                ptt<-paste("Spearman\'s rank correlation coefficient for Exp.", i, ":", sep=" ")
                print(ptt, quote=F)
                print("", quote=F)

                if(sum(!is.na(CVmatrix[,1])) > 0 & sum(!is.na(CVmatrix[,2])) > 0){
                    result1<-cor(CVmatrix[, 1:2], use="complete.obs", method="spearman")
                    
                    print(paste("Replicates 1 and 2:", result1[1, 2], sep=" "), quote=F)
                    
                    write.table(paste(i, 1, 2, result1[1, 2], sep="\t"), file=filnam, quote=F,
                    col.names=F, row.names=F, append=T)
                }else{
                    print(paste("Replicates 1 and 2:", NA, sep=" "), quote=F)
                    
                    write.table(paste(i, 1, 2, NA, sep="\t"), file=filnam, quote=F,
                    col.names=F, row.names=F, append=T)
                }
                
                if(sum(!is.na(CVmatrix[,2])) > 0 & sum(!is.na(CVmatrix[,3])) > 0){
                    result2<-cor(CVmatrix[, 2:3], use="complete.obs", method="spearman")
                    
                    print(paste("Replicates 2 and 3:", result2[1, 2], sep=" "), quote=F)
                    
                    write.table(paste(i, 2, 3, result2[1, 2], sep="\t"), file=filnam, quote=F, 
                    col.names=F, row.names=F, append=T)
                }else{
                    print(paste("Replicates 2 and 3:", NA, sep=" "), quote=F)
                    
                    write.table(paste(i, 2, 3, NA, sep="\t"), file=filnam, quote=F, 
                    col.names=F, row.names=F, append=T)
                }
                
                if(sum(!is.na(CVmatrix[,1])) > 0 & sum(!is.na(CVmatrix[,3])) > 0){
                    result3<-cor(CVmatrix[, 1:3], use="complete.obs", method="spearman")
                    
                    print(paste("Replicates 1 and 3:", result3[1, 2], sep=" "), quote=F)
                    
                    write.table(paste(i, 1, 3, result3[1, 2], sep="\t"), file=filnam, quote=F, 
                    col.names=F, row.names=F, append=T)
                }else{
                    print(paste("Replicates 1 and 3:", NA, sep=" "), quote=F)
                    
                    write.table(paste(i, 1, 3, NA, sep="\t"), file=filnam, quote=F, 
                    col.names=F, row.names=F, append=T)
                }


            }else{
                write.table(paste(i, 1, 2, NA, sep="\t"), file=filnam, quote=F, col.names=F, 
                row.names=F, append=T)
                write.table(paste(i, 2, 3, NA, sep="\t"), file=filnam, quote=F, col.names=F, 
                row.names=F, append=T)
                write.table(paste(i, 1, 3, NA, sep="\t"), file=filnam, quote=F, col.names=F, 
                row.names=F, append=T)
            }
        }
    }
    
    if (flag == 2){
    
        numOfScreens<-max(dataset$ScreenNb)
        minOfScreens<-min(dataset$ScreenNb)

        tempSubset<-createSubset(dataset, dataset$ScreenNb, 1)
        vek<-c("SigIntensity", "Background", "NbCells", "PcCells")
        tempVar<-summarizeRepsNoFiltering(tempSubset, rms, vek, "GeneName", 
        NA_character_)
        matrixSummarizedReplicates<-matrix(0, nrow(tempVar), 
        length(unique(dataset$ScreenNb)))

        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                tempSubset<-createSubset(dataset, dataset$ScreenNb, i)
                tempVar<-summarizeRepsNoFiltering(tempSubset, rms, vek, "GeneName", 
                NA_character_)
                matrixSummarizedReplicates[, i]<-as.vector(tempVar[[get("col4val")]])
                
            }
        }
        result<-cor(matrixSummarizedReplicates, use="complete.obs", method="spearman")

        headerTemp<-strsplit(header[1], ",")
        filnam<-paste(headerTemp[[1]][2], "Spearmancor_AllExp.txt", sep="_")
        write.table(paste("Experiment", "Experiment", "Correlation_coeff", sep="\t"), 
        file=filnam, quote=F, col.names=F, row.names=F)

        print("", quote=F)
        print("Spearman\'s rank correlation coefficient:", quote=F)
        print("", quote=F)

        for (j in 1:(numOfScreens-minOfScreens)){
            for (k in (j+1):(numOfScreens-minOfScreens+1)){
                if (length(which(dataset$ScreenNb == (j+minOfScreens-1)))>0 
                & length(which(dataset$ScreenNb == (k+minOfScreens-1)))>0){
                
                    print(paste("Experiments", (j+minOfScreens-1), "and", (k+minOfScreens-1), 
                    ":", result[j, k], sep=" "), quote=F)
                    
                    tow<-paste((j+minOfScreens-1), (k+minOfScreens-1), result[j, k], sep="\t")
                    write.table(tow, file=filnam, quote=F, col.names=F, row.names=F, append=T)            
                }else{
                    tow<-paste((j+minOfScreens-1), (k+minOfScreens-1), NA, sep="\t")
                    write.table(tow, file=filnam, quote=F, col.names=F, row.names=F, append=T)
                }
            }
        }
    }
    Tabelle<-read.table(filnam, header=T)
}


replicatesCV<-function(header, dataset, PlotTitle, col4val, col4anno, 
plotDesign, showPlot){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    
    if (showPlot == 1){
        if(interactive()){
    
            if (plotDesign == 1){
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
                c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
                par(mfrow=c(c1, c2))
            }

            for (m in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == m))>0){

                    subset<-dataset[which(dataset$ScreenNb == m), ]

                    temp1<-generateReplicateMat(subset, 2, "Intensities", col4val, col4anno)
                    CVmatrix<-temp1[[1]]
                    indexPosControls<-temp1[[2]]
                    indexNegControls<-temp1[[3]]


                    CVvec<-rep(0, nrow(CVmatrix))
                    names(CVvec)<-rownames(CVmatrix)
                    meanVec<-rep(0, nrow(CVmatrix))
                    names(meanVec)<-rownames(CVmatrix)

                    for (i in 1:length(CVvec)){
    
                        temp<-CVmatrix[i, !is.na(CVmatrix[i, ])]
                        CVvec[i]<-sd(temp, na.rm=T)/mean(temp, na.rm=T)
                        meanVec[i]<-mean(temp, na.rm=T)
        
                    }
    
                    xxlim<-c(min(meanVec)-min(meanVec)/20, max(meanVec)+max(meanVec)/20)
                    yylim<-c(min(CVvec)-min(CVvec)/20, max(CVvec)+max(CVvec)/20)
                    plot(meanVec, CVvec, main=paste(PlotTitle, ", Exp. ", m, sep=""), 
                    xlab="Mean Intensity", ylab="CV", cex.main=0.6, xlim=xxlim, ylim=yylim)
        
                    points(meanVec[indexPosControls], CVvec[indexPosControls], col="green")
                    points(meanVec[indexNegControls], CVvec[indexNegControls], col="red")
                    identify(meanVec, CVvec, labels=names(meanVec))
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], PlotTitle, sep="_")
    
    if (plotDesign == 1){
        pdf(paste(plotName, ".pdf", sep=""))

            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
            for (m in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == m))>0){

                    subset<-dataset[which(dataset$ScreenNb == m), ]

                    temp1<-generateReplicateMat(subset, 2, "Intensities", col4val, col4anno)
                    CVmatrix<-temp1[[1]]
                    indexPosControls<-temp1[[2]]
                    indexNegControls<-temp1[[3]]

                    CVvec<-rep(0, nrow(CVmatrix))
                    names(CVvec)<-rownames(CVmatrix)
                    meanVec<-rep(0, nrow(CVmatrix))
                    names(meanVec)<-rownames(CVmatrix)

                    for (i in 1:length(CVvec)){
        
                        temp<-CVmatrix[i, !is.na(CVmatrix[i, ])]
                        CVvec[i]<-sd(temp, na.rm=T)/mean(temp, na.rm=T)
                        meanVec[i]<-mean(temp, na.rm=T)
                    }
    
                    plot(meanVec, CVvec, main=paste(PlotTitle, ", Exp. ", m, sep=""), 
                    xlab="Mean Intensity", ylab="CV", cex.main=0.6)
                    points(meanVec[indexPosControls], CVvec[indexPosControls], col="green")
                    points(meanVec[indexNegControls], CVvec[indexNegControls], col="red")
                    identify(meanVec, CVvec, labels=names(meanVec))
                }
            }
        dev.off()
        
        png(paste(plotName, ".png", sep=""))

            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
            for (m in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == m))>0){

                    subset<-dataset[which(dataset$ScreenNb == m), ]

                    temp1<-generateReplicateMat(subset, 2, "Intensities", col4val, col4anno)
                    CVmatrix<-temp1[[1]]
                    indexPosControls<-temp1[[2]]
                    indexNegControls<-temp1[[3]]

                    CVvec<-rep(0, nrow(CVmatrix))
                    names(CVvec)<-rownames(CVmatrix)
                    meanVec<-rep(0, nrow(CVmatrix))
                    names(meanVec)<-rownames(CVmatrix)

                    for (i in 1:length(CVvec)){
        
                        temp<-CVmatrix[i, !is.na(CVmatrix[i, ])]
                        CVvec[i]<-sd(temp, na.rm=T)/mean(temp, na.rm=T)
                        meanVec[i]<-mean(temp, na.rm=T)
                    }
    
                    plot(meanVec, CVvec, main=paste(PlotTitle, ", Exp. ", m, sep=""), 
                    xlab="Mean Intensity", ylab="CV", cex.main=0.6)
                    points(meanVec[indexPosControls], CVvec[indexPosControls], col="green")
                    points(meanVec[indexNegControls], CVvec[indexNegControls], col="red")
                    identify(meanVec, CVvec, labels=names(meanVec))
                }
            }
        dev.off()
    }else{
    
        for (m in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == m))>0){

                subset<-dataset[which(dataset$ScreenNb == m), ]

                temp1<-generateReplicateMat(subset, 2, "Intensities", col4val, col4anno)
                CVmatrix<-temp1[[1]]
                indexPosControls<-temp1[[2]]
                indexNegControls<-temp1[[3]]

                CVvec<-rep(0, nrow(CVmatrix))
                names(CVvec)<-rownames(CVmatrix)
                meanVec<-rep(0, nrow(CVmatrix))
                names(meanVec)<-rownames(CVmatrix)

                for (i in 1:length(CVvec)){
        
                    temp<-CVmatrix[i, !is.na(CVmatrix[i, ])]
                    CVvec[i]<-sd(temp, na.rm=T)/mean(temp, na.rm=T)
                    meanVec[i]<-mean(temp, na.rm=T)
                }
    
                pdf(paste(plotName, " (Exp.", m, ").pdf", sep=""))
                    plot(meanVec, CVvec, main=paste(PlotTitle, ", Exp. ", m, sep=""), 
                    xlab="Mean Intensity", ylab="CV", cex.main=0.6)
                    points(meanVec[indexPosControls], CVvec[indexPosControls], col="green")
                    points(meanVec[indexNegControls], CVvec[indexNegControls], col="red")
                    identify(meanVec, CVvec, labels=names(meanVec))
                dev.off()
                
                png(paste(plotName, " (Exp.", m, ").png", sep=""), width=300, height=300)
                    plot(meanVec, CVvec, main=paste(PlotTitle, ", Exp. ", m, sep=""), 
                    xlab="Mean Intensity", ylab="CV", cex.main=0.6)
                    points(meanVec[indexPosControls], CVvec[indexPosControls], col="green")
                    points(meanVec[indexNegControls], CVvec[indexNegControls], col="red")
                    identify(meanVec, CVvec, labels=names(meanVec))
                dev.off()
            }
        }
    }
    invisible(list(plotName, minOfScreens, numOfScreens))
}

makeBoxplotControls<-function(header, dataset, channel, plotTitle, showPlot){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]
    tempVec<-rep(0, length(dataset$SpotType))
    tempVec[dataset$SpotType == 0]<-"Neg. contr."
    tempVec[dataset$SpotType == 1]<-"Pos. contr."
    tempVec[dataset$SpotType == 2]<-"Exp. data"

    if (showPlot == 1){
        if(interactive()){
            boxplot(dataset[[get("channel")]]~tempVec, main=plotTitle, ylab=channel, 
            cex=0.8)
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(plotName, ".pdf", sep=""))
        boxplot(dataset[[get("channel")]]~tempVec, main=plotTitle, ylab=channel, 
        cex.main=0.8)
    dev.off()
    
    png(paste(plotName, ".png", sep=""), width=300, height=300)
        boxplot(dataset[[get("channel")]]~tempVec, main=plotTitle, ylab=channel, 
        cex.main=0.8, cex.lab=0.7, cex.axis=0.6)
    dev.off()
    invisible(plotName)
}


makeBoxplotControlsPerPlate<-function(header, dataset, channel, plotTitle, 
plotDesign, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    if (showPlot == 1){
        if(interactive()){
            for (j in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == j))>0){
    
                    subset<-dataset[which(dataset$ScreenNb == j), ]
                    minOfPlates<-min(subset$LabtekNb)
                    numOfPlates<-max(subset$LabtekNb)
            
                    if (plotDesign == 1){
                        c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                        c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                        par(mfrow=c(c1, c2))
                    }

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                
                            subsubset<-subset[which(subset$LabtekNb == i), ]
        
                            tempVec<-rep(0, length(subsubset$SpotType))
                            tempVec[subsubset$SpotType == 0]<-"Neg. controls"
                            tempVec[subsubset$SpotType == 1]<-"Pos. controls"
                            tempVec[subsubset$SpotType == 2]<-"Exp. data"

                            boxplot(subsubset[[get("channel")]]~tempVec, 
                            main=paste(plotTitle, "for plate", i, sep=" "), ylab=channel, cex=0.8)
                        }
                    }
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){

            subset<-dataset[which(dataset$ScreenNb == j), ]
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)

            plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
            
            if (plotDesign == 1){
                pdf(paste(plotName, "_Exp_", j, ".pdf", sep=""))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2))

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]

                            tempVec<-rep(0, length(subsubset$SpotType))
                            tempVec[subsubset$SpotType == 0]<-"Neg. controls"
                            tempVec[subsubset$SpotType == 1]<-"Pos. controls"
                            tempVec[subsubset$SpotType == 2]<-"Exp. data"

                            boxplot(subsubset[[get("channel")]]~tempVec, 
                            main=paste(plotTitle, "for plate", i, sep=" "), ylab=channel, cex=0.8)
                        }
                    }
                dev.off()
                
                png(paste(plotName, "_Exp_", j, ".png", sep=""))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2))

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]

                            tempVec<-rep(0, length(subsubset$SpotType))
                            tempVec[subsubset$SpotType == 0]<-"Neg. controls"
                            tempVec[subsubset$SpotType == 1]<-"Pos. controls"
                            tempVec[subsubset$SpotType == 2]<-"Exp. data"

                            boxplot(subsubset[[get("channel")]]~tempVec, 
                            main=paste(plotTitle, "for plate", i, sep=" "), ylab=channel, cex=0.8)
                        }
                    }
                dev.off()
            }else{
            
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                    
                        subsubset<-subset[which(subset$LabtekNb == i), ]

                        tempVec<-rep(0, length(subsubset$SpotType))
                        tempVec[subsubset$SpotType == 0]<-"Neg. contr."
                        tempVec[subsubset$SpotType == 1]<-"Pos. contr."
                        tempVec[subsubset$SpotType == 2]<-"Exp. data"
                        pdf(paste(plotName, "_Exp_", j, "_PerPlate", i, ".pdf", sep=""))
                            boxplot(subsubset[[get("channel")]]~tempVec, 
                            main=paste(plotTitle, "for plate", i, sep=" "), ylab=channel, cex=0.8)
                        dev.off()
                        
                        png(paste(plotName, "_Exp_", j, "_PerPlate", i, ".png", sep=""), width=300, height=300)
                            tit<-paste(plotTitle, "for plate", i, sep=" ")
                            boxplot(subsubset[[get("channel")]]~tempVec, main=tit, ylab=channel, 
                            cex.main=0.6, cex.lab=0.7, cex.axis=0.6)
                        dev.off()
                    }
                }            
            }
        }
    }
    invisible(list(plotName, c(minOfScreens, numOfScreens), 
    c(minOfPlates, numOfPlates)))
}


makeBoxplotControlsPerScreen<-function(header, dataset, channel, plotTitle, 
plotDesign, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    
    if (showPlot == 1){
        if(interactive()){
    
            if (plotDesign == 1){
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
                c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
                par(mfrow=c(c1, c2))
            }

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
    
                    subset<-dataset[which(dataset$ScreenNb == i), ]
    
                    tempVec<-rep(0, length(subset$SpotType))
                    tempVec[subset$SpotType == 0]<-"Neg. contr."
                    tempVec[subset$SpotType == 1]<-"Pos. contr."
                    tempVec[subset$SpotType == 2]<-"Exp. data"

                    boxplot(subset[[get("channel")]]~tempVec, 
                    main=paste(plotTitle, "for Exp.", i, sep=" "), ylab=channel, cex=0.8)
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    if (plotDesign == 1){
        pdf(paste(plotName, ".pdf", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]

                    tempVec<-rep(0, length(subset$SpotType))
                    tempVec[subset$SpotType == 0]<-"Neg. controls"
                    tempVec[subset$SpotType == 1]<-"Pos. controls"
                    tempVec[subset$SpotType == 2]<-"Exp. data"

                    boxplot(subset[[get("channel")]]~tempVec, 
                    main=paste(plotTitle, "for Exp.", i, sep=" "), ylab=channel, cex=0.8)
                }
            }
        dev.off()
        
        png(paste(plotName, ".png", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]

                    tempVec<-rep(0, length(subset$SpotType))
                    tempVec[subset$SpotType == 0]<-"Neg. controls"
                    tempVec[subset$SpotType == 1]<-"Pos. controls"
                    tempVec[subset$SpotType == 2]<-"Exp. data"

                    boxplot(subset[[get("channel")]]~tempVec, 
                    main=paste(plotTitle, "for Exp.", i, sep=" "), ylab=channel, cex=0.8)
                }
            }
        dev.off()
    }else{
    
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-dataset[which(dataset$ScreenNb == i), ]

                tempVec<-rep(0, length(subset$SpotType))
                tempVec[subset$SpotType == 0]<-"Neg. contr."
                tempVec[subset$SpotType == 1]<-"Pos. contr."
                tempVec[subset$SpotType == 2]<-"Exp. data"
                    
                pdf(paste(plotName, " (Exp.", i, ").pdf", sep=""))
                    boxplot(subset[[get("channel")]]~tempVec, 
                    main=paste(plotTitle, "for Exp.", i, sep=" "), ylab=channel, cex.main=0.8)
                dev.off()
                
                png(paste(plotName, " (Exp.", i, ").png", sep=""), width=300, height=300)
                    boxplot(subset[[get("channel")]]~tempVec, 
                    main=paste(plotTitle, "for Exp.", i, sep=" "), ylab=channel, cex.main=0.6, cex.lab=0.7, cex.axis=0.6)
                dev.off()
            }
        }    
    }
    invisible(list(plotName, minOfScreens, numOfScreens))
}


makeBoxplotPerScreen<-function(header, dataset, channel, plotTitle, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]
    
    if (showPlot == 1){
        if(interactive()){
            boxplot(dataset[[get("channel")]]~dataset$ScreenNb, 
            main=plotTitle, xlab="Exp. ##", ylab=channel, cex.main=0.8)
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(plotName, ".pdf", sep=""))
        boxplot(dataset[[get("channel")]]~dataset$ScreenNb, 
        main=plotTitle, xlab="Exp. ##", ylab=channel, cex.main=0.8)
    dev.off()
    png(paste(plotName, ".png", sep=""), width=300, height=300)
        boxplot(dataset[[get("channel")]]~dataset$ScreenNb, 
        main=plotTitle, xlab="Exp. ##", ylab=channel, cex.main=0.8)
    dev.off()
    invisible(plotName)
}

makeBoxplotPerPlate<-function(header, dataset, channel, plotTitle, plotDesign, showPlot){

    dataset<-dataset[which(dataset$SpotType!=-1), ]

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    
    
    if(plotDesign == 1){
    
        if (showPlot == 1){
            if(interactive()){
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
                c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
                par(mfrow=c(c1, c2))

                for (i in minOfScreens:numOfScreens){
                    if (length(which(dataset$ScreenNb == i))>0){

                        subset<-dataset[which(dataset$ScreenNb == i), ]
                        tit<-paste(plotTitle, "for Exp.", i, sep=" ")
                        boxplot(subset[[get("channel")]]~subset$LabtekNb, main=tit, xlab="Plate", 
                        ylab=channel, cex=0.7, cex.main=0.7)
                    }
                }
            }
        }

        headerTemp<-strsplit(header[1], ",")
        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
        pdf(paste(plotName, ".pdf", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]
                    tit<-paste(plotTitle, "for Exp.", i, sep=" ")
                    boxplot(subset[[get("channel")]]~subset$LabtekNb, main=tit, xlab="Plate", 
                    ylab=channel, cex=0.8, cex.main=0.7)
                }
            }    
        dev.off()
        png(paste(plotName, ".png", sep=""))
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))
    
                for (i in minOfScreens:numOfScreens){
                    if (length(which(dataset$ScreenNb == i))>0){
            
                        subset<-dataset[which(dataset$ScreenNb == i), ]
                        tit<-paste(plotTitle, "for Exp.", i, sep=" ")
                        boxplot(subset[[get("channel")]]~subset$LabtekNb, main=tit, xlab="Plate", 
                        ylab=channel, cex=0.8, cex.main=0.7)
                    }
                }    
        dev.off()
    }
    
    
    if(plotDesign == 2){

        if (showPlot == 1){
            if(interactive()){

                for (i in minOfScreens:numOfScreens){
                    if (length(which(dataset$ScreenNb == i))>0){

                        subset<-dataset[which(dataset$ScreenNb == i), ]
                        tit<-paste(plotTitle, "for Exp.", i, sep=" ")
                        boxplot(subset[[get("channel")]]~subset$LabtekNb, main=tit, xlab="Plate", 
                        ylab=channel, cex=0.7, cex.main=0.8)
                    }
                }
            }
        }

        headerTemp<-strsplit(header[1], ",")
        plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
        
        for (i in minOfScreens:numOfScreens){
            if (length(which(dataset$ScreenNb == i))>0){
        
                subset<-dataset[which(dataset$ScreenNb == i), ]
                
                pdf(paste(plotName, "(Exp. ", i, ").pdf", sep=""))
                    tit<-paste(plotTitle, "for Exp.", i, sep=" ")
                    boxplot(subset[[get("channel")]]~subset$LabtekNb, main=tit, xlab="Plate", ylab=channel, cex=0.8, cex.main=0.8)
                dev.off()
                
                png(paste(plotName, "(Exp. ", i, ").png", sep=""), width=300, height=300)    
                    tit<-paste(plotTitle, "for Exp.", i, sep=" ")
                    boxplot(subset[[get("channel")]]~subset$LabtekNb, main=tit, xlab="Plate", ylab=channel, cex=0.8, cex.main=0.8)
                dev.off()
            }
        }    
    }
    invisible(list(plotName, minOfScreens, numOfScreens))
}

makeBoxplot4PlateType<-function(header, dataset, channel, plotTitle, showPlot){
    
    dataset<-dataset[which(dataset$SpotType!=-1), ]
    numOfPlates<-max(dataset$LabtekNb)
    minOfPlates<-min(dataset$LabtekNb)
    
    if (showPlot == 1){
        if(interactive()){
            for (i in minOfPlates:numOfPlates){
                if (length(which(dataset$LabtekNb == i))>0){
                    subsetPlateType<-createSubset(dataset, dataset$LabtekNb, i)
                    x11()
                    tit<-paste(plotTitle, "for Plate", i, sep=" ")
                    boxplot(subsetPlateType[[get("channel")]]~subsetPlateType$ScreenNb, main=tit, 
                    xlab="Exp. ##", ylab=channel, cex=0.8, cex.main=0.8)    
                }
            }
        }
    }
    
    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, "", sep="_")
    
    for (i in minOfPlates:numOfPlates){
        if (length(which(dataset$LabtekNb == i))>0){
            subsetPlateType<-createSubset(dataset, dataset$LabtekNb, i)
                
            pdf(paste(plotName, " for Plate", i, ".pdf", sep=""))
                tit<-paste(plotTitle, "for Plate", i, sep=" ")
                boxplot(subsetPlateType[[get("channel")]]~subsetPlateType$ScreenNb, main=tit, xlab="Exp. ##", ylab=channel, cex=0.8, cex.main=0.8)    
            dev.off()
            
            png(paste(plotName, " for Plate", i, ".png", sep=""), width=300, height=300)
                tit<-paste(plotTitle, "for Plate", i, sep=" ")
                boxplot(subsetPlateType[[get("channel")]]~subsetPlateType$ScreenNb, main=tit, xlab="Exp. ##", ylab=channel, cex=0.8, cex.main=0.8)
            dev.off()
        }
    }    
    invisible(list(plotName, minOfPlates, numOfPlates))
}



spatialDistrib<-function(header, dataset, plotTitle, col4plot, col4anno, 
showPlot){
    
    dataset[[get("col4plot")]][which(dataset$SpotType == -1)]<-NA_integer_

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    headerTemp<-strsplit(header[1], ",")
    basicPlotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
    
            subset<-dataset[which(dataset$ScreenNb == j), ]
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
            
            for (i in minOfPlates:numOfPlates){
                if (length(which(subset$LabtekNb == i))>0){
                
                    subsubset<-subset[which(subset$LabtekNb == i), ]
                
                    if (sum(!is.na(subsubset[[get("col4plot")]])) > 0){
                
                        posNegLabelVec<-rep(NA_character_, nrow(subsubset))
                        posNegLabelVec[which(subsubset$SpotType == 0)]<-"N"
                        posNegLabelVec[which(subsubset$SpotType == 1)]<-"P"

                        if (showPlot == 1){
                            if(interactive()){
                                x11()
                            }
                        }
                        ColNo<-max(subsubset$ColNb)
                        RowNo<-max(subsubset$RowNb)
                    
                        dings<-plotPlate((as.numeric(t(subsubset[[get("col4plot")]]))), ncol=ColNo, 
                        nrow=RowNo, char=posNegLabelVec, cex.char=0.8, 
                        main=paste(plotTitle, "plate", i, "Exp.", j, sep=" "), 
                        col=brewer.pal(9, "YlOrBr"), cex.main=0.8, na.action='xout')
                    
                        if (ColNo == 12 & RowNo == 8){
                        ##plotPlate function coordinates are wrong for these parameters...
                            dings$coord[, 2]<-dings$coord[, 2]-100
                            dings$coord[, 4]<-dings$coord[, 4]-100
                        }

                        png(paste(basicPlotName, "Exp", j, "Plate", i, "4html.png", sep="_"), width=300, height=300)
                            plotPlate((as.numeric(t(subsubset[[get("col4plot")]]))), ncol=ColNo, 
                            nrow=RowNo, char=posNegLabelVec, cex.char=0.8, 
                            main=paste(plotTitle, "plate", i, "Exp.", j, sep=" "), 
                            col=brewer.pal(9, "YlOrBr"), cex.main=0.8, na.action='xout')
                        dev.off()
                    
                        png(paste(basicPlotName, "Exp", j, "Plate", i, ".png", sep="_"), 
                        width=dings$width, height=dings$height)
                    
                            plotPlate((as.numeric(t(subsubset[[get("col4plot")]]))), ncol=ColNo, 
                            nrow=RowNo, char=posNegLabelVec, cex.char=0.8, 
                            main=paste(plotTitle, "plate", i, "Exp.", j, sep=" "), 
                            col=brewer.pal(9, "YlOrBr"), cex.main=0.8, na.action='xout')
                        dev.off()

                        title<-subsubset[[get("col4anno")]]
                        con<-file(paste(basicPlotName, "Exp", j, "Plate", i, ".html", sep="_"), open="w")
                    
                        imageMap(dings$coord, con, list(TITLE=title, 
                        TARGET=rep("main", nrow(dings$coord))), paste(basicPlotName, "Exp", j, 
                        "Plate", i, ".png", sep="_"))
                        close(con)
                    }
                }
            }
        }
    }
    invisible(list(basicPlotName, c(minOfScreens, numOfScreens), c(minOfPlates, numOfPlates)))
}

compareReplicates<-function(header, dataset, plotTitle, col4val, col4anno, 
plotDesign, showPlot){

    headerTemp<-strsplit(header[1], ",")

    dataset<-dataset[which(dataset$SpotType!=-1), ]

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    
    maxCombinationNum<-0
    
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
    
            subset<-dataset[which(dataset$ScreenNb == i), ]


            temp1<-generateReplicateMat(subset, 2, "Intensities", col4val, col4anno)

            replicaMatrix<-temp1[[1]]
            indexPosControls<-temp1[[2]]
            indexNegControls<-temp1[[3]]

            ColNo<-0
            for (m in 1:ncol(replicaMatrix)){
                if (sum(is.na(replicaMatrix[, m]))>0.4*nrow(replicaMatrix)){
                    ColNo<-m-1
                    break
                }
            }
            if (ColNo == 0){
                ColNo<-ncol(replicaMatrix)
            }
            
            combinationNum<-choose(ColNo, 2)
            if (maxCombinationNum<combinationNum){
                maxCombinationNum<-combinationNum
            }
            
            if (showPlot == 1){
                if(interactive()){
                    x11()
                    if (plotDesign == 1){
                        c1<-ceiling(sqrt(combinationNum))
                        c2<-ceiling(combinationNum/ceiling(sqrt(combinationNum)))
                        par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
                    }
    
                    for (j in 1:(ColNo-1)){
                        for (k in (j+1):ColNo){

                            max1<-max(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
                            min1<-min(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)

                            xxlim<-c(min1-(min1)/20, max1+(max1)/20)
                            yylim<-c(min1-(min1)/20, max1+(max1)/20)
                            plot(replicaMatrix[, j], replicaMatrix[, k], xlab=paste("replicate", j, sep=" "), 
                            ylab=paste("replicate", k, sep=" "), xlim=xxlim, ylim=yylim)
                        
                            points(replicaMatrix[indexPosControls, j], replicaMatrix[indexPosControls, k], 
                            col="green")
                            points(replicaMatrix[indexNegControls, j], replicaMatrix[indexNegControls, k], 
                            col="red")
                            lines(c(min1-(min1)/20, max1+(max1)/20), c(min1-(min1)/20, max1+(max1)/20))
                            mtext(paste(plotTitle, "for Exp.", i, sep=" "), side=3, outer=T, cex=0.8)
                            identify(replicaMatrix[, j], replicaMatrix[, k], labels=rownames(replicaMatrix))
                        }
                    }
                }
            }

            plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
            
            if (plotDesign == 1){
                pdf(paste(plotName, "Exp", i, ".pdf", sep=""))
                    c1<-ceiling(sqrt(combinationNum))
                    c2<-ceiling(combinationNum/ceiling(sqrt(combinationNum)))
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))

                    for (j in 1:(ColNo-1)){
                        for (k in (j+1):ColNo){

                            max1<-max(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
                            min1<-min(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
    
                            xxlab<-paste("replicate", j, sep=" ")
                            yylab<-paste("replicate", k, sep=" ")
                            plot(replicaMatrix[, j], replicaMatrix[, k], xlab=xxlab, ylab=yylab, 
                            xlim=c(min1, max1), ylim=c(min1, max1))
                            
                            points(replicaMatrix[indexPosControls, j], replicaMatrix[indexPosControls, k], 
                            col="green")
                            points(replicaMatrix[indexNegControls, j], replicaMatrix[indexNegControls, k], 
                            col="red")
                            
                            lines(c(min1, max1), c(min1, max1))
                            mtext(paste(plotTitle, "for Exp", i, sep=" "), side=3, outer=T, cex=0.8)
                            identify(replicaMatrix[, j], replicaMatrix[, k], labels=rownames(replicaMatrix))

                        }
                    }                
                dev.off()
                
                png(paste(plotName, "Exp", i, ".png", sep=""))
                    c1<-ceiling(sqrt(combinationNum))
                    c2<-ceiling(combinationNum/ceiling(sqrt(combinationNum)))
                    par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))

                    for (j in 1:(ColNo-1)){
                        for (k in (j+1):ColNo){

                            max1<-max(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
                            min1<-min(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
    
                            xxlab<-paste("replicate", j, sep=" ")
                            yylab<-paste("replicate", k, sep=" ")
                            plot(replicaMatrix[, j], replicaMatrix[, k], xlab=xxlab, ylab=yylab, 
                            xlim=c(min1, max1), ylim=c(min1, max1))
                            
                            points(replicaMatrix[indexPosControls, j], replicaMatrix[indexPosControls, k], 
                            col="green")
                            points(replicaMatrix[indexNegControls, j], replicaMatrix[indexNegControls, k], 
                            col="red")
                            
                            lines(c(min1, max1), c(min1, max1))
                            mtext(paste(plotTitle, "for Exp", i, sep=" "), side=3, outer=T, cex=0.8)
                            identify(replicaMatrix[, j], replicaMatrix[, k], labels=rownames(replicaMatrix))

                        }
                    }                  
                dev.off()
            }else{
                counter<-0
                for (j in 1:(ColNo-1)){
                    for (k in (j+1):ColNo){

                        max1<-max(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
                        min1<-min(rbind(replicaMatrix[, j], replicaMatrix[, k]), na.rm=T)
                        counter<-counter+1
    
                        pdf(paste(plotName, "_Exp_", i, "(", counter, ").pdf", sep=""))
                            xxlab<-paste("replicate", j, sep=" ")
                            yylab<-paste("replicate", k, sep=" ")
                            plot(replicaMatrix[, j], replicaMatrix[, k], xlab=xxlab, ylab=yylab, 
                            xlim=c(min1, max1), ylim=c(min1, max1))
                            
                            points(replicaMatrix[indexPosControls, j], replicaMatrix[indexPosControls, k], 
                            col="green")
                            points(replicaMatrix[indexNegControls, j], replicaMatrix[indexNegControls, k], 
                            col="red")
                            
                            lines(c(min1, max1), c(min1, max1))
                            mtext(paste(plotTitle, "for Exp", i, sep=" "), side=3, outer=T, cex=0.8)
                            identify(replicaMatrix[, j], replicaMatrix[, k], labels=rownames(replicaMatrix))
                        dev.off()
                        
                        png(paste(plotName, "_Exp_", i, "(", counter, ").png", sep=""), width=300, height=300)
                            xxlab<-paste("replicate", j, sep=" ")
                            yylab<-paste("replicate", k, sep=" ")
                            plot(replicaMatrix[, j], replicaMatrix[, k], xlab=xxlab, ylab=yylab, 
                            xlim=c(min1, max1), ylim=c(min1, max1))
                            
                            points(replicaMatrix[indexPosControls, j], replicaMatrix[indexPosControls, k], 
                            col="green")
                            points(replicaMatrix[indexNegControls, j], replicaMatrix[indexNegControls, k], 
                            col="red")
                            
                            lines(c(min1, max1), c(min1, max1))
                            mtext(paste(plotTitle, "for Exp", i, sep=" "), side=3, outer=T, cex=0.8)
                            identify(replicaMatrix[, j], replicaMatrix[, k], labels=rownames(replicaMatrix))
                        dev.off()
                    }
                }
            }
        }
    }
    invisible(list(plotName, minOfScreens, numOfScreens, maxCombinationNum))
}

compareReplicaPlates<-function(header, dataset, plotTitle, col4val){

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, ".pdf", sep="_")

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    numOfPlates<-max(dataset$LabtekNb)

    if (length(unique(dataset$ScreenNb))<2){
        stop("Comparison not possible - Only one experiment in dataset")
    }

    dataset[[get("col4val")]][which(dataset$SpotType == -1)]<-NA_integer_

    
    for (i in 1:numOfPlates){
        count<-0
        for (j in minOfScreens:(numOfScreens-1)){
        
            count<-count+1
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))-count))
            c2<-ceiling((length(unique(dataset$ScreenNb))-count)/c1)
            
            if(interactive()){
                x11()
                par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
            }
            
            for (k in (j+1):numOfScreens){
            
                if (length(which(dataset$ScreenNb == j))>0 
                & length(which(dataset$ScreenNb == k))>0){
                
                    testSub1<-createSubset(dataset, dataset$ScreenNb, j)
                    testSub2<-createSubset(dataset, dataset$ScreenNb, k)
                    
                    if (length(which(testSub1$LabtekNb == i))>0 
                    & length(which(testSub2$LabtekNb == i))>0){
                    
                        tempSubset<-createSubset(dataset, dataset$LabtekNb, i)
                        firstPlate<-createSubset(tempSubset, tempSubset$ScreenNb, j)
                        secondPlate<-createSubset(tempSubset, tempSubset$ScreenNb, k)
                        max1<-max(rbind(firstPlate[[get("col4val")]], secondPlate[[get("col4val")]]), 
                        na.rm=T)
                        min1<-min(rbind(firstPlate[[get("col4val")]], secondPlate[[get("col4val")]]), 
                        na.rm=T)

                        xxlab<-paste("Plate", i, "Exp.", j, sep=" ")
                        yylab<-paste("Plate", i, "Exp.", k, sep=" ")
                        if(interactive()){
                            plot(firstPlate[[get("col4val")]], secondPlate[[get("col4val")]], 
                            xlab=xxlab, ylab=yylab, xlim=c(min1, max1), ylim=c(min1, max1))
                            lines(c(min1, max1), c(min1, max1))
                        
                            textt<-paste(plotTitle, "plate", i, "from Exp.", j, "/ plate", i, 
                            "from other exps", sep=" ")
                            mtext(textt, side=3, outer=T, cex=0.6)
                        }
                    }
                }
            }
        }
    }

    pdf(plotName)

        for (i in 1:numOfPlates){
            count<-0
            for (j in minOfScreens:(numOfScreens-1)){
                count<-count+1
                c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))-count))
                c2<-ceiling((length(unique(dataset$ScreenNb))-count)/c1)
                par(mfrow=c(c1, c2), oma=c(0, 0, 2, 0))
            
                for (k in (j+1):numOfScreens){
                
                    if (length(which(dataset$ScreenNb == j))>0 
                    & length(which(dataset$ScreenNb == k))>0){

                        testSub1<-createSubset(dataset, dataset$ScreenNb, j)
                        testSub2<-createSubset(dataset, dataset$ScreenNb, k)
                        
                        if (length(which(testSub1$LabtekNb == i))>0 
                        & length(which(testSub2$LabtekNb == i))>0){

                            tempSubset<-createSubset(dataset, dataset$LabtekNb, i)
                            firstPlate<-createSubset(tempSubset, tempSubset$ScreenNb, j)
                            secondPlate<-createSubset(tempSubset, tempSubset$ScreenNb, k)
                            max1<-max(rbind(firstPlate[[get("col4val")]], secondPlate[[get("col4val")]]), 
                            na.rm=T)
                            min1<-min(rbind(firstPlate[[get("col4val")]], secondPlate[[get("col4val")]]), 
                            na.rm=T)

                            xxlab<-paste("Plate", i, "Exp.", j, sep=" ")
                            yylab<-paste("Plate", i, "Exp.", k, sep=" ")
                            plot(firstPlate[[get("col4val")]], secondPlate[[get("col4val")]], 
                            xlab=xxlab, ylab=yylab, xlim=c(min1, max1), ylim=c(min1, max1))
                        
                            lines(c(min1, max1), c(min1, max1))
                            textt<-paste(plotTitle, "plate", i, "from Exp.", j, "/ plate", i, 
                            "from other exps", sep=" ")
                            mtext(textt, side=3, outer=T, cex=0.6)
                        }
                    }
                }
            }
        }    
    dev.off()
}


compareReplicateSD<-function(header, dataset, plotTitle, colname4SD, col4anno){

    temp1<-generateReplicateMat(dataset, 2, "Intensities", colname4SD, col4anno)

    replicaMatrix<-temp1[[1]]
    indexPosControls<-temp1[[2]]
    indexNegControls<-temp1[[3]]

    sdVec<-rep(0, nrow(replicaMatrix))
    for (i in 1:nrow(replicaMatrix)){
        sdVec[i]<-sd(replicaMatrix[i, ], na.rm=T)
    }
    names(sdVec)<-rownames(replicaMatrix)

    RowNo<-c(ceiling(sqrt(length(sdVec))))
    ColNo<-ceiling(length(sdVec)/ceiling(sqrt(length(sdVec))))
    
    if (RowNo*ColNo > length(sdVec)){
        nbEmptyWells<-(RowNo*ColNo)-length(sdVec)
        sdVec<-c(sdVec, rep(as.double(NA_character_), nbEmptyWells))
    }

    posNegLabelVec<-rep(NA_character_, length(sdVec))
    posNegLabelVec[indexNegControls]<-"N"
    posNegLabelVec[indexPosControls]<-"P"

    dings<-plotPlate(sdVec, ncol=ColNo, nrow=RowNo, char=posNegLabelVec, 
    cex.char=0.8, main=plotTitle, col=brewer.pal(9, "YlOrBr"), cex.main=0.8, 
    na.action='xout')

    headerTemp<-strsplit(header[1], ",")
    basicPlotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

    png(paste(basicPlotName, "4html.png", sep=" "), width=300, height=300)
        plotPlate(sdVec, ncol=ColNo, nrow=RowNo, char=posNegLabelVec, 
        cex.char=0.8, main=plotTitle, col=brewer.pal(9, "YlOrBr"), cex.main=0.8, 
        na.action='xout')
    dev.off()

    png(paste(basicPlotName, ".png", sep=" "), width=dings$width, height=dings$height)
        plotPlate(sdVec, ncol=ColNo, nrow=RowNo, char=posNegLabelVec, 
        cex.char=0.8, main=plotTitle, col=brewer.pal(9, "YlOrBr"), cex.main=0.8, 
        na.action='xout')
    dev.off()

    title<-names(sdVec)
    con<-file(paste(basicPlotName, ".html", sep=" "), open="w")
    imageMap(dings$coord, con, list(TITLE=title, 
    TARGET=rep("main", nrow(dings$coord))), paste(basicPlotName, ".png", sep=" "))
    
    close(con)
    invisible(basicPlotName)
}



compareReplicateSDPerScreen<-function(header, dataset, plotTitle, colname4SD, 
col4anno, showPlot){

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
    
        subset<-dataset[which(dataset$ScreenNb == i), ]

            temp1<-generateReplicateMat(subset, 2, "Intensities", colname4SD, col4anno)

            replicaMatrix<-temp1[[1]]
            indexPosControls<-temp1[[2]]
            indexNegControls<-temp1[[3]]

            sdVec<-rep(0, nrow(replicaMatrix))
            for (j in 1:nrow(replicaMatrix)){
                sdVec[j]<-sd(replicaMatrix[j, ], na.rm=T)
            }
            names(sdVec)<-rownames(replicaMatrix)

            RowNo<-c(ceiling(sqrt(length(sdVec))))
            ColNo<-ceiling(length(sdVec)/ceiling(sqrt(length(sdVec))))
    
            if (RowNo*ColNo > length(sdVec)){
                nbEmptyWells<-(RowNo*ColNo)-length(sdVec)
                sdVec<-c(sdVec, rep(as.double(NA_character_), nbEmptyWells))
            }
        
            posNegLabelVec<-rep(NA_character_, length(sdVec))
            posNegLabelVec[indexNegControls]<-"N"
            posNegLabelVec[indexPosControls]<-"P"

            if (showPlot == 1){
                if(interactive()){
                    x11()
                }
            }
            dings<-plotPlate(sdVec, ncol=ColNo, nrow=RowNo, char=posNegLabelVec, 
            cex.char=0.8, main=paste(plotTitle, ", Exp", i, sep=""), 
            col=brewer.pal(9, "YlOrBr"), cex.main=0.8, na.action='xout')

            headerTemp<-strsplit(header[1], ",")
            basicPlotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")

            png(paste(basicPlotName, "Exp", i, "4html.png", sep="_"), width=300, height=300)
                dings<-plotPlate(sdVec, ncol=ColNo, nrow=RowNo, char=posNegLabelVec, 
                cex.char=0.8, main=paste(plotTitle, ", Exp", i, sep=""), 
                col=brewer.pal(9, "YlOrBr"), cex.main=0.8, na.action='xout')
            dev.off()

            png(paste(basicPlotName, "Exp", i, ".png", sep="_"), width=dings$width, height=dings$height)
                dings<-plotPlate(sdVec, ncol=ColNo, nrow=RowNo, char=posNegLabelVec, 
                cex.char=0.8, main=paste(plotTitle, ", Exp", i, sep=""), 
                col=brewer.pal(9, "YlOrBr"), cex.main=0.8, na.action='xout')
            dev.off()

            title<-names(sdVec)
            con<-file(paste(basicPlotName, "Exp", i, ".html", sep="_"), open="w")
            imageMap(dings$coord, con, list(TITLE=title, TARGET=rep("main", 
            nrow(dings$coord))), paste(basicPlotName, "Exp", i, ".png", sep="_"))
            
            close(con)
        }
    }
    invisible(list(basicPlotName, minOfScreens, numOfScreens))
}




controlDensity<-function(header, dataset, channel, plotTitle, showPlot){

    subsetPosControls<-dataset[which(dataset$SpotType == 1), ]
    subsetNegControls<-dataset[which(dataset$SpotType == 0), ]

    posDens<-density(subsetPosControls[[get("channel")]])
    negDens<-density(subsetNegControls[[get("channel")]])

    minxlim=min(c(posDens$x, negDens$x))
    maxxlim=max(c(posDens$x, negDens$x))
    
    minylim=min(c(posDens$y, negDens$y))
    maxylim=max(c(posDens$y, negDens$y))
    
    if (showPlot == 1){
        if(interactive()){
            plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
            col="green", main=plotTitle)
        
            lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
            
            legend("topleft", c("Positive Controls", "Negative Controls"), 
            col=c("green", "red"), lty=1)
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    pdf(paste(plotName, ".pdf", sep=""))
        plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
        col="green", main=plotTitle)
        
        lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
        
        legend("topleft", c("Positive Controls", "Negative Controls"), 
        col=c("green", "red"), lty=1)
    dev.off()
    
    png(paste(plotName, ".png", sep=""), width=300, height=300)
        plot(posDens, xlim=c(minxlim,maxxlim), ylim=c(minylim,maxylim), 
        col="green", main=plotTitle, cex.lab=0.7, cex.axis=0.7)
        
        lines(negDens, xlim=c(minxlim,maxxlim), ylim=c(minylim,maxylim), col="red")
        
        legend("topleft", c("Positive Controls", "Negative Controls"), 
        col=c("green","red"), lty=1, bty="n", cex=0.4)
    dev.off()
    invisible(plotName)
}



controlDensityPerScreen<-function(header, dataset, channel, plotTitle, 
showPlot){

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)
    
    if (showPlot == 1){
        if(interactive()){
            c1<-ceiling(sqrt(length(unique(dataset$ScreenNb))))
            c2<-ceiling((length(unique(dataset$ScreenNb)))/c1)
            par(mfrow=c(c1, c2))

            for (i in minOfScreens:numOfScreens){
                if (length(which(dataset$ScreenNb == i))>0){
        
                    subset<-dataset[which(dataset$ScreenNb == i), ]

                    subsetPosControls<-subset[which(subset$SpotType == 1), ]
                    subsetNegControls<-subset[which(subset$SpotType == 0), ]

                    posDens<-density(subsetPosControls[[get("channel")]])
                    negDens<-density(subsetNegControls[[get("channel")]])

                    minxlim=min(c(posDens$x, negDens$x))
                    maxxlim=max(c(posDens$x, negDens$x))
    
                    minylim=min(c(posDens$y, negDens$y))
                    maxylim=max(c(posDens$y, negDens$y))
    
                    plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
                    col="green", main=plotTitle)
                
                    lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                
                    legend("topleft", c("Positive Controls", "Negative Controls"), 
                    col=c("green", "red"), lty=1)
                }
            }
        }
    }

    headerTemp<-strsplit(header[1], ",")
    plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
    
    for (i in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == i))>0){
            
            subset<-dataset[which(dataset$ScreenNb == i), ]

            subsetPosControls<-subset[which(subset$SpotType == 1), ]
            subsetNegControls<-subset[which(subset$SpotType == 0), ]

            posDens<-density(subsetPosControls[[get("channel")]])
            negDens<-density(subsetNegControls[[get("channel")]])

            minxlim=min(c(posDens$x, negDens$x))
            maxxlim=max(c(posDens$x, negDens$x))
    
            minylim=min(c(posDens$y, negDens$y))
            maxylim=max(c(posDens$y, negDens$y))
    
            pdf(paste(plotName, "(Exp. ", i, ").pdf", sep=""))
                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="green", 
                main=plotTitle)
                
                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                
                legend("topleft", c("Positive Controls", "Negative Controls"), 
                col=c("green", "red"), lty=1)
            dev.off()
            
            png(paste(plotName, "(Exp. ", i, ").png", sep=""), width=300, height=300)
                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="green", 
                main=plotTitle, cex.lab=0.7, cex.axis=0.7)
                
                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                
                legend("topleft", c("Positive Controls", "Negative Controls"), 
                col=c("green", "red"), lty=1, bty="n", cex=0.4)
            dev.off()            
        }
    }
    invisible(list(plotName, minOfScreens, numOfScreens))
}


controlDensityPerPlate<-function(header, dataset, channel, plotTitle, 
plotDesign, showPlot){

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    headerTemp<-strsplit(header[1], ",")

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
    
            subset<-dataset[which(dataset$ScreenNb == j), ]
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
            
            if (showPlot == 1){
                if(interactive()){
            
                    if (plotDesign == 1){
                        x11()
                        c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                        c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                        par(mfrow=c(c1, c2))
                    }

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                
                            subsubset<-subset[which(subset$LabtekNb == i), ]

                            subsetPosControls<-subsubset[which(subsubset$SpotType == 1), ]
                            subsetNegControls<-subsubset[which(subsubset$SpotType == 0), ]

                            if(length(subsetPosControls[[get("channel")]])>1 
                            & length(subsetNegControls[[get("channel")]])>1){
                        
                                posDens<-density(subsetPosControls[[get("channel")]])
                                negDens<-density(subsetNegControls[[get("channel")]])
    
                                minxlim=min(c(posDens$x, negDens$x))
                                maxxlim=max(c(posDens$x, negDens$x))
    
                                minylim=min(c(posDens$y, negDens$y))
                                maxylim=max(c(posDens$y, negDens$y))
    
                                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
                                col="green", main=plotTitle)
                            
                                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                            
                                legend("topleft", c("Positive Controls", "Negative Controls"), 
                                col=c("green", "red"), lty=1)
                            }else{
                                print(paste("Too little controls on plate", i, "to compute density"), sep=" ")
                            }
                        }
                    }
                }
            }
        }
    }

    for (j in minOfScreens:numOfScreens){
        if (length(which(dataset$ScreenNb == j))>0){
    
            subset<-dataset[which(dataset$ScreenNb == j), ]
            minOfPlates<-min(subset$LabtekNb)
            numOfPlates<-max(subset$LabtekNb)
            
            plotName<-paste(headerTemp[[1]][2], plotTitle, sep="_")
            
            if (plotDesign == 1){
            
                pdf(paste(plotName, "_Exp_", j, "_.pdf", sep=""))
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2))

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){    ##labtek vorhanden
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]

                            subsetPosControls<-subsubset[which(subsubset$SpotType == 1), ]
                            subsetNegControls<-subsubset[which(subsubset$SpotType == 0), ]

                            if(length(subsetPosControls[[get("channel")]])>1 
                            & length(subsetNegControls[[get("channel")]])>1){
                            
                                posDens<-density(subsetPosControls[[get("channel")]])
                                negDens<-density(subsetNegControls[[get("channel")]])

                                minxlim=min(c(posDens$x, negDens$x))
                                maxxlim=max(c(posDens$x, negDens$x))
    
                                minylim=min(c(posDens$y, negDens$y))
                                maxylim=max(c(posDens$y, negDens$y))
    
                                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
                                col="green", main=plotTitle)
                                
                                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                                
                                legend("topleft", c("Positive Controls", "Negative Controls"), 
                                col=c("green", "red"), lty=1)
                            }else{
                                print(paste("Too little controls on plate", i, "to compute density"), sep=" ")
                            }
                        }
                    }
                dev.off()
                
                png(paste(plotName, "_Exp_", j, "_.png", sep=""))
                    
                    c1<-ceiling(sqrt(length(unique(subset$LabtekNb))))
                    c2<-ceiling(length(unique(subset$LabtekNb))/c1)
                    par(mfrow=c(c1, c2))

                    for (i in minOfPlates:numOfPlates){
                        if (length(which(subset$LabtekNb == i))>0){
                    
                            subsubset<-subset[which(subset$LabtekNb == i), ]

                            subsetPosControls<-subsubset[which(subsubset$SpotType == 1), ]
                            subsetNegControls<-subsubset[which(subsubset$SpotType == 0), ]

                            if(length(subsetPosControls[[get("channel")]])>1 
                            & length(subsetNegControls[[get("channel")]])>1){
                            
                                posDens<-density(subsetPosControls[[get("channel")]])
                                negDens<-density(subsetNegControls[[get("channel")]])

                                minxlim=min(c(posDens$x, negDens$x))
                                maxxlim=max(c(posDens$x, negDens$x))
    
                                minylim=min(c(posDens$y, negDens$y))
                                maxylim=max(c(posDens$y, negDens$y))
    
                                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
                                col="green", main=plotTitle)
                                
                                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                                
                                legend("topleft", c("Positive Controls", "Negative Controls"), 
                                col=c("green", "red"), lty=1)
                            }else{
                                print(paste("Too little controls on plate", i, "to compute density"), sep=" ")
                            }
                        }
                    }
                dev.off()

            }else{
            
                for (i in minOfPlates:numOfPlates){
                    if (length(which(subset$LabtekNb == i))>0){
                    
                        subsubset<-subset[which(subset$LabtekNb == i), ]

                        subsetPosControls<-subsubset[which(subsubset$SpotType == 1), ]
                        subsetNegControls<-subsubset[which(subsubset$SpotType == 0), ]

                        if(length(subsetPosControls[[get("channel")]])>1 
                        & length(subsetNegControls[[get("channel")]])>1){
                        
                            posDens<-density(subsetPosControls[[get("channel")]])
                            negDens<-density(subsetNegControls[[get("channel")]])

                            minxlim=min(c(posDens$x, negDens$x))
                            maxxlim=max(c(posDens$x, negDens$x))
    
                            minylim=min(c(posDens$y, negDens$y))
                            maxylim=max(c(posDens$y, negDens$y))
    
                            pdf(paste(plotName, "_Exp_", j, "_PerPlate_", i, "_.pdf", sep=""))
                                
                                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
                                col="green", main=plotTitle)
                                
                                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                                
                                legend("topleft", c("Positive Controls", "Negative Controls"), 
                                col=c("green", "red"), lty=1)
                            dev.off()
                            
                            png(paste(plotName, "_Exp_", j, "_PerPlate_", i, "_.png", sep=""), width=300, height=300)
                                
                                plot(posDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), 
                                col="green", main=plotTitle, cex.lab=0.7, cex.axis=0.7)
                                
                                lines(negDens, xlim=c(minxlim, maxxlim), ylim=c(minylim, maxylim), col="red")
                                
                                legend("topleft", c("Positive Controls", "Negative Controls"), 
                                col=c("green", "red"), lty=1, bty="n", cex=0.4)
                            dev.off()
                            
                        }else{
                            print(paste("Too little controls on plate", i, "to compute density"), sep=" ")
                        }
                    }
                }
            }
        }
    }
    invisible(list(plotName, c(minOfScreens, numOfScreens), 
    c(minOfPlates, numOfPlates)))
}
