mainAnalysis <- function(header, dataset, flagForSameExp, listOfNormalizations, 
listOfArgs4norm, listOfStatTests, listOfArgs4stat, multTestAdj, hitScoringVec1, 
hitScoringVec2, posNegFlag, flag4Gsea, vecOfChannels, whichOnto){

    firstLine<-strsplit(header[[1]], ",")
    firstLineBis<-firstLine[[1]][2]
    secondLine<-strsplit(header[[2]], ",")
    secondLineBis<-secondLine[[1]][2]
    thirdLine<-strsplit(header[[3]], "comments,")
    thirdLineBis<-thirdLine[[1]][2]
    

    #################################
    ##########QUALITY PLOTS##########
    #################################
    
    print(paste("Computing quality plots on raw data...", date(), sep=": "))
    
    ##MAINPAGE titles:
    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file="index.html", quote=F, row.names=F, col.names=F)
    
    write.table(paste("<BODY><CENTER><H1>", firstLineBis, "</H1></CENTER>", sep=""), 
    file="index.html", append=T, quote=F, row.names=F, col.names=F)


    ##SUBPAGE1 titles:
    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file="sub1.html", quote=F, row.names=F, col.names=F)
    
    write.table(paste("<BODY><CENTER><H2>Boxplots per Plate</H2></CENTER>", sep=""), 
    file="sub1.html", append=T, quote=F, row.names=F, col.names=F)

    print(paste("Generating boxplots of signals...", date(), sep=": "))
    function1(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
    "index.html", "sub1.html", 0)
    
    print(paste("Generating spatial distribution of signals...", date(), sep=": "))
    function2(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
    "index.html", "sub2.html", 0)
    
    print(paste("Plotting signals per well...", date(), sep=": "))
    function2b(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
    "index.html", "sub2b.html", 0)
    
    if (length(vecOfChannels)>1){
        print(paste("Comparing channels...", date(), sep=": "))
        function2c(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
        "index.html", "sub2c.html", 0)
    }
    
    print(paste("Plotting histogram distribution of signals...", date(), sep=": "))
    function3(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
    secondLineBis, "index.html", "sub3.html", 0)
    
    print(paste("Plotting control densities...", date(), sep=": "))
    function4a(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
    secondLineBis, "index.html", "sub4a.html", 0)
    
    print(paste("Computing Z'factors...", date(), sep=": "))
    function4b(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
    secondLineBis, "index.html", "sub4b.html", 0)
    
    print(paste("Generating boxplots of control signals...", date(), sep=": "))    
    function5(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
    secondLineBis, "index.html", "sub5.html", 0)
    
    print(paste("Generating QQ plots...", date(), sep=": "))
    function6(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
    "index.html", "sub6.html", 0)


    write.table("</BODY>", file="index.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub1.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub2.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub2b.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub2c.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub3.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub4a.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub4b.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub5.html", append=T, quote=F, row.names=F, 
    col.names=F)
    write.table("</BODY>", file="sub6.html", append=T, quote=F, row.names=F, 
    col.names=F)


    ################################
    ##########NORMALIZATION#########
    ################################
    
    print(paste("Normalizing data...", date(), sep=": "))
    
    for (i in 1:length(listOfNormalizations)){
    
        result<-listOfNormalizations[[i]](header, dataset, listOfArgs4norm[[i]])
        header<-result[[1]]
        dataset<-result[[2]]
        
        ##change Inf and -Inf values that potentially arise during normalization to NA
        ##it is IMPORTANT for this that the argument order for normalization functions stays the same!
        col4plot<-listOfArgs4norm[[i]][[1]][1]
        if(length(which(dataset[[get("col4plot")]]==Inf))>0){
            places<-which(dataset[[get("col4plot")]]==Inf)
            dataset[[get("col4plot")]][places]=NA
        }
                        
        if(length(which(dataset[[get("col4plot")]]==-Inf))>0){
	    places<-which(dataset[[get("col4plot")]]==-Inf)
	    dataset[[get("col4plot")]][places]=NA
        }
    
        saveDataset(header, dataset, paste(firstLineBis, "_Fulldataset.txt", sep=""))
    
        #################################
        ##########QUALITY PLOTS##########
        #################################
    
        print("Computing quality plots on normalized data...")
    
        thirdLine<-strsplit(header[[3]], "comments,")
        thirdLineBis<-thirdLine[[1]][2]
    
        ##MAINPAGE titles:
        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=paste("indexnorm", i, ".html", sep=""), 
        quote=F, row.names=F, col.names=F)
        
        write.table(paste("<BODY><CENTER><H1>", firstLineBis, " after normalization", 
        "</H1></CENTER>", sep=""), file=paste("indexnorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        
        write.table(paste("<CENTER>(", thirdLineBis, ")</CENTER>", sep=""), 
        file=paste("indexnorm", i, ".html", sep=""), append=T, quote=F, 
        row.names=F, col.names=F)
    
    
        ##SUBPAGE1 titles:
        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=paste("sub1norm", i, ".html", sep=""), 
        quote=F, row.names=F, col.names=F)
        
        write.table(paste("<BODY><CENTER><H1>Boxplots per Plate</H1></CENTER>", sep=""), 
        file=paste("sub1norm", i, ".html", sep=""), append=T, quote=F, row.names=F, 
        col.names=F)

    
        print(paste("Generating boxplots of signals...", date(), sep=": "))
        function1(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
        paste("indexnorm", i, ".html", sep=""), paste("sub1norm", i, ".html", sep=""), 
        i)
        
        print(paste("Generating spatial distribution of signals...", date(), sep=": "))
        function2(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
        paste("indexnorm", i, ".html", sep=""), paste("sub2norm", i, ".html", sep=""), 
        i)
        
        print(paste("Plotting signals per well...", date(), sep=": "))
        function2b(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
        paste("indexnorm", i, ".html", sep=""), paste("sub2bnorm", i, ".html", sep=""), 
        i)
        
        if (length(vecOfChannels)>1){
            print(paste("Comparing channels...", date(), sep=": "))
            function2c(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
            paste("indexnorm", i, ".html", sep=""), paste("sub2cnorm", i, ".html", sep=""), 
            i)
        }
        
        print(paste("Plotting histogram distribution of signals...", date(), sep=": "))
        function3(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
        secondLineBis, paste("indexnorm", i, ".html", sep=""), 
        paste("sub3norm", i, ".html", sep=""), i)
        
        print(paste("Plotting control densities...", date(), sep=": "))
        function4a(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
        secondLineBis, paste("indexnorm", i, ".html", sep=""), 
        paste("sub4anorm", i, ".html", sep=""), i)
        
        print(paste("Computing Z'factors...", date(), sep=": "))
        function4b(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
        secondLineBis, paste("indexnorm", i, ".html", sep=""), 
        paste("sub4bnorm", i, ".html", sep=""), i)    
        
        print(paste("Generating boxplots of control signals...", date(), sep=": "))    
        function5(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
        secondLineBis, paste("indexnorm", i, ".html", sep=""), 
        paste("sub5norm", i, ".html", sep=""), i)
        
        print(paste("Generating QQ plots...", date(), sep=": "))
        function6(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
        paste("indexnorm", i, ".html", sep=""), paste("sub6norm", i, ".html", sep=""), 
        i)
        
        print(paste("Computing Spearman's correlation coefficients...", date(), 
        sep=": "))
        function7(header, dataset, vecOfChannels, firstLineBis, secondLineBis, 
        paste("indexnorm", i, ".html", sep=""), paste("sub7anorm", i, ".html", sep=""), 
        paste("sub7bnorm", i, ".html", sep=""), flagForSameExp, i)
    
        write.table("</BODY>", file=paste("indexnorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub1norm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub2norm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub2bnorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub2cnorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub3norm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub4anorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub4bnorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub5norm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub6norm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub7anorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        write.table("</BODY>", file=paste("sub7bnorm", i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)

    }

    ##############################
    ##########STAT TESTS##########
    ##############################
    
    print(paste("Performing statistical tests...", date(), sep=": "))
    
    ##MAINPAGE titles:
    write.table(paste("<HTML><HEAD><TITLE>Statistical analysis report, ", 
    secondLineBis, ", ", firstLineBis, "</TITLE></HEAD>", sep=""), 
    file="stats.html", quote=F, row.names=F, col.names=F)

    ##MAINPAGE Seitentitel:
    write.table(paste("<BODY><CENTER><H1>Statistical analysis report, ", 
    firstLineBis, "</H1></CENTER>", sep=""), file="stats.html", append=T, 
    quote=F, row.names=F, col.names=F)

    write.table(paste("<CENTER>(Normalizations performed: ", thirdLineBis, 
    ")</CENTER>", sep=""), file="stats.html", append=T, quote=F, 
    row.names=F, col.names=F)

    ##Apply statistical tests and score hits:
    hitVectors<-list("platzhalter")
    hitVectorsTestNames<-list("platzhalter")
    counterHitVectorsTestnames<-0
    
    for (i in 1:length(listOfStatTests)){


        testResult<-listOfStatTests[[i]](dataset, listOfArgs4stat[[i]])
        pValVec<-testResult[[1]]
        dataset<-testResult[[2]]
        namepValVec<-testResult[[3]]
        
        saveDataset(header, dataset, paste(firstLineBis, "_Fulldataset.txt", sep=""))
        
        ##adjust for multiple testing and reinsert pValVec into dataset:
        pValVec<-multTestAdjust(pValVec, multTestAdj)
        temp<-generateReplicateMat(dataset, 1, "Intensities", listOfArgs4stat[[i]][[3]], 
        listOfArgs4stat[[i]][[4]])        
        
        dataset<-incorporatepValVec(dataset, pValVec, temp[[1]], 
        listOfArgs4stat[[i]][[4]], paste(namepValVec, "_adjusted_", multTestAdj, 
        sep=""))
        
        print(paste("Scoring hits...", date(), sep=": "))

        result<-hitselectionPval(dataset, pValVec, listOfArgs4stat[[i]][[3]], 
        paste("hits_", namepValVec, "_adjusted_", multTestAdj, sep=""), 
        hitScoringVec2[1], listOfArgs4stat[[i]][[4]], paste("pval_hits_", 
        testResult[[4]], " (", i, ").txt", sep=""))

        dataset<-result[[1]]
        hitVectors[[i]]<-result[[2]]
        hitVectorsTestNames[[i]]<-testResult[[4]]
        hitThreshold<-round(result[[4]], digits=2)
        

        if (hitScoringVec1[1] == 1){
        ##scoring according to only p-value

            genesnPvals<-read.table(paste("pval_hits_", testResult[[4]], " (", i, 
            ").txt", sep=""))

            if (hitVectorsTestNames[[i]] == "Mann-Whitney test" 
            | hitVectorsTestNames[[i]] == "t test"){    
            
                if (listOfArgs4stat[[i]][[1]] == "g"){
                    upOrDownReg<-"Upregulated"
                }
                if (listOfArgs4stat[[i]][[1]] == "l"){
                    upOrDownReg<-"Downregulated"
                }
                if (listOfArgs4stat[[i]][[1]] == "two.sided"){
                    upOrDownReg<-"Up and downregulated"
                }
            }
            if (hitVectorsTestNames[[i]] == "Rank product test"){
                if (listOfArgs4stat[[i]][[2]] == 1){
                    upOrDownReg<-"Downregulated"
                }
                if (listOfArgs4stat[[i]][[2]] == 2){
                    upOrDownReg<-"Upregulated"
                }
            }


            write.table(paste("<br><br><br><CENTER><H2>", upOrDownReg, 
            " genes according to ", hitVectorsTestNames[[i]], 
            "</H2> (multiple testing adjument: ", multTestAdj, ") (p-value <  ", 
            hitThreshold, ") - <a href=\"", "pval_hits_", testResult[[4]], " (", i, ").txt", 
            "\">Textfile</a> -</CENTER>", sep=""), file="stats.html", append=T, quote=F, 
            row.names=F, col.names=F)


            s1<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">Gene name</TH>"
            s2<-"<TH BGCOLOR=\"#e0e0ff\">p-value</TH><TH BGCOLOR=\"#d0d0f0\">Median</TH>"
            s3<-paste(s1, s2, sep="")
            write.table(s3, file="stats.html", append=T, quote=F, row.names=F, col.names=F)

            for (m in 1:(ncol(genesnPvals)-2)){
                write.table("<TH BGCOLOR=\"#d0d0f0\">Norm val</TH>", file="stats.html", 
                append=T, quote=F, row.names=F, col.names=F)
            }
            
            write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
            col.names=F)

            for (n in 1:nrow(genesnPvals)){

                write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", rownames(genesnPvals)[n], 
                "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)

                for (p in 1:ncol(genesnPvals)){

                    write.table(paste("<TD BGCOLOR=\"#f0f0ff\" align=center>", genesnPvals[n, p], 
                    "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                    col.names=F)
                }
                write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)
            }
            write.table("</CENTER></TABLE>", file="stats.html", append=T, quote=F, 
            row.names=F, col.names=F)


            ##VOLCANO PLOT###
            volcPlotName<-volcanoPlot(header, dataset, listOfArgs4stat[[i]][[3]], 
            paste(namepValVec, "_adjusted_", multTestAdj, sep=""), 
            listOfArgs4stat[[i]][[4]], paste("Volcano plot (", i, ")_p", sep=""), 
            hitScoringVec2[1], 0)

            write.table("<br><br><CENTER>Volcano plot</CENTER><br>", file="stats.html", 
            append=T, quote=F, row.names=F, col.names=F)
            
            write.table(paste("<CENTER><A HREF=\"", volcPlotName, ".pdf\"><IMG SRC=\"", 
            volcPlotName, ".png\"/></A></CENTER>", sep=""), file="stats.html", append=T, 
            quote=F, row.names=F, col.names=F)


            ##SPATIAL DISTRIBUTION OF HITS
            res<-spatialDistribHits(header, dataset, paste("Spatial hits distribution (", 
            i, ")_p", sep=""), paste("hits_", namepValVec, "_adjusted_", multTestAdj, 
            sep=""), listOfArgs4stat[[i]][[4]], 0)


            write.table("<BR><BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\"></TH>", 
            file="stats.html", append=T, quote=F, row.names=F, col.names=F)

            for (j in res[[3]][1]:res[[3]][2]){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)
            }
            write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, col.names=F)

            for (k in res[[2]][1]:res[[2]][2]){
                write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Exp. ", k, "</TH>", sep=""), 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                for (m in res[[3]][1]:res[[3]][2]){

                    write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", res[[1]], 
                    "_Exp_", k, "_Plate_", m, "_.html", "\"><IMG SRC=\"", res[[1]], "_Exp_", k, 
                    "_Plate_", m, "_4html.png\"/></A></CENTER><BR></TD>", sep=""), 
                    file="stats.html", append=T, quote=F, row.names=F, col.names=F)
                }
            }
            write.table("</TR></CENTER></TABLE>", file="stats.html", append=T, quote=F, 
            row.names=F, col.names=F)


        }

        if (hitScoringVec1[2] == 2 ){
        ##hitscoring according to zSc < thresh & pval < thresh in hitScoringVec2[1]
        
            if (((hitVectorsTestNames[[i]] == "Mann-Whitney test" 
            | hitVectorsTestNames[[i]] == "t test") 
            & (listOfArgs4stat[[i]][[1]] == "two.sided" 
            | listOfArgs4stat[[i]][[1]] == "l")) 
            | (hitVectorsTestNames[[i]] == "Rank product test" 
            & listOfArgs4stat[[i]][[2]] == 1)){
        
                result<-hitselectionZscorePval(dataset, pValVec, listOfArgs4stat[[i]][[3]], 
                paste("hits_", namepValVec, "_adjusted_", multTestAdj, "_ZScore<", 
                hitScoringVec2[2], sep=""), hitScoringVec2[2], hitScoringVec2[1], -2, 
                listOfArgs4stat[[i]][[4]], median, paste("pval_hits_", testResult[[4]], 
                "combined_ZSc_sm_", hitScoringVec2[2], ".txt", sep=""))
                
                
                dataset<-result[[1]]
                counterHitVectorsTestnames<-counterHitVectorsTestnames+1
                hitVectors[[length(listOfStatTests)+counterHitVectorsTestnames]]<-result[[2]]

                hitThresholdZsc<-round(result[[3]], digits=3)
                hitThresholdPval<-round(result[[4]], digits=2)
                
                s3<-paste(testResult[[4]], " & ZScore sm ", hitThresholdZsc, sep="")
                hitVectorsTestNames[[length(listOfStatTests)+counterHitVectorsTestnames]]<-s3

                genesnPvals<-read.table(paste("pval_hits_", testResult[[4]], "combined_ZSc_sm_", 
                hitScoringVec2[2], ".txt", sep=""))
            
                upOrDownReg<-"Downregulated"

            
            
                write.table(paste("<br><br><br><CENTER><H2>", upOrDownReg, 
                " genes according to ", hitVectorsTestNames[[length(listOfStatTests)
                +counterHitVectorsTestnames]], "</H2> (multiple testing adjument: ", 
                multTestAdj, ") (Thresholds: pval < ", hitThresholdPval, ", Z-score < ", 
                hitThresholdZsc, ") - <a href=\"", "pval_hits_", testResult[[4]], 
                "combined_ZSc_sm_", hitScoringVec2[2], ".txt", "\">Textfile</a> -</CENTER>", 
                sep=""), file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                s5<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">Gene name</TH>"
                s6<-"<TH BGCOLOR=\"#e0e0ff\">sumZScore</TH><TH BGCOLOR=\"#e0e0ff\">p-value</TH>"
                s7<-paste(s5, s6, sep="")
                write.table(s7, file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                for (m in 1:(ncol(genesnPvals)-2)){
                    write.table("<TH BGCOLOR=\"#d0d0f0\">Norm val</TH>", file="stats.html", 
                    append=T, quote=F, row.names=F, col.names=F)
                }
                
                write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)
                
                
                for (n in 1:nrow(genesnPvals)){

                    write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", rownames(genesnPvals)[n], 
                    "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                    col.names=F)

                    for (p in 1:ncol(genesnPvals)){

                        write.table(paste("<TD BGCOLOR=\"#f0f0ff\" align=center>", 
                        genesnPvals[n, p], "</TD>", sep=""), file="stats.html", append=T, 
                        quote=F, row.names=F, col.names=F)
                    }
                    write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                    col.names=F)
                }
                write.table("</CENTER></TABLE>", file="stats.html", append=T, 
                quote=F, row.names=F, col.names=F)
                
                ##VOLCANO PLOT###

                if (hitScoringVec1[3] == 2){
                    volcPlotName<-volcanoPlot(header, dataset, listOfArgs4stat[[i]][[3]], 
                    paste(namepValVec, "_adjusted_", multTestAdj, sep=""), 
                    listOfArgs4stat[[i]][[4]], paste("Volcano plot (", i, ")_zl", sep=""), 
                    c(hitScoringVec2[1], hitScoringVec2[2], hitScoringVec2[3]), 0)
                    
                }else{
                    volcPlotName<-volcanoPlot(header, dataset, listOfArgs4stat[[i]][[3]], 
                    paste(namepValVec, "_adjusted_", multTestAdj, sep=""), 
                    listOfArgs4stat[[i]][[4]], paste("Volcano plot (", i, ")_zl", sep=""), 
                    c(hitScoringVec2[1], hitScoringVec2[2]), 0)
                }

                write.table("<br><br><CENTER>Volcano plot</CENTER><br>", file="stats.html", 
                append=T, quote=F, row.names=F, col.names=F)
                
                write.table(paste("<CENTER><A HREF=\"", volcPlotName, ".pdf\"><IMG SRC=\"", 
                volcPlotName, ".png\"/></A></CENTER>", sep=""), file="stats.html", 
                append=T, quote=F, row.names=F, col.names=F)
                
                res<-spatialDistribHits(header, dataset, paste("Spatial hits distribution (", 
                i, ")_zl", sep=""), paste("hits_", namepValVec, "_adjusted_", multTestAdj, 
                "_ZScore<", hitScoringVec2[2], sep=""), listOfArgs4stat[[i]][[4]], 0)
                
                
                write.table("<BR><BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\"></TH>", 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                for (j in res[[3]][1]:res[[3]][2]){
                    write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                    file="stats.html", append=T, quote=F, row.names=F, col.names=F)
                }
                
                write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)

                for (k in res[[2]][1]:res[[2]][2]){
                    write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Exp. ", k, "</TH>", sep=""), 
                    file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                    for (m in res[[3]][1]:res[[3]][2]){

                        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", res[[1]], 
                        "_Exp_", k, "_Plate_", m, "_.html", "\"><IMG SRC=\"", res[[1]], "_Exp_", k, 
                        "_Plate_", m, "_4html.png\"/></A></CENTER><BR></TD>", sep=""), 
                        file="stats.html", append=T, quote=F, row.names=F, col.names=F)
                    }
                }
                write.table("</TR></CENTER></TABLE>", file="stats.html", append=T, 
                quote=F, row.names=F, col.names=F)
            }
        }

        if (hitScoringVec1[3] == 2){    
        ##hitscoring according to zSc > thresh & pval < thresh in hitScoringVec2[1]

            if (((hitVectorsTestNames[[i]] == "Mann-Whitney test" 
            | hitVectorsTestNames[[i]] == "t test") 
            & (listOfArgs4stat[[i]][[1]] == "two.sided" 
            | listOfArgs4stat[[i]][[1]] == "g")) 
            | (hitVectorsTestNames[[i]] == "Rank product test" 
            & listOfArgs4stat[[i]][[2]] == 2)){
            
                result<-hitselectionZscorePval(dataset, pValVec, listOfArgs4stat[[i]][[3]], 
                paste("hits_", namepValVec, "_adjusted_", multTestAdj, "_ZScore>", 
                hitScoringVec2[3], sep=""), hitScoringVec2[3], hitScoringVec2[1], 2, 
                listOfArgs4stat[[i]][[4]], median, paste("pval_hits_", testResult[[4]], 
                "combined_ZSc_gr_", hitScoringVec2[3], ".txt", sep=""))

                dataset<-result[[1]]
                counterHitVectorsTestnames<-counterHitVectorsTestnames+1
                hitVectors[[length(listOfStatTests)+counterHitVectorsTestnames]]<-result[[2]]

                hitThresholdZsc<-round(result[[3]], digits=3)
                hitThresholdPval<-round(result[[4]], digits=2)

                hitVectorsTestNames[[length(listOfStatTests)+counterHitVectorsTestnames]]<-
                paste(testResult[[4]], " & ZScore gr ", hitThresholdZsc, sep="")

                genesnPvals<-read.table(paste("pval_hits_", testResult[[4]], "combined_ZSc_gr_", 
                hitScoringVec2[3], ".txt", sep=""))

                upOrDownReg<-"Upregulated"
                
                write.table(paste("<br><br><br><CENTER><H2>", upOrDownReg, 
                " genes according to ", hitVectorsTestNames[[length(listOfStatTests)
                +counterHitVectorsTestnames]], "</H2> (multiple testing adjument: ", 
                multTestAdj, ") (Thresholds: p-value < ", hitThresholdPval, ",  Z-score > ", 
                hitThresholdZsc, ") - <a href=\"", "pval_hits_", testResult[[4]], 
                "combined_ZSc_gr_", hitScoringVec2[3], ".txt", "\">Textfile</a> -</CENTER>", 
                sep=""), file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                s8<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">Gene name</TH>"
                s9<-"<TH BGCOLOR=\"#e0e0ff\">sumZScore</TH><TH BGCOLOR=\"#e0e0ff\">p-value</TH>"
                s1<-paste(s8, s9, sep="")
                write.table(s1, file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                for (m in 1:(ncol(genesnPvals)-2)){
                    write.table("<TH BGCOLOR=\"#d0d0f0\">Norm val</TH>", file="stats.html", 
                    append=T, quote=F, row.names=F, col.names=F)
                }
                
                write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)
                
                for (n in 1:nrow(genesnPvals)){

                    write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", rownames(genesnPvals)[n], 
                    "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                    col.names=F)

                    for (p in 1:ncol(genesnPvals)){

                        write.table(paste("<TD BGCOLOR=\"#f0f0ff\" align=center>", genesnPvals[n, p], 
                        "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                        col.names=F)
                    }
                    write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                    col.names=F)
                }
                write.table("</CENTER></TABLE>", file="stats.html", append=T, quote=F, 
                row.names=F, col.names=F)
                
                ##VOLCANO PLOT###

                if (hitScoringVec1[3] == 2){
                    volcPlotName<-volcanoPlot(header, dataset, listOfArgs4stat[[i]][[3]], 
                    paste(namepValVec, "_adjusted_", multTestAdj, sep=""), 
                    listOfArgs4stat[[i]][[4]], paste("Volcano plot (", i, ")_zg", sep=""), 
                    c(hitScoringVec2[1], hitScoringVec2[2], hitScoringVec2[3]), 0)
                    
                }else{
                    volcPlotName<-volcanoPlot(header, dataset, listOfArgs4stat[[i]][[3]], 
                    paste(namepValVec, "_adjusted_", multTestAdj, sep=""), 
                    listOfArgs4stat[[i]][[4]], paste("Volcano plot (", i, ")_zg", sep=""), 
                    c(hitScoringVec2[1], hitScoringVec2[2]), 0)
                }

                write.table("<br><br><CENTER>Volcano plot</CENTER><br>", file="stats.html",
                append=T, quote=F, row.names=F, col.names=F)
                
                write.table(paste("<CENTER><A HREF=\"", volcPlotName, ".pdf\"><IMG SRC=\"", 
                volcPlotName, ".png\"/></A></CENTER>", sep=""), file="stats.html", append=T, 
                quote=F, row.names=F, col.names=F)
                
                res<-spatialDistribHits(header, dataset, paste("Spatial hits distribution (", 
                i, ")_zg", sep=""), paste("hits_", namepValVec, "_adjusted_", multTestAdj, 
                "_ZScore>", hitScoringVec2[3], sep=""), listOfArgs4stat[[i]][[4]], 0)
                
                write.table("<BR><BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\"></TH>", 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                for (j in res[[3]][1]:res[[3]][2]){
                    write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                    file="stats.html", append=T, quote=F, row.names=F, col.names=F)
                }
                
                write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)

                for (k in res[[2]][1]:res[[2]][2]){
                    write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Exp. ", k, "</TH>", sep=""), 
                    file="stats.html", append=T, quote=F, row.names=F, col.names=F)

                    for (m in res[[3]][1]:res[[3]][2]){

                        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", res[[1]], 
                        "_Exp_", k, "_Plate_", m, "_.html", "\"><IMG SRC=\"", res[[1]], "_Exp_", k, 
                        "_Plate_", m, "_4html.png\"/></A></CENTER><BR></TD>", sep=""), 
                        file="stats.html", append=T, quote=F, row.names=F, col.names=F)
                    }
                }
                write.table("</TR></CENTER></TABLE>", file="stats.html", append=T, 
                quote=F, row.names=F, col.names=F)
            }
        }        
    }



    if (hitScoringVec1[2] == 1){    
    ##scoring according to ZScore with ZScore <= threshold
    
        result<-hitselectionZscore(dataset, "SigIntensity", "hits_ZScore_do", 
        hitScoringVec2[2], 2, -2, "GeneName", median, "ZScore_hits_smallerthresh.txt")
        
        dataset<-result[[1]]

        count<-length(hitVectors)
        hitVectors[[count+1]]<-result[[2]]
        hitVectorsTestNames[[count+1]]<-"ZScore_do"

        hitThreshold<-round(result[[3]], digits=3)
        
        res<-spatialDistribHits(header, dataset, "Spatial hits distribution (ZScore)", 
        "hits_ZScore_do", "GeneName", 0)
        
        genesnPvals<-read.table("ZScore_hits_smallerthresh.txt")

        write.table(paste("<br><br><br><CENTER><H2>Downregulated genes according to ", 
        hitVectorsTestNames[[count+1]], "</H2> (Threshold: Z-score < ", hitThreshold, 
        ") - <a href=\"ZScore_hits_smallerthresh.txt\">Textfile</a> -</CENTER>", 
        sep=""), file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        
        s2<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">Gene name</TH>"
        s3<-"<TH BGCOLOR=\"#e0e0ff\">ZScore (median)</TH>"
        s4<-paste(s2, s3, sep="")
        write.table(s4, file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        
        for (m in 1:(ncol(genesnPvals)-1)){
            write.table("<TH BGCOLOR=\"#d0d0f0\">Norm val</TH>", file="stats.html", 
            append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (n in 1:nrow(genesnPvals)){
        
            write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", rownames(genesnPvals)[n], 
            "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
            col.names=F)
            
            for (p in 1:ncol(genesnPvals)){
            
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\" align=center>", genesnPvals[n, p], 
                "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)
            }
            write.table("</TR>", file="stats.html", append=T, quote=F, 
            row.names=F, col.names=F)
        }
        write.table("</CENTER></TABLE>", file="stats.html", append=T, quote=F, 
        row.names=F, col.names=F)
        
        
        ##SPATIAL DISTRIBUTION OF HITS - output
        write.table("<BR><BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\"></TH>", 
        file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        
        for (j in res[[3]][1]:res[[3]][2]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
            file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (k in res[[2]][1]:res[[2]][2]){
            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Exp. ", k, "</TH>", sep=""), 
            file="stats.html", append=T, quote=F, row.names=F, col.names=F)

            for (m in res[[3]][1]:res[[3]][2]){
            
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", res[[1]], 
                "_Exp_", k, "_Plate_", m, "_.html", "\"><IMG SRC=\"", res[[1]], "_Exp_", k, 
                "_Plate_", m, "_4html.png\"/></A></CENTER><BR></TD>", sep=""), 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)
            }
        }
        write.table("</TR></CENTER></TABLE>", file="stats.html", append=T, quote=F, 
        row.names=F, col.names=F)
    
    }
       
    
    if (hitScoringVec1[3] == 1){
    
        result<-hitselectionZscore(dataset, "SigIntensity", "hits_ZScore_up", 
        hitScoringVec2[3], 2, 2, "GeneName", median, "ZScore_hits_higherthresh.txt")
        
        dataset<-result[[1]]
        
        count<-length(hitVectors)
        hitVectors[[count+1]]<-result[[2]]
        hitVectorsTestNames[[count+1]]<-"ZScore_up"
        
        hitThreshold<-round(result[[3]], digits=3)
        
        res<-spatialDistribHits(header, dataset, "Spatial hits distribution (ZScore)", 
        "hits_ZScore_up", "GeneName", 0)

        genesnPvals<-read.table("ZScore_hits_higherthresh.txt")

        write.table(paste("<br><br><br><CENTER><H2>Upregulated genes according to ", 
        hitVectorsTestNames[[count+1]], "</H2> (Threshold: Z-score >", hitThreshold, 
        ") - <a href=\"ZScore_hits_higherthresh.txt\">Textfile</a> -</CENTER>", sep=""), 
        file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        
        s6<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">Gene name</TH>"
        s7<-"<TH BGCOLOR=\"#e0e0ff\">ZScore (median)</TH>"
        s8<-paste(s6, s7, sep="")
        write.table(s8, file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        
        for (m in 1:(ncol(genesnPvals)-1)){
            write.table("<TH BGCOLOR=\"#d0d0f0\">Norm val</TH>", file="stats.html", 
            append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (n in 1:nrow(genesnPvals)){
        
            write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", rownames(genesnPvals)[n], 
            "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
            col.names=F)
            
            for (p in 1:ncol(genesnPvals)){
            
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\" align=center>", genesnPvals[n, p], 
                "</TD>", sep=""), file="stats.html", append=T, quote=F, row.names=F, 
                col.names=F)
            }
            write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
            col.names=F)
        }
        write.table("</CENTER></TABLE>", file="stats.html", append=T, quote=F, 
        row.names=F, col.names=F)
        

        ##SPATIAL DISTRIBUTION OF HITS - output
        write.table("<BR><BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\"></TH>", 
        file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        
        for (j in res[[3]][1]:res[[3]][2]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
            file="stats.html", append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file="stats.html", append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (k in res[[2]][1]:res[[2]][2]){
            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Exp. ", k, "</TH>", sep=""), 
            file="stats.html", append=T, quote=F, row.names=F, col.names=F)

            for (m in res[[3]][1]:res[[3]][2]){
            
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", res[[1]], 
                "_Exp_", k, "_Plate_", m, "_.html", "\"><IMG SRC=\"", res[[1]], "_Exp_", k, 
                "_Plate_", m, "_4html.png\"/></A></CENTER><BR></TD>", sep=""), 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)
            }
        }
        write.table("</TR></CENTER></TABLE>", file="stats.html", append=T, quote=F, 
        row.names=F, col.names=F)
    }
    
    
    blob<-vennDiagfunc(hitVectors, hitVectorsTestNames, hitScoringVec1, 
    listOfArgs4stat)
    
    typeOfScoringVec<-blob[[1]]
    typeOfComparison<-blob[[2]]
    
    write.table("<br><br><br><CENTER><H2>Venn diagrams:</H2></CENTER><br>", 
    file="stats.html", append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><CENTER><TABLE border=0><TR>", file="stats.html", append=T, 
    quote=F, row.names=F, col.names=F)
    
    for (element in c("l", "g", "two.sided")){
        if (length(which(typeOfScoringVec == element))>1){
        
            if (element == "l"){
                write.table("<TH BGCOLOR=\"#e0e0ff\">Downregulated genes</TH>", 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)
            }
            if (element == "g"){
                write.table("<TH BGCOLOR=\"#e0e0ff\">Upregulated genes</TH>", 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)
            }
            if (element == "two.sided"){
                write.table("<TH BGCOLOR=\"#e0e0ff\">Up- and downregulated genes</TH>", 
                file="stats.html", append=T, quote=F, row.names=F, col.names=F)
            }
        }
    }
    write.table("</TR><TR>", file="stats.html", append=T, quote=F, row.names=F, 
    col.names=F)
    
    
    ##make the venn diagram plots:
    for (element in c("l", "g", "two.sided")){
        if (length(which(typeOfScoringVec == element))>1){

            write.table("<TD BGCOLOR=\"#f0f0ff\"><CENTER><TABLE border=0><TR>", 
            file="stats.html", append=T, quote=F, row.names=F, col.names=F)

            for (compType in c("P", "Z", "PZ")){
                if (length(which(typeOfComparison == compType 
                & typeOfScoringVec == element))>1){

                    arg4venn2<-as.list(hitVectorsTestNames[which(typeOfComparison == compType 
                    & typeOfScoringVec == element)])

                    chosenScoringVecIndexes<-which((typeOfScoringVec == element) 
                    & (typeOfComparison == compType))
                    
                    arg4venn1<-listChooseIndex(hitVectors, chosenScoringVecIndexes)

                    if (length(arg4venn1) == 2 | length(arg4venn1) == 3){
                        vennPlotName<-vennDiag(header, arg4venn1, arg4venn2, 
                        paste("Venn Diagram (", element, ") (", compType, ")", sep=""), 0)
                        
                        write.table(paste("<TH><A HREF=\"", vennPlotName, ".pdf\"><IMG SRC=\"", 
                        vennPlotName, ".png\"/></A></TH>", sep=""), file="stats.html", append=T, 
                        quote=F, row.names=F, col.names=F)
                    }
                }
            }
            write.table("</TD></TR></CENTER></TABLE>", file="stats.html", append=T, 
            quote=F, row.names=F, col.names=F)
        }
    }
    write.table("</TR></CENTER></TABLE>", file="stats.html", append=T, 
    quote=F, row.names=F, col.names=F)
        
    ##LINK TO FULL DATASET:
    saveDataset(header, dataset, paste(firstLineBis, "_Fulldataset.txt", sep=""))
    
    write.table(paste("<br><br><br><CENTER><a href=\"", firstLineBis, 
    "_Fulldataset.txt\">Full scored dataset</a></CENTER><br>", sep=""), 
    file="stats.html", append=T, quote=F, row.names=F, col.names=F)
    
    ##GSEA analysis
    
    if(length(flag4Gsea)==1){
        if(flag4Gsea==1){
            print(paste("Performing GSEA analysis...", date(), sep=": "))
            
            GOlist<-gseaAnalysisPt1(hitVectors[[1]], whichOnto)
    
            for (g in 1:length(hitVectors)){

                if (sum(hitVectors[[g]])>2 & sum(hitVectors[[g]])<length(hitVectors[[g]])){
                    geneTable<-gseaAnalysisPt2(hitVectors[[g]], GOlist, whichOnto)
    
                    write.table(geneTable, file=paste(firstLineBis, "_GSEAanalysis_", 
                    hitVectorsTestNames[[g]], " (", typeOfScoringVec[[g]], ").txt", sep=""))

                    gseaTab<-read.table(paste(firstLineBis, "_GSEAanalysis_", 
                    hitVectorsTestNames[[g]], " (", typeOfScoringVec[[g]], ").txt", sep=""))
            
                    gseaTableFunc(gseaTab, paste(firstLineBis, "_GSEAanalysis_", 
                    hitVectorsTestNames[[g]], " (", typeOfScoringVec[[g]], ").txt", sep=""), 
                    paste("Results of GSEA analysis for ", hitVectorsTestNames[[g]], " (", 
                    typeOfScoringVec[[g]], ")", sep=""), paste("subGSEA", g, ".html", sep=""))
    
                    write.table(paste("<br><CENTER><a href=\"subGSEA", g, 
                    ".html\">Results of GSEA analysis for ", hitVectorsTestNames[[g]], " (", 
                    typeOfScoringVec[[g]], ")</a></CENTER>", sep=""), file="stats.html", 
                    append=T, quote=F, row.names=F, col.names=F)
                }
            }
        }
    }else{
        print(paste("Performing GSEA analysis...", date(), sep=": "))
        GOlist<-gseaAnalysisPt1(hitVectors[[1]], whichOnto)
        
        for(i in 1:length(listOfStatTests)){
            for (j in 1:length(hitScoringVec1)){
                if(hitScoringVec1[j]>0){
                
                
                    indexNum<-(3-length(which(hitScoringVec1==0)))
                    
                    g<-(j+(i-1)*indexNum)
                                        
                    if (flag4Gsea[g]==1 & sum(hitVectors[[g]])>2 & sum(hitVectors[[g]])<length(hitVectors[[g]])){

                        geneTable<-gseaAnalysisPt2(hitVectors[[g]], GOlist, whichOnto)
    
                        write.table(geneTable, file=paste(firstLineBis, "_GSEAanalysis_", 
                        hitVectorsTestNames[[g]], " (", typeOfScoringVec[[g]], ").txt", sep=""))

                        gseaTab<-read.table(paste(firstLineBis, "_GSEAanalysis_", 
                        hitVectorsTestNames[[g]], " (", typeOfScoringVec[[g]], ").txt", sep=""))
            
                        gseaTableFunc(gseaTab, paste(firstLineBis, "_GSEAanalysis_", 
                        hitVectorsTestNames[[g]], " (", typeOfScoringVec[[g]], ").txt", sep=""), 
                        paste("Results of GSEA analysis for ", hitVectorsTestNames[[g]], " (", 
                        typeOfScoringVec[[g]], ")", sep=""), paste("subGSEA", g, ".html", sep=""))
    
                        write.table(paste("<br><CENTER><a href=\"subGSEA", g, 
                        ".html\">Results of GSEA analysis for ", hitVectorsTestNames[[g]], " (", 
                        typeOfScoringVec[[g]], ")</a></CENTER>", sep=""), file="stats.html", 
                        append=T, quote=F, row.names=F, col.names=F)
                    }

                }
            }
        }
    }
    
    write.table("</BODY>", file="stats.html", append=T, quote=F, row.names=F, 
    col.names=F)
    
    s1<-"Analysis finished. Quality analysis of raw data can be found in index.html"
    s2<-", quality analysis of normalized data can be found in indexnorm.html. "
    s3<-"Statistical analysis can be found in stats.html. "
    s4<-date()
    print(paste(s1, s2, s3, s4, sep=""))

}



gseaTableFunc<-function(gseaTab, gseaTabfile, gseaTitle, nameOfHtmlPage){

    write.table(paste("<BR><BR><CENTER><H2>", gseaTitle, "</H2><a href=\"", 
    gseaTabfile, "\">Textfile</a></CENTER><br><br><CENTER><TABLE border=0><TR>", 
    sep=""), file=nameOfHtmlPage, quote=F, row.names=F, col.names=F)

    for (i in 1:length(colnames(gseaTab))){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">", colnames(gseaTab)[i], "</TH>", 
        sep=""), file=nameOfHtmlPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=nameOfHtmlPage, append=T, quote=F, row.names=F, 
    col.names=F)

    for (j in 1:nrow(gseaTab)){
        write.table("<TR>", file=nameOfHtmlPage, append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (k in 1:ncol(gseaTab)){
            write.table(paste("<TD BGCOLOR=\"#e0e0f0\">", gseaTab[j, k], "</TD>", sep=""), 
            file=nameOfHtmlPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=nameOfHtmlPage, append=T, quote=F, row.names=F, 
        col.names=F)
    }
    write.table("</TABLE></CENTER>", file=nameOfHtmlPage, append=T, quote=F, 
    row.names=F, col.names=F)
}



function1<-function(header, dataset, vecOfChannels, firstLineBis, secondLineBis, indexOutput, subPage, flag){

    #########
    ####1####
    #########
    
    write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
    quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }    
    write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, 
    col.names=F)
    
    write.table("<TR>", file=subPage, append=T, quote=F, row.names=F, 
    col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        if (flag == 0){
            namePlot1<-makeBoxplotPerScreen(header, dataset, vecOfChannels[i], 
            paste("Data (", i, ") per Exp.", sep=""), 0)
        }else{
            namePlot1<-makeBoxplotPerScreen(header, dataset, vecOfChannels[i], 
            paste("norm. Data (", flag, ") (", i, ") per Exp.", sep=""), 0)
        }
        
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", 
        namePlot1, ".pdf\"><IMG SRC=\"", namePlot1, ".png\"/></A></CENTER><BR></TD>", 
        sep=""), file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, 
    col.names=F)


    for (i in 1:length(vecOfChannels)){
    
        if (flag == 0){
            stuff<-makeBoxplotPerPlate(header, dataset, vecOfChannels[i], 
            paste("Data (", i, ") per Plate", sep=""), 2, 0)
        }else{
            stuff<-makeBoxplotPerPlate(header, dataset, vecOfChannels[i], 
            paste("norm. Data (", flag, ") (", i, ") per Plate", sep=""), 2, 0)
        }
        
        s1<-"<BR><BR><CENTER><TABLE border=0><TR><TD BGCOLOR=\"#e0e0f0\"><H3 align="
        s2<-"left>Channel"
        s3<-"</H3></TD>"
        write.table(paste(s1, s2, i, s3, sep=""), file=subPage, append=T, 
        quote=F, row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
        
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", 
            stuff[[1]], "(Exp. ", j, ").pdf\"><IMG SRC=\"", stuff[[1]], "(Exp. ", j, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, 
            quote=F, row.names=F, col.names=F)
        }
    }

    
    for (i in 1:length(vecOfChannels)){
    
        if (flag == 0){
            stuff<-makeBoxplot4PlateType(header, dataset, vecOfChannels[i], 
            paste("Data (", i, ") per Exp.", sep=""), 0)
        }else{
            stuff<-makeBoxplot4PlateType(header, dataset, vecOfChannels[i], 
            paste("norm. Data (", flag, ") (", i, ") per Exp.", sep=""), 0)
        }
        
        s4<-"<BR><BR><CENTER><TABLE border=0><TR><TD BGCOLOR=\"#e0e0f0\"><H3 align=left"
        s5<-">Channel"
        s6<-"</H3></TD>"
        write.table(paste(s4, s5, i, s6, sep=""), file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){    #for each plate type / for each plot
        
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", 
            stuff[[1]], " for Plate", j, ".pdf\"><IMG SRC=\"", stuff[[1]], 
            " for Plate", j, ".png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, 
            append=T, quote=F, row.names=F, col.names=F)
        }
    }

    write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=indexOutput, 
    append=T, quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }    
    write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)
    
    write.table("<TR>", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        if (flag == 0){
            namePlot1<-makeBoxplotPerScreen(header, dataset, vecOfChannels[i], 
            paste("Data (", i, ") per Exp.", sep=""), 0)
        }else{
            namePlot1<-makeBoxplotPerScreen(header, dataset, vecOfChannels[i], 
            paste("norm. Data (", flag, ") (", i, ") per Exp.", sep=""), 0)
        }
        
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
        "\"><IMG SRC=\"", namePlot1, ".png\"/></A></CENTER><BR></TD>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)
}


function2<-function(header, dataset, vecOfChannels, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    #########
    ####2####
    #########

    s1<-"<br><br><br><CENTER><H2>Spatial distribution of Intensities</H2></CENTER>"
    write.table(s1, file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        
    for (i in 1:length(vecOfChannels)){
        if (flag == 0){
            stuff2<-spatialDistrib(header, dataset, paste("Spatial I distribution (", i, 
            ")", sep=""), vecOfChannels[i], "Internal_GeneID", 0)
        }else{
            stuff2<-spatialDistrib(header, dataset, paste("norm. Spatial I distribution (", 
            flag, ") (", i, ")", sep=""), vecOfChannels[i], "Internal_GeneID", 0)
        }
        
        if (i>1){
            write.table("<BR>", file=indexOutput, append=T, quote=F, row.names=F, 
            col.names=F)
        }
        
        s2<-"<BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\">Channel "
        write.table(paste(s2, i, "</TH>", sep=""), file=indexOutput, append=T, 
        quote=F, row.names=F, col.names=F)
        
        for (j in stuff2[[3]][1]:stuff2[[3]][2]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (k in stuff2[[2]][1]:stuff2[[2]][2]){
            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Exp. ", k, "</TH>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)

            for (m in stuff2[[3]][1]:stuff2[[3]][2]){
            
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                "_Exp_", k, "_Plate_", m, "_.html", "\"><IMG SRC=\"", stuff2[[1]], "_Exp_", k, 
                "_Plate_", m, "_4html.png\"/></A></CENTER><BR></TD>", sep=""), file=indexOutput, 
                append=T, quote=F, row.names=F, col.names=F)
            }
        }
        write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
        row.names=F, col.names=F)
        
    }
    
    
    s4<-"<br><br><br><CENTER><H2>Standard deviation of Intensities</H2></CENTER>"
    write.table(s4, file=indexOutput, append=T, quote=F, row.names=F, col.names=F)

    write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
    
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        
        if (flag == 0){
            zeugs<-compareReplicateSD(header, dataset, paste("SD (", i, ") of replicates", 
            sep=""), vecOfChannels[i], "GeneName")
        }else{
            zeugs<-compareReplicateSD(header, dataset, paste("SD (", i, 
            ") of norm. replicates (", flag, ")", sep=""), vecOfChannels[i], "GeneName")
        }
        
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
        "\"><IMG SRC=\"", zeugs, " 4html.png\"/></A></CENTER><BR></TD>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)

    }
    write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)
    

    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)
    
    s1<-"<CENTER><H2>Standard deviation of Intensities per experiment</H2></CENTER>"
    write.table(s1, file=subPage, append=T, quote=F, row.names=F, col.names=F)

    write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
    quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        if (flag == 0){
            zeugs<-compareReplicateSD(header, dataset, paste("SD (", i, ") of replicates", 
            sep=""), vecOfChannels[i], "GeneName")
        }else{
            zeugs<-compareReplicateSD(header, dataset, paste("SD (", i, 
            ") of norm. replicates (", flag, ")", sep=""), vecOfChannels[i], "GeneName")
        }
        
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs, 
        " .html\"><IMG SRC=\"", zeugs, " 4html.png\"/></A></CENTER><BR></TD>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
    row.names=F, col.names=F)


    write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), file=subPage, 
    append=T, quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){        
    
        write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=subPage, append=T, 
        quote=F, row.names=F, col.names=F)
        
        if (flag == 0){
            stuff<-compareReplicateSDPerScreen(header, dataset, paste("SD (", i, 
            ") of replicates", sep=""), vecOfChannels[i], "GeneName", 0)
        }else{
            stuff<-compareReplicateSDPerScreen(header, dataset, paste("SD (", i, 
            ") of norm. replicates (", flag, ")", sep=""), vecOfChannels[i], "GeneName", 0)
        }
        
        for (j in stuff[[2]]:stuff[[3]]){
        
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
        
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff[[1]], 
            "_Exp_", j, "_.html\"><IMG SRC=\"", stuff[[1]], "_Exp_", j, 
            "_4html.png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, 
            quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</CENTER></TABLE>", file=subPage, append=T, quote=F, 
    row.names=F, col.names=F)
}


function2b <-function (header, dataset, vecOfChannels, 
firstLineBis, secondLineBis, indexOutput, subPage, flag){

    ##########
    ####2b####
    ##########

    subPageClip1<-strsplit(subPage, ".html")
    subPageClip<-subPageClip1[[1]][1]
    subPageClip<-paste(subPageClip, "_", sep="")
    
    s1<-"<BR><BR><BR><CENTER><H2>Distribution of data per well</H2></CENTER>"
    write.table(s1, sep="", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)

    numOfScreens<-max(dataset$ScreenNb)
    minOfScreens<-min(dataset$ScreenNb)

    write.table(paste("<BR><BR><CENTER>", sep=""), file=indexOutput, 
    append=T, quote=F, row.names=F, col.names=F)

    for (i in minOfScreens:numOfScreens){

        write.table(paste("<A HREF=\"", subPageClip, i, ".html\">Experiment ", 
        i, "</A><br>", sep=""), file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)

        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
        quote=F, row.names=F, col.names=F)

        for (o in 1:length(vecOfChannels)){
    
            if (flag == 0){
                stuff2<-ZScorePlotTwo(header, dataset, 2, 1, vecOfChannels[o], 
                "GeneName", paste("Data per well (", o, ")", sep=""), 0)
            }else{
                stuff2<-ZScorePlotTwo(header, dataset, 2, 1, vecOfChannels[o], "GeneName", 
                paste("norm. Data per well (", flag, ") (", o, ")", sep=""), 0)
            }
        
            write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
            row.names=F, col.names=F)
            
            write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=paste(subPageClip, i, 
            ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
            for (j in stuff2[[3]][1]:stuff2[[3]][2]){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)
        
            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", o, "</TH>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, row.names=F, 
            col.names=F)
        
            for (j in stuff2[[3]][1]:stuff2[[3]][2]){
        
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                "_Exp_", i, "_Plate_", j, "_.pdf\"><IMG SRC=\"", stuff2[[1]], "_Exp_", i, 
                "_Plate_", j, "_.png\"/></A></CENTER><BR></TD>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
    }
    write.table(paste("</CENTER>", sep=""), file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)
}




function2c <-function (header, dataset, vecOfChannels, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    ##########
    ####2c####
    ##########

    subPageClip1<-strsplit(subPage, ".html")
    subPageClip<-subPageClip1[[1]][1]
    subPageClip<-paste(subPageClip, "_", sep="")
    
    combinationNum<-choose(length(vecOfChannels), 2)

    write.table("<BR><BR><BR><CENTER><H2>Comparison of channels</H2></CENTER>", 
    file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)
    
    if (flag == 0){
        zeugs<-channelPlot(header, dataset, vecOfChannels, 0, "Channel comparison", 0)
    }else{
        zeugs<-channelPlot(header, dataset, vecOfChannels, 0, 
        paste("Channel comparison (", flag, ")", sep=""), 0)
    }

    for (i in 1:length(combinationNum)){
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
        "\"><IMG SRC=\"", zeugs, "(", i, ").png\"/></A></CENTER><BR></TD>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)

    
    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)
    
    write.table("<CENTER><H2>Comparison of channels - per experiment</H2></CENTER>", 
    file=subPage, append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
    quote=F, row.names=F, col.names=F)
    
    if (flag == 0){
        zeugs<-channelPlot(header, dataset, vecOfChannels, 0, "Channel comparison", 0)
    }else{
        zeugs<-channelPlot(header, dataset, vecOfChannels, 0, 
        paste("Channel comparison (", flag, ")", sep=""), 0)
    }

    for (i in 1:length(combinationNum)){
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs, "(", i, 
        ").pdf\"><IMG SRC=\"", zeugs, "(", i, ").png\"/></A></CENTER><BR></TD>", 
        sep=""), file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
    row.names=F, col.names=F)

    
    if (flag == 0){
        stuff<-channelPlot(header, dataset, vecOfChannels, 1, "Channel comparison", 0)
    }else{
        stuff<-channelPlot(header, dataset, vecOfChannels, 1, 
        paste("Channel comparison (", flag, ")", sep=""), 0)
    }

    write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
    file=subPage, append=T, quote=F, row.names=F, col.names=F)

    for (j in stuff[[2]]:stuff[[3]]){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, 
    col.names=F)
        

    for (i in 1:length(combinationNum)){
        write.table(paste("<TR>", sep=""), file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPageClip, j, 
            ".html\"><IMG SRC=\"", stuff[[1]], "(", i, ") (Exp. ", j, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, quote=F, 
            row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
    }


    for (i in stuff[[2]]:stuff[[3]]){

        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
        quote=F, row.names=F, col.names=F)

        s1<-"<CENTER><H2>Comparison of channels - per plate for experiment "
        s2<-"</H2></CENTER>"
        write.table(paste(s1, i, s2, sep=""), file=paste(subPageClip, i, ".html", 
        sep=""), append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=paste(subPageClip, i, 
        ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
        
        if (flag == 0){
            zeugs<-channelPlot(header, dataset, vecOfChannels, 1, "Channel comparison", 0)
        }else{
            zeugs<-channelPlot(header, dataset, vecOfChannels, 1, 
            paste("Channel comparison (", flag, ")", sep=""), 0)
        }
        
        for (j in 1:length(combinationNum)){
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs[[1]], "(", 
            j, ") (Exp. ", i, ").pdf\"><IMG SRC=\"", zeugs[[1]], "(", j, ") (Exp. ", i, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
        sep=""), append=T, quote=F, row.names=F, col.names=F)

    
        if (flag == 0){
            stuff2<-channelPlot(header, dataset, vecOfChannels, 2, "Channel comparison", 0)
        }else{
            stuff2<-channelPlot(header, dataset, vecOfChannels, 2, 
            paste("Channel comparison (", flag, ")", sep=""), 0)
        }
        
        write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
        file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
        row.names=F, col.names=F)

        for (j in stuff2[[3]][1]:stuff2[[3]][2]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate", j, "</TH>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
            row.names=F, col.names=F)
        }
        write.table(paste("</TR>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        
        for (o in 1:length(combinationNum)){
            write.table(paste("<TR>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
            append=T, quote=F, row.names=F, col.names=F)

            for (j in stuff2[[3]][1]:stuff2[[3]][2]){

                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                "(", o, ")_Exp_", i, "_Plate_", j, "_.pdf\"><IMG SRC=\"", stuff2[[1]], "(", o, 
                ")_Exp_", i, "_Plate_", j, "_.png\"/></A></CENTER><BR></TD>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, row.names=F, 
                col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
    }
}






function3<-function(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    #########
    ####3####
    #########
    
    subPageClip1<-strsplit(subPage, ".html")
    subPageClip<-subPageClip1[[1]][1]
    subPageClip<-paste(subPageClip, "_", sep="")
    
    write.table("<BR><BR><BR><CENTER><H2>Distribution of data</H2></CENTER>", 
    file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)

    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)

    write.table("<TR>", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)
    
    for (i in 1:length(vecOfChannels)){
    
        if (posNegFlag == 0){
            if (flag == 0){
                zeugs<-plotHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of Data (", i, ")", sep=""), 0)
            }else{
                zeugs<-plotHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of norm. Data (", flag, ")  (", i, ")", sep=""), 0)
            }
        }else{
            if (flag == 0){
                zeugs<-plotControlHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of Data and Controls (", i, ")", sep=""), 0)
            }else{
                zeugs<-plotControlHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of norm. Data (", flag, ") and Controls (", i, ")", 
                sep=""), 0)
            }
        }
        
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
        "\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)

    
    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)
    
    write.table("<CENTER><H2>Distribution of data per experiment</H2></CENTER>", 
    file=subPage, append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
    quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)

    write.table("<TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){    
        if (posNegFlag == 0){
            if (flag == 0){
                zeugs<-plotHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of Data (", i, ")", sep=""), 0)
            }else{
                zeugs<-plotHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of norm. Data (", flag, ") (", i, ")", sep=""), 0)
            }
        }else{
            if (flag == 0){
                zeugs<-plotControlHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of Data and Controls (", i, ")", sep=""), 0)
            }else{
                zeugs<-plotControlHisto(header, dataset, vecOfChannels[i], 
                paste("Distribution of norm. Data (", flag, ") and Controls (", i, 
                ")", sep=""), 0)
            }
        }
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs, 
        ".pdf\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
    row.names=F, col.names=F)

    
    for (i in 1:length(vecOfChannels)){
    
        if (posNegFlag == 0){
            if (flag == 0){
                stuff<-plotHistoPerscreen(header, dataset, vecOfChannels[i], 
                paste("Distribution of Data (", i, ")", sep=""), 2, 0)
            }else{
                stuff<-plotHistoPerscreen(header, dataset, vecOfChannels[i], 
                paste("Distribution of norm. Data (", flag, ") (", i, ")", sep=""), 2, 0)
            }
        }else{
            if (flag == 0){
                stuff<-plotControlHistoPerscreen(header, dataset, vecOfChannels[i], 
                paste("Distribution of Data and Controls (", i, ")", sep=""), 2, 0)
            }else{
                stuff<-plotControlHistoPerscreen(header, dataset, vecOfChannels[i], 
                paste("Distribution of norm. Data (", flag, ") and Controls (", i, ")", 
                sep=""), 2, 0)
            }
        }
        
        write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), file=subPage, 
        append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
        
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPageClip, j, 
            ".html\"><IMG SRC=\"", stuff[[1]], "(Exp. ", j, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, 
            quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
    }


    for (i in stuff[[2]]:stuff[[3]]){

        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
        quote=F, row.names=F, col.names=F)

        write.table(paste("<CENTER><H2>Distribution of data per plate for experiment ",
        i, "</H2></CENTER>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
        append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=paste(subPageClip, i, 
        ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
        
        for (j in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", j, "</TH>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
            row.names=F, col.names=F)
        }
        write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
        quote=F, row.names=F, col.names=F)

        write.table("<TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
        quote=F, row.names=F, col.names=F)
        
        for (j in 1:length(vecOfChannels)){
        
            if (posNegFlag == 0){
                if (flag == 0){
                    zeugs<-plotHistoPerscreen(header, dataset, vecOfChannels[j], 
                    paste("Distribution of Data (", j, ")", sep=""), 2, 0)
                }else{
                    zeugs<-plotHistoPerscreen(header, dataset, vecOfChannels[j], 
                    paste("Distribution of norm. Data (", flag, ") (", j, ")", sep=""), 2, 0)
                }
            }else{
                if (flag == 0){
                    zeugs<-plotControlHistoPerscreen(header, dataset, vecOfChannels[j], 
                    paste("Distribution of Data and Controls (", j, ")", sep=""), 2, 0)
                }else{
                    zeugs<-plotControlHistoPerscreen(header, dataset, vecOfChannels[j], 
                    paste("Distribution of norm. Data (", flag, ") and Controls (", j, ")", 
                    sep=""), 2, 0)
                }
            }
        
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs[[1]], 
            "(Exp. ", i, ").pdf\"><IMG SRC=\"", zeugs[[1]], "(Exp. ", i, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
        sep=""), append=T, quote=F, row.names=F, col.names=F)


        for (o in 1:length(vecOfChannels)){
    
            if (posNegFlag == 0){
                if (flag == 0){
                    stuff2<-plotHistoPerplate(header, dataset, vecOfChannels[o], 
                    paste("Distribution of Data (", o, ")", sep=""), 2, 0)
                }else{
                    stuff2<-plotHistoPerplate(header, dataset, vecOfChannels[o], 
                    paste("Distribution of norm. Data (", flag, ") (", o, ")", sep=""), 2, 0)
                }
            }else{
                if (flag == 0){
                    stuff2<-plotControlHistoPerplate(header, dataset, vecOfChannels[o], 
                    paste("Distribution of Data and Controls (", o, ")", sep=""), 2, 0)
                }else{
                    stuff2<-plotControlHistoPerplate(header, dataset, vecOfChannels[o], 
                    paste("Distribution of norm. Data (", flag, ") and Controls (", o, ")", 
                    sep=""), 2, 0)
                }
            }
        
            write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
            row.names=F, col.names=F)
            
            write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=paste(subPageClip, i, 
            ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
            
            for (j in stuff2[[3]][1]:stuff2[[3]][2]){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, row.names=F, 
                col.names=F)
            }
            write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), 
            append=T, quote=F, row.names=F, col.names=F)
        
            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", o, "</TH>", 
            sep=""), file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)
        
            for (j in stuff2[[3]][1]:stuff2[[3]][2]){
        
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                "_Exp_", i, "_PerPlate_", j, "_.pdf\"><IMG SRC=\"", stuff2[[1]], "_Exp_", i, 
                "_PerPlate_", j, "_.png\"/></A></CENTER><BR></TD>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
    }
}



function4a<-function(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    ##########
    ####4a####
    ##########
    
    subPageClip1<-strsplit(subPage, ".html")
    subPageClip<-subPageClip1[[1]][1]
    subPageClip<-paste(subPageClip, "_", sep="")

    if (posNegFlag == 1){

        write.table("<BR><BR><BR><CENTER><H2>Distribution of controls</H2></CENTER>", 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
        quote=F, row.names=F, col.names=F)

        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)
        write.table("<TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)

        for (i in 1:length(vecOfChannels)){

            if (flag == 0){
                zeugs<-controlDensity(header, dataset, vecOfChannels[i], 
                paste("Control density (", i, ")", sep=""), 0, 1)
            }else{
                zeugs<-controlDensity(header, dataset, vecOfChannels[i], 
                paste("norm. Control density (", flag, ")  (", i, ")", sep=""), 0, 1)
            }

            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
            "\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
        row.names=F, col.names=F)

        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)

        write.table("<CENTER><H2>Control density per experiment</H2></CENTER>", 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, 
        append=T, quote=F, row.names=F, col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, 
        col.names=F)

        write.table("<TR>", file=subPage, append=T, quote=F, row.names=F, 
        col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            if (flag == 0){
                zeugs<-controlDensity(header, dataset, vecOfChannels[i], 
                paste("Control density (", i, ")", sep=""), 0, 1)
            }else{
                zeugs<-controlDensity(header, dataset, vecOfChannels[i], 
                paste("norm. Control density (", flag, ")  (", i, ")", sep=""), 0, 1)
            }
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs, 
            ".pdf\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)


        for (i in 1:length(vecOfChannels)){

            if (flag == 0){
                stuff<-controlDensityPerScreen(header, dataset, vecOfChannels[i], 
                paste("Control density (", i, ")", sep=""), 0, 1)
            }else{
                stuff<-controlDensityPerScreen(header, dataset, vecOfChannels[i], 
                paste("norm. Control density (", flag, ") (", i, ")", sep=""), 0, 1)
            }


            write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)

            write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=subPage, append=T, 
            quote=F, row.names=F, col.names=F)
            
            for (j in stuff[[2]]:stuff[[3]]){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), 
                file=subPage, append=T, quote=F, row.names=F, col.names=F)
            }
            write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)

            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)

            for (j in stuff[[2]]:stuff[[3]]){

                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPageClip, j, 
                ".html\"><IMG SRC=\"", stuff[[1]], "(Exp. ", j, 
                ").png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, 
                quote=F, row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
            row.names=F, col.names=F)
        }


        for (i in stuff[[2]]:stuff[[3]]){

            write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
            "</TITLE></HEAD>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
            quote=F, row.names=F, col.names=F)

            write.table(paste("<CENTER><H2>Control density per plate for experiment ", i, 
            "</H2></CENTER>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
            append=T, quote=F, row.names=F, col.names=F)
            
            write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=paste(subPageClip, i, 
            ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
            
            for (j in 1:length(vecOfChannels)){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", j, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)

            write.table("<TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)
            
            for (j in 1:length(vecOfChannels)){
                if (flag == 0){
                    zeugs<-controlDensityPerScreen(header, dataset, vecOfChannels[j], 
                    paste("Control density (", j, ")", sep=""), 0, 1)
                }else{
                    zeugs<-controlDensityPerScreen(header, dataset, vecOfChannels[j], 
                    paste("norm. Control density (", flag, ") (", j, ")", sep=""), 0, 1)
                }

                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs[[1]], 
                "(Exp. ", i, ").pdf\"><IMG SRC=\"", zeugs[[1]], "(Exp. ", i, 
                ").png\"/></A></CENTER><BR></TD>", sep=""), file=paste(subPageClip, i, 
                ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)

            for (o in 1:length(vecOfChannels)){

                if (flag == 0){
                    stuff2<-controlDensityPerPlate(header, dataset, vecOfChannels[o], 
                    paste("Control density (", o, ")", sep=""), 2, 0, 1)
                }else{
                    stuff2<-controlDensityPerPlate(header, dataset, vecOfChannels[o], 
                    paste("norm. Control density (, ", flag, ") (", o, ")", sep=""), 2, 0, 1)
                }

                write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)

                write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
                
                for (j in stuff2[[3]][1]:stuff2[[3]][2]){
                    write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                    file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                    row.names=F, col.names=F)
                }
                write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
                quote=F, row.names=F, col.names=F)

                write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", o, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
                

                for (j in stuff2[[3]][1]:stuff2[[3]][2]){

                    write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                    "_Exp_", i, "_PerPlate_", j, "_.pdf\"><IMG SRC=\"", stuff2[[1]], "_Exp_", i, 
                    "_PerPlate_", j, "_.png\"/></A></CENTER><BR></TD>", sep=""), 
                    file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                    row.names=F, col.names=F)
                }
                write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
                sep=""), append=T, quote=F, row.names=F, col.names=F)
            }
        }    
    }
}





function4b<-function(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    ##########
    ####4b####
    ##########

    if (posNegFlag == 1){

        write.table("<BR><BR><BR><CENTER><H2>Z' Factors</H2></CENTER>", 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, 
        append=T, quote=F, row.names=F, col.names=F)

        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)
        write.table("<TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)

        for (i in 1:length(vecOfChannels)){

            if (flag == 0){
                zeugs<-ZPRIMEQualControl(header, dataset, vecOfChannels[i], 
                paste("Z' Factors (", i, ")", sep=""), 0)
            }else{
                zeugs<-ZPRIMEQualControl(header, dataset, vecOfChannels[i], 
                paste("Z' Factors (", flag, ") (", i, ") (norm.)", sep=""), 0)
            }

            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
            "\"><IMG SRC=\"", zeugs[[1]], ".png\"/></A></CENTER><BR></TD>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
        row.names=F, col.names=F)
        
        
        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)
        
        write.table("<CENTER><H2>Z' Factors</H2></CENTER>", file=subPage, append=T, 
        quote=F, row.names=F, col.names=F)
        write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
        quote=F, row.names=F, col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), file=subPage, 
            append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)

        write.table("<TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            if (flag == 0){
                zeugs<-ZPRIMEQualControl(header, dataset, vecOfChannels[i], 
                paste("Z' Factors (", i, ")", sep=""), 0)
            }else{
                zeugs<-ZPRIMEQualControl(header, dataset, vecOfChannels[i], 
                paste("Z' Factors (", flag, ") (", i, ") (norm.)", sep=""), 0)
            }
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs[[1]], 
            ".pdf\"><IMG SRC=\"", zeugs[[1]], ".png\"/></A></CENTER><BR></TD>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)

        s1<-"<BR><BR><CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">"
        s2<-"Exp. #</TH><TH BGCOLOR=\"#e0e0ff\">Plate #</TH>"
        s3<-paste(s1, s2, sep="")
        write.table(s3, file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#d0d0f0\">Z' Score Ch. ", i, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        
        for (i in 2:nrow(zeugs[[2]])){
            
            write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", zeugs[[2]][i, 1], 
            "</TD><TD BGCOLOR=\"#f0f0ff\" align=center>", zeugs[[2]][i, 2], "</TD>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
            
            for (j in 1:length(vecOfChannels)){
                if (flag == 0){
                    zeugs<-ZPRIMEQualControl(header, dataset, vecOfChannels[j], 
                    paste("Z' Factors (", j, ")", sep=""), 0)
                }else{
                    zeugs<-ZPRIMEQualControl(header, dataset, vecOfChannels[j], 
                    paste("Z' Factors (", flag, ") (", j, ") (norm.)", sep=""), 0)
                }
                write.table(paste("<TD BGCOLOR=\"#e0e0f0\" align=center>", zeugs[[2]][i, 3], 
                "</TD>", sep=""), file=subPage, append=T, quote=F, row.names=F, col.names=F)
                
            }
            write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, 
            col.names=F)
        }
        
        write.table("</TABLE></CENTER>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
    }
}




function5<-function(header, dataset, vecOfChannels, posNegFlag, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    #########
    ####5####
    #########
    
    subPageClip1<-strsplit(subPage, ".html")
    subPageClip<-subPageClip1[[1]][1]
    subPageClip<-paste(subPageClip, "_", sep="")
    
    if (posNegFlag == 1){
    
        write.table("<BR><BR><BR><CENTER><H2>Comparison of controls</H2></CENTER>", 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
        quote=F, row.names=F, col.names=F)

        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)

        write.table("<TR>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)
    
        for (i in 1:length(vecOfChannels)){
    
            if (flag == 0){
                zeugs<-makeBoxplotControls(header, dataset, vecOfChannels[i], 
                paste("Comparison of controls (", i, ")", sep=""), 0)
            }else{
                zeugs<-makeBoxplotControls(header, dataset, vecOfChannels[i], 
                paste("Comparison of controls (", flag, ") (", i, ") (norm.)", sep=""), 0)
            }
        
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
            "\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
            file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
        row.names=F, col.names=F)
    
        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)
    

        write.table("<CENTER><H2>Comparison of controls per experiment<H2></CENTER>", 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
        quote=F, row.names=F, col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)

        write.table("<TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
        for (i in 1:length(vecOfChannels)){
            if (flag == 0){
                zeugs<-makeBoxplotControls(header, dataset, vecOfChannels[i], 
                paste("Comparison of controls (", i, ")", sep=""), 0)
            }else{
                zeugs<-makeBoxplotControls(header, dataset, vecOfChannels[i], 
                paste("Comparison of controls (", flag, ") (", i, ") (norm.)", sep=""), 0)
            }
            
            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs, 
            ".pdf\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)


        for (i in 1:length(vecOfChannels)){
    
            if (flag == 0){
                stuff<-makeBoxplotControlsPerScreen(header, dataset, vecOfChannels[i], 
                paste("Comparison of controls (", i, ")", sep=""), 2, 0)
            }else{
                stuff<-makeBoxplotControlsPerScreen(header, dataset, vecOfChannels[i], 
                paste("Comparison of controls (", flag, ") (", i, ") (norm.)", sep=""), 2, 0)
            }
        
            write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), file=subPage, 
            append=T, quote=F, row.names=F, col.names=F)
        
            write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=subPage, append=T, quote=F, 
            row.names=F, col.names=F)
            
            for (j in stuff[[2]]:stuff[[3]]){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), file=subPage, 
                append=T, quote=F, row.names=F, col.names=F)
            }
            write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)
        
            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)

            for (j in stuff[[2]]:stuff[[3]]){
        
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPageClip, j, 
                ".html\"><IMG SRC=\"", stuff[[1]], " (Exp.", j, 
                ").png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
            row.names=F, col.names=F)
        }



        for (i in stuff[[2]]:stuff[[3]]){

            write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
            "</TITLE></HEAD>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
            quote=F, row.names=F, col.names=F)

            s1<-"<CENTER><H2>Comparison of controls per plate for experiment "
            s2<-"</H2></CENTER>"
            write.table(paste(s1, i, s2, sep=""), file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
            
            write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=paste(subPageClip, i, 
            ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
            
            for (j in 1:length(vecOfChannels)){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", j, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)

            write.table("<TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)
            
            for (j in 1:length(vecOfChannels)){
                if (flag == 0){
                    zeugs<-makeBoxplotControlsPerScreen(header, dataset, vecOfChannels[j], 
                    paste("Comparison of controls (", j, ")", sep=""), 2, 0)
                }else{
                    zeugs<-makeBoxplotControlsPerScreen(header, dataset, vecOfChannels[j], 
                    paste("Comparison of controls (", flag, ") (", j, ") (norm.)", sep=""), 2, 0)
                }
        
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs[[1]], 
                " (Exp.", i, ").pdf\"><IMG SRC=\"", zeugs[[1]], " (Exp.", i, 
                ").png\"/></A></CENTER><BR></TD>", sep=""), file=paste(subPageClip, i, ".html", 
                sep=""), append=T, quote=F, row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)


            for (o in 1:length(vecOfChannels)){

                if (flag == 0){
                    stuff2<-makeBoxplotControlsPerPlate(header, dataset, vecOfChannels[o], 
                    paste("Comparison of controls (", o, ")", sep=""), 2, 0)
                }else{
                    stuff2<-makeBoxplotControlsPerPlate(header, dataset, vecOfChannels[o], 
                    paste("Comparison of controls (", flag, ") (", o, ") (norm.)", sep=""), 2, 0)
                }

        
                write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            

                write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=paste(subPageClip, i, 
                ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
                
                for (j in stuff2[[3]][1]:stuff2[[3]][2]){
                    write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                    file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                    row.names=F, col.names=F)
                }
                write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
                quote=F, row.names=F, col.names=F)
        
                write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", o, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
        
                for (j in stuff2[[3]][1]:stuff2[[3]][2]){
        
                    write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                    "_Exp_", i, "_PerPlate", j, ".pdf\"><IMG SRC=\"", stuff2[[1]], "_Exp_", i, 
                    "_PerPlate", j, ".png\"/></A></CENTER><BR></TD>", sep=""), 
                    file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                    row.names=F, col.names=F)
                }
                write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
                sep=""), append=T, quote=F, row.names=F, col.names=F)
            }
        }
    }
}




function6<-function(header, dataset, vecOfChannels, firstLineBis, 
secondLineBis, indexOutput, subPage, flag){

    #########
    ####6####
    #########

    subPageClip1<-strsplit(subPage, ".html")
    subPageClip<-subPageClip1[[1]][1]
    subPageClip<-paste(subPageClip, "_", sep="")    
    

    write.table("<BR><BR><BR><CENTER><H2>QQ plots</H2></CENTER>", file=indexOutput, 
    append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)

    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)

    write.table("<TR>", file=indexOutput, append=T, quote=F, row.names=F, 
    col.names=F)

    for (i in 1:length(vecOfChannels)){
        if (flag == 0){
            zeugs<-plotQQ(header, dataset, vecOfChannels[i], paste("QQ plot (", i, ")", 
            sep=""), 0)
        }else{
            zeugs<-plotQQ(header, dataset, vecOfChannels[i], paste("QQ plot (", flag, 
            ") (", i, ") (norm.)", sep=""), 0)
        }

        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPage, 
        "\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
        file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)


    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file=subPage, quote=F, row.names=F, col.names=F)

    write.table("<CENTER><H2>QQ plots per experiment</H2></CENTER>", file=subPage, 
    append=T, quote=F, row.names=F, col.names=F)
    
    write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=subPage, append=T, 
    quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, 
    col.names=F)

    write.table("<TR>", file=subPage, append=T, quote=F, row.names=F, 
    col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        if (flag == 0){
            zeugs<-plotQQ(header, dataset, vecOfChannels[i], paste("QQ plot (", 
            i, ")", sep=""), 0)
        }else{
            zeugs<-plotQQ(header, dataset, vecOfChannels[i], paste("QQ plot (", 
            flag, ") (", i, ") (norm.)", sep=""), 0)
        }
        
        write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs, 
        ".pdf\"><IMG SRC=\"", zeugs, ".png\"/></A></CENTER><BR></TD>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)
    }
    write.table("</TR></CENTER></TABLE>", file=subPage, append=T, quote=F, 
    row.names=F, col.names=F)


    for (i in 1:length(vecOfChannels)){

        if (flag == 0){
            stuff<-plotQQperscreen(header, dataset, vecOfChannels[i], paste("QQ plot (", 
            i, ")", sep=""), 2, 0)
        }else{
            stuff<-plotQQperscreen(header, dataset, vecOfChannels[i], paste("QQ plot (", 
            flag, ") (", i, ") (norm.)", sep=""), 2, 0)
        }

        write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), file=subPage, 
        append=T, quote=F, row.names=F, col.names=F)

        write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=subPage, append=T, quote=F, 
        row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), 
            file=subPage, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subPage, append=T, quote=F, row.names=F, col.names=F)

        write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subPage, append=T, quote=F, row.names=F, col.names=F)

        for (j in stuff[[2]]:stuff[[3]]){

            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", subPageClip, j, 
            ".html\"><IMG SRC=\"", stuff[[1]], " (Exp.", j, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=subPage, append=T, quote=F, 
            row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subPage, append=T, 
        quote=F, row.names=F, col.names=F)
    }


    for (i in stuff[[2]]:stuff[[3]]){

        write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
        "</TITLE></HEAD>", sep=""), file=paste(subPageClip, i, ".html", sep=""), 
        quote=F, row.names=F, col.names=F)

        s1<-"<CENTER><H2>Comparison of controls per plate for experiment "
        s2<-"</H2></CENTER>"
        write.table(paste(s1, i, s2, sep=""), file=paste(subPageClip, i, ".html", 
        sep=""), append=T, quote=F, row.names=F, col.names=F)
        
        write.table("<BR><BR><CENTER><TABLE border=0><TR>", file=paste(subPageClip, i, 
        ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
        
        for (j in 1:length(vecOfChannels)){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Channel ", j, "</TH>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, row.names=F, 
            col.names=F)
        }
        write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
        quote=F, row.names=F, col.names=F)

        write.table("<TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
        quote=F, row.names=F, col.names=F)
        
        for (j in 1:length(vecOfChannels)){

            if (flag == 0){
                zeugs<-plotQQperscreen(header, dataset, vecOfChannels[j], paste("QQ plot (", j, 
                ")", sep=""), 2, 0)
            }else{
                zeugs<-plotQQperscreen(header, dataset, vecOfChannels[j], paste("QQ plot (", 
                flag, ") (", j, ") (norm.)", sep=""), 2, 0)
            }

            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", zeugs[[1]], 
            " (Exp.", i, ").pdf\"><IMG SRC=\"", zeugs[[1]], " (Exp.", i, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
        sep=""), append=T, quote=F, row.names=F, col.names=F)


        for (o in 1:length(vecOfChannels)){

            if (flag == 0){
                stuff2<-plotQQperplate(header, dataset, vecOfChannels[o], paste("QQ plot (", o, 
                ")", sep=""), 2, 0)
            }else{
                stuff2<-plotQQperplate(header, dataset, vecOfChannels[o], paste("QQ plot (", 
                flag, ") (", o, ") (norm.)", sep=""), 2, 0)
            }

            write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
            row.names=F, col.names=F)

            write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=paste(subPageClip, i, 
            ".html", sep=""), append=T, quote=F, row.names=F, col.names=F)
            
            for (j in stuff2[[3]][1]:stuff2[[3]][2]){
                write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Plate ", j, "</TH>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR>", file=paste(subPageClip, i, ".html", sep=""), append=T, 
            quote=F, row.names=F, col.names=F)

            write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", o, "</TH>", sep=""), 
            file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
            row.names=F, col.names=F)

            for (j in stuff2[[3]][1]:stuff2[[3]][2]){

                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff2[[1]], 
                "Exp", i, "_PerPlate", j, ".pdf\"><IMG SRC=\"", stuff2[[1]], "Exp", i, 
                "_PerPlate", j, ".png\"/></A></CENTER><BR></TD>", sep=""), 
                file=paste(subPageClip, i, ".html", sep=""), append=T, quote=F, 
                row.names=F, col.names=F)
            }
            write.table("</TR></CENTER></TABLE>", file=paste(subPageClip, i, ".html", 
            sep=""), append=T, quote=F, row.names=F, col.names=F)
        }
    }
}



function7<-function(header, dataset, vecOfChannels, firstLineBis, 
secondLineBis, indexOutput, subpageA, subpageB, flag, flagForNorm){

    #########
    ####7####
    #########
    
    s1<-"<BR><BR><BR><CENTER><H2>Spearmans correlation coefficients</H2></CENTER>"
    write.table(s1, file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    
    
    write.table("<BR><CENTER><TABLE border=0><TR>", file=indexOutput, append=T, 
    quote=F, row.names=F, col.names=F)
    
    if (flag == 1){
        for (j in 1:2){
            for (i in 1:length(vecOfChannels)){
                write.table(paste("<TH>Channel", i, "</TH>", sep=""), file=indexOutput, 
                append=T, quote=F, row.names=F, col.names=F)
            }
        }
    }else{
        for (i in 1:length(vecOfChannels)){
            write.table(paste("<TH>Channel", i, "</TH>", sep=""), file=indexOutput, 
            append=T, quote=F, row.names=F, col.names=F)
        }
    }
    write.table("</TR><TR>", file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)


    for (i in 1:length(vecOfChannels)){

        write.table("<TD>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)

        tableOfCoeffs<-replicatesSpearmancor(header, dataset, 1, vecOfChannels[i], 
        "GeneName", "per_experiment")
        
        s2<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\"></TH><TH BGCOLOR="
        s3<-"\"#e0e0ff\">Rep. 1 vs. 2</TH><TH BGCOLOR=\"#d0d0f0\">Rep. 2 vs. 3</TH>"
        s4<-"<TH BGCOLOR=\"#e0e0ff\">Rep. 1 vs. 3</TH></TR>"
        s5<-paste(s2, s3, s4, sep="")
        write.table(s5, file=indexOutput, append=T, quote=F, row.names=F, col.names=F)
    
        counter<-1
        for (j in 1:max(tableOfCoeffs[, 1])){

            write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">Exp. ", j, 
            "</TD><TD BGCOLOR=\"#f0f0ff\" align=center>", tableOfCoeffs[counter, 4], 
            "</TD><TD BGCOLOR=\"#e0e0f0\" align=center>", tableOfCoeffs[(counter+1), 4], 
            "</TD><TD BGCOLOR=\"#f0f0ff\" align=center>", tableOfCoeffs[(counter+2), 4], 
            "</TD></TR>", sep=""), file=indexOutput, append=T, quote=F, row.names=F, 
            col.names=F)
    
            counter<-counter+3
        }
        write.table("</CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
        row.names=F, col.names=F)

        write.table("</TD>", file=indexOutput, append=T, quote=F, row.names=F, 
        col.names=F)
    
    }
    

    if (flag == 1){
        for (i in 1:length(vecOfChannels)){

            write.table("<TD>", file=indexOutput, append=T, quote=F, row.names=F, 
            col.names=F)
    
            tableOfCoeffs<-replicatesSpearmancor(header, dataset, 2, vecOfChannels[i], 
            "GeneName", "between_experiments")
            
            s6<-"<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#d0d0f0\">Exp.</TH><TH BGCOLOR="
            s7<-"\"#e0e0ff\">Exp.</TH><TH BGCOLOR=\"#d0d0f0\">Correlation coeff.</TH></TR>"
            s8<-paste(s6, s7, sep="")
            write.table(s8, file=indexOutput, append=T, quote=F, row.names=F, col.names=F)    
    
            for (j in 1:nrow(tableOfCoeffs)){
                write.table(paste("<TR><TD BGCOLOR=\"#e0e0f0\">", tableOfCoeffs[j, 1], 
                "</TD><TD BGCOLOR=\"#f0f0ff\" align=center>", tableOfCoeffs[j, 2], 
                "</TD><TD BGCOLOR=\"#e0e0f0\" align=center>", tableOfCoeffs[j, 3], 
                "</TD></TR>", sep=""), file=indexOutput, append=T, quote=F, 
                row.names=F, col.names=F)
            }        
    
            write.table("</CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
            row.names=F, col.names=F)
            
            write.table("</TD>", file=indexOutput, append=T, quote=F, row.names=F, 
            col.names=F)
        
        }
    }


    write.table("</TR></CENTER></TABLE>", file=indexOutput, append=T, quote=F, 
    row.names=F, col.names=F)
    
    s9<-"<A HREF=\""
    s1<-"\">Correlation of variation plots</A>, <A HREF=\""
    s2<-"\">Plots for replicate comparison</A>"
    write.table(paste(s9, subpageA, s1, subpageB, s2, sep=""), file=indexOutput, 
    append=T, quote=F, row.names=F, col.names=F)

    
    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file=subpageA, quote=F, row.names=F, col.names=F)
    
    write.table("<CENTER><H2>Correlation of variation</H2></CENTER>", file=subpageA, 
    append=T, quote=F, row.names=F, col.names=F)

    for (i in 1:length(vecOfChannels)){

        stuff<-replicatesCV(header, dataset, paste("Correlation of variation (", 
        flagForNorm, ") (", i, ")", sep=""), vecOfChannels[i], "GeneName", 2, 0)

        write.table(paste("<BR><BR><CENTER><TABLE border=0><TR>", sep=""), 
        file=subpageA, append=T, quote=F, row.names=F, col.names=F)

        write.table("<TH BGCOLOR=\"#e0e0ff\"></TH>", file=subpageA, append=T, quote=F, 
        row.names=F, col.names=F)
        
        for (j in stuff[[2]]:stuff[[3]]){
            write.table(paste("<TH BGCOLOR=\"#e0e0ff\">Exp. ", j, "</TH>", sep=""), 
            file=subpageA, append=T, quote=F, row.names=F, col.names=F)
        }
        write.table("</TR>", file=subpageA, append=T, quote=F, row.names=F, col.names=F)

        write.table(paste("<TR><TH BGCOLOR=\"#e0e0ff\">Channel ", i, "</TH>", sep=""), 
        file=subpageA, append=T, quote=F, row.names=F, col.names=F)

        for (j in stuff[[2]]:stuff[[3]]){

            write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff[[1]], 
            " (Exp.", j, ").pdf\"><IMG SRC=\"", stuff[[1]], " (Exp.", j, 
            ").png\"/></A></CENTER><BR></TD>", sep=""), file=subpageA, append=T, 
            quote=F, row.names=F, col.names=F)
        }
        write.table("</TR></CENTER></TABLE>", file=subpageA, append=T, quote=F, 
        row.names=F, col.names=F)
    }


    write.table(paste("<HTML><HEAD><TITLE>", secondLineBis, ", ", firstLineBis, 
    "</TITLE></HEAD>", sep=""), file=subpageB, quote=F, row.names=F, col.names=F)
    
    write.table("<CENTER><H2>Comparison of replicates</H2></CENTER>", 
    file=subpageB, append=T, quote=F, row.names=F, col.names=F)
    
    for (i in 1:length(vecOfChannels)){
        
        if (i>1){
            write.table("<BR>", file=subpageB, append=T, quote=F, row.names=F, col.names=F)
        }
    
        write.table(paste("<BR><CENTER>Channel ", i, "</CENTER>", sep=""), 
        file=subpageB, append=T, quote=F, row.names=F, col.names=F)
    
        stuff<-compareReplicates(header, dataset, paste("Comparison of replicates (", 
        flagForNorm, ") (", i, ")", sep=""), vecOfChannels[i], "GeneName", 2, 0)
        
        for (j in stuff[[2]]:stuff[[3]]){
        
            write.table(paste("<CENTER><TABLE border=0><TR><TH BGCOLOR=\"#e0e0ff\">Exp.", 
            j, "</TH>", sep=""), file=subpageB, append=T, quote=F, row.names=F, col.names=F)
            
            for (k in 1:stuff[[4]]){
                write.table(paste("<TD BGCOLOR=\"#f0f0ff\"><CENTER><A HREF=\"", stuff[[1]], 
                "_Exp_", j, "(", k, ").pdf\"><IMG SRC=\"", stuff[[1]], "_Exp_", j, "(", k, 
                ").png\"/></A></CENTER><BR></TD>", sep=""), file=subpageB, append=T, quote=F, 
                row.names=F, col.names=F)
            }
        }
        write.table("</TR></CENTER></TABLE>", file=subpageB, append=T, quote=F, 
        row.names=F, col.names=F)
    }
}



vennDiagfunc<-function(hitVectors, hitVectorsTestNames, hitScoringVec1, 
listOfArgs4stat){

    typeOfScoringVec<-rep(NA_character_, length(hitVectors))
    typeOfComparison<-rep(NA_character_, length(hitVectors))
    saveif<-0
    
    if (hitScoringVec1[1] == 1){
    
        for (i in 1:length(hitVectorsTestNames)){

            if (hitVectorsTestNames[[i]] == "t test" 
            | hitVectorsTestNames[[i]] == "Mann-Whitney test"){

                if (listOfArgs4stat[[i]][1] == "l"){
                    typeOfScoringVec[i]<-"l"
                    typeOfComparison[i]<-"P"
                }

                if (listOfArgs4stat[[i]][1] == "g"){
                    typeOfScoringVec[i]<-"g"
                    typeOfComparison[i]<-"P"
                }

                if (listOfArgs4stat[[i]][1] == "two.sided"){
                    typeOfScoringVec[i]<-"two.sided"
                    typeOfComparison[i]<-"P"
                }
            }

            if (hitVectorsTestNames[[i]] == "Rank product test"){

                if (listOfArgs4stat[[i]][2] == 1){
                    typeOfScoringVec[i]<-"l"
                    typeOfComparison[i]<-"P"
                }else{
                    typeOfScoringVec[i]<-"g"
                    typeOfComparison[i]<-"P"
                }
            }
        }
    }
    
    if (hitScoringVec1[2] == 1 | hitScoringVec1[2] == 2){
    
        saveif<-length(which(typeOfScoringVec == "l"))
        for (g in 1:saveif){
        
            if(sum(is.na(typeOfScoringVec))!=0){
            ##if there are NAs in typeOfScoringVec
        
                ind<-min(which(is.na(typeOfScoringVec)))
                typeOfScoringVec[ind]<-"l"
    
                if (hitScoringVec1[2] == 1){
                    typeOfComparison[ind]<-"Z"
                }
                if (hitScoringVec1[2] == 2){
                    typeOfComparison[ind]<-"PZ"
                }
            }
        }
    }    
    
    if (hitScoringVec1[3] == 1 | hitScoringVec1[3] == 2){
    
        saveif<-length(which(typeOfScoringVec == "g"))
        for (g in 1:saveif){
        
            if(sum(is.na(typeOfScoringVec))!=0){
            ##if there are NAs in typeOfScoringVec
        
                ind<-min(which(is.na(typeOfScoringVec)))
                typeOfScoringVec[ind]<-"g"

                if (hitScoringVec1[3] == 1){
                    typeOfComparison[ind]<-"Z"
                }
                if (hitScoringVec1[3] == 2){
                    typeOfComparison[ind]<-"PZ"
                }
            }
        }
    }
    
    invisible(list(typeOfScoringVec, typeOfComparison))
}


listChooseIndex<-function(fullList, chooseVec){

    for (g in 1:length(chooseVec)){
        if (g == 1){
            partList<-list(fullList[[chooseVec[g]]])
        }else{
            partList[[g]]<-fullList[[chooseVec[g]]]
        }
    }
    invisible(partList)    
}
