# RNAither - statistical analysis of RNAi primary screens
# This file contains extensions added to RNAither in version 2,
# most notably the "rnaither" function to carry out a full analysis,
# as well as a few additional plotting functions
# rnaither should be called with
# the data read in an RNAither data frame, and with options specifying what
# normalization to perform. Data analysis will then be carried out fully automatical.

#library("geneplotter")
#library("car")
#library("splots")

#makeallplots <<- TRUE

  setClass("rnaistat_data",
           representation(expname="character",
                          excludeCellcounts="character",
                          logtransform="logical",
                          normalization="character",
			  test="character",
                          scorethresh="numeric",
                          pvalthresh="numeric",
                          dogo="logical",
                          thedata="data.frame",
			  rawdata="data.frame",
			  normalizeddata="data.frame",
                          status="list",
                          header="character",
                          layoutnames="character"),
           prototype(expname="No Name",
                     excludeCellcounts="none",
                     logtransform=FALSE,
                     normalization=c("lowess","bscore"),
                     test="ttest",
                     scorethresh=2.0,
                     pvalthresh=0.05,
                     dogo=FALSE,
                     thedata=data.frame(NA),
                     rawdata=data.frame(NA),
                     normalizeddata=data.frame(NA),
                     status=list(inited=FALSE,
                                 ckexcluded=FALSE,
                                 logtransformed=FALSE,
                                 normalized=FALSE,
                                 scored=FALSE,
                                 annotated=FALSE),
                     header="empty",
                     layoutnames="NA"))


#################################################################
# html writing

whtml <- function(s,f,app=T){
  cat(s,"\n",sep="",file=f,append=app)
}


#################################################################
# Additional new plotting functions

## plot histogram of data
plothist <- function(thedata,pos,neg,title,xl,fname){
  jpeg(fname)
  if (sum(!is.na(thedata))>0){
    r <- hist(thedata,freq=T,breaks=100,add=F,main=title,xlab=xl,ylab="Frequency")
    if (length(pos)>0)
      hist(thedata[pos],freq=T,add=T,col="green",breaks=r$breaks,density=50)
    if (length(neg)>0)
      hist(thedata[neg],freq=T,add=T,col="red",breaks=r$breaks,density=50)
  }
  dev.off()
}

## make qqplot of data
plotqq <- function(thedata,pos,neg,title,xl,fname){
  jpeg(fname)
  if (sum(!is.na(thedata))>0){
    thecol <- rep("black",length(thedata))
    thecol[pos] <- "green"
    thecol[neg] <- "red"
    ix <- sort(thedata,index.return=T)$ix
    qqPlot(thedata,ylab=xl,main=title,col=thecol[ix])
  }
  dev.off()
}

## make overview boxplot
makeBoxplotOverview <- function(data,channel,i,fn)
{
  jpeg(fn)
  itsdata <- data@thedata
  itsdata <- itsdata[itsdata$ScreenNb == i,]
  boxplot(itsdata[[get("channel")]] ~ itsdata$LabtekNb,xlab="Plate",ylab=channel)
  dev.off()
}

## scatter plot signal vs cell counts
plotscatter <-function(ck,sig,pos,neg,fn,width=480,height=480){
  jpeg(fn,width,height)
  use <- !is.na(sig)
  theck <- ck[use]
  thesig <- sig[use]
  ctrl <- rep(0,length(sig))
  ctrl[pos] <- 1
  ctrl[neg] <- -1
  thepos <- which(ctrl[use] == 1)
  theneg <- which(ctrl[use] == -1)
  if (sum(!is.na(sig))>0){
    plot(theck,thesig,main="Signal vs. Cell Counts",xlab="Cellcounts",ylab="Signal")
    lines(theck[thepos],thesig[thepos],col="green",type="p")
    lines(theck[theneg],thesig[theneg],col="red",type="p")
    lines(lowess(theck,thesig,f=2/3),col="orange")
  }
  dev.off()
}

## plot of control distribution for experiment i, per plate
makeCtrlBoxPlot <- function(data,channel,i,fn)
{
  jpeg(fn)
  itsdata <- data@thedata
  itsdata <- itsdata[itsdata$ScreenNb == i,]
  itsdatas <- itsdata[[get("channel")]]
  posctrl <- which(itsdata$SpotType==1)
  negctrl <- which(itsdata$SpotType==0)
  if (length(!is.na(itsdatas))>0){
    boxplot(itsdatas~itsdata$LabtekNb)
    lines(itsdatas[posctrl]~itsdata$LabtekNb[posctrl],col="green",pch="o",type="p",)
    lines(itsdatas[negctrl]~itsdata$LabtekNb[negctrl],col="red",pch="o",type="p")
  }
  dev.off()
}


## plot of control density for experiment i
makeCtrlDensityPlot <- function(data,channel,i,fn)
{
  jpeg(fn)
  itsdata <- data@thedata
  itsdata <- itsdata[itsdata$ScreenNb == i,]
  itsdatas <- itsdata[[get("channel")]]
  posctrl <- itsdatas[itsdata$SpotType==1]
  negctrl <- itsdatas[itsdata$SpotType==0]
  if (length(!is.na(itsdatas))>0){
    mx <- max(density(itsdatas,na.rm=T)$y)
    hist(itsdatas,freq=F,breaks=50,col="grey",ylim=c(0,mx*1.5),main="Histogram",xlab=channel)
    if (sum(!is.na(posctrl))>1)
      lines(density(posctrl,na.rm=T),col="green",lwd=3)
    if (sum(!is.na(negctrl))>1)
      lines(density(negctrl,na.rm=T),col="red",lwd=3)
  }
  dev.off()
}

## plot image overview of entire screen
# makeScreenOverview <- function (data,fn) 
# {
#   plateno <- max(data@thedata$LabtekNb)
#   repno <- max(data@thedata$ScreenNb)
#   wellno <- max(data@thedata$Spotnumber)
#   xmax <- max(data@thedata$ColNb)
#   ymax <- max(data@thedata$RowNb)
#   xdim <- (xmax+5)*plateno
#   ydim <- (ymax+5)*repno
#   mat <- matrix(as.numeric(NA),nrow=xdim,ncol=ydim)
#   zrange <- range(data@thedata$SigIntensity,na.rm=TRUE)
#   rdbu = rev(brewer.pal(11, "RdBu"))[c(1:5, rep(6, 3),7:11)]                                                                                
#   colrs = colorRampPalette(rdbu)(256) 
#   # fill matrix with numbers
#   zr <- 256/(zrange[2]-zrange[1]);
#   for (i in 1:length(data@thedata$SigIntensity))
#   {
#     x <- (xmax+5)*(data@thedata$LabtekNb[i]-1)+data@thedata$ColNb[i]-1
#     y <- (ymax+5)*(data@thedata$ScreenNb[i]-1)+data@thedata$RowNb[i]-1
#     mat[x,y] <- round((data@thedata$SigIntensity[i]-zrange[1])*zr);
#   }
#   jpeg(fn,width=600,height=600)
#   image(mat,col=colrs);
#   dev.off()
# }

makeScreenOverview <- function (data,fn)
{
  repno <- max(data@thedata$ScreenNb)
  plateno <- max(data@thedata$LabtekNb)
  wellno <- max(data@thedata$Spotnumber)
  colno <- max(data@thedata$ColNb)
  rowno <- max(data@thedata$RowNb)
  wd <- plateno*(10*colno+10)
  ht <- repno*(10*rowno+10)
  jpeg(fn,width=wd,height=ht)
  z <- lapply(1:(plateno*repno), 
    function(i)
    {
      theplate = ((i-1) %% plateno) + 1
      therep = ((i-1) %/% plateno) + 1
      itsdata <- data@thedata[(data@thedata$LabtekNb==theplate)&(data@thedata$ScreenNb==therep),]
      itsdata$SigIntensity[itsdata$SpotType == -1] <- NA
      val <- rep(NA,wellno)
      for (j in 1:wellno)
        if (length(itsdata$SigIntensity[itsdata$Spotnumber==j])==1)
          val[j] <- (itsdata$SigIntensity[itsdata$Spotnumber==j])
      return(val)
    })
  plotScreen(z,ncol=plateno,nx=colno,ny=rowno)
  dev.off()
}

createReplicateMatrix <- function(data,channel)
{
  replicateno <- max(data@thedata$ScreenNb)
  spotno <- max(data@thedata$Spotnumber)
  mat <- matrix(NA,nrow=spotno,ncol=replicateno+1) # last column is SpotType!
  chdata <- data@thedata[[get("channel")]]
  for (j in 1:spotno){
    mat[j,replicateno+1] <- mean(data@thedata$SpotType[data@thedata$Spotnumber==j],na.rm=T)
    for (i in 1:replicateno){
      xd <- chdata[(data@thedata$Spotnumber==j)&(data@thedata$ScreenNb==i)]
      if (length(xd)==1)
        mat[j,i] <- xd
      else if (length(xd)>1)
        stop("**ERROR: Replacement length > 1 in createReplicateMatrix**")
    }
  }
  return(mat)
}

makeCorrPlots <- function(data,channel,fnprefix,htmlf)
{
  if (length(data@thedata$Spotnumber)==0)
    return();
  replicateno <- max(data@thedata$ScreenNb)
  spotno <- max(data@thedata$Spotnumber)
  repmat <- createReplicateMatrix(data,channel)
  whtml("<table>",htmlf)
  col <- rep("black",spotno)
  col[repmat[,replicateno+1]==1] <- "green"
  col[repmat[,replicateno+1]==0] <- "red"
  for (i in 1:replicateno){
    whtml("<tr>",htmlf)
    for (j in 1:replicateno){
      fn <- paste(fnprefix,"_",i,"x",j,".jpg",sep="")
      jpeg(fn,width=300,height=300)
      if (sum(!is.na(repmat[,i]))>0)
        if (sum(!is.na(repmat[,j]))>0)
          plot(repmat[,i],repmat[,j],col=col,pch="o",main=paste("Corr=",cor(repmat[,i],repmat[,j],use="c"),sep=""),
               xlab=paste("Replicate ",i,sep=""),ylab=paste("Replicate ",j,sep=""))
      dev.off()
      whtml(paste("<td><img src=\"",fn,"\"></img></td>",sep=""),htmlf)
    }
    whtml("</tr>",htmlf)
  }
  whtml("</table><br/>",htmlf)
}

PlotSpatialDistrib <- function(data,exp,plt,fn,channel)
{
  thedata <- data@thedata
  mCol <- max(thedata$ColNb)                           
  mRow <- max(thedata$RowNb)
  jpeg(fn,width=350,height=round(350*mRow/mCol));
  thedata[[get("channel")]][which(thedata$SpotType == -1)] <- NA_integer_
  subset <- thedata[which(thedata$ScreenNb == exp),]
  if (length(subset)>0){
    subsubset <- subset[which(subset$LabtekNb == plt),]
    if (sum(!is.na(subsubset[[get("channel")]])) >  0){ 
      posNegLabelVec <- rep(NA_character_, nrow(subsubset))   
      posNegLabelVec[which(subsubset$SpotType == 0)] <- "N"                                            
      posNegLabelVec[which(subsubset$SpotType == 1)] <- "P"                                            
      ColNo <- max(subsubset$ColNb)                           
      RowNo <- max(subsubset$RowNb)
      WellNo <- ColNo*RowNo
      plotPlate((as.numeric(t(subsubset[[get("channel")]]))),ncol = ColNo, nrow = RowNo, 
                char = posNegLabelVec, ind=subsubset$Spotnumber,
                main = paste(channel," Plate ",plt," Replicate ",exp, sep = " "), 
                col = brewer.pal(9,"YlOrBr"), na.action = "xout")
    }
  }
  dev.off()
}

plotmybar <- function(data,pos,neg,thefile1,thefile2,iszscore=FALSE,width=1600,height=600){
  if (length(data$Spotnumber)==0)
    len <- 0
  else
    len <- max(data$Spotnumber,na.rm=T)
  colors <- rep("black",len)
  colors[pos] <- "green"
  colors[neg] <- "red"
  med = median(data$SigIntensity,na.rm=T)
  medck = median(data$NbCells,na.rm=T)
  dev = mad(data$SigIntensity,na.rm=T)
  dev2 = 2*dev
  devck = mad(data$NbCells,na.rm=T)
  if (iszscore){
    med <- 0
    dev <- 1.65
    dev2 <- 2.3
  }
  if (is.na(med))
    med <- 0
  if (is.na(medck))
    medck <-0
  if (is.na(dev))
    dev <- 1
  if (is.na(devck))
    devck <- 1
  minim = max(min(as.numeric(data$SigIntensity),na.rm=TRUE),med-8*dev)
  maxim = min(max(as.numeric(data$SigIntensity),na.rm=TRUE),med+8*dev)
  minimck = max(min(as.numeric(data$NbCells),na.rm=TRUE),med-8*devck)
  maximck = min(max(as.numeric(data$NbCells),na.rm=TRUE),med+8*devck)
  if (minim == Inf)
    minim <- -5
  if (maxim == -Inf)
    maxim <- 5
  if (minimck == Inf)
    minimck <- -5
  if (maximck == -Inf)
    maximck <- 5
  jpeg(thefile1,width=width,height=height)
  if (sum(!is.na(data$SigIntensity))>0){
    boxplot(SigIntensity~Spotnumber,data=data,border=colors,col=colors,na.action=function(x){return(x)}) #,ylim=c(minim,maxim))
    abline(med,0,col="blue")
    abline(med+dev,0,col="green")
    abline(med-dev,0,col="green")
    abline(med+dev2,0,col="red")
    abline(med-dev2,0,col="red")
  }
  dev.off()
  jpeg(thefile2,width=width,height=height)  
  if (sum(!is.na(data$SigIntensity))>0){
    boxplot(NbCells~Spotnumber,data=data,border=colors,col=colors,na.action=function(x){return(x)}) #,ylim=c(minimck,maximck))
    abline(medck,0,col="blue")
    abline(medck+devck,0,col="green")
    abline(medck-devck,0,col="green")
    abline(medck+2*devck,0,col="red")
    abline(medck-2*devck,0,col="red")
  }
  dev.off()
}


##############################
# Vulcano Plot of Screen
makevolcanoplot <- function(data,fn)
{
  jpeg(fn,width=1600,height=600)
  col <- rep("black",length(data@thedata$SpotType))
  col[data@thedata$SpotType == -1] <- "grey"
  col[data@thedata$SpotType == 0]<-"red"                                  
  col[data@thedata$SpotType == 1]<-"green"                                  
  pch <- rep(".",length(data@thedata$SpotType))                                  
  pch[col!="black"] <- "o"
  pch[col=="grey"] <- "."
#  pch[col=="yellow"] <- "."
  ctrl <- which(data@thedata$SpotType %in% c(0,1))
  if (sum(is.na(data@thedata$pvalue))<length(data@thedata$pvalue)){
    plot(data@thedata$SigIntensity,data@thedata$pvalue,pch=pch,col=col,xlab="Score",ylab="p-Value")
    lines(data@thedata$SigIntensity[ctrl],data@thedata$pvalue[ctrl],pch=pch[ctrl],col=col[ctrl],type="p")
    abline(0.05,0,col="orange")
    abline(v=1.65,col="blue")
    abline(v=-1.65,col="blue")
    abline(v=2.3,col="magenta")
    abline(v=-2.3,col="magenta")
  }
  dev.off()
  return();
}

makefullbarplot <- function(data,fn)
{
  jpeg(fn,width=1600,height=600)
  idx <- (data@thedata$LabtekNb-1)*max(data@thedata$Spotnumber)+data@thedata$Spotnumber
  col <- rep("black",length(data@thedata$SpotType))
  col[data@thedata$SpotType == -1] <- "yellow"
  col[data@thedata$SpotType == 0]<-"red"                                  
  col[data@thedata$SpotType == 1]<-"green"                                  
  pch <- rep(".",length(data@thedata$SpotType))                                  
  pch[col!="black"] <- "o"
  pch[col=="yellow"] <- "x"
  ctrl <- which(data@thedata$SpotType %in% c(0,1))
  plot(idx,data@thedata$SigIntensity,pch=pch,col=col,xlab="Index",ylab="Score")
  lines(idx[ctrl],data@thedata$SigIntensity[ctrl],pch=pch[ctrl],col=col[ctrl],type="p")
  abline(1.65,0,col="blue")
  abline(-1.65,0,col="blue")
  abline(2.3,0,col="magenta")
  abline(-2.3,0,col="magenta")
  abline(0,0)
  lines(sort(data@thedata$SigIntensity),lw=2,col="black")
  dev.off()
 return();
}


##############################
# Summarize Replicates
summarizeReplicates <- function(mydata,NoOfLayouts,NoOfWells,test)
{
  itssize <- NoOfWells*NoOfLayouts
  maxcol <- max(mydata@thedata$ColNb)
  # first, create an empty RNAither data from of size NoOfLayouts x NoOfWells
  newdata <- data.frame(Spotnumber=rep(NA,itssize),
                               SpotType=rep(NA,itssize),
                               Internal_GeneID=rep(NA,itssize),
                               GeneName=rep(NA,itssize),
                               SigIntensity=rep(NA,itssize),
                               LabtekNb=rep(NA,itssize),
                               RowNb=rep(NA,itssize),
                               ColNb=rep(NA,itssize),
                               ScreenNb=rep(1,itssize),
                               NbCells=rep(NA,itssize),
                               pvalue=rep(NA,itssize),
                               meanscore=rep(NA,itssize),
                               mad=rep(NA,itssize),
                               sd=rep(NA,itssize))
  # Now fill this thing with the data!
  for (i in 1:NoOfLayouts){
    for (j in 1:NoOfWells){
      whodest <- ((i-1)*NoOfWells)+j
      whosrc <- which((mydata@normalizeddata$LabtekNb==i)&(mydata@normalizeddata$Spotnumber==j))
      ones <- whosrc[1]
      newdata$Spotnumber[whodest] <- j
      newdata$SpotType[whodest] <- max(mydata@thedata$SpotType[whosrc])
      newdata$Internal_GeneID[whodest] <- mydata@thedata$Internal_GeneID[ones]
      newdata$GeneName[whodest] <- mydata@thedata$GeneName[ones]
      xv <- mydata@thedata$SigIntensity[whosrc]
      newdata$SigIntensity[whodest] <- median(xv,na.rm=TRUE)
      newdata$meanscore[whodest] <- mean(xv,na.rm=TRUE)
      newdata$mad[whodest] <- mad(xv,na.rm=TRUE)
      newdata$sd[whodest] <- sd(xv,na.rm=TRUE)
      newdata$LabtekNb[whodest] <- i
      newdata$RowNb[whodest] <- ((j-1) %/% maxcol) + 1
      newdata$ColNb[whodest] <- ((j-1) %% maxcol) + 1
      newdata$NbCells[whodest] <- median(mydata@thedata$NbCells[whosrc],na.rm=TRUE)
      if (sum(!is.na(xv))>1){
         if (min(xv,na.rm=T)==max(xv,na.rm=T))
           newdata$pvalue[whodest] <- NA
         else if (test=="ttest")
           newdata$pvalue[whodest] <- t.test(xv,mu=0,alternative="two.sided")$p.value
         else if (test=="wilcox")
           newdata$pvalue[whodest] <- wilcox.test(xv,mu=0,alternative="two.sided")$p.value
      }  
    }
  }
  mydata@thedata <- newdata
  return(mydata)
}

#################################################################
# Create plots for individual Layout!
makePlateDetailPlots <- function(data,layoutnumber,directory,iszscore=FALSE)
{
  dir.create(directory)
  currdir <- getwd()
  setwd(directory)
  
  replicateno <- max(data@thedata$ScreenNb)
  thename <- paste("Layout ",layoutnumber,sep="")
  if (length(data@layoutnames)>0)
    thename <- data@layoutnames[layoutnumber]
  # write html file header
  f <- "index.html"
  whtml("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">",f,F)
  whtml("<HTML>",f)
  whtml("<HEAD>",f)
  whtml("  <META HTTP-EQUIV=\"CONTENT-TYPE\" CONTEN=\"text/html; charset=utf-8\">",f)
  whtml(paste("  <TITLE>",data@expname," - Layout ",layoutnumber,"</TITLE>",sep=""),f)
  whtml("</HEAD>",f)
  whtml("<BODY LANG=\"en-US\" DIR=\"LTR\">",f)
#  whtml(paste("  <H1 ALIGN=\"CENTER\">",data@expname," - ",title,"</H1>",sep=""),f)
#  whtml(paste("  <P ALIGN=CENTER STYLE=\"margin-top: 0.17in; page-break-after: avoid\"><FONT FACE=\"Albany, sans-serif\"><FONT SIZE=4>",date(),"</FONT></FONT></P>",sep=""),f)
  whtml(paste("  <H1 ALIGN=CENTER>",data@expname," - Layout ",layoutnumber,"(",thename,")</H1>"),f)
  whtml("  <HR>",f)

  # retrieve data for this layout, and find the controls
  procdata <- data
  procdata@thedata <- data@thedata[data@thedata$LabtekNb == layoutnumber,]
  pos <- which(procdata@thedata$SpotType == 1)
  neg <- which(procdata@thedata$SpotType == 0)

  # data without "bad" values
  cleaneddata <- procdata
  cleaneddata@thedata <- procdata@thedata[procdata@thedata$SpotType != -1,]
  poscl <- which(cleaneddata@thedata$SpotType == 1)
  negcl <- which(cleaneddata@thedata$SpotType == 0)

  # now create the individual plots...
  #  1 full distribution of data (histograms), raw and log-transformed
  message("\n     - Histograms...",appendLF = F)
  plothist(cleaneddata@thedata$NbCells,poscl,negcl,"Cell Counts","Cell Number","1_ck.jpg")
  plothist(log(cleaneddata@thedata$NbCells),poscl,negcl,"Log Cell Counts","Log Cell Number","1_cklog.jpg")
  plothist(cleaneddata@thedata$SigIntensity,poscl,negcl,"Signal","Signal","1_sig.jpg")
  plothist(log(cleaneddata@thedata$SigIntensity),poscl,negcl,"Log Signal","Log Signal","1_siglog.jpg")
  whtml("<h2>Data Overview</h2>",f)
  whtml("<table>",f)
  whtml("<tr><th colspan=\"2\">Cell Counts</th><th colspan=\"2\">Signal</th></tr>",f)
  whtml("<tr><td><img src=\"1_ck.jpg\"></img></td>",f)
  whtml("    <td><img src=\"1_cklog.jpg\"></img></td>",f)
  whtml("    <td><img src=\"1_sig.jpg\"></img></td>",f)
  whtml("    <td><img src=\"1_siglog.jpg\"></img></td></tr>",f)

  #  2 qq-plot of raw and log data
  message("\n     - QQ-Plots...",appendLF = F)
  plotqq(cleaneddata@thedata$NbCells,poscl,negcl,"Cell Counts","Cell Number","2_ck.jpg")
  plotqq(log(cleaneddata@thedata$NbCells),poscl,negcl,"Log Cell Counts","Log Cell Number","2_cklog.jpg")
  plotqq(cleaneddata@thedata$SigIntensity,poscl,negcl,"Signal","Signal","2_sig.jpg")
  plotqq(log(cleaneddata@thedata$SigIntensity),poscl,negcl,"Log Signal","Log Signal","2_siglog.jpg")
  whtml("<tr><td><img src=\"2_ck.jpg\"></img></td>",f)
  whtml("    <td><img src=\"2_cklog.jpg\"></img></td>",f)
  whtml("    <td><img src=\"2_sig.jpg\"></img></td>",f)
  whtml("    <td><img src=\"2_siglog.jpg\"></img></td></tr>",f)
  whtml("</table>",f)  

  #  3 Boxplot of controls
  message("\n     - Boxplots of Controls...",appendLF = F)
  #    A6b Controls per Plate and replicate
  for (i in 1:replicateno){
    fn1 <- paste("3_boxes_sig_exp",i,".jpg",sep="")
    fn2 <- paste("3_density_sig_exp",i,".jpg",sep="")
    fn3 <- paste("3_boxes_ck_exp",i,".jpg",sep="")
    fn4 <- paste("3_density_ck_exp",i,".jpg",sep="")
    makeCtrlBoxPlot(cleaneddata,"SigIntensity",i,fn1);
    makeCtrlDensityPlot(cleaneddata,"SigIntensity",i,fn2);
    makeCtrlBoxPlot(cleaneddata,"NbCells",i,fn3);
    makeCtrlDensityPlot(cleaneddata,"NbCells",i,fn4);
    whtml(paste("  <h3>Replicate ",i,"</h3>",sep=""),f)
    whtml("<table>",f)
    whtml("  <tr><th>Cellcounts</th><th>Signal Intensity</th><th>Cellcounts</th><th>Signal Intensity</th></tr>",f)
    whtml(paste("  <tr><td><img src=\"",fn3,"\"></img></td><td><img src=\"",fn1,"\"></img></td>",sep=""),f)
    whtml(paste("      <td><img src=\"",fn4,"\"></img></td><td><img src=\"",fn2,"\"></img></td></tr>",sep=""),f)
    whtml("</table>",f)
  }
  whtml("<BR><HR><BR>",f)

  #  4 correlation between signal and cell counts
  message("\n     - Correlation between channels...",appendLF = F)
  whtml("<H2>Correlation Between Channels</H2>",f)
  whtml("<table><tr>",f)
  for (i in 1:replicateno)
    whtml(paste("<th>Replicate ",i,"</th>",sep=""),f)
  whtml("</tr><tr>",f)
  for (i in 1:replicateno)
  {
    ck <- cleaneddata@thedata$NbCells[cleaneddata@thedata$ScreenNb==i]
    sig <- cleaneddata@thedata$SigIntensity[cleaneddata@thedata$ScreenNb==i]
    mpos <- which(cleaneddata@thedata$SpotType[cleaneddata@thedata$ScreenNb==i] == 1)
    mneg <- which(cleaneddata@thedata$SpotType[cleaneddata@thedata$ScreenNb==i] == 0)
    fn <- paste("4_scatter_rep",i,".jpg",sep="")
    plotscatter(ck,sig,mpos,mneg,fn,width=300,height=300)
    whtml(paste("<td><img src=\"",fn,"\"></img></td>",sep=""),f)
  }
  whtml("</tr></table>",f)

  #  5 plate plots signal and cell counts
  message("\n     - Plate plots...",appendLF = F)
  whtml("<HR><H2>Spatial Plate Plots</H3>",f)
  whtml("<table>",f)
  whtml("<tr><th> </th>",f)
  for (i in 1:replicateno)
    whtml(paste("<th>Replicate ",i,"</th>",sep=""),f)
  whtml("</tr><tr><th>Cell Counts</th>",f)
  for (i in 1:replicateno){
    fn <- paste("5_ck_",i,".jpg",sep="")
    PlotSpatialDistrib(procdata,i,layoutnumber,fn,"NbCells")
    whtml(paste("<td><img src=\"",fn,"\"></img></td>",sep=""),f)
  }
  whtml("</tr><tr><th>Signal</th>",f)
  for (i in 1:replicateno){
    fn <- paste("5_sig_",i,".jpg",sep="")
    PlotSpatialDistrib(procdata,i,layoutnumber,fn,"SigIntensity")
    whtml(paste("<td><img src=\"",fn,"\"></img></td>",sep=""),f)
  }
  whtml("</tr></table>",f)

  #  6 bar plots signal and cell counts - individual and joining replicates
  message("\n     - Bar plots...",appendLF = F)
  plotmybar(cleaneddata@thedata,poscl,negcl,"6_sig.jpg","6_ck.jpg",iszscore)
  whtml("<br/><hr><h2>Boxplot of Signal and Cell Counts over Replicates</h2>",f)
  whtml("<h3>Cell Counts</h3>",f)
  whtml("<img src=\"6_ck.jpg\"></img>",f)
  whtml("<h3>Signal Intensity</h3>",f)
  whtml("<img src=\"6_sig.jpg\"></img>",f)
  whtml("<br/>",f)
  whtml("<table>",f)
  whtml("<tr><th> </th>",f)
  for (i in 1:replicateno)
    whtml(paste("<th>Replicate ",i,"</th>",sep=""),f)
  whtml("</tr><tr><th>Cell Counts</th>",f)
  for (i in 1:replicateno){
    itsdata <- cleaneddata@thedata[cleaneddata@thedata$ScreenNb==i,]
    f1 <- paste("6_sig_",i,".jpg",sep="")
    f2 <- paste("6_ck_",i,".jpg",sep="")
    itspos <- which(itsdata$SpotType==1)
    itsneg <- which(itsdata$SpotType==0)
    plotmybar(itsdata,itspos,itsneg,f1,f2,iszscore,width=500,height=400)
    whtml(paste("<td><img src=\"",f2,"\"></img></td>",sep=""),f)
  }
  whtml("</tr><tr><th>Signal</th>",f)
  for (i in 1:replicateno){
    f1 <- paste("6_sig_",i,".jpg",sep="")
    whtml(paste("<td><img src=\"",f1,"\"></img></td>",sep=""),f)
  }
  whtml("</tr></table>",f)

  #  7 correlation between replicates (Signal and cell counts)
  message("\n     - Correlation plots...",appendLF = F)
  whtml("<HR><H2>Correlation of Cell Counts between Replicates</H2>",f)
  makeCorrPlots(cleaneddata,"NbCells","7_cells",f)
  whtml("<H2>Correlation of Signal Intensities between Replicates</H2>",f)
  makeCorrPlots(cleaneddata,"SigIntensity","7_sig",f)

  # write html file end
  whtml(" <BR/><HR>",f)
  whtml("  Return to <a href=\"../index.html\">previous page</a>",f)
  whtml("  <HR>\n",f)
  whtml("</BODY>",f)
  whtml("</HTML>",f)
  setwd(currdir)
  return();
}

#################################################################
# Create plots showing data
makeQCPlots <- function(data,directory,title,iszscore=FALSE,makeplots=TRUE)
{
  if (!makeplots)
    return();
  # creates a subdirectory "directory", in which an html file "index.html" is generated
  # containing details about the data passed in "data". The "title" will be used to specify
  # contents of the file further.
  # makeQCPlots then creates several subplots, which are included in the html page:
  #  * for the entire screen:
  #    1 histogram of data (signal and cell count), if data@status$logtransformed=FALSE also
  #      including a histogram after log-transformation.
  #    2 A qq-plot of the full dataset, raw and possibly also log-transformed
  #    3 boxplots including all plates and replicates
  #    4 correlation between signal and cell counts
  #    5 MA-Plot
  #    6 boxplots of controls
  #    7 color-coded overview of entire screen
  #  * for each individual plate:
  #    1 full distribution of data (histograms), raw and log-transformed
  #    2 qq-plot of raw and log data
  #    3 Boxplot of controls
  #    4 correlation between signal and cell counts
  #    5 plate plots signal and cell counts
  #    6 bar plots signal and cell counts - individual and joining replicates
  #    7 correlation between replicates (Signal and cell counts)
  #    8 M-A plots

  dir.create(directory)
  currdir <- getwd()
  setwd(directory)

  plateno <- max(data@thedata$LabtekNb)
  replicateno <- max(data@thedata$ScreenNb)

  # write html file header
  f <- "index.html"
  whtml("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">",f,F)
  whtml("<HTML>",f)
  whtml("<HEAD>",f)
  whtml("  <META HTTP-EQUIV=\"CONTENT-TYPE\" CONTEN=\"text/html; charset=utf-8\">",f)
  whtml(paste("  <TITLE>",data@expname," - ",title,"</TITLE>",sep=""),f)
  whtml("</HEAD>",f)
  whtml("<BODY LANG=\"en-US\" DIR=\"LTR\">",f)
#  whtml(paste("  <H1 ALIGN=\"CENTER\">",data@expname," - ",title,"</H1>",sep=""),f)
#  whtml(paste("  <P ALIGN=CENTER STYLE=\"margin-top: 0.17in; page-break-after: avoid\"><FONT FACE=\"Albany, sans-serif\"><FONT SIZE=4>",date(),"</FONT></FONT></P>",sep=""),f)
  whtml(paste("  <H1 ALIGN=CENTER>",data@expname," - ",title,"</H1>"),f)
  whtml("  <HR>",f)

  # where are the controls?
  procdata <- data
  procdata@thedata <- data@thedata[data@thedata$SpotType != -1,]
  pos <- which(procdata@thedata$SpotType == 1)
  neg <- which(procdata@thedata$SpotType == 0)

  # A1. Create histograms of raw data, signal and cell count, with and without log-transformation
  # A2. Create QQ-Plot of raw data, signal and cell count, with and without log-transformation
  message("\n - Histograms and QQ-Plots...",appendLF = F)
  whtml("<H2>Data Overview</H2>",f)
  plothist(procdata@thedata$NbCells,pos,neg,"Cell Counts","Cell Number","A1_ck.jpg")
  plotqq(procdata@thedata$NbCells,pos,neg,"Cell Counts","Cell Number","A2_ck.jpg")
  plothist(procdata@thedata$SigIntensity,pos,neg,"Signal","Signal Intensity","A1_sig.jpg")
  plotqq(procdata@thedata$SigIntensity,pos,neg,"Signal","Signal Intensity","A2_sig.jpg")
  if (procdata@status$logtransformed == FALSE)
  {
    plothist(log(procdata@thedata$SigIntensity),pos,neg,"Log Signal","ln(Signal Intensity)","A1_logsig.jpg")
    plotqq(log(procdata@thedata$SigIntensity),pos,neg,"Log Signal","ln(Signal Intensity)","A2_logsig.jpg")
    whtml("  <table>",f)
    whtml("    <tr><td><img src=\"A1_ck.jpg\"></img></td>",f)
    whtml("        <td><img src=\"A1_sig.jpg\"></img></td>",f)
    whtml("        <td><img src=\"A1_logsig.jpg\"></img></td>",f)
    whtml("    <tr><td><img src=\"A2_ck.jpg\"></img></td>",f)
    whtml("        <td><img src=\"A2_sig.jpg\"></img></td>",f)
    whtml("        <td><img src=\"A2_logsig.jpg\"></img></td>",f)
    whtml("    </tr>",f)
    whtml("  </table>",f)
  }
  else
  {
    whtml("  <table>",f)
    whtml("    <tr>",f)
    whtml("        <td><img src=\"A1_ck.jpg\"></img></td>",f)
    whtml("        <td><img src=\"A1_sig.jpg\"></img></td>",f)
    whtml("    </tr>",f)
    whtml("    <tr>",f)
    whtml("        <td><img src=\"A2_ck.jpg\"></img></td>",f)
    whtml("        <td><img src=\"A2_sig.jpg\"></img></td>",f)
    whtml("    </tr>",f)
    whtml("  </table>",f)
  }

  #    4 correlation between signal and cell counts
  whtml("<br/>",f)
  message("\n - Scatter Plots...",appendLF = F)
  plotscatter(procdata@thedata$NbCells,procdata@thedata$SigIntensity,pos,neg,"A4_scatter.jpg")
  whtml("<img src=\"A4_scatter.jpg\"></img>",f)

  # A3. Create a Box-Plot of signal and cell counts of the entire screen
  # Loop over all replicates, and show a boxplot showing the different plates within a replicate!
  #whtml("<BR/><HR><BR/>",f)
  #cat("\n - Boxplots of Replica and Plates...")
  #whtml("<H2>Distribution of Replicates and Layouts</H2>",f)
  #for (i in 1:replicateno){
  #  fn1 <- paste("A3_ck_exp",i,".jpg",sep="")
  #  fn2 <- paste("A3_sig_exp",i,".jpg",sep="")
  #  makeBoxplotOverview(procdata,"NbCells",i,fn1)
  #  makeBoxplotOverview(procdata,"SigIntensity",i,fn2)
  #  whtml(paste("  <h3>Replicate ",i,"</h3>",sep=""),f)
  #  whtml("<table>",f)
  #  whtml("  <tr><th>Cellcounts</th><th>Signal</th></tr>",f)
  #  whtml(paste("  <tr><td><img src=\"",fn1,"\"></img></td><td><img src=\"",fn2,"\"></img></td></tr>",sep=""),f)
  #  whtml("</table>",f)
  #}
  #whtml("<HR>",f)

  #    A6 boxplots of controls
  message("\n - Boxplot of Replica and Plates, including Controls...",appendLF = F)
  whtml("<BR/><HR><BR/>",f)
  whtml("<H2>Boxplot of Plates per Replicate, and Controls</h2>",f)
  # first of all, a boxplot of the controls for the entire experiment
  pdata <- procdata@thedata
  tempVec <- rep(0, length(pdata$SpotType))
  tempVec[pdata$SpotType == 0] <- "Neg. contr."
  tempVec[pdata$SpotType == 1] <- "Pos. contr."
  tempVec[pdata$SpotType == 2] <- "Exp. data"
  jpeg("A6_controls_sig.jpg")
  boxplot(pdata$SigIntensity ~ tempVec,
        ylab = "Signal Intensity", cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6,
	col=c("grey","red","green"))
  dev.off()
  jpeg("A6_controls_ck.jpg")
  boxplot(pdata$NbCells ~ tempVec,
        ylab = "Cell Counts", cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.6,
	col=c("grey","red","green"))
  dev.off()
  whtml("<table>",f)
  whtml("  <tr><th>Cellcounts</th><th>Signal Intensity</th></tr>",f)
  whtml("  <tr><td><img src=\"A6_controls_ck.jpg\"></img></td>",f)
  whtml("      <td><img src=\"A6_controls_sig.jpg\"></img></td></tr>",f)
  whtml("</table><br/>",f)

  #    A6b Controls per Plate and replicate
  for (i in 1:replicateno){
    fn1 <- paste("A6b_boxes_sig_exp",i,".jpg",sep="")
    fn2 <- paste("A6b_density_sig_exp",i,".jpg",sep="")
    fn3 <- paste("A6b_boxes_ck_exp",i,".jpg",sep="")
    fn4 <- paste("A6b_density_ck_exp",i,".jpg",sep="")
    makeCtrlBoxPlot(procdata,"SigIntensity",i,fn1);
    makeCtrlDensityPlot(procdata,"SigIntensity",i,fn2);
    makeCtrlBoxPlot(procdata,"NbCells",i,fn3);
    makeCtrlDensityPlot(procdata,"NbCells",i,fn4);
    whtml(paste("  <h3>Replicate ",i,"</h3>",sep=""),f)
    whtml("<table>",f)
    whtml("  <tr><th>Cellcounts</th><th>Signal Intensity</th><th>Cellcounts</th><th>Signal Intensity</th></tr>",f)
    whtml(paste("  <tr><td><img src=\"",fn3,"\"></img></td><td><img src=\"",fn1,"\"></img></td>",sep=""),f)
    whtml(paste("      <td><img src=\"",fn4,"\"></img></td><td><img src=\"",fn2,"\"></img></td></tr>",sep=""),f)
    whtml("</table>",f)
  }
  whtml("<BR><HR><BR>",f)
  
  #    7 color-coded overview of entire screen
  whtml("<H2>Overview over entire screen</H2>",f)
  message("\n - Screen overview...",appendLF = F)
  whtml("<img src=\"A7_overview.jpg\"></img>",f)
  whtml("<BR><HR><BR>",f)
  makeScreenOverview(data,"A7_overview.jpg") 
  
  ###    5 MA-Plot ????

  # create detailed plots for the individual plates
  message("\n - Plate overview...",appendLF = F)
  whtml("<H2>Plate details</H2>",f)
  whtml("<BR/>",f)
  whtml("<p><ul>",f) 
  for (pnumber in 1:plateno){
    message(paste("\n   - Plate ",pnumber,"...",sep=""),appendLF = F)
    fn <- paste("plate_",pnumber,sep="")
    makePlateDetailPlots(data,pnumber,fn,iszscore)
    if (length(procdata@layoutnames==pnumber))
      whtml(paste(" <li> <a href=\"",fn,"/index.html\">Layout ",pnumber," (",procdata@layoutnames[pnumber],")</a></li>",sep=""),f)
    else
      whtml(paste(" <li> <a href=\"",fn,"/index.html\">Layout ",pnumber,"</a></li>",sep=""),f)
  }
  whtml("</ul></p>",f)
  
  # write html file end
  whtml(" <BR/><HR>",f)
  whtml("  Return to <a href=\"../index.html\">main page</a>",f)
  whtml("  <HR>\n",f)
  whtml("</BODY>",f)
  whtml("</HTML>",f)
  setwd(currdir)
  return();
}

#################################################################
# write hit list to html file and to text file
writehitlist <- function(ledata,f)
{
  data <- ledata
  data@thedata <- ledata@thedata[ledata@thedata$SpotType != -1,]
  data@normalizeddata <- ledata@normalizeddata[ledata@normalizeddata$SpotType != -1,] 
  whtml("<h2>Hit lists</h2>",f)
  whtml("<p><ul>",f)
  whtml("  <li><a href=\"sirna.html\">Hits at siRNA level</a></li>",f)
  whtml("  <li><a href=\"gene.html\">Hits at Gene level</a></li>",f)  
  whtml("</ul></p>",f)

  # create hit list at siRNA level
  fullsirna <- data@thedata
  fullsirna$ZScore <- as.numeric(data@thedata$SigIntensity)
  fullsirna <- fullsirna[!is.na(fullsirna$ZScore),]
  sirnahits <- fullsirna[(abs(fullsirna$ZScore) >= data@scorethresh) & (fullsirna$pvalue <= data@pvalthresh),]
  h1 <- "sirna.html"
  whtml("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">",h1,F)
  whtml("<HTML>",h1)
  whtml("<HEAD>",h1)
  whtml("  <META HTTP-EQUIV=\"CONTENT-TYPE\" CONTEN=\"text/html; charset=utf-8\">",h1)
  whtml(paste("  <TITLE>",data@expname," - Hitlist (siRNA)</TITLE>",sep=""),h1)
  whtml("</HEAD>",h1)
  whtml("<BODY LANG=\"en-US\" DIR=\"LTR\">",h1)
  whtml(paste("  <H1 ALIGN=\"CENTER\">",data@expname," - Hitlist (siRNA Level)</H1>",sep=""),h1)
  whtml(paste("  <P ALIGN=CENTER STYLE=\"margin-top: 0.17in; page-break-after: avoid\"><FONT FACE=\"Albany, sans-serif\"><FONT SIZE=4>",date(),"</FONT></FONT></P>",sep=""),h1)
  whtml("  <HR>",h1)
  whtml(paste("  <p>Hits are marked if abs(zscore)>=",data@scorethresh," and pvalue<=",data@pvalthresh,".</p>",sep=""),h1)
  whtml("<table>",h1)
  whtml(" <tr><th>Layout</th><th>Well</th><th>siRNA ID</th><th>Gene Name</th><th>Median Score</th><th>MAD Score</th><th>Mean Score</th><th>SD Score</th><th>p-value</th></tr>",h1)
  ix <- sort(fullsirna$ZScore,index.return=TRUE)$ix
  for (i in ix){
    if (!is.na(fullsirna$ZScore[i])&(!is.na(fullsirna$pvalue[i]))){
      if ((abs(fullsirna$ZScore[i])>data@scorethresh)&(fullsirna$pvalue[i]<data@pvalthresh))
      {
         whtml(paste("<tr><td><font color=\"red\">",fullsirna$LabtekNb[i],"</font></td><td><font color=\"red\">",fullsirna$Spotnumber[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",fullsirna$Internal_GeneID[i],"</font></td><td><font color=\"red\">",fullsirna$GeneName[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",fullsirna$ZScore[i],"</font></td><td><font color=\"red\">",fullsirna$mad[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",fullsirna$meanscore[i],"</font></td><td><font color=\"red\">",fullsirna$sd[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",fullsirna$pvalue[i],"</font></td></tr>",sep=""),h1)
      } else {
         whtml(paste("<tr><td>",fullsirna$LabtekNb[i],"</td><td>",fullsirna$Spotnumber[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$Internal_GeneID[i],"</td><td>",fullsirna$GeneName[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$ZScore[i],"</td><td>",fullsirna$mad[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$meanscore[i],"</td><td>",fullsirna$sd[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$pvalue[i],"</td></tr>",sep=""),h1)
      }
    } else{
         whtml(paste("<tr><td>",fullsirna$LabtekNb[i],"</td><td>",fullsirna$Spotnumber[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$Internal_GeneID[i],"</td><td>",fullsirna$GeneName[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$ZScore[i],"</td><td>",fullsirna$mad[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$meanscore[i],"</td><td>",fullsirna$sd[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",fullsirna$pvalue[i],"</td></tr>",sep=""),h1)
    }
  }
  whtml(" </table><BR/><HR>",h1)
  whtml("  Return to <a href=\"index.html\">main page</a>",h1)
  whtml("  <HR>\n",h1)
  whtml("</BODY>",h1)
  whtml("</HTML>",h1)
  write.table(fullsirna,file="sirna_scores.txt",row.names=F,col.names=T,sep="\t",quote=F)
  
  # create hit list at gene level
  pvals <- Ttest(data@normalizeddata,list("two.sided",0,"SigIntensity","GeneName"))
  pvalvec <- pvals[[1]]
  scoredset <- pvals[[2]]
  upscoredHits <- hitselectionZscorePval(scoredset,pvalvec,"SigIntensity","Zscore_pval_hits",data@scorethresh,data@pvalthresh,2,"GeneName",median,"gene_scores_up.txt")
  dnscoredHits <- hitselectionZscorePval(scoredset,pvalvec,"SigIntensity","Zscore_pval_hits",data@scorethresh,data@pvalthresh,-2,"GeneName",median,"gene_scores_dn.txt")
  write.table(scoredset,file="gene_scores.txt",row.names=F,col.names=T,sep="\t",quote=F)
  genes <- unique(scoredset$GeneName)
  scores <- rep(NA,length(genes))
  meanscore <- scores
  madscore <- scores
  sdscore <- scores
  pval <- rep(NA,length(genes))
  for (i in 1:length(genes)){
    xv <- scoredset$SigIntensity[scoredset$GeneName==genes[i]]
#    scores[i] <- mean(xv,na.rm=T)
    scores[i] <- median(xv,na.rm=T)
    meanscore[i] <- mean(xv,na.rm=T)
    madscore[i] <- mad(xv,na.rm=T)
    sdscore[i] <- sd(xv,na.rm=T)
    if (sum(!is.na(xv))>1){
      if (data@test=="ttest")
        pval[i] <- t.test(xv,mu=0,alternative="two.sided")$p.value
      else if (data@test=="wilcox")
        pval[i] <- wilcox.test(xv,mu=0,alternative="two.sided")$p.value
    }
  }
  wm <- (!is.na(scores))
  hits <- data.frame(gene=genes[wm],score=scores[wm],mad=madscore[wm],meanscore=meanscore[wm],sd=sdscore[wm],pval=pval[wm])
  write.table(hits,file="gene_scores.txt",row.names=F,col.names=T,sep="\t",quote=F)
  ix <- sort(hits$score,index.return=TRUE)$ix
  h1 <- "gene.html"
  whtml("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">",h1,F)
  whtml("<HTML>",h1)
  whtml("<HEAD>",h1)
  whtml("  <META HTTP-EQUIV=\"CONTENT-TYPE\" CONTEN=\"text/html; charset=utf-8\">",h1)
  whtml(paste("  <TITLE>",data@expname," - Hitlist (Gene Name)</TITLE>",sep=""),h1)
  whtml("</HEAD>",h1)
  whtml("<BODY LANG=\"en-US\" DIR=\"LTR\">",h1)
  whtml(paste("  <H1 ALIGN=\"CENTER\">",data@expname," - Hitlist (Gene Level)</H1>",sep=""),h1)
  whtml(paste("  <P ALIGN=CENTER STYLE=\"margin-top: 0.17in; page-break-after: avoid\"><FONT FACE=\"Albany, sans-serif\"><FONT SIZE=4>",date(),"</FONT></FONT></P>",sep=""),h1)
  whtml("  <HR>",h1)
  whtml(paste("  <p>Hits are marked if abs(zscore)>=",data@scorethresh," and pvalue<=",data@pvalthresh,".</p>",sep=""),h1)
  whtml("<table>",h1)
  whtml("  <tr><th>Gene Name</th><th>Median Score</th><th>MAD Score</th><th>Mean Score</th><th>SD Score</th><th>p-value</th></tr>",h1)
  for (i in ix){
    if (!is.na(hits$score[i])&(!is.na(hits$pval[i]))){
      if ((abs(hits$score[i])>data@scorethresh)&(hits$pval[i]<data@pvalthresh))
      {
         whtml(paste("<tr><td><font color=\"red\">",hits$gene[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",hits$score[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",hits$mad[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",hits$meanscore[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",hits$sd[i],"</font></td>",sep=""),h1)
         whtml(paste("    <td><font color=\"red\">",hits$pval[i],"</font></td></tr>",sep=""),h1)
      } else {
         whtml(paste("<tr><td>",hits$gene[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$score[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$mad[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$meanscore[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$sd[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$pval[i],"</td></tr>",sep=""),h1)
      }
    } else{
         whtml(paste("<tr><td>",hits$gene[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$score[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$mad[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$meanscore[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$sd[i],"</td>",sep=""),h1)
         whtml(paste("    <td>",hits$pval[i],"</td></tr>",sep=""),h1)
    }
  }
  whtml(" </table><BR/><HR>",h1)
  whtml("  Return to <a href=\"index.html\">main page</a>",h1)
  whtml("  <HR>\n",h1)
  whtml("</BODY>",h1)
  whtml("</HTML>",h1)
  
  return();
}


writegoresults <- function(data,f)
{
  whtml("<h2>Results - Pathway Analysis</h2>",f)
  whtml("<p>Pathway Analysis has been skipped.</p>",f)
  return();
}

#################################################################
# Create main html file
writeSummaryHtml <- function(data)
{
  # creates a file "index.html" in current directory, summarizing normalization procedure
  # used and linking to all results.
  f <- "index.html"
  whtml("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">",f,F)
  whtml("<HTML>",f)
  whtml("<HEAD>",f)
  whtml("  <META HTTP-EQUIV=\"CONTENT-TYPE\" CONTEN=\"text/html; charset=utf-8\">",f)
  whtml(paste("  <TITLE>",data@expname," - Statistical Analysis</TITLE>",sep=""),f)
  whtml("</HEAD>",f)
  whtml("<BODY LANG=\"en-US\" DIR=\"LTR\">",f)
  whtml(paste("  <H1 ALIGN=\"CENTER\">",data@expname," - Statistical Analysis</H1>",sep=""),f)
  whtml(paste("  <P ALIGN=CENTER STYLE=\"margin-top: 0.17in; page-break-after: avoid\"><FONT FACE=\"Albany, sans-serif\"><FONT SIZE=4>",date(),"</FONT></FONT></P>",sep=""),f)
  whtml("  <HR>",f)
  whtml("<h2>Normalization procedure used:</h2>",f)
  whtml("<p>Analyzed in R using the Bioconductor package RNAither.</p>",f)
  whtml("<p><b>Please cite: </b><br>",f)
  whtml("N. Rieber, B. Knapp, R. Eils, L. Kaderali (2009). <em>RNAither, an automated pipeline for the statistical analysis of high-throughput RNAi screens</em>. <b>Bioinformatics</b> 25, 678-679.</p>",f)
  whtml("<br>",f)
  whtml("<h3>The following steps were done in normalization:</h3>",f)
  whtml("<p><ol>",f)
  whtml("  <li>Read in raw data</li>",f)
  if (data@excludeCellcounts=="none")
    whtml("  <li>No exclusion of wells based on cell counts</li>",f)
  else if (data@excludeCellcounts=="lowest")
    whtml("  <li>Exclusion of 5% of wells with lowest cell counts in entire screen</li>",f)
  else if (data@excludeCellcounts=="both")
    whtml("  <li>Exclusion of 5% of wells with lowest and 5% of highest cell counts in entire screen</li>",f)
  else if (data@excludeCellcounts=="lowestperplate")
    whtml("  <li>Exclusion of 5% of wells with lowest cell counts per plate</li>",f)
  else if (data@excludeCellcounts=="bothperplate")
    whtml("  <li>Exclusion of 5% of wells with lowest and 5% of highest cell counts per plate</li>",f)
  if (data@logtransform)
    whtml("  <li>Log-Transformation of data</li>",f)
  for (i in data@normalization)
    whtml(paste("  <li>",i," normalization</li>",sep=""),f)
  whtml("  <li>Summarizing replicates using the mean</li>",f)
  whtml(paste("  <li>Assess statistical significance using statistical test: ",data@test,"</li>",sep=""),f)
  whtml(paste("  <li>Hit calling based on score threshold ",data@scorethresh," and p-value threshold ",data@pvalthresh,"</li>",sep=""),f)
  whtml("</ol>",f)
  whtml("<hr>",f)
  writehitlist(data,f)
  whtml("<hr>",f)
  writegoresults(data,f)
  whtml("<hr>",f)
  whtml("<h2>Data Overview</h2>",f)
  whtml("<table>",f)
  whtml("<tr><th>Cell Counts</th><th>Signal Raw</th><th>Signal Log</th><th>Signal Normalized</th><th>Signal Summarized</th></tr>",f)
  whtml("<tr><td><img src=\"raw/A1_ck.jpg\"></img></td>",f)
  whtml("    <td><img src=\"raw/A1_sig.jpg\"></img></td>",f)
  whtml("    <td><img src=\"raw/A1_logsig.jpg\"></img></td>",f)
  whtml(paste("    <td><img src=\"",data@normalization[length(data@normalization)],"/A1_sig.jpg\"></img></td>",sep=""),f)
  whtml("    <td><img src=\"summarized/A1_sig.jpg\"></img></td></tr>",f)
  whtml("<tr><td><img src=\"raw/A2_ck.jpg\"></img></td>",f)
  whtml("    <td><img src=\"raw/A2_sig.jpg\"></img></td>",f)
  whtml("    <td><img src=\"raw/A2_logsig.jpg\"></img></td>",f)
  whtml(paste("    <td><img src=\"",data@normalization[length(data@normalization)],"/A2_sig.jpg\"></img></td>",sep=""),f)
  whtml("    <td><img src=\"summarized/A2_sig.jpg\"></img></td></tr>",f)
  whtml("<tr><td></td>",f)
  whtml("    <td><img src=\"raw/A4_scatter.jpg\"></img></td>",f)
  whtml("    <td></img></td>",f)
  whtml(paste("    <td><img src=\"",data@normalization[length(data@normalization)],"/A4_scatter.jpg\"></img></td>",sep=""),f)
  whtml("    <td><img src=\"summarized/A4_scatter.jpg\"></img></td></tr>",f)
  whtml("</table>",f)
  
  ## Volcano PlotSpatialDistrib
  makevolcanoplot(data,fn="volcano.jpg")
  whtml("<hr>",f)
  whtml("<p><img src=\"volcano.jpg\"></img></p>",f)
    
  ## Barplot full screen
  makefullbarplot(data,fn="fullbar.jpg")
  whtml("<p><img src=\"fullbar.jpg\"></img></p>",f)
  whtml("<hr>",f)
  
  ## Links to detailed steps of each normalization
  whtml("<h2>Detailed Results</h2>",f)
  whtml("<p><ol>",f)
  whtml("<li><a href=\"raw/index.html\">Raw Data</a></li>",f)
  if (data@excludeCellcounts!="none")
    whtml("<li><a href=\"ckremoved/index.html\">Wells removed based on Cell Numbers</a></li>",f)
  if (data@logtransform==TRUE)
    whtml("<li><a href=\"log/index.html\">Log Transformation</a></li>",f)
  for (norm in data@normalization)
    whtml(paste("<li><a href=\"",norm,"/index.html\">",norm," normalized</a></li>",sep=""),f)
  whtml("<li><a href=\"summarized/index.html\">Summarized</a></li>",f)
  whtml("</ol>",f)
  whtml(" <BR/><HR>",f)
  whtml("</BODY>",f)
  whtml("</HTML>",f)
  return();
}


#################################################################
# Additional normalization functions

# normalization on negative controls!!

doCtrlNorm <- function (header, dataset) 
{
    dataset$SigIntensity.BeforeCtrlNorm <- dataset$SigIntensity
    numOfScreens <- max(dataset$ScreenNb)
    minOfScreens <- min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens) {
        if (length(which(dataset$ScreenNb == i)) > 0) {
            subset <- createSubset(dataset, dataset$ScreenNb, i)
            ixx <- indexSubset(dataset$ScreenNb, i)
            minOfPlates = min(subset$LabtekNb)
            numOfPlates = max(subset$LabtekNb)
            for (j in minOfPlates:numOfPlates) {
                if (length(which(subset$LabtekNb == j)) > 0) {
                    subsubset <- createSubset(subset, subset$LabtekNb, j)
                    iyy <- indexSubset(subset$LabtekNb, j)
                    normSubset <- createSubset(subsubset, subsubset$SpotType, 0)
                    c1 <- subsubset$SigIntensity
                    med <- median(normSubset$SigIntensity, na.rm = T)
                    themad <- mad(normSubset$SigIntensity, na.rm=T)
                    message(paste("Normalizing Plate ",j,", Replicate ",i,"; median=",med,", mad=",themad,sep=""))
                    subsubset$SigIntensity <- (c1 - med) / themad
                    subset[iyy, ] <- subsubset
                }
            }
            dataset[ixx, ] <- subset
        }
    }
    header[3] <- gsub("NA", "", header[3])
    header[3] <- paste(header[3], "Normalization of SigIntensity on Negative Controls (std. of negative)")
    invisible(list(header, dataset))
}

doPocNorm <- function (header, dataset) 
{
    dataset$SigIntensity.BeforeCtrlNorm <- dataset$SigIntensity
    numOfScreens <- max(dataset$ScreenNb)
    minOfScreens <- min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens) {
        if (length(which(dataset$ScreenNb == i)) > 0) {
            subset <- createSubset(dataset, dataset$ScreenNb, i)
            ixx <- indexSubset(dataset$ScreenNb, i)
            minOfPlates = min(subset$LabtekNb)
            numOfPlates = max(subset$LabtekNb)
            for (j in minOfPlates:numOfPlates) {
                if (length(which(subset$LabtekNb == j)) > 0) {
                    subsubset <- createSubset(subset, subset$LabtekNb, j)
                    iyy <- indexSubset(subset$LabtekNb, j)
                    normSubsetNeg <- createSubset(subsubset, subsubset$SpotType, 0)
		    normSubsetPos <- createSubset(subsubset, subsubset$SpotType, 1)
                    c1 <- subsubset$SigIntensity
                    medpos <- median(normSubsetPos$SigIntensity, na.rm = T)
		    medneg <- median(normSubsetNeg$SigIntensity, na.rm = T)
                    message(paste("Normalizing Plate ",j,", Replicate ",i,"; medpos=",medpos,", medneg=",medneg,sep=""))
                    subsubset$SigIntensity <- 100 * (c1 - medpos) / (medneg - medpos)
                    subset[iyy, ] <- subsubset
                }
            }
            dataset[ixx, ] <- subset
        }
    }
    header[3] <- gsub("NA", "", header[3])
    header[3] <- paste(header[3], "Normalization of SigIntensity on Controls (Percentage of Controls)")
    invisible(list(header, dataset))
}


doPercNegNorm <- function (header, dataset) 
{
    dataset$SigIntensity.BeforeCtrlNorm <- dataset$SigIntensity
    numOfScreens <- max(dataset$ScreenNb)
    minOfScreens <- min(dataset$ScreenNb)
    for (i in minOfScreens:numOfScreens) {
        if (length(which(dataset$ScreenNb == i)) > 0) {
            subset <- createSubset(dataset, dataset$ScreenNb, i)
            ixx <- indexSubset(dataset$ScreenNb, i)
            minOfPlates = min(subset$LabtekNb)
            numOfPlates = max(subset$LabtekNb)
            for (j in minOfPlates:numOfPlates) {
                if (length(which(subset$LabtekNb == j)) > 0) {
                    subsubset <- createSubset(subset, subset$LabtekNb, j)
                    iyy <- indexSubset(subset$LabtekNb, j)
                    normSubsetNeg <- createSubset(subsubset, subsubset$SpotType, 0)
                    c1 <- subsubset$SigIntensity
		    medneg <- median(normSubsetNeg$SigIntensity, na.rm = T)
                    message(paste("Normalizing Plate ",j,", Replicate ",i,"; medneg=",medneg,sep=""))
                    subsubset$SigIntensity <- c1/medneg*100
                    subset[iyy, ] <- subsubset
                }
            }
            dataset[ixx, ] <- subset
        }
    }
    header[3] <- gsub("NA", "", header[3])
    header[3] <- paste(header[3], "Normalization of SigIntensity on Controls (Percentage of Negative Controls)")
    invisible(list(header, dataset))
}


#######################################
# reorder data so that plates are in consecutive frame, and all spots are filled - might otherwise cause problems in RNAither

reorderdata <- function(olddata)
{
  message("\nReordering data set...",appendLF = F)
  maxspot <- max(olddata$Spotnumber)
  maxrow <- max(olddata$RowNb)
  maxcol <- max(olddata$ColNb)
  maxreplicate <- max(olddata$ScreenNb)
  maxplate <- max(olddata$LabtekNb)
  message(paste("\n",maxspot," Wells; ",maxrow," Rows; ",maxcol," Columns; ",maxplate," Layouts; ",maxreplicate," Replicates.",sep=""),appendLF = F)
  # first plausibility test: is row*col = spot?
  if (maxrow * maxcol != maxspot)
    stop(paste("Number of rows * columns (",maxrow," * ",maxcol," = ",maxrow*maxcol," is not equal to number of spots = ",maxspot,"!\n",
               "Please check your input data!",sep=""))
  # now create a new dataframe containing a value for everything!
  framecols <- dim(olddata)[2]
  newdata <- olddata
  newdata$SpotType <- as.numeric(as.character(olddata$SpotType))
  i <- 1
  for (rep in 1:maxreplicate)
    for (plt in 1:maxplate){
      message(paste("\nPlate ",plt,", Replicate ",rep,sep=""),appendLF = F)
      for (wll in 1:maxspot)
      {
        mycol <- which((olddata$ScreenNb==rep)&(olddata$LabtekNb==plt)&(olddata$Spotnumber==wll))
        if (length(mycol)==1){
          message(".",appendLF = F)
          newdata[i,] <- olddata[mycol,]
        }
        else if (length(mycol)==0){
          message("M",appendLF = F)
          newdata[i,] <- rep(NA,framecols)
          newdata$Spotnumber[i] <- wll
          newdata$RowNb[i] <- ((wll-1) %/% maxcol) + 1
          newdata$SpotType[i] <- -1
          newdata$ColNb[i] <- ((wll-1) %% maxcol) + 1
          newdata$LabtekNb[i] <- plt
          newdata$ScreenNb[i] <- rep
        }
        else
          stop(paste("Duplicate entry in input file: Plate ",plt,", Replicate ",rep,", Well ",wll,sep=""))
        i <- i + 1
      }
    }
  return(newdata);
}

attempt_fix <- function(itsdata)
{
  newdata <- itsdata
  NoOfWells <- max(itsdata$Spotnumber)
  NoOfRows <- max(itsdata$RowNb)
  NoOfCols <- max(itsdata$ColNb)
  NoOfReps <- max(itsdata$ScreenNb)
  NoOfLayouts <- max(itsdata$LabtekNb)
  
  # for each layout separately, take the format of the gene ids as factor levels, sort and use this
  # sorting further for all replicates! Just ensure that Spotnumber, row and column are subsequently
  # corrected -> in consequence, this means the layout of the screen will be distorted, hence any
  # spatial normalization is infeasible subsequently.
  
  genenames <- as.numeric(as.factor(itsdata$Internal_GeneID))
  for (i in 1:NoOfLayouts){
    for (j in 1:NoOfReps){
      whois <- (itsdata$LabtekNb==i)&(itsdata$ScreenNb==j)
      gsorted <- sort(genenames[whois],index.return=T)
      thedata <- itsdata[whois,]
      newdata[whois,] <- thedata[gsorted$ix,]
    }
  }
  # and finally fix the Layout
  newdata$Spotnumber <- itsdata$Spotnumber
  newdata$RowNb <- itsdata$RowNb
  newdata$ColNb <- itsdata$ColNb
  return(newdata)
}

  
  
check_layout_consistent <- function(itsdata)
{
  newdata <- itsdata
  NoOfWells <- max(itsdata$Spotnumber)
  NoOfRows <- max(itsdata$RowNb)
  NoOfCols <- max(itsdata$ColNb)
  NoOfReps <- max(itsdata$ScreenNb)
  NoOfLayouts <- max(itsdata$LabtekNb)
  
  # let us first check, if the layouts are identical between replicates
  diffLayout <- FALSE
  for (i in 1:NoOfLayouts){
    for (j in 1:NoOfWells){
      whosrc <- which((itsdata$LabtekNb==i)&(itsdata$Spotnumber==j))
      if (length(unique(itsdata$Internal_GeneID[whosrc]))>1)
        diffLayout <- TRUE
    }
  }
  if (diffLayout){
    warning("*** Layouts differ between replicates. This WILL result in errors. ****")
    message("\n\n***********************************************************************")
    message(  "\n*** Layouts differ between replicates. This WILL result in errors. ****")
    message(  "\n***********************************************************************\n\n")
    message("Attempting fix...")
    # check further if at least the same genes are present in the different replicates...
    # if so, we can attempt a fix...
    fixable <- TRUE
    for (i in 1:NoOfLayouts){
      # for each layout test, if the layout contains the same genes, just in different order,
      # in the different replicates
      genes_l1 <- sort(itsdata$Internal_GeneID[(itsdata$LabtekNb==i)&(itsdata$ScreenNb==1)])
      len1 <- length(genes_l1)
      for (j in 2:NoOfReps){
        genes_lr <- sort(itsdata$Internal_GeneID[((itsdata$LabtekNb==i)&(itsdata$ScreenNb==j))])
        lenr <- length(genes_lr)
        if (len1 != lenr)
          fixable <- FALSE
        else if (sum(genes_l1 == genes_lr) != len1)
          fixable <- FALSE      
      }
    }
    if (!fixable)
      stop("\n** Error: Replica plates contain different genes. Unable to fix this. Aborting **")
    newdata <- attempt_fix(itsdata)
  } else message("\nLayout is identical between replicates.")
    
  return(newdata);
}
  
#################################################################
# rnaither -- main function to carry out analysis of screen.

rnaither <- function(data,                           # rnaither data frame containing data
                     expname,                        # experiment name
                     excludeCellcounts="none",       # one of "none", "lowest","both","lowestperplate","bothperplate"
                     logtransform=FALSE,             # log transformation of data
                     normalization=                  # normalization steps to carry out
                           c("lowess","bscore"),
                     test="ttest",                   # what statistical test to perform
                     scorethresh=2.0  ,              # threshold on score
                     pvalthresh=0.05,                # threshold on p-value
                     dogo=FALSE,                     # carry out gene ontology analysis
                     outdir="results",               # directory for results
                     layoutnames="NA",               # list of strings for layout names
		     makeplots=TRUE,                 # should plots not be made (for debugging)
                     reorder=TRUE                    # reorder dataset before processing?
)
{
  # excludeCellcounts:  none             - no wells will be excluded based on cell counts
  #                     lowest           - wells with lowest 5% of cell counts in entire screen
  #                     both             - wells with lowest and highest 5% of cell counts in entire screen
  #                     lowestperplate   - 5% lowest cell counts per plate are excluded
  #                     bothperplate     - 5% lowest and highest cell counts per plate are excluded
  # normalization:      lowess           - carry out lowess normalization
  #                     liwong           - carry out Li-Wong rank normalization
  #                     varadjust        - adjust variance
  #                     divnorm          - divide by median
  #                     quantile         - quantile normalization
  #                     bscore           - bscore normalization
  #                     zscore           - zscore normalization
  #                     negcontrol       - normalization on negative controls
  #                     percontrol       - normalization on positive controls and negative controls
  # test:               ttest            - carry out a t-test (mu=0)
  #                     wilcox           - carry out a wilcoxon test (mu=0)
  #                     none             - no statistical test

  message("Initializing...",appendLF = F)

  mydata <- new("rnaistat_data")
  if (reorder)
    mydata@thedata <- reorderdata(data)
  else
    mydata@thedata <- data

  # Let us check first if the layout is the same for all plates in the different replicates!
  mydata@thedata <- check_layout_consistent(mydata@thedata)
  
  mydata@expname <- expname
  mydata@excludeCellcounts <- excludeCellcounts
  mydata@logtransform <- logtransform
  mydata@normalization <- normalization
  mydata@test <- test
  mydata@scorethresh <- scorethresh
  mydata@pvalthresh <- pvalthresh
  mydata@dogo <- dogo
  mydata@status$inited=TRUE
  header <- c(as.character(paste("external_experiment_name,",expname,sep="")),
              as.character(paste("type_of_data,",max(data$Spotnumber)," well plate data",sep="")),
              as.character("comments,Analyzed with RNAither"))
  mydata@header <- header
  mydata@layoutnames <- layoutnames
  mydata@rawdata <- data

  if (layoutnames[1]=="NA")
    message("\nNo layout names provided...",appendLF = F)
  else if (length(layoutnames)!=max(data$LabtekNb))
    message("\nLength of LayoutNames and Layout Number differ!",appendLF = F)
  else
    message("\nUsing LayoutNames...",appendLF = F)

  NoOfLayouts = max(mydata@thedata$LabtekNb)
  NoOfReplicates = max(mydata@thedata$ScreenNb)
  NoOfWells = max(mydata@thedata$Spotnumber)

  # ok, we now have the data set up. first thing to do is to create the quality plots!
  message("\nQuality Control Plots...",appendLF = F)
  dir.create(outdir)
  currdir <- getwd()
  setwd(outdir)
  makeQCPlots(mydata,"raw","Raw Data",iszscore=FALSE)

  # next remove cellcounts, and carry out log transformation
  # excludeCellcounts:  none             - no wells will be excluded based on cell counts
  #                     lowest           - wells with lowest 5% of cell counts in entire screen
  #                     both             - wells with lowest and highest 5% of cell counts in entire screen
  #                     lowestperplate   - 5% lowest cell counts per plate are excluded
  #                     bothperplate     - 5% lowest and highest cell counts per plate are excluded
  
  # MODIFIED LK 2010 Nov 03: Do NOT remove Wells annotated as controls!
  myctrl0 <- which(mydata@thedata$SpotType==0)
  myctrl1 <- which(mydata@thedata$SpotType==1)
  if (excludeCellcounts=="lowest")
  {
    message("\nExclusion of 5% of wells with lowest cell count (total)...",appendLF = F)
    ctoff <- quantile(mydata@thedata$NbCells,probs=0.05,na.rm=TRUE)
    mydata@thedata$SpotType[mydata@thedata$NbCells < ctoff] <- -1
    mydata@thedata$SpotType[myctrl0] <- 0
    mydata@thedata$SpotType[myctrl1] <- 1
    makeQCPlots(mydata,"ckremoved","Removed Lowest 5% of Cellcounts",iszscore=FALSE)
  }
  else if (excludeCellcounts=="both")
  {
    message("\nExclusion of 10% of wells with lowest and highest cell count (total)...",appendLF = F)
    ctoff <- quantile(mydata@thedata$NbCells,probs=0.05,na.rm=TRUE)
    ctoff2 <- quantile(mydata@thedata$NbCells,probs=0.95,na.rm=TRUE)
    mydata@thedata$SpotType[mydata@thedata$NbCells < ctoff] <- -1
    mydata@thedata$SpotType[mydata@thedata$NbCells > ctoff2] <- -1
    mydata@thedata$SpotType[myctrl0] <- 0
    mydata@thedata$SpotType[myctrl1] <- 1
    makeQCPlots(mydata,"ckremoved","Removed Extreme 10% of Cellcounts",iszscore=FALSE)
  }
  else if (excludeCellcounts=="lowestperplate")
  {
    message("\nExclusion of 5% of wells with lowest cell count (per plate)...",appendLF = F)
    for (i in 1:NoOfLayouts){
      for (j in 1:NoOfReplicates){
        subset <- which((mydata@thedata$LabtekNb == i)&(mydata@thedata$ScreenNb == j))
        ctoff <- quantile(mydata@thedata$NbCells[subset],probs=0.05,na.rm=TRUE)
        issmaller <- which(mydata@thedata$NbCells[subset]<ctoff)
        mydata@thedata$SpotType[subset][issmaller] <- -1
      }
    }
    mydata@thedata$SpotType[myctrl0] <- 0
    mydata@thedata$SpotType[myctrl1] <- 1
    makeQCPlots(mydata,"ckremoved","Removed Lowest 5% of Cellcounts per Plate",iszscore=FALSE)
  }
  else if (excludeCellcounts=="bothperplate")
  {
    message("\nExclusion of 10% of wells with lowest and highest cell count (per plate)...",appendLF = F)
    for (i in 1:NoOfLayouts){
      for (j in 1:NoOfReplicates){
        subset <- which((mydata@thedata$LabtekNb == i)&(mydata@thedata$ScreenNb == j))
        ctoff <- quantile(mydata@thedata$NbCells[subset],probs=0.05,na.rm=TRUE)
        ctoff2 <- quantile(mydata@thedata$NbCells[subset],probs=0.95,na.rm=TRUE)
        issmaller <- which(mydata@thedata$NbCells[subset]<ctoff)
        islarger <- which(mydata@thedata$NbCells[subset]>ctoff2)
        mydata@thedata$SpotType[subset][issmaller] <- -1
        mydata@thedata$SpotType[subset][islarger] <- -1
      }
    }
    mydata@thedata$SpotType[myctrl0] <- 0
    mydata@thedata$SpotType[myctrl1] <- 1
    makeQCPlots(mydata,"ckremoved","Removed Extreme 10% of Cellcounts per Plate",iszscore=FALSE)
  }
  else if (excludeCellcounts=="none")
    message("\nNo exclusion of wells based on cell count...",appendLF = F)
  else
    stop(paste("Unrecognized option: excludeCellcounts=",excludeCellcounts,sep=""))

  # log-transform data
  if (logtransform){
    message("\nLog transformation...",appendLF = F)
    mydata@thedata$SigIntensity <- log(mydata@thedata$SigIntensity)
    makeQCPlots(mydata,"log","Log Transformed",iszscore=FALSE)
  }

  # now perform normalization
  # normalization:      lowess           - carry out lowess normalization
  #                     liwong           - carry out Li-Wong rank normalization
  #                     varadjust        - adjust variance
  #                     divnorm          - divide by median
  #                     quantile         - quantile normalization
  #                     bscore           - bscore normalization
  #                     zscore           - zscore normalization
  #                     negcontrol       - normalization on negative controls
  #                     percontrol       - normalization on positive controls and negative controls
  #                     percneg          - set median of negative controls to 100%
  # loop over the different normalization procedures!
  iszscore <- FALSE
  for (normtodo in normalization){
    if (normtodo=="lowess"){
      # carry out lowess normalization
      message("\nLowess Normalization...",appendLF = F)
      normres <- lowessNorm(mydata@header,mydata@thedata,list("NbCells","SigIntensity"))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="liwong"){
      message("\nLiWong Rank Normalization...",appendLF = F)
      iszscore <- TRUE
      # carry out liwong rank normalization
      normres <- LiWongRank(mydata@header,mydata@thedata,list("SigIntensity","GeneName"))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="varadjust"){
      # carry out variance adjustment normalization
      message("\nVariance Adjustment...",appendLF = F)
      normres <- varAdjust(mydata@header,mydata@thedata,list(2,"SigIntensity",1))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="divnorm"){
      # carry out division normalization
      message("\nDivision Normalization...",appendLF = F)
      normres <- divNorm(mydata@header,mydata@thedata,list(median,2,1,"SigIntensity",1))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="quantile"){
      # carry out quantile normalization
      message("\nQuantile Normalization...",appendLF = F)
      iszscore <- TRUE
      normres <- quantileNormalization(mydata@header,mydata@thedata,list(2,"SigIntensity"))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="bscore"){
      # carry out bscore normalization
      message("\nBScore Normalization...",appendLF = F)
      iszscore <- TRUE
      normres <- BScore(mydata@header,mydata@thedata,list("SigIntensity",1))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="zscore"){
      # carry out zscore normalization
      message("\nZScore Normalization...",appendLF = F)
      iszscore <- TRUE
      normres <- ZScore(mydata@header,mydata@thedata,list("SigIntensity",1))
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="negcontrol"){
      # carry out negcontrol normalization
      message("\nNormalization on Negative Controls...",appendLF = F)
      normres <- doCtrlNorm(mydata@header,mydata@thedata)
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else if (normtodo=="percontrol"){
      # carry out percentage of controls normalization
      message("\nPercentage of Controls Normalization...",appendLF = F)
      normres <- doPocNorm(mydata@header,mydata@thedata)
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
#      mydata@thedata$SigIntensity = 100*mydata@thedata$SigIntensity
    }
    else if (normtodo=="percneg"){
      # carry out negcontrol normalization
      message("\nPercentage of Negative Controls Normalization...",appendLF = F)
      normres <- doPercNegNorm(mydata@header,mydata@thedata)
      mydata@header <- normres[[1]]
      mydata@thedata <- normres[[2]]
    }
    else{
      stop(paste("Unknown normalization method ",normtodo,sep=""))
    }
    # create plots
    makeQCPlots(mydata,normtodo,paste(normtodo," Normalized",sep=""),iszscore)
  }

  # summarize replicates, carry out statistical test
  message("\nSummarizing replicates, performing statistical test...",appendLF = F) 
  mydata@normalizeddata <- mydata@thedata
  mydata <- summarizeReplicates(mydata,NoOfLayouts,NoOfWells,test)
  makeQCPlots(mydata,"summarized","Summarized over Replicates",iszscore)
   
  # carry out go analysis
  if (dogo)
    warning("**Go analysis is not supported by the rnaither wrapper at present**")

  message("\nWriting summary html file...",appendLF = F)
  writeSummaryHtml(mydata)
  message("\nWriting results.Rda file...",appendLF = F)
  save(mydata,file="results.Rda")
  setwd(currdir)
  message(paste("\nAnalysis is complete, results are in directory ",outdir,"\n",sep=""))
}

