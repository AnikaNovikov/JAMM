########################################################################
# JAMMv1.0.7rev2 is a peak finder for joint analysis of NGS replicates.
# Copyright (C) 2014-2015  Mahmoud Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: mahmoud.ibrahim@mdc-berlin.de
########################################################################





# ======================= 
# User-defined variables
# =======================
presched = TRUE #Set preschedule in mclapply to TRUE or FALSE
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000).
#=======================> DONE! 



# ========================= 
# Required Libraries check 
# ========================= 
if ((is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package'parallel' is required. Please install it!")
}
if ((is.element('Rsamtools', installed.packages()[,1]) == FALSE)) {
  stop("R package'Rsamtools' is required. Please install it!")
}
if ((is.element('GenomicRanges', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicRanges' is required. Please install it!")
}
if ((is.element('GenomicAlignments', installed.packages()[,1]) == FALSE)) {
  stop("R package'GenomicAlignments' is required. Please install it!")
}
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("GenomicAlignments"))
#=======================> DONE! 




# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
ibam = NA #input bam file
sFile = NA # chromosome size
storeFile = NA # file to store result
ReadChromVector = NA

for (each.arg in args) {
	if (grepl('^-s=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			sFile <- arg.split[2] 
		} else {
			stop('No chromosome size file')
		}
		chromosomes = read.table(sFile, header=FALSE)
		chromName = as.character(chromosomes$V1) #which chromosome
		chromSize = as.numeric(chromosomes$V2) #chromosome size
		rm(chromosomes)
	}

	#input bam file
	if (grepl('^-ibam=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			ibam <- arg.split[2]
		} else {
			message(paste0(chromName, ": No reads found in one or more replicates!"))
			quit()
		}
	}
	if (grepl('^-iindex=',each.arg)) {
	  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	  if (! is.na(arg.split[2]) ) {
	    iindex <- arg.split[2]
	  } else {
	    #iindex <- index(BamFile(ibam))
	    stop('No index file for one or multiple bam files given')
	    #quit()
	  }
	}
	#background files directory
	if (grepl('-bkgd=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bkgd <- arg.split[2]
		} else {
			stop('No background file')
		} 
	}	
#	if (grepl('^-rl=',each.arg)) {			
#		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
#		if (! is.na(arg.split[2]) ) {
#			#rl <- arg.split[2] 
#      #print (paste("rl",rl))
#		} else {
#			stop('Read length missing')
#		}
#	}
	if (grepl('^-d=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			storeFile <- arg.split[2] 
		} else {
			stop('No file to store result')
		}
	}
	if (grepl('^-p=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			cornum <- as.numeric(arg.split[2])
      #print(paste("cornum",cornum))
		} else {
			stop('No number of cores given')
		}
	}
	#number of replicates
	if (grepl('-nreps=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			nreps <- as.numeric(arg.split[2])
		} 
	}
  #####
	#if (grepl('^-readwhichchroms=',each.arg)) {  		
	#  arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
	#  if (! is.na(arg.split[2]) ) {
	#    ReadChromVector <- arg.split[2] 
	#  }
	#}
}

#Read in variables
if(is.na(ReadChromVector))
{
  #ReadChromVector <- rep(1, times=length(chromName))
  ReadChromVector <- chromName
}
#ReadChromVector = strsplit(ReadChromVector, ",", fixed = TRUE)[[1]]
ibam = strsplit(ibam, ",", fixed = TRUE)[[1]]
iindex = strsplit(iindex, ",", fixed = TRUE)[[1]]
#print(paste("iindex",iindex))

numdup = length(ibam) #number of replicates
numbkgd = numdup - nreps
if (bkgd != "None") {
	nreps = nreps + 1
}

#if (numdup != nreps) {
#	message(paste0(chromName, ": No reads found in one or more replicates!"))
#	quit()
#}

#####rl = as.numeric(strsplit(rl, ",", fixed = TRUE)[[1]])
#=======================> DONE! 

matchChromNames <- function(actlevels,userchrom,userchromlength)
{
  userchrom <- tolower(userchrom)
  userchrom <- gsub("chr","",userchrom)
  userchrom <- gsub("seq","",userchrom)
  tempactlevels <- tolower(actlevels)
  tempactlevels <- gsub("chr","",tempactlevels)
  tempactlevels <- gsub("seq","",tempactlevels)
  temp <- match(tempactlevels,userchrom)
  if(!any(!is.na(temp)))
  {
    stop("Chromosome names in the bam file do not correspond to information in the size file")
  }
}

makeCountList <- function(actlevels, userlevels, userlevellength, binchromv)
{
  countlist <- list()
  userlevels <- tolower(userlevels)
  userlevels <- gsub("chr","",userlevels)
  userlevels <- gsub("seq","",userlevels)
  if(identical(class(binchromv),"character")) 
  {
    #print("yes")
    chromv <- tolower(binchromv)
    chromv <- gsub("chr","",chromv)
    chromv <- gsub("seq","",chromv)
    st <- match(userlevels,chromv)
    for(i in 1:length(st))
    {
      if(is.na(st[[i]])) {binchromv[[i]] <- 0} 
      else {binchromv[[i]] <- 1}
    }
  }
  #binchromv <- rep(1, times=length(actlevels))
  tempactlevels <- tolower(actlevels)
  tempactlevels <- gsub("chr","",tempactlevels)
  tempactlevels <- gsub("seq","",tempactlevels)
  temp <- match(userlevels,tempactlevels)
  temp2 <- match(tempactlevels,userlevels)
  #false if 1)what's in the file is provided and what's provided is in the file (or more)
  if (any(is.na(temp))||any(is.na(temp2)))
  {
    actlevels <- actlevels[which(!is.na(temp2))]
    #what's in the file is provided, but there's more provided, so remove the unnecessary information
    userlevels <- userlevels[which(!is.na(temp))]
    binchromv <- binchromv[which(!is.na(temp))]
    userlevellength <- userlevellength[which(!is.na(temp))]
    temp <- Filter(Negate(function(x) is.na(unlist(x))), temp)
    if(any(is.na(temp))){
      warning("There are less chromosomes in the bam file than the sizes file indicates")
    }
  }
  newbinchromv <- binchromv[order(temp)]
  newuserlevellength <- userlevellength[order(temp)]
  
  actlevels <- actlevels[which(newbinchromv==1)]
  newuserlevellength <- newuserlevellength[which(newbinchromv==1)]
  #print(actlevels)
  #length(countlist) <<- length(actlevels)
  #names(countlist) <<- actlevels
  length(countlist) <- length(actlevels)
  names(countlist) <- actlevels
  for (i in 1:length(actlevels))
  {
    countlist[[actlevels[[i]]]] <- newuserlevellength[[i]]
  }
  #mapply(defineLength,actlevels,newuserlevellength)
  return(countlist)
  #return()
}

#defineLength <- function(chromname,chromlength)
#{
#  countlist[[chromname]] <<- chromlength
#  return()
#}

counthelp <- function(readstarts, breaks, cornum, presched){
  if (max(readstarts) > max(breaks)) {
    return(rep(0, length = length(breaks)))
  } else {
    return(hist(readstarts, breaks = breaks, plot = FALSE)$counts)
  }
}

#countingreads = function(bedfile, readlen, filelist, gv) {
#  reads = read.table(bedfile, header = FALSE)
#  reads = cbind(as.character(reads[[1]]), as.character(reads[[2]]))
#  o = which(filelist == bedfile)
#  tri = list(as.numeric(reads[reads[,2] == "+",,drop = FALSE][,1]), (as.numeric(reads[reads[,2] == "-",,drop = FALSE][,1])) + readlen[o] - 1)
#  readcounts = mclapply(tri, counting, br = gv, mc.cores = cornum, mc.preschedule = presched)
#  return(readcounts)
#}


normalize <- function (vector)
{
  return(vector/mean(vector))
}

countup <- function (bamfile, indexfile, filelist, cornum, presched, rmduplicates=FALSE, RCV=ReadChromVector, calcreadlen=FALSE) 
{
  #print("rcv")
  #print(RCV)
  #print(paste("mc.cores",cornum))
  alsm <- NA; alsp <- NA
  cat("\n", "xcorr: Counting",bamfile,"\n")
  o = which(filelist == bamfile)
  bamname <- bamfile
  bamfile <- BamFile(bamfile)
  param <- ScanBamParam()
  #countlist <<- list()
  countlist <- list()
  #creates a list for all chromosomes to calculate xcorr for
  #print(seqnames(seqinfo(bamfile)));print(chromName);print(chromSize);print(paste("READCHROMVECTOR",RCV))
  #print("rcv")
  #makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,RCV)
  countlist <- makeCountList(seqnames(seqinfo(bamfile)),chromName,chromSize,RCV)
  #matchChromNames(seqnames(seqinfo(bamfile)),as.character(chromName),as.integer(chromSize))
  #countlist[[as.character(chromName)]] <<- as.integer(chromSize)
  #names(countlist) <<- as.character(chromName)
  chromnames <- names(countlist)
  #print(paste(chromnames,class(chromnames)))
  #open(bamfile)
  curnum <- 1
  store <- 0
  for (element in chromnames)
  {
    if(curnum == 1) {
      cat("\n","Counting:",length(chromnames),"elements","\n")
      cat(paste0(element,", ",curnum,"; "))
    } else {
      cat(paste0(element,", ",curnum,"; "))
    }
    
    chromlength <- as.integer(countlist[[element]])
    #print(paste("chromlength",chromlength))
    #print(bamfile);print(indexfile);print(element);print(class(element));print(as.integer(chromlength))
    alsm <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE),which=GRanges(element,IRanges(1,as.integer(chromlength)))))
    alsp <- readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),which=GRanges(element,IRanges(1,as.integer(chromlength)))))
    #print(paste("length alsm",length(alsm)))
    #print(paste("length alsp",length(alsp)))
    if ((!length(alsm)) && (!length(alsp)))
    { 
      countlist[[element]] <- NULL
      curnum=curnum+1
      next()
    }
    #readstarts <- list()
    readstarts <- c()
    if(!length(alsm)){
      alsm <- GRanges()
    } else {
      alsm <- GRanges(alsm)
      alsm <- unique(alsm)
      #alsm <- resize(alsm,1)
    }
    if(!length(alsp)){
      alsp <- GRanges()
    } else {
      alsp <- GRanges(alsp)
      alsp <- unique(alsp)
      #alsp <- resize(alsp,1)
    }
    
    if(!length(alsm)){
      #readstarts <- list(start(ranges(alsp)))
      if (calcreadlen){
        readstarts <- start(ranges(alsp))
        readends <- end(ranges(alsp))
        len <- length(readstarts)
        result <- readends-readstarts
        store <- store+c(len,sum(result))
      }
      alsp <- resize(alsp,1)
      readstarts <- list(start(ranges(alsp)))
      curvnames <- c("+")
    } else if (!length(alsp)) {
      #readstarts <- list(start(ranges(alsm)))
      if(calcreadlen){
        readstarts <- start(ranges(alsm))
        readends <- end(ranges(alsm))
        len <- length(readstarts)
        result <- readends-readstarts
        store <- store+c(len,sum(result))
      }
      alsm <- resize(alsm,1)
      readstarts <- list(start(ranges(alsm)))
      curvnames <- c("-")
    } else {
      #readstarts <- list(start(ranges(alsp)),start(ranges(alsm)))
      if(calcreadlen){
        readstarts <- c(start(ranges(alsp)),start(ranges(alsm)))
        readends <- c(end(ranges(alsp)),end(ranges(alsm)))
        len <- length(readstarts)
        result <- readends-readstarts
        store <- store+c(len,sum(result))
      }
      alsm <- resize(alsm,1)
      alsp <- resize(alsp,1)
      readstarts <- list(start(ranges(alsp)),start(ranges(alsm)))
      #print(head(readstarts));print(class(readstarts))
      curvnames <- c("+","-")
    }
    #rm(readends)
    #alsp <- GRanges(readGAlignments(bamfile,index=indexfile,param=ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),which=GRanges(element,IRanges(1,as.integer(chromlength))))))
    #if (rmduplicates) {
    #  alsm <- unique(alsm)
    #  alsp <- unique(alsp)
    #}
    #alsm <- unique(alsm)
    #alsp <- unique(alsp)
    
    
    #make the window vector (a vector with start positions of chromosome windows)
    breaks = seq(1, chromlength, by = 20)
    if (max(breaks) < chromlength) {
      breaks = c(breaks, chromlength);
    }
    curvector <- list()
    curvector <- mclapply(readstarts, counthelp, breaks=breaks, cornum = cornum, presched = presched)
    #normalize the read counts
    #curvector <- lapply(curvector,normalize)
    
    #names(curvector) <- c("+","-")
    names(curvector) <- curvnames
    #####in case alsm or alsp are empty... adapt condition like above
    if (all(curvector[[1]] == 0) && all(curvector[[2]] == 0)) {
      message(paste0(element, ", Warning: Read alignments do not match chromosome length, Skipped!"))
      quit()
    }
    
    ################COMMENT FROM MAHMOUD###########
    ##the function "xc" can be run here directly without ever constructing this huge "countlist" variable.
    ##I think this is what's causing the memory issues and not "curvector"
    ##In this case the function "countup" would actually return the results from the function "xc" directly
    ###############################################
    #countlist[[element]] <<- curvector
    countlist[[element]] <- curvector
    #print(head(countlist[[element]][["+"]]))
    #print(paste("head curvector",head(curvector[[1]])))
    alsm <- NA
    alsp <- NA
    curnum=curnum+1
    #gc()
    #bamfile <- NA
  }
  
    ################COMMENT FROM MAHMOUD###########
    ##what does this if statement do?
    ###############################################
  if(calcreadlen){
    rl <- calculateAverageReadLength(store)
    #print(rl)
    filename = strsplit(bamname, "/", fixed = TRUE)[[1]]
    filename = filename[length(filename)]
    filename = substr(filename,1,nchar(filename)-4)
    write(paste(bamname, rl, sep= ","), file = paste0(storeFile, "/xc.rl.", filename, ".tab"), append = TRUE)
    #write(paste(bamname, rl, sep= ","), file = paste0(storeFile, "xc.rl.", bamname, ".tab"), append = TRUE);
  }
  rm(alsm,alsp,curvector,readstarts)
  gc()
  return (countlist)
}



calculateAverageReadLength <- function(store)
{
  return(trunc((store[[2]]/store[[1]])+0.5))
}

xc = function(countlist, cornum, presched) {
  if (sum(unlist(countlist)) != 0) {
    cat("next, ","\n")
		crossCorrelation = ccf(countlist[[2]], countlist[[1]], plot = FALSE); #xcorr
		crossCorrelation$lag = crossCorrelation$lag * 20; #correct lag for counts window
		maxCorr = which.max(crossCorrelation$acf);
		maxCorr = abs(crossCorrelation$lag[maxCorr]);
		return(maxCorr)
	} else {
		return(0)
	}
}

#####what if default fragment length, this should work as well

# ================== 
# Read Start Count
# ==================
#if(is.na(ReadChromVector))
#{
#  ReadChromVector <- rep(1, times=length(chromName))
#}


areinbkgd <- c()
if(numbkgd!=0)
{
  ibkgd <- ibam[(length(ibam)-numbkgd+1):length(ibam)]
  ibkgdindex <- iindex[(length(iindex)-numbkgd+1):length(iindex)]
  datain = mclapply(ibkgd, countup, indexfile=ibkgdindex, filelist = ibam, cornum = cornum, presched = presched, RCV=ReadChromVector, calcreadlen=FALSE)
  print(paste("size allocated datain",object.size(datain)))
  if(!length(datain[[1]][[1]])){
    quit()
  }
  for (i in 1:length(datain)) 
  {
    areinbkgd <- unique(c(areinbkgd,names(datain[[i]])))
  }
  ifrgd <- ibam[1:(length(ibam)-numbkgd)]
  ifrgdindex <- iindex[1:(length(iindex)-numbkgd)]
  datainf = mclapply(ifrgd, countup, indexfile=ifrgdindex, filelist = ibam, cornum = cornum, presched = presched,RCV=as.character(areinbkgd), calcreadlen=TRUE)
  datain = c(datainf,datain)
  #rm(datainf);gc()
} else {
  datain = mclapply(ibam, countup, indexfile=iindex, filelist = ibam, cornum = cornum, presched = presched, RCV=ReadChromVector, calcreadlen=TRUE)
}
#print(length(datain));print(length(datain[[1]]));print(length(datain[[1]][[1]]));print(length(datain[[1]][["chrY"]]))



# ================== 
# Cross Correlation
# ==================

xcorr = list()
if(numdup-nreps > 0)
{
  chrsinbkgd <- names(datain[[nreps]])
  for (i in (nreps+1):numdup)#nreps would be the first background file, the one the others will be added to
  {
    chrsini <- names(datain[[i]])
    for (j in 1:length(datain[[i]]))
    {
      if (chrsini[[j]] %in% chrsinbkgd){
        datain[[nreps]][[chrsini[[j]]]][[1]] <- datain[[nreps]][[chrsini[[j]]]][[1]] + datain[[i]][[chrsini[[j]]]][[1]]
        datain[[nreps]][[chrsini[[j]]]][[2]] <- datain[[nreps]][[chrsini[[j]]]][[2]] + datain[[i]][[chrsini[[j]]]][[2]]
      } else {
        chrsinbkgd <- c(chrsinbkgd,chrsini[[j]])
        #print(paste("chrsinbkgd",chrsinbkgd))
        datain[[nreps]][[chrsini[[j]]]][[1]] <- datain[[i]][[chrsini[[j]]]][[1]]
        datain[[nreps]][[chrsini[[j]]]][[2]] <- datain[[i]][[chrsini[[j]]]][[2]]
      }
    }
  }
}

#print(paste("size allocated datain",object.size(datain)))

cat("\n","Starting xcorr calculation: ",nreps," files")
for (i in 1:nreps)
{
  xcorr[[i]] <- list()
  cnames <- names(datain[[i]])
  cat("\n","File",i,"with",length(cnames),"elements","\n")
  xcorr[[i]] = mclapply(datain[[i]], xc, cornum = cornum, presched = presched)
}
#=======================> DONE! 


# ===================== 
# Write result to File
# =====================
for (i in 1:length(xcorr)) {
  bamname=ibam[i]
  filename = strsplit(ibam[i], "/", fixed = TRUE)[[1]]
  filename = filename[length(filename)]
  filename = substr(filename,1,nchar(filename)-4)
  if (i == length(xcorr)) {
    filename = "ctrl"
  }
  cnames <- names(datain[[i]])
  for (element in cnames)
  {
    correlation <- xcorr[[i]][[element]]
    if ((correlation <= 500) && (correlation >= 50)) { #only write this to file if xcorr value was plausible
      #message(paste0(chromName, ", ", filename, ": Ok!"))
      message(paste0(element, ", ", filename, ": Ok!"))
      #write(paste(chromName, xcorr[i], sep = "\t"), file = paste0(storeFile, "/xc.", filename, ".tab"), append = TRUE)
      write(paste(element, correlation, sep = "\t"), file = paste0(storeFile, "/xc.", filename, ".tab"), append = TRUE)
    } else {
      if (correlation) {
      } else {
        message(paste0(element, ", ", filename, ": Value Not Used!"))
      }
    }
  }
  
}
rm(datain)
rm(xcorr);gc()
#=======================> DONE!


#=======================> DONE! 


