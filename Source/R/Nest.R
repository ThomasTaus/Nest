library(data.table)
library(foreach)
library(stringi)
library(matrixStats)


#-----------------------------------------
# Simulate allele frequency trajectories -
#-----------------------------------------

wf.traj <- function(p0, Ne, t, s=0, h=0.5, haploid=FALSE) {
  
  # initialize trajectory matrix
  traj <- matrix(NA, ncol=length(t), nrow=max(length(p0), length(s), length(h)), dimnames=list(c(), paste0("F", t)))
  if(0 %in% t)
    traj[,"F0"] <- p0
  
  if(!haploid)
    Ne <- 2*Ne
  
  g <- 1
  p <- p0
  q <- 1-p0
  wAA <- 1+s
  wAa <- 1+h*s
  waa <- 1
  
  # simulate allele frequencies across time
  while(g <= max(t)) {
    # compute mean fitness
    w <- if(haploid) p*wAA + q*waa else p^2*wAA + 2*p*q*wAa + q^2*waa
    
    # apply selection and random drift
    p <- if(haploid) p*wAA/w else (wAA*p^2 + wAa*p*q)/w
    if(!is.na(Ne))
      p <- rbinom(length(p), Ne, p) / Ne
    q <- 1-p
    
    # if necessary then save current allele frequency to results
    if(g %in% t)
      traj[,paste0("F", g)] <- p
    
    g <- g+1
  }
  
  if(nrow(traj) == 1)
    return(as.vector(traj))
  
  return(traj)
}
  
sample.alleles <- function(p, size, mode=c("coverage", "individuals"), Ncensus=NA, ploidy=2) {  
  # determin number of return values
  maxlen <- max(length(p), length(size))
  # check length of a, n and size parameter
  if(maxlen %% length(p) != 0 || maxlen %% length(size) != 0)
    warning("Parameters differ in length and are not multiple of one another: length(p)=", length(p), ", length(size)=", length(size))
  
  # check mode parameter
  mode <- match.arg(mode)
  
  # sample allele counts according to the specified method
  return(switch(mode,
                coverage={
                  # if length of 'size' equals '1' then generate target coverage values using the Poisson distribution, otherwise use values of 'size' directly
                  cov <- if(length(size) == 1) rpois(n=maxlen, lambda=size) else size
                  # sample allele frequencies from Binomial distribution based on 'cov' and 'p'
                  p.smpld <- rbinom(n=maxlen, size=cov, prob=p) / cov
                  # return results, including coverage values if they were drawn from a Poisson distribution
                  if(length(size) == 1) data.table(p.smpld=p.smpld, size=cov) else p.smpld
                },
                individuals={
                  # if length of 'size' is larger than 1, then send warning message
                  if(length(size) > 1)
                    warning("Only the first element in 'size' will be used, because sampling mode is 'individuals'")
                  
                  # sample random allele frequencies from Hypergeometric distribution
                  rhyper(nn=maxlen, m=p*Ncensus*ploidy, n=(1-p)*Ncensus*ploidy, k=size[1]*ploidy)/(size[1]*ploidy)
                }))
}


#---------------------------------------
# The Sync-class and related functions -
#---------------------------------------

#----- class definition -----
setClass(Class="sync",
         
         # member variables
         slots=c(
           gen="numeric",
           repl="numeric",
           isAF="logical",
           alleles="data.table"), 
         
         # default constructor
         prototype=list(
           gen=numeric(0),
           repl=numeric(0),
           isAF=logical(0),
           alleles=NULL),
         
         # validation
         validity=function(object) {
           # length of 'gen', 'repl' and 'isAF' do not all match
           if(!(length(object@gen) == length(object@repl) && length(object@gen) == length(object@isAF))) {
             return(paste0("Lengths of 'gen' (", length(object@gen), "), 'repl' (", length(object@repl), ") and 'isAF' (", length(object@isAF), ") do not all match."))
           }
           # first six column names are not set as expected
           if(!is.null(object@alleles) && !isTRUE(all.equal(colnames(object@alleles)[1:6], c("chr", "pos", "ref", "major", "minor", "rising")))) {
             return(paste0("The names of the first five columns in 'alleles' must be 'c(\"chr\", \"pos\", \"ref\", \"major\", \"minor\", \"rising\") in exactly that order"))
           }
           # column number of 'alleles' does not match the length of 'gen'
           if((is.null(object@alleles) && length(object@gen) != 0) || (!is.null(object@alleles) && ncol(object@alleles) != length(object@gen)+6)) {
             return(paste0("Column number of 'alleles' (", ncol(object@alleles), ") has to be equal to the length of 'gen' (", length(object@gen), ") + 6."))
           }
           # neither allele frequencies nor sequence coverages are provided
           if(length(object@gen) == 0) {
             return("Neither allele frequencies nor sequence coverage are provided")
           }
           # data types of columns in 'alleles' do not match
           for(cc in c(2, seq(7, length(object@gen)+6))) {
             if(!is.numeric(object@alleles[[cc]])) {
               return(paste0("Column ", cc, " in 'alleles' should be numeric."))
             }
           }
           for(cc in c(1, 3:6)) {
             if(!is.character(object@alleles[[cc]])) {
               return(paste0("Column ", cc, " in 'alleles' should be character."))
             }
           }
           
           # allele frequency is outside the range [0, 1]
           for(cc in which(object@isAF)+6) {
             if(!all(object@alleles[[cc]] >= 0, object@alleles[[cc]] <= 1, na.rm=TRUE)) {
               return(paste0("Allele frequencies in column ", cc, " are outside the range [0, 1]."))
             }
           }
           
           # coverage is < 0
           for(cc in which(!object@isAF)+6) {
             if(!all(object@alleles[[cc]] >= 0, na.rm=TRUE)) {
               return(paste0("Sequence coverages in column ", cc, " are < 0."))
             }
           }
           
           # object is valid
           return(TRUE)
         }
)

#----- initializer -----
setMethod("initialize",
          signature="sync",
          definition=function(.Object, gen, repl, isAF, alleles) {
            .Object@gen <- gen
            .Object@repl <- repl
            .Object@isAF <- isAF
            .Object@alleles <- alleles
            # call the inspector
            validObject(.Object)
            # if 'alleles' is specified then add column with 'posID' and set key
            if(!is.null(.Object@alleles)) {
              .Object@alleles[,posID:=paste(chr, pos, sep=".")]
              setkey(.Object@alleles, posID)
            }
            return(.Object)
          }
)

#----- related methods -----
is.sync <- function(x) {
  return(inherits(x, "sync"))
}

af.traj <- function(sync, chr, pos, repl) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")
  
  # extract selected SNPs
  sel <- paste(chr <- as.character(chr), pos <- as.numeric(pos), sep=".")
  subAlleles <- sync@alleles[sel]
  
  # extract time series data for specific replicate
  traj.mat <- function(r, include=FALSE) {
    cols <- which(sync@repl == r & sync@isAF)
    res <- subAlleles[,cols[order(sync@gen[cols])]+6,with=FALSE]
    res <- as.matrix(setnames(res, sub("\\.freq", "", colnames(res))))
    rownames(res) <- if(include) paste0(sel, ".R", r) else sel
    colnames(res) <- sub("\\.R[0-9]+", "", colnames(res))
    return(res)
  }
  
  # if only one replicate is specified then create 1xt matrix for trajectory with t time points
  if(length(repl <- as.numeric(repl)) == 1) {
    traj <- traj.mat(repl, include=FALSE)
  # if multiple replicates are specified 
  } else if(length(repl) > 1) {
    # then create list of trajectory matrices, one for each replicate
    traj <- foreach(r=repl) %do% {
      traj.mat(r, include=TRUE)
    }
    # if the number of time points (columns) is equal for all replicates then combine all matrices to one
    if(all((x <- sapply(traj, ncol)) == x[1])) {
      traj <- Reduce(rbind, traj)
    }
  }
  
  return(traj)
}

af <- function(sync, repl, gen) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")
  
  # extract allele frequencies of specified replicates and time points
  cols <- which(sync@gen %in% (gen <- as.numeric(gen)) & sync@repl %in% (repl <- as.numeric(repl)) & sync@isAF)
  if(length(cols) == 0)
    stop("The combination of 'repl' (", paste(repl, sep=","), ") and 'gen' (", paste(gen), ") does not match the data available in 'sync'.")
  subAlleles <- sync@alleles[,cols[order(sync@repl[cols], sync@gen[cols])]+6,with=FALSE]
  afMat <- as.matrix(subAlleles)
  rownames(afMat) <- sync@alleles$posID
  
  return(afMat[order(sync@alleles$chr, sync@alleles$pos),])
}

coverage <- function(sync, repl, gen) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")
  
  # extract allele frequencies of specified replicates and time points
  cols <- which(sync@gen %in% (gen <- as.numeric(gen)) & sync@repl %in% (repl <- as.numeric(repl)) & !sync@isAF)
  if(length(cols) == 0)
    stop("The combination of 'repl' (", paste(repl, sep=","), ") and 'gen' (", paste(gen), ") does not match the data available in 'sync'.")
  subAlleles <- sync@alleles[,cols[order(sync@repl[cols], sync@gen[cols])]+6,with=FALSE]
  afMat <- as.matrix(subAlleles)
  rownames(afMat) <- sync@alleles$posID
  
  return(afMat[order(sync@alleles$chr, sync@alleles$pos),])
}

alleles <- function(sync) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")
  
  return(sync@alleles[order(chr, pos),1:6,with=FALSE])
}

splitLocusID <- function(id, sep=".") {
  splitMat <- stri_split(id <- as.character(id), fixed=(sep <- as.character(sep)), simplify=TRUE)
  
  if(!is.matrix(splitMat) || dim(splitMat)[1] != length(id) || dim(splitMat)[2] != 2)
    stop("Locus IDs could not be split successfully.")
  
  return(data.table(chr=splitMat[,1], pos=as.numeric(splitMat[,2])))
}


#----------------------------------------------------
# Convert Sync-file to allele frequency data object -
#----------------------------------------------------

read.sync <- function(file, gen, repl, rising=FALSE) {
  
  cat("Reading sync file ...\n")
  # load sync-file
  syncDt <- read.delim(file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  setDT(syncDt)
  # if either 'gen' or 'repl' are not of propper length then stop
  if((ncol(syncDt)-3) %% length(gen) != 0 || (ncol(syncDt)-3) %% length(repl) != 0)
    stop("Either 'gen' (", length(gen), ") or 'repl' (", length(repl), ") is not a multiple of the number of populations (", ncol(syncDt)-3, ") specified in the sync-file.")
  gc()
  
  cat("Extracting biallelic counts ...\n")
  # extract numeric allele counts for A:T:C:G
  syncCnts <- lapply(syncDt[,-1:-3,with=FALSE], function(col) {
    matrix(as.numeric(stri_split(col, fixed=":", simplify=TRUE)), ncol=6)[,1:4]
  } ) # <- bottleneck
  
  # get rid of data that is no longer needed
  snpCnt <- nrow(syncDt)
  chr <- syncDt$V1
  pos <- syncDt$V2
  ref <- syncDt$V3
  popCnt <- ncol(syncDt)-3
  rm(syncDt)
  gc()
  
  # sum up counts across populations (time points and replicates) and add e (random uniform >=0 & <= 0.99) to make each count value unique
  sumCnts <- Reduce("+", syncCnts)
  sumCnts <- sumCnts + runif(nrow(sumCnts)*ncol(sumCnts), min=0, max=0.99)
  # deterime allele ranks for each position
  alleleRank <- rowRanks(sumCnts, ties.method="min")
  # extract 2 most common alleles (considering all populations)
  alleleCnts <- lapply(syncCnts, function(pop) {
    cbind(major=t(pop)[t(alleleRank) == 4], minor=t(pop)[t(alleleRank) == 3])
  } ) # <- bottleneck
  
  cat("Creating result object ...\n")
  # compute chromosome IDs (to later replace character vector by numeric one)
  chrNames <- unique(chr)
  chrID <- 1:length(chrNames)
  names(chrID) <- chrNames
  
  # extract major and minor allele for each position
  syncCntCol <- 1:4
  names(syncCntCol) <- c("A", "T", "C", "G")
  alleles <- data.table(chr=chr, pos=pos, ref=ref,
                        major=names(syncCntCol)[which(t(alleleRank) == 4) - 4*seq(0, nrow(alleleRank)-1)],
                        minor=names(syncCntCol)[which(t(alleleRank) == 3) - 4*seq(0, nrow(alleleRank)-1)],
                        rising=NA_character_)
  
  # combine generation and replicate info
  popInfo <- data.table(pop=1:popCnt, gen=gen, repl=repl)
  
  # add allele frequency and sequence coverage for each population
  for(r in unique(repl)) {
    for(i in seq(1:nrow(popInfo))[popInfo$repl == r]) {
      seqCov <- rowSums(alleleCnts[[i]])
      alleles[,paste0("F", popInfo$gen[i], ".R", r, ".freq"):=alleleCnts[[i]][,"minor"]/seqCov]
      alleles[,paste0("F", popInfo$gen[i], ".R", r, ".cov"):=seqCov]
    }
  }
  
  # if required then polarize allele counts for the rising allele
  if(rising && length(unique(popInfo$gen)) > 1) {
    ugens <- unique(popInfo$gen)
    minGen <- min(popInfo$gen)
    
    # if minGen allele frequency column is not available for all replicates then stop execution
    if(sum(grepl(paste0("F", minGen, "\\.[R0-9]+\\.freq"), colnames(alleles))) != length(unique(repl)))
      stop("Not all replicates provide allele frequency estimates at generation F", minGen)
    
    # calculate mean allele frequency change per SNP and replicate
    meanAF <- foreach(r=unique(repl), .combine=cbind, .final=function(x) { if(is.matrix(x)) return(rowMeans(x)) else return(x) }) %do% {
      rowMeans(alleles[,grepl(paste0("F[0-9]+\\.R", r, "\\.freq"), colnames(alleles)),with=FALSE]-alleles[[which(grepl(paste0("F", minGen, "\\.[R0-9]+\\.freq"), colnames(alleles)))[1]]])
    }
    
    # polarize allele frequencies
    needsPolarization <- meanAF < 0
    for(pop in grep("F[0-9]+\\.R[0-9]+\\.freq", colnames(alleles), value=TRUE)) {
      alleles[,eval(pop):=ifelse(needsPolarization, 1-alleles[[pop]], alleles[[pop]])]
    }
    
    # set column with rising allele
    alleles[,rising:=ifelse(needsPolarization, alleles$major, alleles$minor)]
  }
  
  # return sync-object for loaded sync-file
  return(new(Class="sync",
             gen=as.numeric(sub("F([0-9]+)\\.R[0-9]+.*", "\\1", colnames(alleles)[-1:-6])),
             repl=as.numeric(sub(".*\\.R([0-9]+)\\..*", "\\1", colnames(alleles)[-1:-6])),
             isAF=grepl(".*\\.freq$", colnames(alleles)[-1:-6]),
             alleles=alleles))
}


#--------------------------------------
# Estimate effective populations size -
#--------------------------------------

checkSNP <- function(p0, pt, cov0, covt, truncAF=NA) {
  # return false if any of the following conditions is met: xi==0, xi==1
  # mask extreme allele frequencies if truncAF is unequal to 'NA'
  return(p0 != 0 & p0 != 1 & cov0 != 0 & covt != 0 & if(is.na(truncAF)) TRUE else p0 >= truncAF & p0 <= 1-truncAF)
}

estimateNe <- function(p0, pt, cov0, covt, t, ploidy=2, truncAF=NA, method="P.planI", Ncensus=NA, poolSize=rep(Ncensus, times=2), asList=FALSE) { 
  
  # check if parameter 'method' has been set properly - stop execution otherwise
  mm <- match.arg(method, choices=c("P.planI", "P.planII",
                                    "JR.planI", "JR.planII",
                                    "W.planI", "W.planII",
                                    "P.alt.1step.planII", "P.alt.2step.planI", "P.alt.2step.planII"), several.ok=TRUE)
  if(length(mm) != length(method))
    stop("Unable to resolve the following method(s): ", paste("'", method[!method %in% mm], "'", sep="", collapse=", "))
  else
    method <- mm
  
  # check if poolSize parameter is set propperly - stop execution otherwise
  if(length(poolSize <- as.numeric(poolSize)) != 2)
    stop("Length of 'poolSize' parameter has to be equal to 2: length(poolSize) = ", length(poolSize))
  
  # remove SNPs for which Ne cannot/should not be estimated
  keep <- checkSNP(p0, pt, cov0, covt, truncAF=truncAF)
  xi <- p0[keep]
  yi <- pt[keep]
  zi = (xi + yi)/2
  n = length(xi)
  # coverage is divided by 2, because later 2*S0, 2*St correction term will be used
  s0 <- cov0[keep]
  st <- covt[keep]
  s0_mean <- mean(s0)
  st_mean <- mean(st)
  
  # estimate Ne using the specified method(s)
  res <- numeric(length=0)
  
  #--- Waples (1989) ---
  if(any(grepl("^W\\.plan(I|II)$", method))) {
    Fc <- ((xi-yi)^2)/(zi-xi*yi)
    # W_planI
    if(any(grepl("^W\\.planI$", method))) {
      Fc_planI <- (1/n)*sum( Fc - (1/s0 + 1/st) + 1/Ncensus )
      res <- c(res, Nw.planI=-t/(ploidy*log(1-Fc_planI)))
    }
    # W_planII
    if(any(grepl("^W\\.planII$", method))) {
      Fc_planII <- (1/n)*sum( Fc - (1/s0 + 1/st))
      res <- c(res, Nw.planII=-t/(ploidy*log(1-Fc_planII)))
    }
  }
  
  #--- Jorde and Ryman (2007) ---
  if(any(grepl("^JR\\.plan(I|II)$", method))) {
    F_nom <- (xi-yi)^2
    F_denom <- zi*(1-zi)
    n_harmonic <- 1/(1/s0 + 1/st)
    # JR_planI
    if(any(grepl("^JR\\.planI$", method))) {
      Fs_planI_nom <- sum(F_nom * (1 - 1/(4*n_harmonic) + 1/(4*Ncensus)) * n_harmonic * Ncensus + F_denom * (n_harmonic - Ncensus))/sum(F_denom * n_harmonic * Ncensus)
      Fs_planI_denom <- sum((4*F_denom + F_nom) * (1 - 1/st))/sum(4*F_denom)
      Fs_planI <- Fs_planI_nom/Fs_planI_denom
      res <- c(res, Njr.planI=-t/(ploidy*log(1-Fs_planI)))
    }
    # JR_planII
    if(any(grepl("^JR\\.planII$", method))) { 
      Fs_planII_nom = sum(F_nom * (1-1/(4*n_harmonic)) * n_harmonic - F_denom )/sum( F_denom * n_harmonic )
      Fs_planII_denom = sum((4*F_denom + F_nom) * (1 - 1/st))/sum( 4*F_denom )
      Fs_planII = Fs_planII_nom/Fs_planII_denom
      res <- c(res, Njr.planII=-t/(ploidy*log(1-Fs_planII)))
    }
  }
  
  #--- Andreas Futschik ---
  if(any(grepl("^P\\.alt", method))) {
    # set correction term for 1-step sampling (plan II)
    if(any(grepl("^P\\.alt\\.1step\\.planII$", method))) {
      S <- sum(xi*(1-xi)*1/s0 + yi*(1-yi)*1/st)
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np.alt.1step.planII=-t/(ploidy*log(1-Ft)))
    }
    # ... 2-step sampling (plan II)
    if(any(grepl("^P\\.alt\\.2step\\.planII$", method))) {
      S <- sum(xi*(1-xi)*(1/s0+1/(ploidy*poolSize[1])-1/(s0*ploidy*poolSize[1])) + yi*(1-yi)*(1/st+1/(ploidy*poolSize[2])-1/(st*ploidy*poolSize[2])))
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np.alt.2step.planII=-t/(ploidy*log(1-Ft)))
    }
    # ... 2-step sampling (plan I)    ----------------TODO: NCENSUS AND POOLSIZE MUST BE MULTIPLIED WITH PLOIDY------------------------
    if(any(grepl("^P\\.alt\\.2step\\.planI$", method))) {
      S <- sum(xi*(1-xi)*(1/s0+1/(ploidy*poolSize[1])*(Ncensus-poolSize[1])/(Ncensus-1)-1/(s0*ploidy*poolSize[1])) + yi*(1-yi)*(1/st+1/(ploidy*poolSize[2])*(Ncensus-poolSize[2])/(Ncensus-1)-1/(st*ploidy*poolSize[2])))
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np.alt.2step.planI=-t/(ploidy*log(1-Ft)))
    } 
  }
  
  #--- Agnes Jonas ---
  if(any(grepl("^P\\.plan", method))) {
    C0i <- 1/s0 + 1/(ploidy*poolSize[1]) - 1/(s0*ploidy*poolSize[1])
    Cti <- 1/st + 1/(ploidy*poolSize[2]) - 1/(st*ploidy*poolSize[2])
    
    # sampling plan II
    if(any(grepl("^P\\.planII$", method))) {
      Ft <- sum((xi - yi)^2 - (zi-xi*yi)*( C0i + Cti )) / sum( (zi-xi*yi) * (1 - Cti))
      res <- c(res, Np.planII=-t/(ploidy*log(1-Ft)))
    }
    # sampling plan I    ----------------TODO: NCENSUS MUST BE MULTIPLIED WITH PLOIDY------------------------
    if(any(grepl("^P\\.planI$", method))) {
      Ft <- sum((xi - yi)^2 * (1 - 1/(ploidy*Ncensus)) - (zi-xi*yi)*( C0i + Cti - 1/Ncensus)) / sum( (zi-xi*yi) * (1 - Cti))
      res <- c(res, Np.planI=-t/(ploidy*log(1-Ft)))
    }
  }
  
  # return either numeric vector or list
  if(asList)
    return(as.list(res))
  
  return(res)
}

estimateWndNe <- function(chr, pos, wndSize, p0, pt, cov0, covt, t, unit=c("bp", "SNP"), ploidy=2, truncAF=NA, method="P.planI", Ncensus=NA, poolSize=rep(Ncensus, times=2)) {

  # check unit parameter
  unit <- match.arg(unit)
  
  # organize all input vectors into one data table
  dataDt <- data.table(chr=as.character(chr), pos=pos, p0=p0, pt=pt, cov0=cov0, covt=covt)
  setkey(dataDt, chr)
  
  # remove SNPs for which Ne cannot/should not be estimated
  dataDt <- dataDt[checkSNP(p0, pt, cov0, covt),]
  
  # go through each chromosome and return results
  return(foreach(cc=unique(dataDt$chr), .combine=rbind) %do% {
    # get SNPs of current chromosome
    dataSubDt <- dataDt[cc]
    # assign window ID to each SNP based on specified windowSize and unit
    if(is.na(wndSize)) {
      dataSubDt[,wndID:=1]
    } else {
      dataSubDt[,wndID:=switch(unit, bp=floor(dataSubDt$pos/wndSize)+1, SNP=floor((rank(dataSubDt$pos)-1)/wndSize)+1)]
    }
    # estimate Ne for each window
    wndRes <- dataSubDt[,estimateNe(p0=p0, pt=pt, cov0=cov0, covt=covt, t=t, ploidy=ploidy, truncAF=truncAF, method=method, poolSize=poolSize, Ncensus=Ncensus, asList=TRUE), by=wndID]
    # add columh for chromosome name, window start/stop pos and number of SNPs per window
    wndRes <- wndRes[order(wndID)]
    wndMinMax <- switch(unit, bp=NA, SNP=dataSubDt[,list(min=min(pos), max=max(pos)),by=wndID])
    wndRes[,c("chr", "wndStart", "wndStop", "SNPs"):=list(cc,
                                                          switch(unit, bp=1+(wndID-1)*wndSize, SNP=wndMinMax$min[wndRes$wndID]),
                                                          switch(unit, bp=wndID*wndSize, SNP=wndMinMax$max[wndRes$wndID]),
                                                          table(dataSubDt$wndID)[as.character(wndRes$wndID)])]
    
    return(wndRes[,!grepl("^wndID$", colnames(wndRes)),with=FALSE])
  } )
}
