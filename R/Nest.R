library(data.table)
library(foreach)
library(stringi)
library(matrixStats)


#-----------------------------------------
# Simulate allele frequency trajectories -
#-----------------------------------------

wf.traj <- function(p0, Ne, t, s=0, h=0.5) {
  
  # initialize trajectory matrix
  traj <- matrix(NA, ncol=length(t), nrow=max(length(p0), length(s), length(h)), dimnames=list(c(), paste0("F", t)))
  if(0 %in% t)
    traj[,"F0"] <- p0
  
  g <- 1
  p <- p0
  q <- 1-p0
  wAA <- 1+s
  wAa <- 1+h*s
  waa <- 1
  # simulate allele frequencies across time
  while(g <= max(t)) {
    # apply selection and random drift
    p <- (wAA*p^2 + wAa*p*q) / (wAA*p^2 + wAa*2*p*q + waa*q^2)
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
                  if(length(size) == 1) list(p.smpld=p.smpld, size=cov) else p.smpld
                },
                individuals={
                  # if length of 'size' is larger than 1, then send warning message
                  if(length(size) > 1)
                    warning("Only the first element in 'size' will be used, because sampling mode is 'individuals'")
                  
                  # sample random allele frequencies from Hypergeometric distribution
                  rhyper(nn=maxlen, m=p*census, n=(1-p)*census, k=size[1]*ploidy)
                }))
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
  alleles <- data.table(chr=chr, pos=pos,
                        major=names(syncCntCol)[which(t(alleleRank) == 4) - 4*seq(0, nrow(alleleRank)-1)],
                        minor=names(syncCntCol)[which(t(alleleRank) == 3) - 4*seq(0, nrow(alleleRank)-1)])
  
  # combine generation and replicate info
  popInfo <- data.table(pop=1:popCnt, gen=gen, repl=repl)
  
  # combine allele frequency data into result matrix for each replicate separately
  cntRes <- foreach(r=unique(repl), .final=function(x) { setNames(x, paste("R", unique(repl), sep="")) }) %do% {
    # create result matrix for current replicate
    g <- popInfo$gen[popInfo$repl == r]
    mat <- matrix(NA, nrow=snpCnt, ncol=sum(popInfo$repl == r)*2+2, dimnames=list(c(), c("chrID", "pos", paste("F", g, "_cnt", sep=""), paste("F", g, "_cov", sep=""))))
    mat[,"chrID"] <- chrID[chr]
    mat[,"pos"] <- pos
  
    # insert minor allele count and sequence coverage for one population after the other
    for(i in seq(1:nrow(popInfo))[popInfo$repl == r]) {
      mat[,paste("F", popInfo[i,"gen",with=FALSE], "_cnt", sep="")] <- alleleCnts[[i]][,"minor"]
      mat[,paste("F", popInfo[i,"gen", with=FALSE], "_cov", sep="")] <- rowSums(alleleCnts[[i]])
    }
    
    return(mat)
  }
  
  # if required then polarize allele counts for the rising allele
  if(rising && length(unique(popInfo$gen)) > 1) {
    # calculate mean allele frequency change per SNP and replicate
    ugens <- unique(popInfo$gen)
    minGen <- min(popInfo$gen)
    meanAf <- foreach(r=names(cntRes), .combine=cbind, .final=rowMeans) %do% {
      foreach(t = ugens[ugens != minGen], .combine=cbind, .final=rowMeans) %do% {
        cntRes[[r]][,paste0("F", t, "_cnt")] / cntRes[[r]][,paste0("F", t, "_cov")] - cntRes[[r]][,paste0("F", minGen, "_cnt")] / cntRes[[r]][,paste0("F", minGen, "_cov")]
      }
    }
    
    # polarize allele counts
    for(r in names(cntRes)) {
      for(t in ugens) {
        cntRes[[r]][,paste0("F", t, "_cnt")] <- ifelse(meanAf >= 0, cntRes[[r]][,paste0("F", t, "_cnt")], cntRes[[r]][,paste0("F", t, "_cov")]-cntRes[[r]][,paste0("F", t, "_cnt")])
      }
    }
    
    # add rising allele to 'allele' info
    alleles[,rising:=ifelse(meanAf >= 0, alleles$minor, alleles$major)]
  }
  
  # return list containing the chromosome information (ID and name) and allele count data (list of matrices, one matrix for each replicate)
  return(list(chr=chrID, alleles=alleles, cnts=cntRes))
}


#--------------------------------------
# Estimate effective populations size -
#--------------------------------------

checkSNP <- function(p0, pt, cov0, covt, truncAF=NA) {
  # return false if any of the following conditions is met xi=yi, xi=0, xi=1, (xi+yi)/2=0, (xi+yi)/2=1
  return(p0 != pt & p0 != 0 & p0 != 1 & (p0+pt)/2 != 0 & (p0+pt)/2 != 1 & cov0 != 0 & covt != 0 & if(is.na(truncAF)) TRUE else p0 >= truncAF & p0 <= 1-truncAF) #& pt != 0 & pt != 1
}

estimateNe <- function(p0, pt, cov0, covt, t, ploidy=2, truncAF=NA, method="P_AF2_planI", poolSize=NA, Ncensus=NA, asList=FALSE) { 
  
  # check if parameter 'method' has been set properly - stop execution otherwise
  mm <- match.arg(method, choices=c("P_AF", "P_AF2_planI", "P_AF2_planII", "P_AJ_planI", "P_AJ_planII", "JR_planI", "JR_planII", "W_planI", "W_planII"), several.ok=TRUE)
  if(length(mm) != length(method))
    stop("Unable to resolve the following method(s): ", paste("'", method[!method %in% mm], "'", sep="", collapse=", "))
  else
    method <- mm
  
  # check if poolSize parameter is set propperly - stop execution otherwise
  if(!all(is.na(poolSize)) && length(poolSize) != 2)
    stop("If not NA, length of 'poolSize' parameter has to be equal to 2: length(poolSize) = ", length(poolSize))
  
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

  #--- Nei and Tajima (1981) ---
  #if(any(grepl("^(a1|a2|b|c)$", method))) {
  #  corre <- 1/s0_mean + 1/st_mean
  #  # Na1 or Na2
  #  if(any(grepl("^a[12]$", method))) {
  #    Fa <- (1/n)*sum(((xi-yi)^2)/(xi*(1-xi)))
  #    # Na1
  #    if(any(grepl("^a1$", method))) { res <- c(res, Na1=(t-2-Fa) / (2*(Fa*(1-1/(s0_mean)) - corre))) }
  #    # Na2
  #    if(any(grepl("^a2$", method))) { res <- c(res, Na2=(t-2) / (2*(Fa - corre))) }
  #  }
  #  # Nb
  #  if(any(grepl("^b$", method))) { Fb <- (1/n)*sum(((xi-yi)^2)/(zi*(1-zi))); res <- c(res, Nb=(t-2)*(1+Fb/4) / (2*(Fb*(1-1/4*corre) - corre))) }
  #  # Nc
  #  if(any(grepl("^c$", method))) { Fc <- (1/n)*sum(((xi-yi)^2)/(zi-xi*yi)); res <- c(res, Nc=(t-2) / (2*(Fc - corre))) }
  #}
  
  #--- Waples (1989) ---
  if(any(grepl("^W_plan(I|II)$", method))) {
    Fc <- ((xi-yi)^2)/(zi-xi*yi)
    # W_planI
    if(any(grepl("^W_planI$", method))) {
      Fc_planI <- (1/n)*sum( Fc - (1/s0 + 1/st) + 1/Ncensus )
      res <- c(res, Nw_planI=-t/(ploidy*log(1-Fc_planI)))
    }
    # W_planII
    if(any(grepl("^W_planII$", method))) {
      Fc_planII <- (1/n)*sum( Fc - (1/s0 + 1/st))
      res <- c(res, Nw_planII=-t/(ploidy*log(1-Fc_planII)))
    }
  }
  
  #--- Jorde and Ryman (2007) ---
  if(any(grepl("^JR_plan(I|II)$", method))) {
    F_nom <- (xi-yi)^2
    F_denom <- zi*(1-zi)
    n_harmonic <- 1/(1/s0 + 1/st)
    # JR_planI
    if(any(grepl("^JR_planI$", method))) {
      Fs_planI_nom <- sum(F_nom * (1 - 1/(4*n_harmonic) + 1/(4*Ncensus)) * n_harmonic * Ncensus + F_denom * (n_harmonic - Ncensus))/sum(F_denom * n_harmonic * Ncensus)
      Fs_planI_denom <- sum((4*F_denom + F_nom) * (1 - 1/st))/sum(4*F_denom)
      Fs_planI <- Fs_planI_nom/Fs_planI_denom
      res <- c(res, Njr_planI=-t/(ploidy*log(1-Fs_planI)))
    }
    # JR_planII
    if(any(grepl("^JR_planII$", method))) { 
      Fs_planII_nom = sum(F_nom * (1-1/(4*n_harmonic)) * n_harmonic - F_denom )/sum( F_denom * n_harmonic )
      Fs_planII_denom = sum((4*F_denom + F_nom) * (1 - 1/st))/sum( 4*F_denom )
      Fs_planII = Fs_planII_nom/Fs_planII_denom
      res <- c(res, Njr_planII=-t/(ploidy*log(1-Fs_planII)))
    }
  }
  
  #--- Andreas Futschik ---
  if(any(grepl("^P_AF", method))) {
    # set correction term for 1-step sampling (plan II)
    if(any(grepl("^P_AF$", method))) {
      S <- sum(xi*(1-xi)*1/s0 + yi*(1-yi)*1/st)
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np_af=-t/(ploidy*log(1-Ft)))
    }
    # ... 2-step sampling (plan II)
    if(any(grepl("^P_AF2_planII$", method))) {
      S <- sum(xi*(1-xi)*(1/s0+1/(ploidy*poolSize[1])-1/(s0*ploidy*poolSize[1])) + yi*(1-yi)*(1/st+1/(ploidy*poolSize[2])-1/(st*ploidy*poolSize[2])))
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np_af2_planII=-t/(ploidy*log(1-Ft)))
    }
    # ... 2-step sampling (plan I)    ----------------TODO: NCENSUS AND POOLSIZE MUST BE MULTIPLIED WITH PLOIDY------------------------
    if(any(grepl("^P_AF2_planI$", method))) {
      S <- sum(xi*(1-xi)*(1/s0+1/(ploidy*poolSize[1])*(Ncensus-poolSize[1])/(Ncensus-1)-1/(s0*ploidy*poolSize[1])) + yi*(1-yi)*(1/st+1/(ploidy*poolSize[2])*(Ncensus-poolSize[2])/(Ncensus-1)-1/(st*ploidy*poolSize[2])))
      Ft <- (sum((xi-yi)^2) - S) / sum(xi*(1-xi))
      res <- c(res, Np_af2_planI=-t/(ploidy*log(1-Ft)))
    } 
  }
  
  #--- Agnes Jonas ---
  if(any(grepl("^P_AJ", method))) {
    C0i <- 1/s0 + 1/(ploidy*poolSize[1]) - 1/(s0*ploidy*poolSize[1])
    Cti <- 1/st + 1/(ploidy*poolSize[2]) - 1/(st*ploidy*poolSize[2])
    
    # sampling plan II
    if(any(grepl("^P_AJ_planII$", method))) {
      Ft <- sum((xi - yi)^2 - (zi-xi*yi)*( C0i + Cti )) / sum( (zi-xi*yi) * (1 - Cti))
      res <- c(res, Np_aj2_planII=-t/(ploidy*log(1-Ft)))
    }
    # sampling plan I    ----------------TODO: NCENSUS MUST BE MULTIPLIED WITH PLOIDY------------------------
    if(any(grepl("^P_AJ_planI$", method))) {
      Ft <- sum((xi - yi)^2 * (1 - 1/(ploidy*Ncensus)) - (zi-xi*yi)*( C0i + Cti - 1/Ncensus)) / sum( (zi-xi*yi) * (1 - Cti))
      res <- c(res, Np_aj2_planI=-t/(ploidy*log(1-Ft)))
    }
  }
  
  # return either numeric vector or list
  if(asList)
    return(as.list(res))
  
  return(res)
}

estimateWndNe <- function(chr, pos, wndSize, p0, pt, cov0, covt, t, unit=c("bp", "SNP"), ploidy=2, truncAF=NA, method="P_AF2_planI", poolSize=NA, Ncensus=NA) {

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
