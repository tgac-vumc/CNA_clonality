clonality.analysis.mod <- function (data, ptlist, pfreq = NULL, refdata = NULL, nmad = 1.25, 
                                    reference = TRUE, allpairs = TRUE, segmethod = "oneseg", 
                                    segpar = NULL) 
{
  if (!inherits(data, "CNA")) 
    stop("First arg must be a copy number array(CNA) object\n")
  if (segmethod == "oneseg") {
    if (missing(segpar)) {
      segpar$alpha <- 0.01
      segpar$nperm <- 2000
      segpar$sbdry <- get.oneseg.bdry(alpha = 0.01, nperm = 2000)
    }
    else {
      if (exists("alpha", segpar) & exists("nperm", segpar)) {
        segpar$sbdry <- get.oneseg.bdry(alpha = segpar$alpha, 
                                        nperm = segpar$nperm)
      }
      else {
        stop("segpar missing alpha and/or nperm")
      }
    }
  }
  else {
    warning(" you have given an alternate segmentation function\n ensure that this function and segpar meet the requirements\n")
  }
  npts <- length(unique(ptlist))
  if (!(any(table(ptlist) >= 2))) 
    stop("No pairs of tumors from the same patient. Check  ptlist")
  if (npts == 1 & reference) {
    warning("Just one patient - can't evaluate reference distribution\n")
    reference <- FALSE
  }
  data <- data[data$chrom != "chrX" & data$chrom != "chrY" & 
                 data$chrom != "ChrX" & data$chrom != "ChrY" & data$chrom != 
                 23 & data$chrom != 24, ]
  if (all(is.numeric(data$chrom))) 
    warning("Chromosomes should be split into p and q arms to increase power")
  if (any(table(ptlist) < 2)) 
    warning("Warning: some patients have only 1 tumor")
  if (nrow(data) >= 15000) 
    warning("Averaging is highly recommended; use ave.adj.probes() function first")
  samnms <- names(data)[-c(1, 2)]
  chrlist <- unique(data$chrom)
  if (any(table(data$chrom) < 5)) {
    cat("Removing the following chromosomes since they have fewer than 5 markers\n")
    cat(paste(names(table(data$chrom))[table(data$chrom) < 
                                         5], "\n"))
  }
  data <- data[!(data$chrom %in% chrlist[table(data$chrom) < 
                                           5]), ]
  chrlist <- unique(data$chrom)
  nchr <- length(chrlist)
  data.seg1 <- segment1(data, segmethod = segmethod, segpar = segpar)
  data.seg1GL <- GL(data.seg1, nmad)
  data.seg1 <- data.seg1GL[[1]]
  classall <- class.all(data.seg1)
  pfreq <- calc.freq(pfreq, refdata, classall, nmad, segmethod = segmethod, 
                     segpar = segpar)
  ptLR <- calculateLR(data.seg1, classall, ptlist, pfreq, reference = FALSE, 
                      allpairs = allpairs, gainthres = data.seg1GL[[2]], lossthres = data.seg1GL[[3]], 
                      segmethod = segmethod, segpar = segpar)
  if (reference) {
    refLR <- calculateLR(data.seg1, classall, ptlist, pfreq, 
                         reference = TRUE, allpairs = allpairs, gainthres = data.seg1GL[[2]], 
                         lossthres = data.seg1GL[[3]], segmethod = segmethod, 
                         segpar = segpar)
    LR2pvalue <- NULL
    for (i in 1:nrow(ptLR)) LR2pvalue <- c(LR2pvalue, mean(ptLR[i, 
                                                                4] <= refLR[, 4], na.rm = TRUE))
    ptLR <- cbind(ptLR, LR2pvalue)
  }
  if (reference) 
    return(list(LR = ptLR, OneStepSeg = data.seg1, ChromClass = classall, 
                refLR = refLR))
  else return(list(LR = ptLR, OneStepSeg = data.seg1, ChromClass = classall))
}

segment1 <-
  function(xcna, segmethod, segpar) {
    ######One step CBS - finds at most one (most prominent) copy number change 
    #on each  chromosome arm of every sample
    if (!inherits(xcna, "CNA"))
      stop("First arg must be a copy number array object")
    
    segres <- list()
    segres$data <- xcna
    outp<-NULL
    for (sam in c(1:(ncol(xcna)-2))) {
      for (chr in unique(xcna$chrom)) {
        if (!all(is.na(xcna[xcna$chrom==chr,2+sam]))) {
          mapl <- xcna$maploc[xcna$chrom==chr]
          x <- xcna[xcna$chrom==chr,2+sam]
          if (any(is.na(x))) {mapl<-mapl[!is.na(x)];  x<-x[!is.na(x)]}
          n <- length(mapl)
          
          segdata <- c(list(x=x), segpar)
          seg <- do.call(segmethod, segdata)
          if (seg[1]==0) 
            outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[1],mapl[n],n,mean(x,na.rm=TRUE)))
          if (seg[1]==1) {
            outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[1],
                               mapl[seg[3]],seg[3],mean(x[1:seg[3]],na.rm=TRUE)))
            outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[seg[3]+1],
                               mapl[n],n-seg[3],mean(x[(seg[3]+1):n],na.rm=TRUE)))
          }
          if (seg[1]==2) {
            outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[1],mapl[seg[2]-1],
                               seg[2]-1,mean(x[c(1:(seg[2]-1),(seg[3]+1):n)],na.rm=TRUE)))
            outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[seg[2]],mapl[seg[3]],
                               seg[3]-seg[2]+1,mean(x[seg[2]:seg[3]],na.rm=TRUE)))
            outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[seg[3]+1],mapl[n],
                               n-seg[3],mean(x[c(1:(seg[2]-1),(seg[3]+1):n)],na.rm=TRUE)))
          }
        }
      }
    }
    outp<-as.data.frame(outp)
    names(outp)<-c( "ID",  "chrom", "loc.start", "loc.end" ,"num.mark", "seg.mean")
    outp$seg.mean<-round(as.numeric(as.character(outp$seg.mean)),5)
    outp$num.mark<-as.numeric(as.character(outp$num.mark))
    outp$loc.start<-as.numeric(as.character(outp$loc.start))
    outp$loc.end<-as.numeric(as.character(outp$loc.end))
    segres$output <- outp
    class(segres) <- "DNAcopy"
    segres
  }

get.oneseg.bdry <-
  function(alpha, nperm) {
    max.ones <- floor(nperm*alpha) + 1
    sbdry <- DNAcopy::getbdry(eta=0.05, nperm, max.ones)
  }


oneseg <-
  function(x, alpha, nperm, sbdry) {
    # setup data and parameters for fortran code
    data.type <- "logratio"
    current.genomdat <- x - mean(x)
    current.tss <- sum(current.genomdat^2)
    n <- current.n <- length(x)
    sbn <- length(sbdry)
    
    # there are used only by hybrid code but need to set it
    kmax <- 25
    nmin <- 200
    ngrid <- 100
    tol <- 1e-6
    
    if (n>200) {
      hybrid <- TRUE
      delta<-(kmax+1)/current.n
    } else {
      hybrid <- FALSE
      delta <- 0
    }
    
    # call the fortran segmentation cose
    zzz <- .Fortran("fndcpt",
                    n = as.integer(current.n),
                    x = as.double(current.genomdat),
                    tss = as.double(current.tss),
                    px = double(current.n),
                    sx = double(n),
                    nperm = as.integer(nperm),
                    cpval = as.double(alpha),
                    ncpt = integer(1),
                    icpt = integer(2),
                    ibin = as.logical(data.type=="binary"),
                    hybrid = as.logical(hybrid),
                    al0 = as.integer(2),
                    hk = as.integer(kmax),
                    delta = as.double(delta),
                    ngrid = as.integer(ngrid),
                    sbn = as.integer(sbn),
                    sbdry = as.integer(sbdry),
                    tol = as.double(tol),
                    PACKAGE = "DNAcopy")
    
    # number of changes detected at alpha level (could be 0, 1 or 2)
    switch(1+zzz$ncpt, c(zzz$ncpt, 0, 0),
           c(zzz$ncpt, 1, zzz$icpt[1]),
           c(zzz$ncpt, c(1,0)+zzz$icpt))
  }

GL <-
  function(seg,nmad)
  {
    ###### Assignment of Gain/Loss/Normal status to each segment based on MAD criteria. 
    #Segment if Gain/Loss of its segment mean is above/below nmad MADS of noise level from the median
    ns<-dim(seg$data)[2]-2
    samnms<-names(seg$data)[-c(1,2)]
    gainthres<-rep(NA,ns)
    lossthres<-rep(NA,ns)
    seg$output[,7]<- "Normal"
    for (i in c(1:ns))
    {resid<-seg$data[!is.na(seg$data[,2+i]),2+i]-
      rep(seg$output$seg.mean[seg$output[,1]==samnms[i]],
          seg$output$num.mark[seg$output[,1]==samnms[i]])
    madd<-mad(resid,na.rm=TRUE)
    seg$output[seg$output[,1]==samnms[i]& seg$output$seg.mean>
                 median(seg$data[,2+i],na.rm=TRUE)+nmad*madd,7]<-"Gain"
    seg$output[seg$output[,1]==samnms[i] & seg$output$seg.mean<
                 median(seg$data[,2+i],na.rm=TRUE)-nmad*madd,7]<-"Loss"
    gainthres[i] <-median(seg$data[,2+i],na.rm=TRUE)+nmad*madd
    lossthres[i] <-median(seg$data[,2+i],na.rm=TRUE)-nmad*madd
    }
    names(seg$output)[7]<-"state"
    return(list(seg,gainthres,lossthres))
  }

class.all <-
  function(seg1)
  {
    ############ classifies all chomosome in 2 tumors based on one-step CBS pattern: by middle segment if there
    ## are 3 segments, by most outstanding segment if there are 2 segments.
    chrlist<-sort(unique(seg1$output$chrom))
    samnms<- names(seg1$data)[-c(1,2)]
    ns<-length(samnms)
    tum<-matrix("Normal",nrow=length(chrlist),ncol=ns)
    for (chr in c(1:length(chrlist)))
      for (pt in 1:ns)
      {
        maplocs<-seg1$data$maploc[seg1$data$chrom==chrlist[chr]]
        ss1<-subset(seg1,chrom=chrlist[chr],sample=samnms[pt])
        s1<-ss1$output
        n1<-nrow(s1)
        w1<-which(s1[,7]!="Normal")
        if (n1==1) tum[chr,pt]<-s1[,7]
        else if (length(w1)==0 | (n1==3 & s1[2,7]=="Normal")) tum[chr,pt]<-"Normal"
        else
        {if (n1==2 & length(w1)==1) ind1<-w1
        if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])>=abs(s1$seg.mean[w1[2]])) 
          ind1<-w1[1]
        if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])<abs(s1$seg.mean[w1[2]]))
          ind1<-w1[2]
        if (n1==3) ind1<-2
        tum[chr,pt]<-s1[ind1,7]
        }
      }
    colnames(tum)<-samnms
    rownames(tum)<- chrlist
    return(tum)
  }

calc.freq <-
  function(pfreq, refdata, classall, nmad, segmethod, segpar) {
    chrlist<-rownames(classall)
    if (!is.null(pfreq)) {
      if (ncol(pfreq)!=4) 
        stop("pfreq should have 4 columns: chromosome arm name, frequency of gain, loss, and normal.") 

      pfreq<-pfreq[match(chrlist,as.character(pfreq[,1])),2:4]
      pfreq<-as.data.frame(pfreq)
      pfreq[,1]<-as.numeric(as.character(pfreq[,1]))
      pfreq[,2]<-as.numeric(as.character(pfreq[,2]))
      pfreq[,3]<-as.numeric(as.character(pfreq[,3]))
      if (!is.null(refdata)) cat("Refdata is not used since pfreq is specified")
      return(pfreq)
    }
    
    if (is.null(pfreq)) {
      if (!is.null(refdata)) {
        cat("Resolution of the reference cohort arrays should be similar to the resolution of the averaged data. Please make sure. Proceeding...\n")
        refseg1<-segment1(refdata[refdata$chrom %in% chrlist,],segmethod=segmethod,segpar=segpar)
        classall<-class.all(GL(refseg1,nmad)[[1]])
      }
      
      nppl<-ncol(classall)
      if (nppl<10) 
        warning("too few patients to estimate the frequencies reliably!! Proceeding anyway...\n")
      
      pg<-apply(classall=="Gain",1,mean,na.rm=TRUE)
      pl<-apply(classall=="Loss",1,mean,na.rm=TRUE)
      
      pg[pg<0.05]<-0.05
      pl[pl<0.05]<-0.05
      pg[pg>0.9]<-0.9
      pl[pl>0.9]<-0.9
      
      w<-which(pg+pl>=0.95)    
      diff<-pg[w]+pl[w]-0.95
      pg[w]<-pg[w]- diff/2
      pl[w]<-pl[w]- diff/2
      
      pfreq<-cbind(pg,pl,1-pg-pl)
      
      if (any(pfreq<0)) {stop("negative frequencies")}
      return(pfreq)
    }
  }

calculateLR <-
  function(data.seg1, classall, ptlist, pfreq, reference,
           allpairs=TRUE, gainthres, lossthres, segmethod,
           segpar) {
    npts<-length(unique(ptlist))
    samnms<-names(data.seg1$data)[-c(1,2)]
    chrlist<-unique(data.seg1$data$chrom)
    nchr<-length(chrlist)
    testset<-NULL
    if (reference) {
      for (i in (1:(npts-1)))
        for (j in ((i+1):npts)) {
          w1<-which(ptlist==unique(ptlist)[i])
          w2<-which(ptlist==unique(ptlist)[j])
          for (p1 in 1:length(w1))
            for (p2 in 1:length(w2))
              if (allpairs | (!allpairs & p1!=p2))   testset<-rbind(testset,c(w1[p1],w2[p2]))
        }
      if (nrow(testset)<20) warning("too few patients to estimate the reference distribution reliably!! Proceeding anyway...\n")
    } else {       
      for (i in unique(ptlist)) {
        w<-which(ptlist==i) 
        ns<- length(w)
        if (ns>1) {
          for (p1 in c(1:(ns-1)))
            for (p2 in c((p1+1):ns))
              testset<-rbind(testset,c(w[p1],w[p2]))
        }
      }
    }
    if (!reference) cat("Calculating LR")
    else cat("Calculating reference LR: %completed ")
    ncomp<-nrow(testset)
    iLRs<-NULL
    for (i in (1:ncomp)) {
      if (reference & (i %in% round(c(1:10)*ncomp/10))) cat(paste(round(100*i/ncomp),", ",sep=""))
      if (!reference) cat(".")
      x<-rep(NA,nchr)
      y<-rep(NA,nchr)
      for (chr in c(1:nchr)) {
        b<-indiv.test(subset(data.seg1,chrom=chrlist[chr],sample=samnms[testset[i,1]]),
                      subset(data.seg1,chrom=chrlist[chr],sample=samnms[testset[i,2]]),func,
                      gainthres[testset[i,]],lossthres[testset[i,]],segmethod=segmethod,segpar=segpar)
        if  (!is.na(b[1]) )
        {x[chr]<-b[[3]]
        y[chr]<-b[[4]]
        }
      }
      tum1<-classall[,testset[i,1]]
      tum2<-classall[,testset[i,2]]
      a<-grantLR(tum1,tum2,pfreq,x,y,cvalue=0.5,chrlist=chrlist)
      
      pattern<-c(sum(tum1==tum2 & tum1!="Normal"),sum(tum1==tum2 & tum1=="Normal"),
                 sum(tum1!=tum2 & tum1!="Normal" & tum2!="Normal") )
      pattern<-c(pattern,length(chrlist)-sum(pattern))
      pattern<-c(pattern,paste(names(a[-c(1,2)]),round(a[-c(1,2)],2),collapse="; "))
      iLRs<-rbind(iLRs,c(samnms[testset[i,]],a[c(1,2)],pattern))
    }
    
    iLRs<-as.data.frame(iLRs)  
    names(iLRs)<-c("Sample1",	"Sample2",	"LR1"	,"LR2",	"GGorLL",	"NN",	"GL",	"GNorLN",	"IndividualComparisons")
    iLRs[,4]<-as.numeric(as.character(iLRs[,4]))
    iLRs[,3]<-as.numeric(as.character(iLRs[,3]))
    cat("\n")
    return(iLRs)
  }

indiv.test <-
  function(ss1, ss2, func, gainthres, lossthres, Nsim=100,
           segmethod, segpar) {
    ###### tests whether two chromosomes have identical changes: returns distribution of the test statistic t# under hypotheses of clonality and independence, NA if discordant changes
    
    s1<-ss1$output
    s2<-ss2$output
    
    s1c<-cumsum(s1$num.mark)
    n1<-nrow(s1)
    s2c<-cumsum(s2$num.mark)
    n2<-nrow(s2)
    
    xx1<-ss1$data[,3]
    xx2<-ss2$data[,3]
    resid1<-xx1
    resid1[!is.na(xx1)]<-
      resid1[!is.na(xx1)]-rep(ss1$output$seg.mean,ss1$output$num.mark)
    resid2<-xx2
    resid2[!is.na(xx2)]<-
      resid2[!is.na(xx2)]-rep(ss2$output$seg.mean,ss2$output$num.mark)
    chrlen<-length(xx1)
    
    ts<-func(ss1,ss2)
    if (is.na(ts[1]))
      return(NA)
    else {
      b1<-ts[2]
      b2<-ts[3]
      mn1main<-s1$seg.mean[b1]
      mn1other<-s1$seg.mean[-b1][1]
      mn2main<-s2$seg.mean[b2]
      mn2other<-s2$seg.mean[-b2][1]
      res2<-rep(NA,Nsim)
      ### indep
      resi<-0
      res<-rep(NA,Nsim)
      ct<-0
      while (resi<Nsim & ct<500) {
        ct<-ct+1
        if (ts[4]==1) # there is overlap
        {hotspot<-sample(max(1,s1c[b1-1],s2c[b2-1]):min(s1c[b1],s2c[b2]),1)
        breaks1<-c(sort(sample(hotspot,max(0,b1-1))),
                   sort(sample(c((hotspot+1):chrlen),n1-b1)))
        breaks2<-c(sort(sample(hotspot,max(0,b2-1))),
                   sort(sample(c((hotspot+1):chrlen),n2-b2)))
        }
        else
        {hotspot<-sample(max(1,s1c[b1-1]):min(s1c[b1]),1)
        breaks1<-c(sort(sample(hotspot,max(0,b1-1))),
                   sort(sample(c((hotspot+1):chrlen),n1-b1)))
        hotspot<-sample(max(1,s2c[b2-1]):min(s2c[b2]),1)
        breaks2<-c(sort(sample(hotspot,max(0,b2-1))),
                   sort(sample(c((hotspot+1):chrlen),n2-b2)))
        }
        
        mns1<-rep((mn1other+mn2other)/2,n1)
        mns1[b1]<-(mn1main+mn2main)/2
        
        mns2<-rep((mn1other+mn2other)/2,n2)
        mns2[b2]<-(mn1main+mn2main)/2
        
        x1<-sample(resid1)+rep(mns1,c(breaks1,chrlen)-c(0,breaks1))
        x2<-sample(resid2)+rep(mns2,c(breaks2,chrlen)-c(0,breaks2))
        
        sseg1<-segment1(CNA(x1,rep(1,length(x1)),ss1$data$maploc),segmethod=segmethod,segpar=segpar)
        sseg2<-segment1(CNA(x2,rep(1,length(x2)),ss2$data$maploc),segmethod=segmethod,segpar=segpar)
        
        sseg1$output[,7]<-"Normal"
        sseg2$output[,7]<-"Normal"
        
        sseg1$output[sseg1$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
        sseg1$output[sseg1$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"
        
        sseg2$output[sseg2$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
        sseg2$output[sseg2$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"
        
        names(sseg1$output)[7]<-"state"
        names(sseg2$output)[7]<-"state"
        
        tsp<-func(sseg1,sseg2)
        
        if (!is.na(tsp[1])) 
        {resi<-resi+1
        res[resi]<-tsp[1]
        }
      }
      
      if (sum(!is.na(res))<10) return(NA)
      
      ### clonal
      
      res2<-rep(NA,Nsim)
      resi<-0
      ct<-0
      while (resi<Nsim & ct<500)
      {ct<-ct+1
      
      if (ts[4]==1) hotspot<-sample(max(1,s1c[b1-1],s2c[b2-1]):min(s1c[b1],s2c[b2]),1)
      else          hotspot<-sample(min(s1c[b1],s2c[b2]):max(s1c[b1-1],s2c[b2-1]),1)
      
      if (n1==n2 & b1==b2) {
        breaks1<-breaks2<-c(sort(sample(hotspot,max(0,b1-1))),
                            sort(sample(c((hotspot+1):s1c[n1]),n1-b1)))
        mns1<-mns2<-(s1$seg.mean+s2$seg.mean)/2
      }
      else if (n1<n2) 
      {breaks1<-breaks2<-c(sort(sample(hotspot,max(0,b1-1))),
                           sort(sample(c((hotspot+1):s1c[n1]),n1-b1)))
      if (b1==1) mns1<-mns2<-c((mn1main+mn2main)/2,(mn1other+mn2other)/2)
      else mns1<-mns2<-c((mn1other+mn2other)/2,(mn1main+mn2main)/2)
      }
      else if (n2<n1)
      {breaks1<-breaks2<-c(sort(sample(hotspot,max(0,b2-1))),
                           sort(sample(c((hotspot+1):s2c[n2]),n2-b2)))
      if (b2==1) mns1<-mns2<-c((mn1main+mn2main)/2,(mn1other+mn2other)/2)
      else mns1<-mns2<-c((mn1other+mn2other)/2,(mn1main+mn2main)/2)
      }
      
      lens1<-c(breaks1,chrlen)-c(0,breaks1)
      lens2<-c(breaks2,chrlen)-c(0,breaks2)
      
      x1<-sample(resid1)+rep(mns1,lens1)
      x2<-sample(resid2)+rep(mns2,lens2)
      
      sseg1<-segment1(CNA(x1,rep(1,length(x1)),ss1$data$maploc),segmethod=segmethod,segpar=segpar)
      sseg2<-segment1(CNA(x2,rep(1,length(x2)),ss2$data$maploc),segmethod=segmethod,segpar=segpar)
      
      sseg1$output[,7]<-"Normal"
      sseg2$output[,7]<-"Normal"
      
      sseg1$output[sseg1$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
      sseg1$output[sseg1$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"
      
      sseg2$output[sseg2$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
      sseg2$output[sseg2$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"
      
      names(sseg1$output)[7]<-"state"
      names(sseg2$output)[7]<-"state"
      
      tsp<-func(sseg1,sseg2)
      
      if (!is.na(tsp[1])) 
      {resi<-resi+1
      res2[resi]<-tsp[1]
      }
      }
      if (sum(!is.na(res2))<10) return(NA)
      
      res<-res[!is.na(res)]
      res2<-res2[!is.na(res2)]
      a<-density(res,from=0,to=length(xx1))
      p1<-a$y[sort.list(abs(a$x-ts[1]))[1]]
      a<-density(res2,from=0,to=length(xx1))
      p2<-a$y[sort.list(abs(a$x-ts[1]))[1]]
      pvalue<-mean(res<=ts[1])
      
      return(list(ts[1],pvalue,p1,p2))
    }
  }

func <-
  function(ss1,ss2){complexts(ss1$output,ss2$output,ss1$data$maploc)}

complexts <-
  function(s1,s2,maplocs)
  {
    ###### function calculating test statistic t - closeness between two concordant segments of gain/loss
    ###### input - two one-step CBS segmentations
    ###### maplocs are used to take into account missing values
    
    maplocs<-as.numeric(as.character(maplocs))
    s1$loc.start<-as.numeric(as.character(s1$loc.start))
    s1$loc.end<-as.numeric(as.character(s1$loc.end))
    s2$loc.start<-as.numeric(as.character(s2$loc.start))
    s2$loc.end<-as.numeric(as.character(s2$loc.end))
    
    i1<-match(round(s1$loc.start,2),round(maplocs,2))
    j1<-match(round(s1$loc.end,2),round(maplocs,2))
    
    i2<-match(round(s2$loc.start,2),round(maplocs,2))
    j2<-match(round(s2$loc.end,2),round(maplocs,2))
    
    n1<-nrow(s1)
    n2<-nrow(s2)
    w1<-which(s1[,7]!="Normal")
    w2<-which(s2[,7]!="Normal")
    
    if (n1==1 | n2==1 | length(w1)==0 | length(w2)==0 | all(s1[,7]=="Gain") | 
        all(s1[,7]=="Loss") |  all(s2[,7]=="Gain") | all(s2[,7]=="Loss") ) return(NA)
    else  
      if ((n1==3 & s1[2,7]=="Normal") | (n2==3 & s2[2,7]=="Normal")) return(NA)
    else 
      if (n1==2 & n2==2 & 
          sign(s1$seg.mean[1]-s1$seg.mean[2])!=sign(s2$seg.mean[1]-s2$seg.mean[2])) 
        return(NA)
    else
    {
      if (n1==2 & length(w1)==1) ind1<-w1
      if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])>abs(s1$seg.mean[w1[2]])) 
        ind1<-w1[1]
      if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])<abs(s1$seg.mean[w1[2]])) 
        ind1<-w1[2]
      if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])==abs(s1$seg.mean[w1[2]]) &
          s1$num.mark[w1[1]]>=s1$num.mark[w1[2]]) ind1<-w1[1]
      if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])==abs(s1$seg.mean[w1[2]]) & 
          s1$num.mark[w1[1]]<s1$num.mark[w1[2]]) ind1<-w1[2]
      if (n1==3) ind1<-2
      
      if (n2==2 & length(w2)==1) ind2<-w2
      if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])>abs(s2$seg.mean[w2[2]]))
        ind2<-w2[1]
      if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])<abs(s2$seg.mean[w2[2]])) 
        ind2<-w2[2]
      if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])==abs(s2$seg.mean[w2[2]]) &
          s2$num.mark[w2[1]]>=s2$num.mark[w2[2]]) ind2<-w2[1]
      if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])==abs(s2$seg.mean[w2[2]]) & 
          s2$num.mark[w2[1]]<s2$num.mark[w2[2]]) ind2<-w2[2]
      if (n2==3) ind2<-2
      
      if (n1==2 & length(w1)==2 & n2==2 & length(w2)==2 ) if ( all(s1[,7]==s2[,7])) 
        ind1<-ind2<-1
      if (s1[ind1,7]!=s2[ind2,7] ) return(NA)   # discordnat pairs
      else{
        overl<-abs(i1[ind1]-i2[ind2])+abs(j1[ind1]-j2[ind2])
        nono<-!(i1[ind1]>j2[ind2] | j1[ind1]<i2[ind2])  # is there overlap
        b1<-ind1; b2<-ind2
        return(c(overl,b1,b2,nono))
      }
    }
  }

grantLR <-
  function(tum1,tum2,p,pi,pc,cvalue=0.5,rescale=TRUE,prnfile=NULL,nm=NULL,chrlist)
  {
    ######## calculates likelihood ratio as in the paper
    ######## input is G/L/N status for each chromosome and both tumors, marginal frequencies,density estimates###for individual comparisons
    unlogit<- function(x){exp(x)/(1+exp(x))}
    logit<- function(x){log(x/(1-x))}
    
    
    pi[is.na(pi)]<-1
    pc[is.na(pc)]<-1
    pi<-as.numeric(pi)
    pc<-as.numeric(pc)
    
    tum1<-factor(tum1,levels=c("Gain","Loss","Normal"))
    tum2<-factor(tum2,levels=c("Gain","Loss","Normal"))
    
    nchr<-length(tum1)
    
    PG<-p[,1]
    PL<-p[,2]
    PN<-p[,3]
    
    if (rescale) #rescales the frequencies to conform to number of gains and losses at particular pair of tumors
    {
      currfr<-(table(tum1)+table(tum2))/(2*nchr)
      w<-currfr==0
      if (any(w)) currfr<-(currfr+0.05)/sum(currfr+0.05)
      
      if (currfr[1]!=0.5) PGp<-unlogit(logit(PG)*logit( currfr[1])/( mean(logit(PG)))) else PGp<-PG* currfr[1]/ mean(PG)
      if (currfr[2]!=0.5) PLp<-unlogit(logit(PL)*logit(currfr[2])/( mean(logit(PL)))) else PLp<-PL* currfr[2]/ mean(PL)
      
      
      w<-PLp+PGp
      PLp[w>1]<-PLp[w>1]/(1.05*w[w>1])
      PGp[w>1]<-PGp[w>1]/(1.05*w[w>1])
      PNp<-1-PGp-PLp
      PG<-PGp
      PL<-PLp
      PN<-pmin(1,PNp)
      
    }
    
    
    if (any(round(PG+PL+PN,8)!=1 | round(PG,8)<=0 |round(PL,8)<=0 |round(PN,8)<=0))
      stop("trouble with p's")
    
    pc[tum1!=tum2 | tum1=="Normal" | tum2=="Normal"]<-1
    pi[tum1!=tum2 | tum1=="Normal" | tum2=="Normal"]<-1
    
    
    Rgg<-(tum1==tum2 & tum1=="Gain")
    Rll<-(tum1==tum2 & tum1=="Loss")
    Rnn<-(tum1==tum2 & tum1=="Normal")
    Rgl<-(tum1=="Loss" & tum2=="Gain") | (tum2=="Loss" & tum1=="Gain")
    Rgn<-(tum1=="Normal" & tum2=="Gain") | (tum2=="Normal" & tum1=="Gain")
    Rln<-(tum1=="Loss" & tum2=="Normal") | (tum2=="Loss" & tum1=="Normal")
    
    
    cc<-0
    LI<-(cc*PG+(1-cc)^2*PG^2/(1-cc*PG-cc*PL))^Rgg  * 
      (cc*PL+(1-cc)^2*PL^2/(1-cc*PG-cc*PL))^Rll  * 
      (2*(1-cc)^2*PG*PL/(1-cc*PG-cc*PL))^Rgl  *
      (2*(1-cc)*PG*PN/(1-cc*PG-cc*PL))^Rgn *  
      (2*(1-cc)*PL*PN/(1-cc*PG-cc*PL))^Rln * (PN^2/(1-cc*PG-cc*PL))^Rnn
    
    
    cc<-cvalue
    LC<-(cc*PG+(1-cc)^2*PG^2/(1-cc*PG-cc*PL))^Rgg  * 
      (cc*PL+(1-cc)^2*PL^2/(1-cc*PG-cc*PL))^Rll  * 
      (2*(1-cc)^2*PG*PL/(1-cc*PG-cc*PL))^Rgl  *
      (2*(1-cc)*PG*PN/(1-cc*PG-cc*PL))^Rgn *  (2*(1-cc)*PL*PN/(1-cc*PG-cc*PL))^Rln * 
      (PN^2/(1-cc*PG-cc*PL))^Rnn
    
    Bg<-cc*PG/(cc*PG+(1-cc)^2*PG^2/(1-cc*PG-cc*PL))
    Bl<-cc*PL/(cc*PL+(1-cc)^2*PL^2/(1-cc*PG-cc*PL))
    
    Bg[tum1!="Gain" |tum2!="Gain" ]<-0.5
    Bl[tum1!="Loss" |tum2!="Loss" ]<-0.5
    
    
    BBg<-Bg*pc+(1-Bg)*pi
    BBl<-Bl*pc+(1-Bl)*pi
    BBg[tum1!="Gain" |tum2!="Gain" ]<-1
    BBl[tum1!="Loss" |tum2!="Loss" ]<-1
    
    a<-(BBg*BBl/pi)[pi!=1]
    names(a)<-chrlist[pi!=1]
    
    if (!is.null(prnfile)) 
    {write.table(cbind(nm,chrlist,PG,PL,PN,tum1,tum2,LI, LC, BBg*BBl/pi),
                 file=prnfile,append=TRUE)
    }
    return(c(prod(LC/LI),prod(LC*BBg*BBl/(LI*pi)),a))
    
  }