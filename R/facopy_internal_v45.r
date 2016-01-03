if(getRversion() >= "2.15.1") 
  utils::globalVariables(c("bp_len","code","chr_q_arm","type","altered","xmin","xmax",
                           "pos_st","pos_en","freq","value","pos","ymin","ymax",
                           "Frequency","Chromosome","facopy_msigdb","facopy_msigdbNames",
                           "GOBPPARENTS","GOCCPARENTS","GOMFPARENTS"))

#######################################################################
############################################### INTERNAL CODIFICATION #
#######################################################################
.printList = function(ttext, x, header="aliases") {
  cat(ttext,"\n",sep=":")
  for (i in seq(length(x))) {
    if (!is.null(names(x)[[i]]))
      cat(" - ", names(x)[[i]][1])
    else
      cat(" - ", x[[i]][1])
    if (length(x[[i]])>1)
      cat(" (",header,": ", paste(x[[i]],collapse=", "), ")", sep="")
    cat("\n")
  }
}
.getFacopyAnnot = function() {
  dat = data(package="facopy.annot")$results[,3]
  dat = sapply(paste(.genomes(), rep(c("_db","_feature"),each=length(.genomes())), sep=""), grep, dat, value=TRUE)
  dat = lapply(dat, function(i) gsub("^[a-zA-Z0-9]+_[a-zA-Z0-9]+_","",i))
  names(dat) = sapply(names(dat), function(i) gsub("_"," ",i))
  dat
}

getFacopyInfo = function() {
  dat = .getFacopyAnnot()
  .printList("Variable types", .types())
  .printList("Alteration combinations", .alterations())
  .printList("Input data", .methods())
  .printList("Genome builds", .genomes())
  .printList("Genomic features", dat[grep("feature",names(dat))], "bundled")
  .printList("External data sets", dat[grep("db",names(dat))], "available")
}

.genomes = function() c("hg18","hg19","mm8")
.features = function() c("ensembl","cancergene","oncogene","tumorsupressor","lincRNA","mirnas")
.allChrs = function(genome, sex) {
  if (genome%in%c("hg18","hg19")) {
    out = paste(rep(c(1:22,sex),each=2), c("p","q"), sep="")
    out[!out%in%c("13p","14p","15p","22p")]
  } else if (genome%in%c("mm8")) {
    paste(c(1:19,sex), "q", sep="")
  }
}
.methods = function() c("seqcna", "cnanorm", "patchwork", "freec", "oncosnp", "oncosnp-seq", "gap", "exomecnv", "titan", "fastcall", "excavator")
.alterCodes = function() c("1","3","11","13","12") #1:loss, 3:gain; x:noLOH, 1x:LOH
.allCodes = function() c("1","2","3","11","12","13") # in that specific order
.alterations = function() list("amplifications", "deletions", "loh", "cnas",
                               c("any","all"), "onlygain", "someloss")
.getAlterCodes = function(alteration) {
  alterCodes = list(c("3","13"), c("1","11"), c("11","13","12"), c("3","13","1","11"), 
                    c("3","13","1","11","12"), c("3"), c("1","11","12","13"))
  wh = which(!is.na(sapply(.alterations(), function(i) pmatch(tolower(alteration),i))))[1]
  if (is.na(wh))
    stop("Unknown alteration name: ", alteration)
  alterCodes[[wh]]
}
.types = function() list(categorical=c("categorical","discrete","qualitative","enumerative"),
                         quantitative=c("quantitative","continuous"))
.getTypeCode = function(type) {
  wh = which(!is.na(sapply(.types(), function(i) pmatch(tolower(type),i))))[1]
  if (is.na(wh))
    stop("Unknown variable type: ", type)
  names(.types())[[wh]]
}
.isCateg = function(fad, x) fad@vData$varTypes[x]=="categorical"
.n = function(fad, vn=NULL, vv=NULL) {
  if (!is.null(vn)) sel = fad@vData$varTable[,vn]%in%vv
  else sel = TRUE
  length(fad@vData$varTable$code[sel])
}
.unLOH = function(x) { x[x>=10] = x[x>=10] - 10; x }
.getLOH = function(x) floor(x/10)
.ord = function(x, arms, ...) order(match(x, arms), ...)
.refcols = function(n) {
	cc = c(colorRampPalette(c("#F8CA00","#bd1550"))(3),"#e00000","#926239", colorRampPalette(c("#25e361","#59BCE9","#5123c1"))(4), "#ff00cc")
	cc = cc[c(3:5,2:1,6:10)]
	cc = cc[c(2:10,1)]
	colorRampPalette(cc)(n+1)[1:n]
}

#######################################################################
########################### READING OUTPUT FROM CNA-DETECTION METHODS #
#######################################################################

.readFromInt = function(samples, samplenames, f, fargs, FUN=NULL, ...) {
  own = c()
  for (i in 1:length(samples)) {
    data = f(read.delim(file=samples[i], ...), fargs)
    own = rbind(own, cbind(rep(samplenames[i],nrow(data)), data))
    message(paste("(",i,"/",length(samples),") ",nrow(data)," calls in sample '",samplenames[i],"'",sep=""))
  }
  own[,1] = as.vector(as.character(own[,1]))
  if (!is.null(FUN)) own[,1] = sapply(own[,1], FUN)
  own
}
.pack = function(data, sex) {
  chrs = rle(data$chr)$values
  chrs = c(chrs[chrs!=c("X","Y")], sex)
  out = c()
  for (chr in chrs) {
    xx = data[data$chr==chr,]
    en = cumsum(rle(xx$type)$lengths)
    if (length(en)>1) st = c(1, en[seq(length(en)-1)]+1)
    else st = 1
    if (length(en)>0) {
      xx2 = xx[st,]
      xx2$pos_en = xx$pos_en[en]
      xx2$bp_len = sapply(seq(st), function(i) sum(xx$bp_len[st[i]:en[i]]+1))
      out = rbind(out, xx2)
    }
  }
  out
}
.makeInfo = function(cData, sex, pack=FALSE) {
  if (!is.null(sex)) if (sex=="") sex = NULL
  out = c()
  en = cumsum(rle(cData$code)$lengths)
  st = c(1, en[seq(length(en)-1)]+1)
  if (pack) for (s in seq(st)) out = rbind(out, .pack(cData[st[s]:en[s],], sex))
  if (!pack) { # sex chrs selection has already been done if 'pack', no need to redo this
    chrs = rle(cData$chr)$values
    chrs = c(chrs[!chrs%in%c("X","Y","XY","chrX","chrY","chrXY")], sex)
    out = cData[cData$chr%in%chrs,]
  }
  new("facopyInfo", cData=out, sex=sex)
}

readCNData = function(folder, method=NULL, sex=c("X"), FUN=NULL, 
                      version=NULL, window=NULL, rankThr=NULL,
                      pfbFilename, lengthThr=10, ...) {
  tryFUN = function(old, new) {
    if (is.null(old)) {
      message("Note:\nFunction applied to filenames will be: ",as.character(enquote(new)[2]))
      new
    } else { old }
  }
  methodArg = .methods()[pmatch(tolower(method),.methods())]
  if (length(methodArg)==0) {
    .readFromFolder(folder, sex, FUN, ...)
  } else if (is.na(methodArg)) {
    stop("Output data from specified copy number calling method not supported.") 
  } else if (methodArg=="seqcna") {
    .readFromSEQCNA(folder, sex, version, FUN)
  } else if (methodArg=="cnanorm") {
    if (is.null(window)) stop("Window size is required.")
    .readFromCNANORM(folder, sex, version, window, FUN)
  } else if (methodArg=="patchwork") {
    FUN = tryFUN(FUN, function(x)gsub("_Copynumbers.csv$","",x))
    .readFromPATCHWORK(folder, sex, version, FUN)
  } else if (methodArg=="freec") {
    FUN = tryFUN(FUN, function(x)gsub("_ratio.txt$","",x))
    .readFromFREEC(folder, sex, version, FUN)
  } else if (methodArg=="oncosnp") {
    if (is.null(version)) stop("Version is required for OncoSNP.")
    FUN = tryFUN(FUN, function(x)gsub(".cnvs$","",x))
    .readFromONCOSNP(folder, sex, version, rankThr, FUN)
  } else if (methodArg=="oncosnp-seq") {
    FUN = tryFUN(FUN, function(x)gsub(".cnvs$","",x))
    .readFromONCOSNPSEQ(folder, sex, version, rankThr, FUN)
  } else if (methodArg=="gap") {
    if (is.null(pfbFilename)) stop("PFB filename is required.")
    .readFromGAP(folder, pfbFilename, sex, lengthThr, FUN)
  } else if (methodArg=="exomecnv") {
    FUN = tryFUN(FUN, function(x)gsub(".segment.copynumber.txt$","",x))
    .readFromEXOMECNV(folder, sex, version, FUN)
  } else if (methodArg=="titan") {
    FUN = tryFUN(FUN, function(x)gsub(".txt$","",x))
    .readFromTITAN(folder, sex, version, FUN)
  } else if (methodArg%in%c("fastcall","excavator")) {
    FUN = tryFUN(function(x)gsub("^FastCallResults_([^.]*).*", "\\1", x))
    .readFromFASTCALL(folder, sex, version, FUN)
  } else {
    stop("Output data from specified copy number calling method not supported.") 
  }
}

.readFromFolder = function(folder, sex=c("X"), FUN=NULL, ...) {
  own = .readFromInt(list.files(folder,full.names=TRUE), list.files(folder,full.names=FALSE), function(x,fargs)I(x), NULL, FUN, ...)
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "bp_len", "type")
  .makeInfo(own, sex)
}

.readFromSEQCNA = function(folder, sex=c("X"), version=NULL, FUN=NULL) {
  f = function(data, fargs) { # the output is NOT allele-specific
    data = data[,c("chrom","win.start","CN")]
    data = data[!is.na(data$CN),]
    data$CN[data$CN<1] = 1
    data$CN[data$CN>3] = 3
    l = diff(data$win.start[1:2])-1
    cbind(data, data$win.start+l, rep(l,nrow(data)))
  }
  getSample = function(x, f) grep("seqCNA_out.txt", list.files(x, full.names=f, recursive=TRUE), value=TRUE)
  own = .readFromInt(getSample(folder,TRUE), dirname(getSample(folder,FALSE)), f, NULL, FUN, header=TRUE, sep="\t")
  own[,2] = as.numeric(gsub("chr","",as.character(own[,2])))
  own = own[,c(1,2,3,5,6,4)]
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "bp_len", "type")
  .makeInfo(own, sex, TRUE)
}

.readFromCNANORM = function(folder, sex=c("X"), version=NULL, window, FUN=NULL) {
  f = function(data, fargs) { # the output is NOT allele-specific
    data$Cn = round(data$SegMean.n)
    data$Cn[data$Cn<1] = 1
    data$Cn[data$Cn>3] = 3
    l = fargs$w
    data$Start = data$Pos*(l-1)+1
    cbind(data[,c("Chr","Start","Cn")], data$Start+l-1, rep(l,nrow(data)))
  }
  own = .readFromInt(list.files(folder,full.names=TRUE), list.files(folder,full.names=FALSE), f, list(w=window), FUN, header=TRUE, sep="\t")
  own = own[,c(1,2,3,5,6,4)]
  own[,2] = gsub("chr","",as.character(own[,2]))
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "bp_len", "type")
  .makeInfo(own, sex, TRUE)
}

.readFromPATCHWORK = function(folder, sex=c("X"), version=NULL, FUN=function(x)gsub("_Copynumbers.csv$","",x)) {
  f = function(data, fargs) { # the output CAN be allele-specific
    mCn = data[,"mCn"]
    data = data[,c("chr","start","end","Cn")]
    data$Cn[data$Cn<1] = 1
    data$Cn[data$Cn>3] = 3
    if (! all(mCn==0)) # means that allele-specific was used
      data$Cn = data$Cn + 10*as.integer(mCn==0)
    cbind(data, data$end-data$start+1)
  }
  getSample = function(folder,f) grep(".*csv$", list.files(folder,full.names=f), value=TRUE)
  own = .readFromInt(getSample(folder,TRUE), basename(getSample(folder,FALSE)), f, NULL, FUN, header=TRUE, sep=",")
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "type", "bp_len")
  own[,2] = gsub("chr","",as.character(own[,2]))
  .makeInfo(own, sex)
}

.readFromFREEC = function(folder, sex=c("X"), version=NULL, FUN=function(x)gsub("_ratio.txt$","",x)) {
  f = function(data, fargs) { # the output CAN be allele-specific
    hasGeno = FALSE
    if ("Genotype"%in%colnames(data)) { # not verified!!!
      geno = data[,"Genotype"]
      geno[is.na(geno)] = "NA"
      isHom = sapply(strsplit(geno,""), function(i)all(i==i[1]))
      hasGeno = TRUE
    }
    data = data[,c("Chromosome","Start","CopyNumber")]
    data[data[,3]<1,3] = 1
    data[data[,3]>3,3] = 3
    if (hasGeno) data[,3] = data[,3] + as.integer(isHom)
    l = diff(data$Start[1:2])-1
    cbind(data, data$Start+l, rep(l,nrow(data)))
  }
  getSample = function(folder,f) grep(".*_ratio.txt$", list.files(folder,full.names=f), value=TRUE)
  own = .readFromInt(getSample(folder,TRUE), basename(getSample(folder,FALSE)), f, NULL, FUN, header=TRUE, sep="\t")
  own = own[,c(1,2,3,5,6,4)]
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "bp_len", "type")
  .makeInfo(own, sex, TRUE)
}

.readFromTITAN = function(folder, sex=c("X"), version=NULL, FUN=function(x)gsub(".txt$","",x)) {
  f = function(data, fargs) { # the output IS allele-specific
    data = data[,c("Chr","Position","TITANstate")]
    data[data[,3]%in%c(0,1,2),3] = 1
    data[data[,3]%in%c(3,5),3] = 12
    data[data[,3]%in%c(4),3] = 2
    data[data[,3]%in%c(6,9),3] = 13
    data[data[,3]%in%c(7,8),3] = 3
    data[data[,3]%in%c(10,14),3] = 13
    data[data[,3]%in%c(11,12,13),3] = 3
    data[data[,3]%in%c(15,20),3] = 13
    data[data[,3]%in%c(16,17,18,19),3] = 3
    data[data[,3]%in%c(21:100),3] = 3
    data$Start = c(1, data$Position[seq(length(data$Position)-1)]+1)
    cbind(data, data$Position-data$Start+1)
  }
  own = .readFromInt(list.files(folder,full.names=TRUE), list.files(folder,full.names=FALSE), f, NULL, FUN, header=TRUE, sep="\t")
  own = own[,c(1,2,5,3,6,4)]
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "bp_len", "type")
  .makeInfo(own, sex)
}

.readFromONCOSNPSEQ = function(folder, sex=c("X"), version=1.04, rankThr=NULL, FUN=function(x)gsub(".cnvs$","",x)) {
  .readFromONCOSNPInt(folder, sex, version, rankThr, FUN, SEQ=TRUE)
}
.readFromONCOSNP = function(folder, sex=c("X"), version=1.3, rankThr=NULL, FUN=function(x)gsub(".cnvs$","",x)) {
  .readFromONCOSNPInt(folder, sex, version, rankThr, FUN, SEQ=FALSE)
}
.readFromONCOSNPInt = function(folder, sex, version, rankThr, FUN, SEQ) {
  f = function(data, fargs) { # the output IS allele-specific
    if (is.null(data$CopyNumber))
      data$CopyNumber = data$Copy.Number
    data$CopyNumber[data$CopyNumber<1] = 1
    data$CopyNumber[data$CopyNumber>3] = 3
    data$LOH[data$LOH==2] = 0
    data$CopyNumber = data$CopyNumber + 10 * data$LOH
    if (!is.null(data$Rank) && !is.null(fargs$rankThr))
      data = data[data$Rank<=fargs$rankThr,]
    if (!is.null(data$Rank)) cbind(data[,1:4], data$Rank)
    else data[,1:4]
  }
  getSample = function(folder,f) grep("*.cnvs$", list.files(folder,full.names=f), value=TRUE)
  if (SEQ==FALSE && version<1.3 && !is.null(rankThr)) {
    rankThr = NULL
    warning("Ranks are only available in OncoSNP from version 1.3 on.")
  }
  own = .readFromInt(getSample(folder,TRUE), getSample(folder,FALSE), f, list(rankThr=rankThr), FUN, header=TRUE,sep="\t",fill=TRUE)
  own = cbind(own, apply(own, 1, function(i) as.numeric(as.vector(i[4]))-as.numeric(as.vector(i[3]))))
  colnames(own)[1:6] = c("code", "chr", "pos_st", "pos_en", "type", "bp_len")
  if (version>=1.3 && ncol(own)>6) colnames(own)[7] = c("conf")
  .makeInfo(own, sex)
}

.readFromGAP = function(folder, pfbFilename, sex=c("X"), lengthThr=10, FUN=NULL) {
  f = function(data, fargs) { # the output IS allele-specific
    data = data[data$Len>=fargs$lengthThr & !is.na(data$CN) & !is.na(data$Chr), ] 
    # if CN = major allele count and there is no germline LOH (LOH column), there is somatic LOH
    isLOH = 10 * (round(data$CN)==round(data$BA) & !data$LOH)
    data$CN[data$CN<1] = 1
    data$CN[data$CN>3] = 3
    data$CN = round(data$CN) + isLOH
    data = cbind(pfb[data$Ind,c(1,3)], pfb[data$Ind_K,c(1,3)], data)
    data[, c(8,16,1,3,5,6,2,4,18)]
  }
  pfb = read.table(pfbFilename,header=TRUE,sep="\t")
  samples = list.files(folder, full.names=TRUE)
  samplenames = list.files(folder, full.names=FALSE)
  own = .readFromInt(samples, samplenames, f, list(lengthThr=lengthThr), FUN, header=TRUE,sep=" ",fill=TRUE)
  own[,2] = floor(as.numeric(as.character(own[,2]))) # may be problematic if sexual chromosomes used, but floor fun. is necessary
  own = cbind(own, apply(own, 1, function(i) as.numeric(as.vector(i[9]))-as.numeric(as.vector(i[8]))))
  colnames(own) = c("code", "chr", "type", "snp_st", "snp_en", "ind_st", "ind_en", "pos_st", "pos_en", "conf", "bp_len")
  .makeInfo(own, sex)
}

.readFromEXOMECNV = function(folder, sex=c("X"), version=NULL, FUN=function(x)gsub(".segment.copynumber.txt$","",x)) {
  f = function(data, fargs) { # the output is NOT allele-specific
    data[,4][data[,4]<1] = 1
    data[,4][data[,4]>3] = 3
    cbind(data, data[,3]-data[,2]+1)
  }
  getSample = function(folder,f) grep(".segment.copynumber.txt$", list.files(folder,full.names=f), value=TRUE)
  own = .readFromInt(getSample(folder,TRUE), basename(getSample(folder,FALSE)), f, NULL, FUN, header=FALSE, sep="\t")
  own[,2] = gsub("chr","",as.character(own[,2]))
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "type", "bp_len")
  .makeInfo(own, sex)
}

.readFromFASTCALL = function(folder, sex=c("X"), version=NULL, FUN=function(x)gsub("^FastCallResults_([^.]*).*", "\\1", x)) {
  f = function(data, fargs) { # the output is NOT allele-specific
    data = data[,c("Chromosome","Start","End","CN")]
    data = data[!is.na(data$CN),]
    data$CN[data$CN<1] = 1
    data$CN[data$CN>3] = 3
    cbind(data, data[,3]-data[,2]+1)
  }
  getSample = function(x, f) grep("^FastCallResults_", list.files(x, full.names=f, recursive=TRUE), value=TRUE)
  own = .readFromInt(getSample(folder,TRUE), dirname(getSample(folder,FALSE)), f, NULL, FUN, header=TRUE, sep="\t")
  own[,2] = gsub("chr","",as.character(own[,2]))
  colnames(own) = c("code", "chr", "pos_st", "pos_en", "bp_len", "type")
  .makeInfo(own, sex)
}

#######################################################################
########################################## GENERIC INTERNAL FUNCTIONS #
#######################################################################

# necessary because aggregation in data.table package does not fill with not present factor levels
# col: annotation columns, thus 1 to 1 with respect to sample code
# ...: columns from which to produce combinations with themselves and code
.fill = function(varTable, aggr, col, emptyVal, ...) {
  eg = do.call(expand.grid,list(code=varTable$code, ...))
  st = ncol(eg)+1
  if (!is.null(col)) {
    eg = merge(eg, varTable[,c("code",col),drop=FALSE])
    colnames(eg)[st:ncol(eg)] = colnames(aggr)[st:ncol(eg)]
  }
  out = merge(aggr, eg, all=TRUE)
  out[is.na(out)] = emptyVal
  out
}

.calcFreqsInt = function(sel_features, tabs, cn, nums, vari) {
  freqs = c()
  for (i in tabs) {
    tt = i[i[,vari]%in%nums,]
    tab = table(factor(tt[,vari],levels=nums), factor(tt$type,levels=.allCodes()))
    rs = rowSums(tab[,colnames(tab)%in%cn,drop=FALSE])
    rsAll = rowSums(tab[,,drop=FALSE])
    freqs = rbind(freqs, rs / rsAll)
  }
  if (! is.null(dim(freqs))) {
    freqs[is.na(freqs)] = 0
    rownames(freqs) = sel_features$feature
    colnames(freqs)[1:length(nums)] = nums
  }
  data.frame(freqs, check.names=FALSE)
}

.makeArmsInt = function(cData, chrs, genome, varColumn=NULL) {
  if (genome%in%c("hg18","hg19")) {
	sepsText = paste(genome,"_armLimits",sep="")
	data(list=sepsText)
	seps = eval(parse(text=sepsText))
    seps = seps$limit[sapply(seps$chr_q_arm, function(i) length(grep("p",i))==1)]
  } else if (genome%in%c("mm8")) {
    seps = rep(0,length(chrs))
  }
  x = cData[,c("code","type","chr","pos_st","bp_len",varColumn)]
  out = c()
  for (i in 1:length(chrs)) {
    chr = chrs[i]
    xsub = x[x$chr==chr,]
    wh = xsub$pos_st >= seps[i]
    if (length(which(wh==FALSE)) > 0) out = rbind(out, cbind(xsub[!wh,!colnames(xsub)%in%c("chr","pos_st")], chr_q_arm=paste(chr,"p",sep="")))
    if (length(which(wh==TRUE)) > 0) out = rbind(out, cbind(xsub[wh,!colnames(xsub)%in%c("chr","pos_st")], chr_q_arm=paste(chr,"q",sep="")))
  }
  out
}

.aggregateInt = function(x, n, genome, arms, FUN, levels, varTable) {
  x = x[x$chr_q_arm%in%arms,] # disregard data in arms with no features
  aggr = data.frame(data.table(x)[,sum(bp_len), by=list(code, chr_q_arm, match.fun(FUN)(type))])
  aggr = .fill(varTable, aggr, NULL, 0, chr_q_arm=arms, match.fun=levels)
  colnames(aggr) = c("code","chr_q_arm", "type", "altered")
  #thanks to 'n', .fill not necessary
  aggr = data.frame(data.table(aggr)[,sum(altered)/n, by=list(chr_q_arm, type)])
  colnames(aggr) = c("chr_q_arm", "type", "altered")
  armlimitsText = paste(genome,"_armLimits",sep="")
  data(list=armlimitsText)
  armlimits = eval(parse(text=armlimitsText))
  lens = matrix(armlimits$limit, 2)
  lens[2,] = lens[2,]-lens[1,]
  armlengths = data.frame(chr_q_arm=as.character(armlimits$chr_q_arm), length=as.vector(lens))
  aggr = merge(data.frame(aggr), armlengths)
  aggr$freq = as.numeric(as.character(aggr$altered)) / as.numeric(as.character(aggr$length))
  # fill missing rows due to lack of alterations in specific arms
  tab = table(aggr$chr_q_arm, aggr$type)
  wh = which(tab==0, arr.ind=TRUE)
  if (nrow(wh) > 0)
    for (i in 1:nrow(wh)) {
      arm = rownames(tab)[wh[i,1]]
      len = armlengths$length[armlengths$chr_q_arm==arm]
      aggr = rbind(aggr, c(arm, colnames(tab)[wh[i,2]], 0, len, 0))
    }
  aggr[.ord(aggr$chr_q_arm, arms),]
}

#######################################################################
############################## ANNOTATION WITH VARIABLES AND FEATURES #
#######################################################################

# allows to give names to variables and their names
# allows to use only subsets of variable values
.makeVars = function(varTypes, varNames, varColumns, varValuesNames, varValues) {
  varTypes = sapply(varTypes, .getTypeCode)
  if (! length(varTypes)==length(varNames))
    stop("The lengths of the vectors for the variables names and their types do no match.")
	if (! length(varNames)==length(varColumns))
		stop("The lengths of the vectors for the variables names and their column names do no match.")
	if (! length(varNames)==length(varValues))
			stop("The lengths of the vector for the variable names and the list for their values do no match.")
	if (! identical(sapply(varValues,length), sapply(varValuesNames,length)))
		stop("The lengths of the lists for the values and their names do no match.")
	names(varTypes) = names(varColumns) = names(varValues) = names(varValuesNames) = varNames
	list(varTypes=varTypes, varNames=varNames, varColumns=varColumns, varValuesNames=varValuesNames, varValues=varValues)
}

.checkData = function(fad) {
	mandatoryColumns = c("code","chr","type","pos_st","pos_en","bp_len")
	wh = which(! mandatoryColumns%in%colnames(fad@cData))
	if (length(wh) >0)
		stop(paste("The following columns are required in the calls table: ",
				   paste(mandatoryColumns[wh],collapse=", "),".",sep=""))
	if (! all(fad@vData$varColumns%in%colnames(fad@cData)))
		stop("Not all variable column names are present in the calls table.")
	for (i in 1:length(fad@vData$varColumns)) {
		inter = intersect(fad@vData$varValues[[i]], unique(fad@cData[,fad@vData$varColumns[i]]))
		if (length(inter) == 0)
			stop(paste("Variable '",fad@vData$varNames[i],"': its defined values and those in the calls table do not match at all.",sep=""))
	}
}

# variable columns can be taken from file; all but 'code' will be taken
# variable values can be taken from file; values can be taken for one or more variables
# variable names can be taken from variable column names; all or none may be taken
# variable value names can be taken from variable values; can be taken for one or more variables
addVariables = function(fad, varInfo, varTypes, varColumns=NULL, varValues=NULL, varColumnsNames=NULL, varValuesNames=NULL, ...) {
	if (class(varInfo)=="character") vv = read.table(varInfo, header=TRUE, ...)
  else vv = varInfo
	if (! "code"%in%colnames(vv))
		stop("The necessary column 'code' is not present in the variables file.")
	cc = merge(fad@cData, cbind(vv, all=1), by="code")
	fad@cData = cc
  #
	if (is.null(varColumns)) {
	  if (!is.null(varValues))
	    stop("Variables should be specified with 'varColumns' if their values are given with 'varValues'.")
	  varColumns = colnames(vv)[colnames(vv)!="code"]
	}
  .autoComplete = function(x, y, f, ...) {
    wnull = function(o) sapply(o, is.null)
    if (is.null(x)) { x = f(y,...) } else if (any(wnull(x))) { x[wnull(x)] = f(y[wnull(x)],...) }
    x
  }
  f = function(y, vv) {
    splitc = function(x) lapply(seq_len(ncol(x)), function(i) x[,i])
    lapply(splitc(vv[,y,drop=FALSE]), function(i) names(table(i)))
  }
	varValues = .autoComplete(varValues, varColumns, f, vv)
	varColumnsNames = .autoComplete(varColumnsNames, varColumns, I)
	varValuesNames = .autoComplete(varValuesNames, varValues, I)  
  varData = .makeVars(c(varTypes,"categorical"), c(varColumnsNames,"All"), c(varColumns,"all"), 
                      c(varValuesNames,"yes"), c(varValues,1))
	fad@vData = varData
	fad@vData$varTable = cbind(vv, all=1)
	.checkData(fad)
	fad
}

addFeatures = function(fad, what=c("ensembl","cancergene","oncogene","tumorsupressor","lincRNA","mirnas")[1], 
                       genome=c("hg18","hg19","mm8")[1], lMargin=0, rMargin=0, minoverlap=1, data=NULL) {
	if (! genome%in%.genomes())
    stop("Currently, only human genome builds (hg18, hg19) and mouse mm8 build are supported.")
	if (genome%in%c("hg18","hg19")) fad@chrs = c(as.character(1:22),fad@sex)
  if (genome%in%c("mm8")) fad@chrs = c(as.character(1:19),fad@sex)
  if (is.null(data)) {
	xText = paste(genome,"_feature_",what,sep="")
	data(list=xText)
	x = eval(parse(text=xText))
  } else {
    if (class(data)=="character") {
      x = read.table(data, header=TRUE,sep="\t")
    } else if (class(data)=="data.frame") {
      reqCols = c("chr","bp_st","bp_en","feature","chr_q_arm")
      if (! all(reqCols%in%colnames(data)))
        stop("At least the following columns should be present in the external data: ", paste(reqCols,collapse=", "))
      data$chr = gsub("chr","",tolower(data$chr))
      data$chr_q_arm = gsub("chr","",tolower(data$chr_q_arm))
      data$feature = toupper(data$feature)
      x = data
    } else {
      stop("External data with information on genomic features should be a data.frame.")
    }
  }
  x = x[x$chr%in%fad@chrs,]
	if (lMargin > 0) x$bp_st = max(1, x$bp_st-lMargin)
	if (rMargin > 0) x$bp_en = x$bp_en + rMargin
	fad@features = x
	fad@genome = genome
	fad@what = what
	fad@arms = unique(as.character(x$chr_q_arm))
	fad = .calcCallFreqsInt(fad, minoverlap=minoverlap)
	fad
}

# features are assigned the status of the first alteration they overlap, where alterations are ordered in
# genomic order or, if a criterion is available (e.g. OncoSNP version>1.3, GAP), based on such criterion
.calcCallFreqsInt = function(fad, minoverlap=1) {
  # compute sample alteration status for each feature
  if (minoverlap < 1)
    stop("Minimum overlap should be a positive integer number.")
  tabs = list()
  for (chr in fad@chrs) {
    if ("conf" %in% colnames(fad@cData)) {
      cD_chr = fad@cData[fad@cData$chr==chr,c("pos_en","pos_st","code","type","conf",fad@vData$varColumns)]
      cD_chr = cD_chr[order(cD_chr$conf),]
    } else {
      cD_chr = fad@cData[fad@cData$chr==chr,c("pos_en","pos_st","code","type",fad@vData$varColumns)]
    }
    chr_features = fad@features[fad@features$chr==chr,]
    a = IRanges(start=cD_chr$pos_st, end=cD_chr$pos_en)
    b = IRanges(start=chr_features$bp_st, end=chr_features$bp_en)
    overlaps = findOverlaps(a,b, minoverlap=minoverlap)
    for (i in 1:nrow(chr_features)) {
      overTab = cD_chr[overlaps@queryHits[overlaps@subjectHits==i],]
      overTab = overTab[!duplicated(overTab$code),]
      tabs[[length(tabs)+1]] = overTab[,c("code","type",fad@vData$varColumns)]
    }
    message(paste("Alteration status for chr",chr,": done.", sep=""))
  }
  names(tabs) = fad@features$feature
  # compute alteration frequencies per feature by sample and variable
  freqs = list()
  for (i in 1:length(fad@vData$varNames)) {
    freqs[[fad@vData$varNames[i]]] = list()
    for (cn in .alterCodes()) {
      if (i%in%which(fad@vData$varTypes=="categorical")) {
        res = .calcFreqsInt(fad@features, tabs, cn, fad@vData$varValues[[i]], fad@vData$varColumns[i])
        freqs[[fad@vData$varNames[i]]][[cn]] = res
      }
    }
    message(paste(fad@vData$varNames[i], "variable processing: done."))
  }
  fad@fTable = freqs
  fad@gTable = tabs
  fad
}

#######################################################################
############################################### OUTPUT SUMMARY TABLES #
#######################################################################
	
preview = function(fad, folder=NULL) {
	if (is.null(folder))
		folder="."
	s1 = variableSummary(fad, file.path(folder,"variable_summary.txt"))
	s2 = alterationSummary(fad, file.path(folder,"alteration_summary.txt"))
	s3 = variableCor(fad, file.path(folder,"variable_correlations.txt"))
	list(byVar=s1, byAlt=s2, varCor=s3)
}

variableSummary = function(fad, filename=NULL) {
	getName = function(i, x) fad@vData$varValuesNames[i][[1]][which(fad@vData$varValues[i][[1]]==x)]
  ff = function(i, p, a, code) unlist(lapply(p[[1]], function(j) c(getName(i,j), round(mean(a[a[,2]==code&a[,3]==j, 4])))))
	wilcoxtest = function(x) signif(wilcox.test(x[,2] ~ as.integer(x[,1]))$p.value, 3)
	cortest = function(x) signif(unlist(cor.test(x[,2], x[,1], method="kendall")[c("p.value","estimate")]), 3)
	bpsTab = c()
  for (i in fad@vData$varNames[-length(fad@vData$varNames)]) {
    vari_col = fad@vData$varColumns[i]
    aggr = data.frame(data.table(fad@cData)[,sum(bp_len), by=list(code,.unLOH=.unLOH(type),fad@cData[,vari_col])])
    aggr = .fill(fad@vData$varTable, aggr, vari_col, 0, .unLOH=1:3)
    aggrLOH = data.frame(data.table(fad@cData)[,sum(bp_len), by=list(code,.getLOH=.getLOH(type),fad@cData[,vari_col])])
    aggrLOH = .fill(fad@vData$varTable, aggrLOH, vari_col, 0, .getLOH=0:1)
    if (.isCateg(fad,i)) {
      l = apply(combn(fad@vData$varValues[[i]],2), 2, list)
      for (p in l) {
        for (cn in c(1,3)) {
          xx = aggr[aggr[,2]==cn&aggr[,3]%in%p[[1]], 3:4]
          bpsTab = rbind(bpsTab, c(i, ifelse(cn=="3","Amplification","Deletion"), ff(i, p, aggr, cn), wilcoxtest(xx), NA))
        }
        xx = aggrLOH[aggrLOH[,2]==1&aggrLOH[,3]%in%p[[1]], 3:4]
        if (nrow(xx) > 0) bpsTab = rbind(bpsTab, c(i, "LOH", ff(i, p, aggrLOH, 1), wilcoxtest(xx), NA))
      }
    } else {
      vals = fad@vData$varTable[,vari_col]
      for (cn in c(1,3)) {
        xx = aggr[aggr[,2]==cn&aggr[,3]%in%vals, 3:4]
        bpsTab = rbind(bpsTab, c(i, ifelse(cn=="3","Amplification","Deletion"), rep(NA,4), cortest(xx)))
      }
      xx = aggrLOH[aggrLOH[,2]==1&aggrLOH[,3]%in%vals, 3:4]
      if (nrow(xx) > 0) bpsTab = rbind(bpsTab, c(i, "LOH", rep(NA,4), cortest(xx)))
    }
    message(paste(i, "variable processing: done."))
  }
	colnames(bpsTab) = c("Variable", "Alteration", "Value A", "Mean altered A", "Value B", "Mean altered B", "p value", "tau")
	rownames(bpsTab) = NULL
	if (!is.null(filename))
		write.table(bpsTab, filename, sep="\t",quote=FALSE,row.names=FALSE)
	bpsTab
}

alterationSummary = function(fad, filename=NULL) {
	out = .makeArmsInt(fad@cData, fad@chrs, fad@genome)
	if (any(as.integer(names(table(fad@cData$type))) > 10)) {
		aggrLOH = .aggregateInt(out, .n(fad), fad@genome, fad@arms, .getLOH, 0:1, fad@vData$varTable)
		aggrLOH$type = as.numeric(as.character(aggrLOH$type)) + 10
	} else {
		aggrLOH = c()
	}
	aggr = rbind(.aggregateInt(out, .n(fad), fad@genome, fad@arms, .unLOH, 1:3, fad@vData$varTable), aggrLOH)
	aggr = aggr[.ord(aggr$chr_q_arm,fad@arms,as.numeric(aggr$type)),]
	aggr$type = sapply(as.vector(aggr$type), function(i) {
								if (i==1) i="Deletion"
								else if (i==2) i="NormalPloidy"
								else if (i==3) i="Amplification"
								else if (i==10) i="NoLOH"
								else if (i==11) i="LOH"
					})
	aggr$freq = signif(as.numeric(aggr$freq), 3)
  aggr[,3] = as.integer(aggr[,3],2)
	colnames(aggr) = c("Arm", "Alteration", "Altered length", "Arm length", "Alteration frequency")
	if (!is.null(filename))
		write.table(aggr, filename, sep="\t",quote=FALSE,row.names=FALSE)
	aggr	
}

variableCor = function(fad, filename=NULL) {
	if (length(fad@vData$varNames) <= 2) # the last var is 'All'
		stop("Variable correlations can only be calculated if there are more two or more variables.")
	n = length(fad@vData$varNames)
  vc = fad@vData$varColumns[-n]
	aggr = data.frame(data.table(fad@cData)[,sum(bp_len), by=sapply(fad@cData[,c("code",vc)], list)])
	aggr = .fill(fad@vData$varTable, aggr, vc, 0)
	comb = combn(fad@vData$varNames[-n], 2)
	res = c()
	for (i in 1:ncol(comb)) {
	  types = fad@vData$varTypes[comb[1:2,i]]
    aggri = aggr[,fad@vData$varColumns[comb[,i]]]
	  aggri = aggri[apply(sapply(1:2, function(v) aggri[,v]%in%fad@vData$varValues[comb[v,i]][[1]]), 1, all) , ]
    resi = c(comb[1,i], comb[2,i], rep(NA, 6))
    if (all(types=="categorical")) {
      corTab = table(aggri)
      if (any(corTab<5)) chi = NA
      else chi = chisq.test(corTab)$p.value
      resi[3:4] = c(fisher.test(corTab)$p.value, chi)
    } else if (all(types=="quantitative")) {
      resi[7:8] = unlist(cor.test(aggri[,1], aggri[,2], method="k")[c("p.value","estimate")])
    } else {
      categ = which(types=="categorical")
      quant = which(types=="quantitative")
      f = aggri[,quant]~aggri[,categ]
      resi[5:6] = c(oneway.test(f)$p.value, kruskal.test(f)$p.value)
    }
    res = rbind(res, resi)
	}
	colnames(res) = c("Variable A", "Variable B", "Fisher", "Chisquare", "Welch_oneway", "Kruskal-Wallis", "Kendall_p", "Kendall_tau")
	rownames(res) = NULL
	res[,3:ncol(res)] = signif(as.numeric(res[,3:ncol(res)]), 3)
  if (!is.null(filename))
		write.table(res, filename, sep="\t",quote=FALSE,row.names=FALSE)
	res
}

#######################################################################
#################################### ARM-WISE BARPLOTS AND HISTOGRAMS #
#######################################################################

.checkSel = function(sel, arms) {
  wh = which(! sel%in%arms)
  if (length(wh) > 0)
    stop(paste("The following selected arms do not contain features and thus cannot be selected: ",paste(sel[wh],collapse=", "),".",sep=""))
}
.checkSelCols = function(sel, selColors) {
  if (length(sel) != length(selColors))
    stop("The selection and the slection colors vectors should be of the same length.")
}
.checkYlim = function(ylim) {
  if (length(ylim) != 2)
    stop("'ylim' should be a vector of length two.")
  if (ylim[1]>ylim[2])
    ylim = rev(ylim)
  if (ylim[1] < -1)
    ylim[1] = -1
  if (ylim[2] > 1)
    ylim[2] = 1
}
.checkAlter = function(alteration) {
  if (is.null(.getAlterCodes(alteration)))
    stop("The alteration type is not recognized.")
}
.checkVar = function(varName, fad) {
  if (length(varName)!=1)
    stop("No variables defined.")
  if (! varName%in%fad@vData$varNames)
    stop(paste("The variable name is not among those in the data: ",paste(fad@vData$varNames,collapse=", "),".",sep=""))
}
.checkVarValue = function(varName, varValue, fad) {
  values = fad@vData$varValuesNames[[varName]]
  if (! varValue%in%values)
    stop(paste("The variable value is not among those in the data: ",paste(values,collapse=", "),".",sep=""))
}

.plotBarInt = function(aggr, sel, selColors, ylim=c(-1,1), ylab, baseColor="black", genome, sex, arms, features) {
  aggr$freq = as.numeric(aggr$freq)
	aggr$freq[aggr$type==1] = aggr$freq[aggr$type==1] + aggr$freq[aggr$type==11]
	aggr$freq[aggr$type==3] = aggr$freq[aggr$type==3] + aggr$freq[aggr$type==13]
	aggr$freq[aggr$type%in%c(1,11)] = -aggr$freq[aggr$type%in%c(1,11)]
	aggr$negLOH = aggr$posLOH = aggr$neg = aggr$pos = NA
	aggr$negLOH[aggr$type==11] = aggr$freq[aggr$type==11]
	aggr$posLOH[aggr$type==13] = aggr$freq[aggr$type==13]
	aggr$neg[aggr$type==1] = aggr$freq[aggr$type==1]
	aggr$pos[aggr$type==3] = aggr$freq[aggr$type==3]

	cols = rep(baseColor, length(arms))
	if (!is.null(sel) & !is.null(selColors)) {
		pos_marked = which(arms%in%sel)
		cols[pos_marked] = selColors
	} else {
		pos_marked = 1:length(arms)
	}

  # alternate backgrounds for chromosomes
  xx = as.factor(gsub("[p|q]","",arms))
  odds = which(as.integer(as.character(xx))%%2==1)
  st = intersect(odds, which(!duplicated(xx)))-0.5
  wi = matrix(rle(as.vector(xx))$l,2)[1,]
	rect = data.frame(xmin=st, xmax=st+wi, ymin=-Inf, ymax=Inf)

	p = qplot(as.factor(chr_q_arm), pos, data=aggr, ylab=paste(ylab,"deletion/amplification frequency"), xlab="Chromosome arm") +
		geom_rect(data=rect, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="white", alpha=0.5, inherit.aes=FALSE) +
		geom_rect(aes(xmin=0,xmax=Inf,ymin=-diff(ylim)*0.005,ymax=diff(ylim)*0.005), fill="grey") +
		geom_point(aes(x=factor(chr_q_arm), y=aggr$pos, fill=factor(chr_q_arm)), shape=24, size=5, inherit.aes=FALSE) + 
		geom_point(aes(x=factor(chr_q_arm), y=aggr$posLOH, fill=factor(chr_q_arm)), shape=24, size=3, inherit.aes=FALSE) + 
		geom_point(aes(x=factor(chr_q_arm), y=aggr$neg, fill=factor(chr_q_arm)), shape=25, size=5, inherit.aes=FALSE) + 
		geom_point(aes(x=factor(chr_q_arm), y=aggr$negLOH, fill=factor(chr_q_arm)), shape=25, size=3, inherit.aes=FALSE) + 
		scale_fill_manual(values=cols) + theme(legend.position="none") + 
		scale_x_discrete(labels=arms) + scale_y_continuous(expand=c(0,0), limits=ylim)
	if (!is.null(sel)) {
		marked_aggr = aggr[aggr$chr_q_arm%in%sel&aggr$type%in%c(1,3),]
		marked_aggr = marked_aggr[.ord(marked_aggr$chr_q_arm,arms,marked_aggr$type),]
		marked_aggr = matrix(marked_aggr$freq[marked_aggr$chr_q_arm%in%sel], 2, byrow=FALSE)
		vals = t(rbind(marked_aggr,pos_marked-0.05,pos_marked+0.05))
		vals = data.frame(xmin=vals[,3],xmax=vals[,4],ymin=vals[,1],ymax=vals[,2])
		p = p + geom_rect(data=vals, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="black", inherit.aes=FALSE)
	}
	suppressWarnings(print(p))
}

plotBar = function(fad, byFeature=TRUE, sel=NULL, selColors=NULL, ylim=c(-1,1), baseColor="black", varName=NULL, value=NULL) {
  .checkSel(sel, fad@arms)
  .checkSelCols(sel, selColors)
  .checkYlim(ylim)
  if (is.null(varName)+is.null(value)==1) # only one of them is null
    stop("You should provide both a variable and its value or neither.")
  if (!is.null(varName)) {
    if (! .isCateg(fad,varName))
      stop("This type of plot is only applicable to categorical variables.")
    .checkVar(varName, fad)
    .checkVarValue(varName, value, fad)
    vn = fad@vData$varColumns[[varName]]
    vv = fad@vData$varValues[[varName]][fad@vData$varValuesNames[[varName]]%in%value]
  } else {
    vn = vv = NULL
  }
  selColors = selColors[order(factor(sel,levels=fad@arms))]
  if (byFeature) .plotBarFeatures(fad, sel, selColors, ylim, baseColor, vn, vv)
  else .plotBarBps(fad, sel, selColors, ylim, baseColor, vn, vv)
}

.plotBarFeatures = function(fad, sel, selColors, ylim, baseColor, vn, vv) {
	.checkFeatures(fad)
	gt = fad@gTable
	if (!is.null(vv)) gt = lapply(gt, function(i) i[which(i[,vn]%in%vv),])
	xx = sapply(gt, function(i) factor(i$type,levels=.allCodes()))
	xx2 = t(sapply(xx, table))
	xx3 = t(apply(xx2,1,function(i)i/sum(i)))
	rownames(xx3) = fad@features$chr_q_arm
	xx4 = melt(xx3)
	xx4 = xx4[xx4[,2]%in%.alterCodes(),]
	xx4[,1] = factor(xx4[,1], levels=fad@arms)
	xx4[,2] = factor(xx4[,2])
	aggr = data.frame(data.table(xx4)[, mean(value,na.rm=TRUE), by=list(xx4[,1],xx4[,2])])
	colnames(aggr) = c("chr_q_arm", "type", "freq")
	eg = do.call(expand.grid,list(chr_q_arm=fad@arms, type=.alterCodes()))
	aggr = merge(aggr, eg, all=TRUE)
	.plotBarInt(aggr, sel, selColors, ylim, "Mean feature", baseColor=baseColor, 
              fad@genome, fad@sex, fad@arms, fad@features)
}

.plotBarBps = function(fad, sel, selColors, ylim, baseColor, vn, vv) {
	if (!is.null(vv)) {
		out = .makeArmsInt(fad@cData, fad@chrs, fad@genome, vn)
		out = out[out[,vn]%in%vv, colnames(out)!=vn]
	} else {
		out = .makeArmsInt(fad@cData, fad@chrs, fad@genome)
	}
	out$chr_q_arm = factor(out$chr_q_arm, levels=fad@arms)
	out$type = factor(out$type, levels=.allCodes())
	aggr = .aggregateInt(out, .n(fad,vn,vv), fad@genome, fad@arms, identity, .allCodes(), fad@vData$varTable)
	.plotBarInt(aggr, sel, selColors, ylim, "Overall", baseColor=baseColor, fad@genome, fad@sex, fad@arms, fad@features)
}

plotHist = function(fad, alteration, varName, sel=NULL, selColors=NULL, selOnly=FALSE, 
                    baseColor="black", bin=0.05, xmax=1, ymax) {
	.checkSel(sel, fad@arms)
	.checkSelCols(sel, selColors)
	if (xmax <= 0)
		stop("'xmax', defining the limit in the X axis should be a positive number between 0 and 1.")
	if (xmax > 1)
		xmax = 1+bin
	if (ymax <= 0)
		stop("'ymax', defining the limit in the Y axis should be a positive number.")
	if (bin <= 0)
		stop("The bin width should be a positive number between 0 and 1.")
	if (bin > 1)
	  bin = 1
	.checkAlter(alteration)
	.checkVar(varName, fad)
  if (! .isCateg(fad, varName))
    stop("Only supported for categorical variables.")
	x = Reduce('+', lapply(.getAlterCodes(alteration), function(i) fad@fTable[[varName]][[i]]))
	x[x==1] = 0.9999999
	cc = as.factor(fad@features$chr_q_arm)
	cols = rep(baseColor, length(fad@arms))
	ord = order(factor(sel,levels=fad@arms))
	selColors = selColors[ord]
	sel = sel[ord]
	if (!is.null(sel) & !is.null(selColors)) {
		num_marked = which(levels(cc)%in%sel)
		cols[1:length(num_marked)] = selColors
	} else {
		num_marked = 1:length(cc)
	}
	cc_ord = cc[c(which(cc%in%levels(cc)[num_marked]),which(!cc%in%levels(cc)[num_marked]))]
	levs = cc_ord[!duplicated(cc_ord)]
	p = list()
	for (i in ncol(x):1) {
		dat = data.frame(Chromosome=cc, Frequency=x[,i])
		dat[,1] = factor(dat[,1], levels=levs)
		if (selOnly)
			dat = dat[dat$Chromosome%in%levels(cc)[num_marked],]
    lab = paste(varName, fad@vData$varValuesNames[fad@vData$varNames==varName][[1]][i], sep=": ")
		p[[i]] = qplot(Frequency, data=dat, geom="histogram", fill=as.factor(Chromosome), binwidth=bin) +
			theme(legend.position="none") + scale_fill_manual(values=cols) + ylab("Genes") +
			coord_cartesian(ylim=c(0,ymax), xlim=c(0,xmax))	 +
		  annotate("text", label=lab, x=xmax*0.8,y=ymax*0.9, size=6, color="white", fontface="bold", family="mono") +
			annotate("text", label=lab, x=xmax*0.8,y=ymax*0.9, size=6, family="mono")		
		if (i==ncol(x)) p[[i]] = p[[i]] + xlab(paste(alteration,"frequency"))  
		else p[[i]] = p[[i]] + xlab("")
	}
	g_legend<-function(){
		a.gplot = qplot(Frequency, data=dat[dat$Chromosome%in%levels(cc)[num_marked],], geom="histogram", fill=Chromosome, binwidth=bin) + 
			scale_fill_manual(values=selColors, labels=sel) + theme(legend.position="bottom")
		tmp = ggplot_gtable(ggplot_build(a.gplot))
		tmp$grobs[[which(sapply(tmp$grobs, function(x) x$name)=="guide-box")]]
	}
	if (!is.null(sel) & !is.null(selColors)) 
		grid.arrange(do.call(arrangeGrob, c(p,list(ncol=1))), g_legend(), nrow=2, heights=c(length(p)*5, 1))
	else
		grid.arrange(do.call(arrangeGrob, c(p,list(ncol=1))))
}

#######################################################################
################################################ ASSOCIATION ANALYSIS #
#######################################################################

.checkFeatures = function(fad) if (length(fad@gTable)==0) stop("Please run 'addFeatures' first.")

.checkModels = function(vars.m, areCateg, vars.0, argModel, toOrdered, toIntervals, strata) {
  if (length(vars.m)==0)
    stop("There are no variables in the model.")
  if (! all(vars.0%in%vars.m))
    stop("The null model is not nested within the model to test.")
  if ( !is.null(toOrdered)) {
    if (class(toOrdered)!="list") 
      stop("The ordered values should be a list whose names are within the variables used in the model.")
    if (! all(names(toOrdered)%in%names(areCateg[areCateg])))
      stop("Unknown categorical variables in the list of variables to convert to ordered.")
  }
  if ( !is.null(toIntervals)) {
    if (class(toIntervals)!="list") 
      stop("The interval breaks should be a list whose names are within the variables used in the model.")
    if (! all(names(toIntervals)%in%names(areCateg[!areCateg])))
      stop("Unknown continuous variables in the list of variables to convert to intervals.")
  }
  if (argModel=="unknown") {
    message(paste("Note: when the dependence between variables and CN is not defined,",
            "the model is limited to x1 + ... + xn | strata. See ?independence_test for more."))
  } else if (argModel%in%c("predictor","whole")) {
    if (! is.null(strata))
      warning("The 'strata' parameter is not used when CN depends on the variables. Use stratification within 'model' and 'nullModel'.")
  } else if (argModel=="response") {
    if (length(vars.m)!=1)
      stop("When variables depend on CN, association is only implemented for the case of a single variable.")
  }
}

.checkFullModel = function(m) {
  #contains the @ and ~ symbols
  if (length(grep("@",m))==0) stop("The symbol @, representing the copy number variable, must be in the model.")
  if (length(grep("~",m))==0) stop("The symbol ~, separating the response and the predictor, must be in the model.")
  m = try(formula(gsub("@","type",m)), silent=TRUE)
  if (class(m=="try-error")) stop("The model is not a valid formula. Please check it and try again.")
}

facopy = function(fad, alteration, model, nullModel=NULL, modelPart=c("response","predictor","unknown","whole")[1],
                       strata=NULL, toOrdered=NULL, toIntervals=NULL, sel=NULL, plot=FALSE, pvalThr=0.05, db=NULL,
                       link=c("logit","probit")[1], parametric=FALSE, design=c("binary","versus","lvog")[1],
                       FUN=NULL, ...) {
  same = function(a, b) all(a%in%b) & all(b%in%a)
  lwh = function(cc, vv, doFactor) {
    r = rep(0.5, length(tt[,"type"]))
    for (n in seq(length(cc))) r[which(tt[,"type"]%in%.getAlterCodes(cc[n]))] = vv[n]
    factor(r, levels=sort(unique(c(0.5,vv)))) # 0.5 is the middle value
  }
  doMultinom = function(tt, model) {
    x = formula(paste(model,"~type"))
    glm.model = try (multinom(x, data=tt, trace=FALSE), silent=TRUE)
    if (class(glm.model)[1]=="try-error") {
      NA
    } else {
      glm.0 = multinom(formula(paste(model,"~1")), data=tt, trace=FALSE)
      anova(glm.model, glm.0, test="Chisq")$P[2]
    }
  }
  ### checks
  argsModelNames = c("response","predictor","unknown","whole")
  argModel = pmatch(tolower(modelPart), argsModelNames)
  if (length(argModel)!=1 || is.na(argModel))
    stop("The 'modelPart' parameter should match one of the following: ", paste(argsModelNames,collapse=", "))
  argModel = argsModelNames[argModel]
  if (argModel=="whole" && is.null(FUN))
    stop("When providiing a full model, please also provide a function through the parameter 'FUN'.")
  .checkFeatures(fad)
  .checkAlter(alteration)
  varNames = all.vars(parse(text=model))
  dummy = sapply(varNames, .checkVar, fad)
  varNames.0 = all.vars(parse(text=ifelse(is.null(nullModel),"1",nullModel)))
  .checkModels(varNames, .isCateg(fad, varNames), varNames.0, argModel, toOrdered, toIntervals, strata)
  if (plot==TRUE && !is.null(sel)) {
    wwarning = "Plotting only available for genome-wide association (i.e. when 'sel' parameter not used)."
    warning(paste(wwarning, "Plotting disabled."))
    plot = FALSE
  }
  if (argModel=="whole") {
    .checkFullModel(model)
    if (!is.null(nullModel)) .checkFullModel(nullModel)
  }
  ### defs
  alterCode = .getAlterCodes(alteration)
  if (same(.getAlterCodes("cna"),alterCode) & design=="versus") {
    cc = c("amp","del"); vv = c(1,0)
  } else if (same(.getAlterCodes("all"),alterCode) & design=="lvog") {
    cc = c("someloss","onlygain"); vv = c(1,0)
  } else {
    cc = alteration; vv =  1
  }
  if (is.null(nullModel)) {
    nullModel = "1"
    if (argModel=="whole") {
      resp = gsub("~.*","",model)
      nullModel = paste(resp,"~",nullModel)
    }
  }
  if (argModel=="whole") {
    model = gsub("@","type",model)
    nullModel = gsub("@","type",nullModel)
  }
  tabs = fad@gTable
  if (! is.null(sel)) {
    sel = which(fad@features$feature%in%sel)
    tabs = tabs[sel]
  }
  ind = sapply(varNames, function(i) which(fad@vData$varNames==i))
  vars = fad@vData$varColumns[ind]
  ### assoc
  p1 = p2 = p3 = c()
  ps= c()
  for (i in tabs) {
    tt = i[,c(vars,"type",strata)]
    colnames(tt)[seq(ncol(tt)-1)] = fad@vData$varNames[ind]
    if (nrow(tt)==0) { ps = c(ps, NA); next }
    for (n in colnames(tt)[seq(ncol(tt)-1)]) {
      varVals = fad@vData$varValues[[n]]
      tt = tt[tt[,n]%in%varVals,]
      if (.isCateg(fad, n)) {
        if (n%in%names(toOrdered)) {
          for (o in seq(length(varVals)))
            tt[,n][tt[,n]==varVals[o]]=toOrdered[[n]][o]
        } else {
          tt[,n] = factor(tt[,n], levels=varVals)
        }
      } else if (n%in%names(toIntervals)) {
          tt[,n] = cut(tt[,n], c(-Inf,toIntervals[[n]],Inf), labels=seq(length(toIntervals[[n]])+1))
          tt[,n] = as.numeric(tt[,n])
      }
    }
    if (! is.null(strata)) {
      strataVals = fad@vData$varValues[[strata]]
      tt = tt[tt[,strata]%in%strataVals,]
      tt[,ncol(tt)] = factor(tt[,ncol(tt)], levels=strataVals)
    }
    tt[,"type"] = lwh(cc, vv)
    if (length(table(as.vector(tt[,"type"]))) < 2) {
      pvalue = NA
    } else {
      if (argModel=="unknown") {
        x = paste("type~",model)
        if (!is.null(strata)) x = paste(x,"|",strata)
        pvalue = pvalue(independence_test(formula(x), data=tt))
      } else if (argModel=="predictor") {
        glm.model = try (glm(paste("type~",model), data=tt, family=binomial(link)), silent=TRUE)
        if (class(glm.model)[1]=="try-error") {
          pvalue = NA
        } else {
          glm.0 = glm(paste("type~",nullModel),  data=tt, family=binomial(link))
          pvalue = anova(glm.model, glm.0, test="Chisq")$P[2]          
        }
      } else if (argModel=="response") { # model == varNames; length(varNames) == 1
        x = formula(paste(model,"~type"))
        if (.isCateg(fad, model)) {
          if (is.factor(tt[,model])) {
            # nnet::multinom(). FOR CATEGORICAL VALUES.
            # "The variables on the rhs of the formula should be roughly scaled to 
            # [0,1] or the fit will be slow or may not converge at all." THIS IS DONE.
            pvalue = doMultinom(tt, model)
          } else { # ordered
            if (length(table(toOrdered[[model]]))==2) { 
              # nnet::multinom() and stats::glm() do the same with 2 levels, but glm() requires 0 <= y <= 1
              pvalue = doMultinom(tt, model)
            } else {
              #MASS::polr(). FOR ORDERED VALUES. ONLY 3 OR MORE LEVELS.
              tt[,model] = factor(tt[,model], levels=toOrdered[[model]])
              method = ifelse(link=="logit","logistic","probit")
              glm.model = polr(x, data=tt, method=method)
              glm.0 = polr(formula(paste(model,"~1")), data=tt, method=method)
              pvalue = anova(glm.model, glm.0, test="Chisq")$P[2]
            }
          }
        } else {
          if (parametric) FUN = oneway_test
          else FUN = kruskal_test
          pvalue = pvalue(FUN(formula(x), data=tt))
        }
      } else { # argModel=="whole"
        res = FUN(formula(model), formula(nullModel), data=tt, ...)
        if ("htest"%in%class(res)) pvalue = res$p.value
        else if ("lm"%in%class(res)) pvalue = summary(res)$coef[2,4]
        else if (attr(class(res),"package")=="coin") pvalue = pvalue(res)
        else if ("numeric"%in%class(res)) pvalue = res
        else stop("The provided function FUN is unknown. Please provide a function whose output is a pvalue.")
      }
    }
    ps = c(ps, pvalue)
    if (length(ps)%%500==0)
      message(paste(round(length(ps)*100/length(tabs),1), "% done...",sep=""))
  }
  if (length(varNames)!=1)
    varNames = "All"
  if (plot)
    .facopyPlotInt(fad, alteration, varNames, db, ps, pvalThr)
  if (is.null(sel))
    sel = seq(nrow(fad@features))
  cbind(fad@features[sel,c("feature","bp_st","bp_en","chr_q_arm")], p_value=round(ps,4))[,c(1,5,4,2,3)]
}

facopyPlot = function(fad, alteration, varName, db=NULL) .facopyPlotInt(fad, alteration, varName, db)

.facopyPlotInt = function(fad, alteration, varName, db=NULL, pvals=NULL, pvalThr=0.05) {
  ### PLOTTING FUNCTIONS
  ######################################
  emptyPlot = function(xlab=NULL, ylab=NULL) {
    ep = ggplot(data.frame()) + geom_point() + scale_x_continuous(expand=c(0,0)) + ylab("") +
      scale_y_continuous(expand=c(0,0), limits=c(0,1)) + theme_classic() + 
      theme(axis.line=element_line(color='white'), axis.ticks=element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(),
            text=element_text(size=ifelse(is.null(ylab),12,20)))
    if (!is.null(xlab)) ep = ep + xlab(xlab)
    if (!is.null(ylab)) ep = ep + ylab(ylab)
    ep
  }
  printPlot = function(p) {
    g_legend = function(){
      a.gplot = p[[1]] + theme(legend.position="bottom", text=element_text(size=20))
      tmp = ggplot_gtable(ggplot_build(a.gplot))
      tmp$grobs[[which(sapply(tmp$grobs, function(x) x$name)=="guide-box")]]
    }
    if (fad@genome%in%c("hg18","hg19")) ncol = 8
    if (fad@genome%in%c("mm8")) ncol = 7
    grid.arrange(emptyPlot(NULL, paste("\n",alteration,"frequency [0,1]")), 
                 do.call(arrangeGrob, c(p,list(ncol=ncol))), emptyPlot(), g_legend(), 
                 nrow=2, ncol=2, widths=c(0.25, ncol), heights=c(12, 1))
  }
  ### DATABASE PROFILES
  ######################################
	if (! is.null(db)) {
		datText = paste(fad@genome,"_db_",db,sep="")
		data(list=datText)
		dat = eval(parse(text=datText))
		a = unlist(.alterations())[pmatch(tolower(alteration), unlist(.alterations()))]
    if (is.na(a)) stop("Unknown alteration type: ", alteration)
		if (a=="amplifications") { 
      sel = "amp"
		} else if (a=="deletions") {
      sel = "del"
		} else if (a%in%c("cnas","any","loh")) {
		  sel = c("amp", "del")
		} else {
		  sel = c("amp", "del")
      warning("Alterations from databases in parameter 'db', including '",db,"', only involve amplifications and deletions.")
		}
		dat = dat[dat$type%in%sel,]
	}
  ### PLOTTING
  ######################################
  message("Processing plot...")
	p = list()
	for (arm in .allChrs(fad@genome,fad@sex)) {
		i = length(p)+1
    sel = which(fad@features$chr_q_arm==arm)
		featnames = as.vector(fad@features$feature[sel])
    if (length(sel)==0) next
    ### frequency calculation
    varVals = fad@vData$varValuesNames[[varName]]
		subxx = list()
		for (ac in .getAlterCodes(alteration)) {
			a = fad@fTable[[varName]][[ac]][featnames,,drop=FALSE]
			if (!is.null(a)) rownames(a) = featnames
			a[is.na(a)] = 0
			subxx[[length(subxx)+1]] = a
		}
		xx = as.matrix(Reduce(`+`, subxx))
		if (nrow(xx)==0) { p[[i]] = emptyPlot(paste("No visualizable data for arm",arm));	next }
		rownames(xx) = fad@features$bp_st[sel]
		colnames(xx) = varVals
		xx2 = melt(xx)
    xx2[,2] = factor(xx2[,2], levels=varVals)
		colnames(xx2) = c("bp",varName,"value")
		xlim = c(min(xx2$bp), max(xx2$bp))
    ### plot initialization
    p[[i]] = ggplot(xx2, aes_string(x="bp",y="value",colour=varName,group=varName)) + 
		  scale_y_continuous(breaks=seq(0.25,0.75,0.25), labels=rep("",3)) + 
      # scale_x_continuous(breaks=pretty_breaks(n=3)) +
      ylab(NULL) +
		  xlab(paste(arm, " (chr",gsub("[p|q]","",arm)," position)",sep="")) +
		  coord_cartesian(ylim=c(0,1), xlim=xlim)
		### add shading for significant features 
    if (! is.null(pvals)) {
      signif = which(pvals[sel]<pvalThr & !is.na(pvals[sel]))
  		if (length(signif) > 0) {
  		  pos1 = xx2$bp[sapply(signif-1, max,1)]
  		  pos2 = xx2$bp[signif]
  		  pos3 = xx2$bp[sapply(signif+1, min,length(sel))]
  		  irang = IRanges(start=pos2-(pos2-pos1)/2, end=pos2+(pos3-pos2)/2-1)
  		  cov = coverage(irang)
  		  csum = cumsum(cov@lengths)
  		  csum[cov@values==0] = csum[cov@values==0]+1 # csum[cov@values==0] are starts of intervals
  		  if (cov@values[1]==1) csum = c(1, csum)
  		  vals = data.frame(matrix(csum,ncol=2,byrow=TRUE))
  		  colnames(vals) = c("xmin","xmax")
  		  p[[i]] = p[[i]] + geom_rect(data=vals, aes(xmin=xmin,xmax=xmax,ymin=0,ymax=1), fill="black",alpha=0.5,inherit.aes=FALSE)
  		}
    }
    ### database visualization
		if (!is.null(db)) {
			dat_chr = dat[dat$chr==gsub("[p|q]","",arm),]
			dat_chr$freq = dat_chr$freq/100
			if (nrow(dat_chr) > 0)
				p[[i]] = p[[i]] + geom_rect(data=dat_chr, aes(xmin=pos_st,xmax=pos_en,ymin=0,ymax=freq),
											inherit.aes=FALSE, fill="black",alpha=0.75)
		}
    ### points, axes and legends
    cols = .refcols(length(varVals))
    p[[i]] = p[[i]] + geom_line(size=1.5,alpha=0.5) + 
      geom_point(colour=cols[as.numeric(xx2[,2])],size=3,alpha=0.5) +
      theme(legend.position="bottom") + scale_color_manual(values=cols) 
	}
	for (i in 1:length(p))
		p[[i]] = p[[i]] + theme(legend.position="none", plot.margin=unit(c(0.02,0.02,0,0),"npc"))
  message("Generating plot, please wait while it appears on the graphics device...")
  printPlot(p)
}

#######################################################################
############################################ GENE ENRICHMENT ANALYSIS #
#######################################################################

.graphPathway = function(p, pval, genes) {
  g <- pathwayGraph(p)
  gg <- setting.graph.attributes(igraph.from.graphNEL(g))
  V(gg)$color = "#00000022"
  wh = toupper(V(gg)$label)%in%toupper(genes)
  V(gg)$color[wh] = "red"
  gg = delete.edges(gg, which(is.loop(gg, eids=E(gg))==TRUE))
  gg = as.undirected(gg)
  degrees = igraph::degree(gg,v=V(gg)) + 1
  degrees = degrees / min(degrees)
  degrees = sapply(degrees, min, 10)
  plot.igraph(gg, vertex.label.font = 2, vertex.frame.color = V(gg)$color,
              vertex.label.color = "#000000cc", vertex.label.cex = 5/sqrt(length(p@nodes)),
              vertex.size = 50/sqrt(length(p@nodes)) * sqrt(degrees),
              edge.width=10/sqrt(length(p@nodes)), edge.curved=TRUE, edge.arrow.mode=0,
              layout = layout.fruchterman.reingold(graph=gg,repulserad=vcount(gg)^3*0.005))
  mtext(paste(p@title,"\np-value: ",pval,"",sep=""))
}

.enrichment = function(dbName, db, dbCategNames, folder, allsymbols, ress, fad, plotThr) {
  # N  	black+white total	score_tg
  # K		white total			score_all
  # k		white drawn			score_eg
  # n		black+white drawn	score_sg
  allsymbols = unique(allsymbols)
  sigGeneSymbol = ress$feature
  totalgenes = length(allsymbols)
  signifgenes = length(sigGeneSymbol)
  score_sg = length(table(as.vector(ress$chr_q_arm)))
  score_tg = length(table(fad@features$chr_q_arm[sapply(fad@gTable,nrow)>0 & fad@features$feature%in%allsymbols]))
  tab = c()
  for (j in 1:length(db)) {
    wh = which(sigGeneSymbol%in%db[[j]])
    enrichmentGenes = sigGeneSymbol[wh]
    count = length(wh)
    wh_all = which(db[[j]]%in%allsymbols)
    size = length(wh_all)
    score_eg = length(table(as.vector(ress$chr_q_arm[wh])))
    score_all = length(table(as.vector(fad@features$chr_q_arm[which(fad@features$feature%in%db[[j]][wh_all])])))
    score_pvalue = phyper(score_eg, score_all, score_tg-score_all, score_sg, lower.tail=TRUE)
    expcount = signifgenes*size/totalgenes
    pvalue = phyper(count-1L, size, totalgenes-size, signifgenes, lower.tail=FALSE)
    oddsratio = count*(totalgenes-size-signifgenes+count) / ((signifgenes-count) * (size-count))
    tab = rbind(tab, c(dbCategNames[j], size, signif(expcount), count, signif(oddsratio), signif(pvalue),
                       paste(score_eg,round(score_pvalue,3),sep="_"), paste(enrichmentGenes,collapse=", ")))
  }
  tab = cbind(tab[,1:6], signif(p.adjust(tab[,6],'fdr')), tab[,7:8])
  colnames(tab) = c("Category","Size","ExpCount","Count","OddsRatio","Pvalue","PvalueAdj","Score","Genes")
  tab = tab[order(as.numeric(tab[,"Pvalue"])),]
  write.table(tab, file.path(folder,paste(dbName,".txt",sep="")), row.names=FALSE,quote=FALSE, sep="\t")
  
  canonical = c("kegg","reactome","biocarta")
  go = c("gobp","gocc","gomf")
  if (dbName%in%c(canonical,go)) {
    bottom.r = tail(which(tab[,"Pvalue"]<plotThr),1)
    if (length(bottom.r) > 0) {
      if (dbName%in%canonical) {
        dir.create(file.path(folder,dbName))
        for (r in 1:bottom.r) {
          graphObjsText = paste("facopy_",dbName,sep="")
          data(list=graphObjsText)
          graphObjs = eval(parse(text=graphObjsText))
          p = graphObjs[[tab[r,1]]]
          filename = paste(r,"_",gsub("[:|/]","_",substr(p@title,1,32)),".pdf", sep="")
          pdf(file.path(folder, dbName, filename), 12, 9)
          .graphPathway(p, tab[r,"Pvalue"], strsplit(tab[r,"Genes"],", ")[[1]])
          dev.off()
        }
      } else if (dbName%in%go) {
        onto = toupper(substr(dbName, 3,4))
        goparents = eval(parse(text=paste("GO",onto,"PARENTS",sep="")))
        sigGO.fdr = tab[seq(bottom.r),"Pvalue"]
        names(sigGO.fdr) = tab[seq(bottom.r),1]
        sigGO.term = sapply(names(sigGO.fdr), function(j) {
          name = try(getGOTerm(j)[[onto]], silent=TRUE)
          name = ifelse(class(name)=="try-error", "", name)
          paste(j,paste(strwrap(name,24),collapse="\\\n"),sigGO.fdr[j],sep="\\\n")
        })
        colors = sapply(as.numeric(sigGO.fdr), function(j) heat.colors(80)[min(80,max(1,round(1/log10(1/j)*100)))])
        names(colors) = names(sigGO.fdr)
        grph = GOGraph(names(sigGO.fdr),goparents)
        sizeFactor = round(sqrt(length(grph@nodes)))
        isSignifNode = sapply(grph@nodes, function(n) ifelse(n%in%names(sigGO.fdr),TRUE,FALSE))
        shapes = sapply(isSignifNode, ifelse, "rectangle", "circle")
        widths = sapply(isSignifNode, ifelse, 0.23, 0.1)
        heights = sapply(isSignifNode, ifelse, 0.13, 0.1)
        fontsizes = sapply(isSignifNode, ifelse, 20, 14)
        nodeAttrs = list(shape=shapes, label=sigGO.term, fillcolor=colors, 
                         width=widths*sizeFactor, height=heights*sizeFactor, fontsize=fontsizes)
        pdf(file.path(folder, paste(dbName,".pdf",sep="")), 2.0*sizeFactor, 1.4*sizeFactor)
        Rgraphviz::plot(grph, nodeAttrs=nodeAttrs)
        dev.off()
      }
    }
  }
  
  tab
}

calculateCor = function(fad, exprProfile, db=NULL) {
  if (! is.null(db)) {
    if (class(exprProfile)!="character" || class(db)!="character") 
      stop("'exprProfile' and 'db' should be characters indicating database and name of expression profile in cBio portal.")
    cor = .corInt(db, exprProfile, fad@features$feature)
  } else {
    if (class(exprProfile)=="character") {
      exprProfile = read.table(exprProfile, header=TRUE,sep="\t")
    } else if (class(exprProfile)=="data.frame") {
      reqCols = c("chr","bp_st","bp_en","feature","chr_q_arm")
      if (! all(reqCols%in%colnames(exprProfile)))
        stop("At least the following columns should be present in the external data: ", paste(reqCols,collapse=", "))
    } else {
      stop("External data with information on expression should be a data.frame.")
    }
    sampleNames = as.character(fad@vData$varTable$code)
    cnProfile = t(sapply(fad@gTable, function(i) { x = i$type; names(x) = i$code; x[sampleNames] }))
    colnames(cnProfile) = sampleNames
    ## next 2 represent a random exprProfile
    #exprProfile = cnProfile[,ncol(cnProfile)]
    #exprProfile[,] = rnorm(length(exprProfile))
    # exprProfile contains the same genes (rows) and samples (columns) as the CN data
    whrow = intersect(rownames(exprProfile),rownames(cnProfile))
    whcol = intersect(colnames(exprProfile),colnames(cnProfile))
    exprProfile = exprProfile[whrow,whcol]
    cnProfile = cnProfile[whrow,whcol]
    cor = .makeCor(cnProfile, exprProfile)
  }  
  list(cor=cor, db=db, profile=exprProfile)
}

facopyEnrichment = function(fad, geneTable, cor, outFolder, pvalThr=0.05, corThr=0.1, plotThr=0.001) {
  data(facopy_msigdb)
  data(facopy_msigdbNames)
  if (! fad@genome%in%c("hg18","hg19"))
    stop("Gene-set enrichment analysis is only available for human genomes.")
  geneTable = geneTable[order(geneTable$p_value),]
  cors = cor$cor[as.character(geneTable$feature)]
  cors_m = cbind(names(cors), round(cors,3))
  colnames(cors_m) = c("feature", "cor_expr")
  geneTable = merge(geneTable, cors_m, all.x=TRUE)
  geneTable = geneTable[order(geneTable$p_value, -as.numeric(as.character(geneTable$cor_expr))),]
  if (!is.null(outFolder)) {
    dir.create(outFolder)
    wh = which(as.vector(geneTable$cor_expr)>=corThr & as.vector(geneTable$p_value)<=pvalThr)
    univ = intersect(names(cor$cor)[cor$cor>corThr], fad@features$feature)
    for (x in names(facopy_msigdbNames))
      dummy = .enrichment(x, facopy_msigdb[[x]], facopy_msigdbNames[[x]], outFolder, univ, geneTable[wh,], fad, plotThr)
  }
  geneTable
}

.makeCor = function(m, n) {
  out = sapply(seq(nrow(m)), function(i) {
    cn = as.numeric(as.vector(m[i,]))
    expr = as.numeric(as.vector(n[i,]))
    if (!all(is.na(cn)) && !all(is.na(expr)))
      summary(lm(expr~cn))$r.squared
    else NA 
  })
  names(out) = rownames(m)
  out
}

.corInt = function(dbGeneExpr, exprProfile, x) {
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  # test(mycgds)
  #getCancerStudies(mycgds)
  mycancerstudy = getCancerStudies(mycgds)[getCancerStudies(mycgds)[,1]==dbGeneExpr,1]
  mycaselist = grep("_all$", getCaseLists(mycgds,mycancerstudy)[,1], value=TRUE)[1]
  mygeneticprofiles = tolower(getGeneticProfiles(mycgds,mycancerstudy)[,1])
  cn_profile = grep("_gistic|_cna", mygeneticprofiles, value=TRUE)
  expr_profile = paste(mycancerstudy,"_",exprProfile,sep="")	
  n = length(x)
  chunk = 500
  cgds_gene = list()
  for (i in c(cn_profile, expr_profile)) {
    cgds_gene[[i]] = NULL
    for (j in 1:ceiling(n/chunk)) {
      st = 1 + chunk*(j-1)
      en = min(chunk*j,n)
      genes_chunk = toupper(x[st:en])
      genes_chunk = gsub("^MIRN", "MIR", genes_chunk) # match miRNA nomenclature
      trials = 0
      repeat{
        prof = try( t(getProfileData(mycgds,genes_chunk,i,mycaselist)) ,silent=TRUE)
        if (! class(prof)=="try-error") break
        trials = trials+1
        if (trials==3) message("Warning: connection to cBio portal is not stable. You might eventually get an error.")
        if (trials==5) stop("Cannot connect to cBio portal. Maybe connection is down?")
        Sys.sleep(0.5)
      }
      res = matrix(NA, ncol=ncol(prof),nrow=length(genes_chunk))
      rownames(res) = genes_chunk
      wh = genes_chunk%in%rownames(prof)
      res[wh,] = prof[genes_chunk[wh],]
      cgds_gene[[i]] = rbind(cgds_gene[[i]], res)
      message(paste(i," ",round(j*100/(ceiling(n/chunk)),2),"% done...",sep=""))
    }
  }
  .makeCor(cgds_gene[[1]], cgds_gene[[2]])
}

#######################################################################
######################################## PLOT ARM OR FEATURE VICINITY #
#######################################################################
#.plotVicinityInt(fad, NULL, gsub("[p|q]","",arm), st, en, alteration, 0, varName, arm)
.plotVicinityInt = function(fad, name=NULL, chr, st=NULL, en=NULL, 
							alteration=NULL, margin=0, varName=NULL, arm=NULL) {
	oldmargins <- par("mar")
	makeAb = function(col) {
		abline(v=c(st,en), lwd=3, col=col)
		abline(v=c(st-margin,en+margin), lwd=1, col=col)
	}
	alterCode = .getAlterCodes(alteration)
	colname = fad@vData$varColumns[[varName]]
  dfsub = fad@cData[ fad@cData[,colname]%in%fad@vData$varValues[[varName]] , ]
	dfsub$code = as.factor(dfsub$code)
	dfsub = dfsub[dfsub$chr==chr,]
	if (!is.null(st) & !(is.null(en)))
		dfsub = dfsub[dfsub$pos_st<en+margin & dfsub$pos_en>st-margin,]
	dfsub = dfsub[dfsub$type%in%alterCode,]
	if (nrow(dfsub)==0)
	  stop("There are no alterations of the selected type in the vicinity of the feature for that variable.")
	dfsub = dfsub[order(dfsub[,colname]),]
	xlim = c(min(dfsub$pos_st), max(dfsub$pos_en))
	if (is.null(name))
		for (i in levels(dfsub$code)[!levels(dfsub$code)%in%dfsub$code]) {
			ro = rep(NA, ncol(dfsub))
			names(ro) = colnames(dfsub)
			ro = as.data.frame(t(as.matrix(ro)))
			ro$code = i
			ro$pos_st = ro$pos_en = -1
			ro[,colname] = fad@cData[which(fad@cData$code==i)[1],colname]
			dfsub = rbind(dfsub, ro)
		}
	vals = - as.numeric(as.factor(dfsub[,colname])) # negative in order to plot from top to bottom
	aggr = data.frame(data.table(dfsub)[,sum(bp_len), by=list(code)]) #.fill not necessary
	aggr[,2] = rank(aggr[,2], na.last=FALSE)
	dfsub$alterRank = sapply(dfsub$code, function(i) aggr[aggr[,1]==i,2])
	dfsub = dfsub[order(vals,dfsub$alterRank),]
	vals = sort(vals)
	height = 4
	bins = dfsub$alterRank
	bins = vals*max(bins) + bins
	bins = as.numeric(as.factor(bins))
	ylim = c(0, max(bins)*(height + 1))
	layout(matrix(c(rep(1,15),2)))
	plot.new()
	plot.window(xlim, ylim)
	grid(ny=0, col="grey", lty=1)
	if (!is.null(name))
		makeAb("#000000")
	ybottom = bins * (1 + height) - height
  if (.isCateg(fad,varName)) {
	  factors = fad@vData$varValues[[varName]]
  } else {
    qu = quantile(as.numeric(names(table(fad@cData[,fad@vData$varColumns[[varName]]]))),c(0,0.5,1))
    factors = seq(qu[1],qu[3])
  }
	cols = .refcols(length(factors))
	allcols = cols[as.numeric(factor(dfsub[,colname], factors))]
	rect(dfsub$pos_st-0.5, ybottom, dfsub$pos_en+0.5, ybottom+height, col=allcols, border=0)
	grid(ny=0, col="white")
	if (!is.null(name))
		makeAb("#00000033")
	if (!is.null(name))
		title(paste(varName, ", alterations (",alteration,") near ", name, sep=""))
	else if (!is.null(arm))
		title(paste(varName, ", alterations (",alteration,") in chr", arm, sep=""))
	axis(1)
  mtext(text=paste("Chromosome",chr,"position"), cex=1, side=1, line=2.5)
  par(mar=c(0,0,0,0))
	plot.new()
	if (.isCateg(fad,varName)) {
	  legend("top", legend=fad@vData$varValuesNames[[varName]], col=cols, horiz=TRUE, lwd=10, box.col="white", cex=1.3)
	  uniq = fad@vData$varTable[,c("code",colname)]
	  uniq = uniq[uniq[,colname]%in%factors,]
	  uniq[,2] = factor(uniq[,2], levels=factors)
    tab = table(merge(aggr, uniq)[,3])
    if (is.null(arm)) tab = paste(tab, "of", table(uniq[,2]))
	  legend("bottom", legend=tab, col=cols, horiz=TRUE, lwd=10, bty="n", cex=1.3)
  } else {
    image(y=c(0.2,0.6), z=matrix(seq(length(allcols)),ncol=1), col=allcols, axes=FALSE,xlab="",ylab="", add=TRUE)
    legend("topleft", legend=qu[1], cex=1.3, bty="n", adj=1)
    legend("top", legend=qu[2], cex=1.3, bty="n")
    legend("topright", legend=qu[3], cex=1.3, bty="n")
  }
	par(mar=oldmargins)
}

plotZoom = function(fad, what=c("feature","arm"), name, alteration, varName, margin=0) {
  arg = function(x) match.arg(tolower(x),c("feature","arm"))
  if (arg(what)=="feature") .plotVicinity(fad, name, alteration, varName, margin)
  if (arg(what)=="arm") .plotArm(fad, name, alteration, varName)
}

.plotVicinity = function(fad, name, alteration, varName, margin=0) {
	info = fad@features[fad@features$feature==name,]
	.plotVicinityInt(fad, name, info$chr, info$bp_st, info$bp_en, alteration, margin, varName, NULL)
}

.plotArm = function(fad, arm, alteration, varName) {
	sepsText = paste(fad@genome,"_armLimits",sep="")
	data(list=sepsText)
	seps = eval(parse(text=sepsText))
	wh = which(seps$chr_q_arm==arm)
	en = seps$limit[wh]
	if (fad@genome%in%c("hg18","hg19")) st = ifelse(length(grep("p",arm))==1, 0, seps$limit[wh-1])
  if (fad@genome%in%c("mm8")) st = 0
	.plotVicinityInt(fad, NULL, gsub("[p|q]","",arm), st, en, alteration, 0, varName, arm)
}

#######################################################################
######################################## PRINCIPAL COMPONENT ANALYSIS #
#######################################################################

.makeMatrixInt = function(gtab, alterCode, samples, design)	{
	same = function(a, b) all(a%in%b) & all(b%in%a)
	lwh = function(cc, vv) lapply(gtab, function(g) {
		r = rep(0, length(samples))
		for (i in seq(length(cc))) r[which(samples%in%g$code[g$type%in%cc[[i]]])] = vv[i]
		r
	})
	if (same(.getAlterCodes("cna"),alterCode) & design=="versus") {
		xx = lwh(list(.getAlterCodes("amp"),.getAlterCodes("del")), c(1,-1))
	} else if (same(.getAlterCodes("all"),alterCode) & design=="lvog") {
		xx = lwh(list(.getAlterCodes("someloss"),.getAlterCodes("onlygain")), c(-1,1))
	} else {
		xx = lwh(list(alterCode), 1)
	}
	do.call(rbind, xx)
}

plotPCA = function(fad, alteration, varName, sel=NULL, 
					design=c("binary","versus","lvog")[1], do.plot=TRUE, by.size=TRUE, cex=4) {
					# versus design requires alteration="CNA"
					# lvog (loss versus only gain) design requires alteration="any","all"
	.checkFeatures(fad)
	.checkAlter(alteration)
	.checkVar(varName, fad)
	.checkSel(sel, fad@arms)
	if (cex < 1)
		stop("The size of the points 'cex' should be a positive number.")
	alterCode = .getAlterCodes(alteration)
	samples = names(table(fad@cData$code))
	rang = 1:length(fad@gTable)
	if (!is.null(sel)) rang = rang[fad@features$chr_q_arm%in%sel]
	m = .makeMatrixInt(fad@gTable[rang], alterCode, samples, design)
	vc = fad@vData$varColumns[[varName]]
	vv = fad@vData$varValues[[varName]]
	dat = fad@vData$varTable[,c("code",vc)]
	dat = dat[dat[,vc]%in%vv,]
	wh = which(fad@cData[,vc]%in%vv & fad@cData$type%in%alterCode)
	aggr = data.frame(data.table(fad@cData[wh,])[, sum(bp_len), by=list(code)])
	aggr[,2] = rank(aggr[,2], na.last=FALSE)
  alterRank = sapply(dat$code, function(i) aggr[aggr[,1]==i,2])
  alterRank[sapply(alterRank, length)==0] = 0
	dat$alterRank = unlist(alterRank)
	colnames(m) = samples
	m = m[,colnames(m)%in%dat$code]
	rownames(m) = fad@features$feature[rang]
	pc = PCA(t(m), graph=FALSE)
	if (do.plot) {
		ord = unlist(sapply(rownames(pc$ind$coord), function(i) which(dat$code==i)))
		cols = .refcols(length(vv))
		cols_t = paste(cols,"44",sep="")
		allcols = cols[as.numeric(as.factor(dat[ord,vc]))]
		allcols_t = cols_t[as.numeric(as.factor(dat[ord,vc]))]
		s = 2
		if (by.size)
			s = seq(1,6,length=nrow(dat))[dat$alterRank[ord]]
		labs = paste("PC",1:2, " (", round(pc$eig[1:2,2],2), "%)", sep="")
		ymin = min(pc$ind$coord[,2])
		ymax = max(pc$ind$coord[,2])
    ytop = ymax+(ymax-ymin)*0.12
		plot(pc$ind$coord[,1:2], xlab=labs[1], ylab=labs[2], cex=cex, pch=21, col=allcols, lwd=s, ylim=c(ymin,ytop), 
			main=paste("PCA clustering (",alteration,", ",varName,"), ",design," design",sep=""),	frame.plot=FALSE)
		points(pc$ind$coord[,1:2], xlab=labs[1], ylab=labs[2], cex=cex, pch=19, col=allcols_t)
    if (.isCateg(fad,varName)) {
      legend("top", legend=fad@vData$varValuesNames[[varName]], col=cols, horiz=TRUE, lwd=10, box.col="white", cex=1.3)
		} else {
		  xmin = min(pc$ind$coord[,1])
		  xmax = max(pc$ind$coord[,1])
  		image(x=seq(xmin,xmax,l=length(cols)), y=c(ymax+(ytop-ymax)*0.6,ytop-(ytop-ymax)*0.2),
            z=matrix(rep(seq(length(cols)),2),ncol=2), col=cols, axes=FALSE,xlab="",ylab="", add=TRUE)
      qu = quantile(dat[,vc],c(0,0.5,1))
  		legend("topleft", legend=qu[1], cex=1.2, bty="n", adj=1)
  		legend("top", legend=qu[2], cex=1.2, bty="n")
  		legend("topright", legend=qu[3], cex=1.2, bty="n")
		}
	}
	pc
}
