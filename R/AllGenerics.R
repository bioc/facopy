.summary.facopyInfo = function(object) {
	cat("Object of type facopyInfo\n")
	nSamples = length(unique(object@cData$code))
  nCalls = round(nrow(object@cData) / nSamples, 1)
	cat(" - Number of samples:", nSamples, "\n")
	cat(" - Mean calls per sample:", nCalls, "\n")
	variables = object@vData$varNames[seq(length(object@vData$varNames)-1)]
	if (length(variables)>0) {	
		cat(" - Variables (",length(variables),"): ", paste(variables, collapse=", "), "\n",sep="")
	}
  if (nrow(object@features)>0) {
    cat(" - Features (",nrow(object@features),"): ", object@what, " from ",object@genome,"\n",sep="")
  }
	cat(" - Sex chromosomes: ", ifelse(length(object@sex)>0,paste(object@sex,collapse=", "), "not included"), sep="")
  cat("\n")
}
