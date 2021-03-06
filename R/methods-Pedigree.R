validPedigree <- function(object){
  msg <- NULL
  if(nrow(trios(object)) > 0){
    if(!identical(colnames(trios(object)), c("F", "M", "O"))){
      msg <- "column names should be 'F', 'M', and 'O'"
      return(msg)
    }
    if(any(duplicated(offspringNames(object)))){
      msg <- "offspring identifiers must uniquely identify a trio"
      return(msg)
    }
    if(any(is.na(unlist(trios(object))))){
      msg <- "Missing values not allowed in pedigree"
      return(msg)
    }
    if(!all(c(is(fatherNames(object), "character"),
              is(motherNames(object), "character"),
              is(sampleNames(object), "character")))){
      msg <- "sample identifiers must be character strings (e.g., not factors)"
      return(msg)
    }
    if(!all(originalNames(allNames(object)) %in% originalNames(unlist(trios(object))))){
      msg <- "all 'individualId' in slot pedigreeIndex must correspond to an id in the trio slot"
      return(msg)
    }
    if(any(fatherNames(object) == motherNames(object))){
      msg <- "fatherNames can not be the same as the motherNames"
      return(msg)
    }
    if(any(fatherNames(object) == sampleNames(object))){
      msg <- "fatherNames can not be the same as the offspringNames"
      return(msg)
    }
    if(any(motherNames(object) == sampleNames(object))){
      msg <- "motherNames can not be the same as the offspringNames"
      return(msg)
    }
  }
  return(msg)
}

setValidity("Pedigree", function(object){
	msg <- validPedigree(object)
	if(is.null(msg)) return(TRUE) else return(msg)
})

setMethod("initialize", signature(.Object="Pedigree"),
	  function(.Object,
		   trios=data.frame(F=character(),
		   M=character(),
		   O=character(), stringsAsFactors=FALSE),
		   trioIndex=data.frame(individualId=character(),
		   memberId=character(), index.in.pedigree=integer(), stringsAsFactors=FALSE),
		   ...){
		  callNextMethod(.Object, trios=trios, trioIndex=trioIndex, ...)
	  })

#' Deprecated function for constructing an instance of class Pedigree
#'
#' This function is deprecated and will be removed in a future release.
#'
#' @param pedigreeInfo a \code{data.frame} with column names 'F' (father), 'M' (mother), and 'O' (offspring). Elements of the \code{data.frame} are the sample names.
#' @param fatherIds character vector of identifiers for the father
#' @param motherIds character vector of identifiers for the mother
#' @param offspringIds character vector of identifiers for the offspring
#' @examples
#' Pedigree()
#' @export
Pedigree <- function(pedigreeInfo,
		     fatherIds=character(),
		     motherIds=character(),
		     offspringIds=character()){
  if(!missing(pedigreeInfo)){
    msg <- "pedigreeInfo must be a data.frame with column names 'F', 'M', and 'O'"
    if(!is(pedigreeInfo, "data.frame"))
      stop(msg)
    trios <- data.frame(F=as.character(pedigreeInfo[[1]]),
                        M=as.character(pedigreeInfo[[2]]),
                        O=as.character(pedigreeInfo[[3]]),
                        stringsAsFactors=FALSE)
    allIds <- as.character(unlist(trios))
  } else {
    fatherIds <- as.character(fatherIds)
    motherIds <- as.character(motherIds)
    offspringIds <- as.character(offspringIds)
    trios <- data.frame(F=fatherIds,
                        M=motherIds,
                        O=offspringIds,
                        stringsAsFactors=FALSE)
    allIds <- c(fatherIds, motherIds, offspringIds)
  }
  trio.index <- as.integer(matrix(seq_len(nrow(trios)), nrow(trios), 3, byrow=FALSE))
  memberId <- rep(c("F", "M", "O"), each=nrow(trios))
  pedigreeIndex <- data.frame(individualId=allIds,
                              memberId=memberId,
                              index.in.pedigree=trio.index,
                              stringsAsFactors=FALSE)
  rownames(pedigreeIndex) <- NULL
  new("Pedigree", trios=trios, trioIndex=pedigreeIndex)
}

#' @param object a \code{Pedigree} object
#' @aliases trios,Pedigree-method
#' @rdname Pedigree-class
setMethod("trios", signature(object="Pedigree"),
	  function(object) {
            ##.Deprecated()
            object@trios
        })
setMethod("trioIndex", signature(object="Pedigree"),
	  function(object) object@trioIndex)

#' @aliases offspringNames,Pedigree-method
#' @rdname Pedigree-class
setMethod("offspringNames", signature(object="Pedigree"), function(object) trios(object)$O)

setMethod("sampleNames", signature(object="Pedigree"), function(object) offspringNames(object))
setMethod("allNames", signature(object="Pedigree"), function(object) unique(trioIndex(object)$individualId))
setMethod("fatherNames", signature(object="Pedigree"), function(object) trios(object)$F)
setMethod("motherNames", signature(object="Pedigree"), function(object) trios(object)$M)

#' @aliases show,Pedigree-method
#' @rdname Pedigree-class
setMethod("show", signature(object="Pedigree"),
	  function(object){
		  ##cat('An object of class "Pedigree": \n')
		  ##cat('Father (F), Mother (M), and Offspring (O) \n')
		  ##cat('Slot "trios":\n')
		  print(head(trios(object), n=2))
		  if(nrow(trios(object)) > 2){
			  cat(".\n.\n.\n")
			  print(tail(trios(object), n=2))
		  }
		  cat('\nSlot "pedigreeIndex":\n')
		  print(head(trioIndex(object), n=2))
		  if(nrow(trioIndex(object)) > 6){
			  cat(".\n.\n.")
			  print(tail(trioIndex(object), n=2))
		  }
		  cat("\n")
                })

#' @param x a \code{Pedigree} object
#' @param i a numeric vector for subsetting  (optional)
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @aliases "[",Pedigree,ANY-method
#' @rdname Pedigree-class
setMethod("[", signature(x="Pedigree"),
	  function(x, i, j, ..., drop=FALSE){
            if(missing(i)) {
              return(x)
            } else {
              x@trios <- trios(x)[i, ]
              sns <- originalNames(unlist(trios(x)))
              x@trioIndex <- trioIndex(x)[match(sns, trioIndex(x)$individualId), ]
            }
            return(x)
	  })

#' @aliases dim,Pedigree-method
#' @rdname Pedigree-class
setMethod("dim", signature(x="Pedigree"), function(x){
	dim(trios(x))
})

setMethod("annotatedDataFrameFrom", signature(object="Pedigree", byrow="logical"),
	  function (object, byrow, sample.sheet, which=c("offspring", "father", "mother"),
		    row.names=NULL, ...){
		  dims <- dim(object)
		  if (is.null(dims) || all(dims == 0)){
			  return(annotatedDataFrameFrom(NULL, byrow = byrow, ...))
		  }
		  which <- match.arg(which)
		  nms <- switch(which,
				offspring=offspringNames(object),
				father=fatherNames(object),
				mother=motherNames(object))
		  nms <- make.unique2(nms)
		  if(missing(sample.sheet)){
			  n <- length(nms)
			  data <- data.frame(numeric(n), row.names = nms)[, FALSE]
			  dimLabels <-  c("sampleNames", "sampleColumns")
			  phenoData <- new("AnnotatedDataFrame", data = data, dimLabels = dimLabels)
		  } else {
			  if(is.null(row.names)){
				  stop("sample.sheet is not missing, but row.names not specified")
			  } else{
				  stopifnot(originalNames(nms) %in% row.names)
				  index <- match(originalNames(nms), row.names)
				  data <- sample.sheet[index, ]
				  rownames(data) <- nms
				  dimLabels <-  c("sampleNames", "sampleColumns")
				  phenoData <- new("AnnotatedDataFrame", data = data, dimLabels = dimLabels)
			  }
		  }
		  return(phenoData)
	  })
