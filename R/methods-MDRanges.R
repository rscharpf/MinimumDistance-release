.data_frame_posteriorSummaries <- function(g){
  L <- length(g)
  x <- matrix(NA, L, 4)
  colnames(x) <-  c("log_RR",
                    "log_odds",
                    "prob_MAP",
                    "prob_222")
  rownames(x) <- names(g)
  as.data.frame(x)
}

.mcols_posteriors <- function(g){
  L <- length(g)
  df <- DataFrame(log_RR=numeric(L),
                  log_odds=numeric(L),
                  prob_MAP=numeric(L),
                  prob_222=numeric(L),
                  prob_221=numeric(L))
}

.mdranges_mcols <- function(g){
  L <- length(g)
  df <- DataFrame(sample=character(),
                  numberProbes=integer(),
                  seg.mean=numeric(),
                  log_RR=numeric(L),
                  log_odds=numeric(L),
                  prob_MAP=numeric(L),
                  prob_222=numeric(L),
                  prob_221=numeric(L))
}

#' @param ... additional arguments to \code{GRanges} constructor
#' @param posteriors a \code{DataFrame}
#' @export
#' @rdname MDRanges-class
MDRanges <- function(..., posteriors){
  g <- GRanges(...)
  if(missing(posteriors)) posteriors <- .mdranges_mcols(g)
  g2 <- as(g, "MDRanges")
  mcols(g2) <- posteriors
  g2
}



setMethod("posteriorLogOdds", "MDRanges", function(object) object$log_odds)

posteriorOdds <- function(object) exp(posteriorLogOdds(object))
posteriorMAP <- function(object) object$prob_MAP
posteriorRR <- function(object) exp(posteriorLogRR(object))
posterior221 <- function(object) object$prob_221

setMethod("posteriorLogRR", "MDRanges", function(object) object$log_RR)

## # @export
isDenovo <- function(states) (states %in% c(duplicationStates(), deletionStates())) & !is.na(states)

setMethod("is221", "MDRanges", function(object)  object$calls=="221" & !is.na(object$calls))

setMethod("is221", "GRangesList", function(object) lapply(object, is221))

setMethod("numberFeatures", "MDRanges", function(object) object$number_probes)

setMethod("state", "MDRanges", function(object) object$calls)


.apply_ped_filters <- function(g, filters){
  if(length(g)==0) return(g)
  keep <- width(g) > width(filters)
  keep <- keep & state(g) %in% state(filters)
  keep <- keep & numberFeatures(g) >= numberFeatures(filters)
  keep <- keep & seqnames(g) %in% seqnames(filters)
  keep <- keep & g$prob_MAP>= probability(filters)
  g[keep]
}


setMethod("cnvFilter", "MDRanges", function(object, filters=FilterParam()){
  .apply_ped_filters(object, filters)
})


#' @param x a \code{MDRanges} object
#' @rdname MDRanges-class
#' @aliases reduce,MDRanges-method
#' @export
setMethod("reduce", "MDRanges", function(x, ...){
  g <- as(x, "GRanges")
  gr <- reduce(g, ...)
  hits <- findOverlaps(gr, g)
  nf <- split(numberFeatures(x), queryHits(hits))
  nf <- sapply(nf, sum)
  st <- split(state(x), queryHits(hits))
  st <- sapply(st, function(x) paste(unique(x), collapse=","))
  seg.mean <- mapply(function(x, w) sum(x)/sum(w), x=split(g$seg.mean * width(g), queryHits(hits)),
                     w=split(width(g), queryHits(hits)))



  ##  prCall <- sapply(split(x$prCall, queryHits(hits)), min)
  id <- sapply(split(x$sample, queryHits(hits)), unique)
  gr$numberFeatures <- nf
  gr$state <- st
  ##  gr$prCall <- prCall
  gr$sample <- id
  gr$seg.mean <- seg.mean
  gr
})
