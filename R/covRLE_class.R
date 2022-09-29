#' Coverage Rle for both strands or single
#'
#' Given a run of coverage(x) where x are reads,
#' this class combines the 2 strands into 1 object
#' @importFrom methods new
#' @export
#' @family ORFik_experiment
setClass("covRle",
         contains="RleList",
         representation(
           forward="RleList",          # of length N, no names
           reverse="RleList",           # of length N, no names
           strandMode="integer"
         ),
         prototype(
           forward="RleList",
           reverse="RleList",
           strandMode=integer()
         )
)
#' Coverage Rlelist for both strands
#' @export
#' @examples
#' covRle(RleList(), RleList())
covRle <- function(forward, reverse = RleList()) {
  strandMode <- ifelse(length(reverse) > 0, 1L, 0L)
  if (strandMode) {
    if (length(forward) != length(reverse)) stop("Length of forward and reverse must match")
  }
  new("covRle", forward = forward, reverse = reverse, strandMode = strandMode)
}

#' @export
setMethod("show", "covRle",
          function(object)
            print(c(forward = object@forward, reverse = object@reverse))
)

#' Seqlevels covRle
#' Extracted from forward RleList
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("seqlevels",
          "covRle",
          function(x) {
            seqlevels(x@forward)
          }
)

#' Seqinfo covRle
#' Extracted from forward RleList
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("seqinfo",
          "covRle",
          function(x) {
            seqinfo(x@forward)
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("strandMode",
          "covRle",
          function(x) {
            x@strandMode
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return the forward RleList
#' @export
setGeneric("f", function(x) standardGeneric("f"))


#' @inherit f
setMethod("f",
          "covRle",
          function(x) {
            x@forward
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return the forward RleList
#' @export
setGeneric("r", function(x) standardGeneric("r"))

#' @inherit r
setMethod("r",
          "covRle",
          function(x) {
            x@reverse
          }
)


