
#' Read multiple DNA alignments
#'
#' These functions read multiple DNA alignments and store the output in a \linkS4class{multidna} object.
#' They are relying on ape's original functions \code{\link[ape]{read.dna}} and \code{\link[ape]{read.FASTA}}.
#'
#' @rdname readfiles
#' @aliases read.multidna
#' @aliases read.multiFASTA
#' @aliases read.multiphyDat
#'
#'
#' @param files a vector of characters indicating the paths to the files to read from
#' @param ... further arguments to be passed to the functions \code{\link[ape]{read.dna}} and \code{\link[ape]{read.FASTA}}.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#'
#' @seealso
#' \itemize{
#' \item \code{\link[ape]{read.dna}}
#' \item  \code{\link[ape]{read.FASTA}}
#' \item \code{\link[phangorn]{read.phyDat}}
#' }
#'
#' @export
#'
#' @examples
#' ## get path to the files
#' files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
#' files
#'
#' ## read files
#' x <- read.multiFASTA(files)
#' x
#' plot(x)
#'
#' y <- read.multiphyDat(files, format="fasta")
#' y
#'
read.multidna <- function(files, ...){
    gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
    dna <- lapply(files, read.dna, ...)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna)
    return(out)
}


#'
#' @rdname readfiles
#' @export
read.multiFASTA <- function(files){
    gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
    dna <- lapply(files, read.FASTA)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna)
    return(out)
}


#'
#' @rdname readfiles
#' @export
read.multiphyDat <- function(files, ...){
  gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
  dna <- lapply(files, read.phyDat, ...)
  names(dna) <- gene.names
  out <- new("multiphyDat", dna=dna)
  return(out)
}