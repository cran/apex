% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multidna.class.R
\docType{class}
\name{multidna-class}
\alias{multidna-class}
\alias{multidna}
\alias{listOrNULL}
\alias{data.frameOrNULL}
\title{multidna: class for multiple gene data}
\description{
This formal (S4) class is used to store multiple DNA alignments.
Sequences are stored as a (possibly named) list, with each element of the list being a separate DNA alignment stored as a DNAbin matrix.
The rows of the separate matrices all correspond to the same individuals, ordered identically.
}
\section{Slots}{

\describe{
\item{\code{dna}}{a list of DNAbin matrices; empty slot should be NULL}

\item{\code{labels}}{a vector of labels of individuals}

\item{\code{n.ind}}{the number of individuals}

\item{\code{n.seq}}{the total number of sequences (pooling all genes), including gap sequences}

\item{\code{n.seq.miss}}{the total number of gap-only sequences}

\item{\code{ind.info}}{a data.frame containing information on the individuals, where individuals are in rows; empty slot should be NULL}

\item{\code{gene.info}}{a data.frame containing information on the genes, where genes are in rows; empty slot should be NULL}
}}

\examples{

## empty object
new("multidna")

## simple conversion with nicely ordered output
data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x
image(woodmouse)
image(x@dna[[1]])
image(x@dna[[2]])

## trickier conversion with missing sequences / wrong order
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
x <- new("multidna", genes)
x
image(x@dna[[1]])
image(x@dna[[2]])

}
\author{
Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
