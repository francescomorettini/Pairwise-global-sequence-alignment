% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nw_alignment.R
\name{nw_alignment}
\alias{nw_alignment}
\title{Needlman-Wunsch global sequence alignment algorithm.}
\usage{
nw_alignment(s1, s2, match, mismatch, gap)
}
\arguments{
\item{s1}{A string representing sequence 1.}

\item{s2}{A string representing sequence 2.}

\item{match}{An integer which indicates the score for a match between
the two sequences.}

\item{mismatch}{An integer which indicates the score for a mismatch
between the two sequences.}

\item{gap}{An integer which indicates the score for a gap inserted
in one of the two sequences.}
}
\value{
A list containing the alignment of \code{s1} and \code{s2}
}
\description{
Computes the global sequence alignment using the Needleman-Wunsh algorithm.
}
\examples{

nw_alignment('ACGGTTGC','AGCGTC', 1, -1, -2)

}
\references{
\url{https://pubmed.ncbi.nlm.nih.gov/5420325/}\cr
}
