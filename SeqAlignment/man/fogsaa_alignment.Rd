% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fogsaa_alignment.R
\name{fogsaa_alignment}
\alias{fogsaa_alignment}
\title{FOGSAA algorithm for global sequence alignment.}
\usage{
fogsaa_alignment(s1, s2, match, mismatch, gap)
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
Computes the optimal global alignment of two sequences using the FOGSAA
algorithm.
}
\examples{
fogsaa_alignment('ACGGTTGC','AGCGTC', 1, -1, -2)
}
\references{
\url{https://pubmed.ncbi.nlm.nih.gov/23624407/}\cr
}
