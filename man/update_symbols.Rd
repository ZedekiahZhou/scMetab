% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getUpdatedSymbols.R
\name{update_symbols}
\alias{update_symbols}
\title{Get Rid of Outdated Gene Symbols}
\usage{
update_symbols(gene_symbols)
}
\arguments{
\item{gene_symbols}{input gene symbols to update, duplicates and NAs are allowed}
}
\value{
a translation table with three variables:
\itemize{
\item original: original symbols of the same length and order as input
\item updated: updated symbols:
\item type: symbols was matched to approved symbols, previous symbols, or deleted
}
}
\description{
Update previous gene symbols to up-to-date approved gene symbols, return a translation table
}
\details{
Gene symbols are updated according to "hgnc_20210601.txt" file download from HGNC database (custom downloads).
Previous symbols are updated, approved symbols are preserved, unmatched symbols are also preserved.

Input symbols are matched against approved symbols first, unmatched ones are then matched against previous
symbols. Symbols left unmatched finally were also preserved, but matched status were recorded.
Alias symbols were not under consideration to reduce multi-mapping which is confused.

If a previous symbol match to multiple approved symbols, or multiple previous symbols match to an approved
symbol, the first match will be kept according to the alphabetical order of approved symbols.
}
\examples{
update_symbols(c("KRAS2", "DUSP2", NA))
}
