% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify_variant.R
\name{simplify_variant}
\alias{simplify_variant}
\title{simplify_variant}
\usage{
simplify_variant(x)
}
\arguments{
\item{x}{A single row of a expanded VCF object (optionally) needing to be
simplified}
}
\value{
A VCF object of length 1, with the REF, ALT and rowname simplified
}
\description{
Removes extraneous bases from variants after mult-allelic variants have been
expanded in a VCF file. This function only fixes the REF and ALT fields of
the @geno slot as well as the rownames (containing the variant name in the
format chr:pos_ref/alt). The position component of the variant name is also
corrected appropriately.
}
\note{
As an example, the variant chr1:123456_GC/GA is simplified to chr1:123457_C/A
}
\references{
\url{http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/}
}
\author{
Alexander Hopkins
}
