# Van Andel Institute Coding Challenge
This package was created to extract the variants from a VCF file, annotate them
with data from the ExAC API, and export an easy to read table for investigators,
as well as a re-annotated VCF for downstream analyses. 

# Dependencies 
I tried to limit the dependencies, but didn't want to reinvent the wheel. This
package requires `VariantAnnotation` and `dplyr` as R dependencies, and needs 
`jq` to be installed on the system (this is in the Debian and Ubuntu Repos). 

# What is included
`scripts/vai.R` provides all the code needed to generate `out/table.csv` which
contains a simple table of relevant info for all variants (for the
investigator), as well as `out/vai_exac_annotated.vcf` which is an annotated VCF
file for use in other analyses which might need a VCF. 

# Future Directions
* Re-write the exac lookup using the API's batch functionality
  * I used the individual lookup methods in the API because I thought it
would be nice to have a function that did a single query
  * This works fine, but is a little slow (I usually complete the ~7000
queries in about 8 minutes)
* Include info from other sources
  * The big thing that is missing is the gene names for variants that are
not annotated in the ExAC database, expecially insertions/deletions which are
unlikely to be shared by other individuals.  
