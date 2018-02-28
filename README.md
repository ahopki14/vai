# Van Andel Institute Coding Challenge
This package was created to extract the variants from a VCF file, annotate them
with data from the ExAC API, and export an easy to read table for investigators,
as well as a re-annotated VCF for downstream analyses. 

# Dependencies 
I tried to limit the dependencies, but didn't want to reinvent the wheel. This
package requires `VariantAnnotation` and `dplyr` as R dependencies, and needs 
`jq` to be installed on the system (this is in the Debian and Ubuntu Repos). 



