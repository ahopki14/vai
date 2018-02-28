# Van Andel Institute Coding Challenge
This package was created to extract the variants from a VCF file, annotate them
with data from the ExAC API, and export an easy to read table for investigators,
as well as a re-annotated VCF for downstream analyses. 

## Dependencies 
I tried to limit the dependencies, but didn't want to reinvent the wheel. This
package requires `VariantAnnotation` and `dplyr` as R dependencies, and needs 
`jq` to be installed on the system (this is in the Debian and Ubuntu Repos). 

## What is included
`scripts/vai.R` provides all the code needed to generate `out/table.csv` which
contains a simple table of relevant info for all variants (for the
investigator), as well as `out/vai_exac_annotated.vcf` which is an annotated VCF
file for use in other analyses which might need a VCF. 
I packaged the input VCF file in `inst/extdata` so that it can be called by the
script portably. I also put the output table and VCF in `out`. 

## Installation and use
I recommend installing with `devtools::install_github('ahopki14/vai')`. Then
`vai.R` can be run and should have everything it needs. 

## Future Directions
* Re-write the exac lookup using the API's batch functionality
  * I used the individual lookup methods in the API because I thought it
would be nice to have a function that did a single query
  * This works fine, but is a little slow (I usually complete the ~7000
queries in about 8 minutes)
* Include info from other sources
  * The big thing that is missing is the gene names for variants that are
not annotated in the ExAC database, especially insertions/deletions which are
unlikely to be shared by other individuals.  
* I tried to write this to be portable, but I have not tested it on Windows
  systems 


## What are the columns in `table.csv`?
Header        | Description
--------------|------------
Type          | Type of variant, i.e. snp, ins, _etc._
AlleleFreq    | Allele Frequency reported in the VCF
ReadDepth     | Total reads covering the variant
GenoType      | 0/1,0/2, _etc._=Hetrozygous, 1/1=Homozygous    
AltReads      | Number of reads supporting the Alternate
RefReads      | Number of reads supporting the Reference
PctAlt        | Percentage of reads supporting Alternate
VariantNames  | chr-pos-ref-alt, Used to query the ExAC database
Consequence   | The reported Major Consequence from ExAC
PopAlleleFreq | The frequency of the variant in the ExAC data
PolyPhen      | Score summarizing the impact of the mutation (0=benign, 1=major)
GeneName      | If the variant is in ExAC, the gene SYMBOL of the gene containing the variant




