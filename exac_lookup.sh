#!/bin/bash
# Script to do single lookups from ExAC

#get the JSON for variant/variant
path=http://exac.hms.harvard.edu/rest/variant/variant/$2
json=`curl $path` 2>/dev/null

#parse the flags and extract the requested info 
while getopts cfp flag; do
	case $flag in
		c)
		#Return the Consequence of the first transcript
		echo $json | jq -r '.vep_annotations | .[0] | .Consequence'
		;;
		f)
		#Return the population allele frequency
		echo $json | jq -r '.allele_freq'
		#this could have been done without jq, but it is already used
		#echo$json | sed 's;.*\"allele_freq\":\([^,]*\).*;\1;'
		;;
		p)
		#The PolyPhen score looks like it may be helpful as well
		echo $json | jq -r '.vep_annotations | .[0] | .PolyPhen' |
	       		sed 's;.*(\(.*\));\1;' #this just reports the score
		;;
esac
done

