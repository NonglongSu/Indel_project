#!/bin/bash
#create a final filter file
ouF=QC/final.filter.txt
cat $1 $2 $3 $4 | sort | uniq > ${ouF}

#delete the bad file
echo "Filter starts!"
cat ${ouF} | xargs -I % sh -c 'find cds_seq/ -name %.fa -delete '
