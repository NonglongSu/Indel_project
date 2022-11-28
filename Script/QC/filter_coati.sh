#!/bin/bash
#filter out the coati_align via id.hex.txt file
ouF=Raw_data/QC/id.hex.txt
echo "Filter starts!"
cat ${ouF} | xargs -I % sh -c 'find Raw_data/coati_align/ -name %.fa -delete '
