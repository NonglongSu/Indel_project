nameList=Prank/nameFilter.txt

###########################################
ls Prank/mapped_cds_Out > $nameList

while read F; do cp Prank/updated_cds/$F Prank/mapped_cds_In_Out/; done < $nameList

cat $nameList
