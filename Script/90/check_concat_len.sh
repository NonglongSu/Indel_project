target=../../test_90_species/Raw_data/concat/*

for i in ${target}
do
	awk '/^>/ {if (seqlen) print seqlen; print; seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ${i}
done

#remove the > title
sed -i '/#>/d' haha.txt 
#divide by 3.
awk -v c=3 '{ print $0%c }' haha.txt
