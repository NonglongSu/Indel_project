
while read -r line
do
  echo $line
  pars=$line
done < tmp.txt

file=Raw_data/cds/01_FcaCaf_aligned_cds/ENSG00000000003_TSPAN6_raw_NT.fasta

haha="~/coati/builddir/src/coati-alignpair ${file} -o haha.fasta -k 3 ${pars}"

echo  "~/coati/builddir/src/coati-alignpair ${file} -o haha.fasta -k 3 ${pars}"
~/coati/builddir/src/coati-alignpair ${file} -o haha.fasta -k 3 ${pars}
