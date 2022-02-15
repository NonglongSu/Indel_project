echo "Remember to switch the target directories!!!!"

#Generate a nameFilter.txt

#ls Mafft/mapped_cds_Out > Mafft/nameFilter.txt
#ls Mafft_Mul/mapped_cds_Out > Mafft_Mul/nameFilter.txt
#ls Mafft/mapped_cds_Out_cD > Mafft/nameFilter.txt


ls ../../Raw_data.2.outgroups/mafft_sim/sw_gaps > ../../Raw_data.2.outgroups/mafft_sim/nameFilter.txt


#########################################################

#Create mappe_cds_IN_Out files
#while read F;do  cp Mafft/updated_cds/$F Mafft/mapped_cds_In_Out/; done < Mafft/nameFilter.txt

#while read F;do  cp Mafft/updated_cds_cD/$F Mafft/mapped_cds_In_Out_cD/; done < Mafft/nameFilter.txt

while read F;do  cp ../../Raw_data.2.outgroups/mafft_sim/updated_gaps/$F ../../Raw_data.2.outgroups/mafft_sim/mafft_In/; done < ../../Raw_data.2.outgroups/mafft_sim/nameFilter.txt

