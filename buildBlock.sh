

#############################Data

dir0=Data_3
dir1=Data_6
dir2=Data_9
dir3=Data_12

sub_a=Mafft
sub_b=Prank
sub_c=Results

sub1=updated_cds
sub2=mapped_cds

mkdir -p $dir0/{Data/{$sub_a/{$sub1,$sub2},$sub_b/{$sub1,$sub2},$sub_c},Figure/{mafft_edge,prank_edge}}
mkdir -p $dir1/{Data/{$sub_a/{$sub1,$sub2},$sub_b/{$sub1,$sub2},$sub_c},Figure/{mafft_edge,prank_edge}}
mkdir -p $dir2/{Data/{$sub_a/{$sub1,$sub2},$sub_b/{$sub1,$sub2},$sub_c},Figure/{mafft_edge,prank_edge}}
mkdir -p $dir3/{Data/{$sub_a/{$sub1,$sub2},$sub_b/{$sub1,$sub2},$sub_c},Figure/{mafft_edge,prank_edge}}







#####################Raw_data 

#dir_a=Raw_data
#dir_b=Raw_data.2.outgroups


#mkdir -p $dir_a/{cds_seq,aa_seq,Alignments/{mafft_out,prank_out},mapped_cds_mafft,mapped_cds_prank,Figure}
#mkdir -p $dir_b/{cds_seq,aa_seq,Alignments/{mafft_out,prank_out},mapped_cds_mafft,mapped_cds_prank,Figure}


######################Script

#mkdir Script/{Test_pseudo_seq,Test_update}

#######################Other



