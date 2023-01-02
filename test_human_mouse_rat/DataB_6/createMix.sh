#Set up vars
From1=Mafft/updated_cds
From2=Mafft/mapped_cds
To=Mafft/mixed_cds

#copy twice 
cp ${From1}/* to ${To}
cp ${From2}/* to ${To}

