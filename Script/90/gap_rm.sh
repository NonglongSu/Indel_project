#Remove all contained gaps from zou's data

Dir=$1

for subdir in ${Dir}/* 
do
	for files in ${subdir}/*
	do
		sed -i 's/-//g' ${files}
	done
done



