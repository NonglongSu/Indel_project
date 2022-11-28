#Remove all contained gaps from zou's data

Dir1=$1

for subdir in ${Dir1}/*
do
	echo $(basename ${subdir})
	for files in ${subdir}/*
	do
		sed -i 's/-//g' ${files}
	done
done
