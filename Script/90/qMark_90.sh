Dir=$1

for subdir in ${Dir}/* 
do
	for files in ${subdir}/*
	do
		sed -i 's/?/N/g' ${files}
	done
done
