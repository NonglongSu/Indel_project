MAFFT="ginsi --preservecase"
inDir=$1
ouDir=$2

for subDir in ${inDir}/*
do
	seed=$(basename "${subDir}")
	[! -d "${ouDir}/${seed}"] ||  mkdir -p ${ouDir}/${seed}
	for file in ${subDir}/*
	do
		filename=$(basename "${file}")
		${MAFFT} ${file} > ${ouDir}/${seed}/${filename}
	done
done


