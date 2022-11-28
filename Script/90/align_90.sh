MAFFT="mafft --globalpair --maxiterate 1000 --preservecase --quiet"
inDir=$1
ouDir=$2

# for subDir in ${inDir}/*
# do
# 	seed=$(basename "${subDir}")
# 	[ ! -d "${ouDir}/${seed}"] ||  mkdir -p ${ouDir}/${seed}
# 	for file in ${subDir}/*
# 	do
# 		filename=$(basename "${file}")
# 		${MAFFT} ${file} > ${ouDir}/${seed}/${filename}
# 	done
# done

[ ! -d "${ouDir}" ] && mkdir -p ${ouDir}
name=$(basename "$ouDir")
echo ${name}

for file in ${inDir}/*
do
	filename=$(basename "${file}")
	${MAFFT} ${file} > ${ouDir}/${filename}
	done
done
