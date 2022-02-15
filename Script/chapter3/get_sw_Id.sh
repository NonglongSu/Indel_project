inDir1=$1
inDir2=$2
ouFile=$3
ouDir=$4

#grab the sw list
ls ${inDir1} > ${ouFile}

#copy the list-files from raw data
rsync -a ${inDir2} --files-from=${ouFile} ${ouDir}
