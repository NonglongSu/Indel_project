Target=$1

for subDir in ${Target}*
do
	echo "$(basename ${subDir})"
	ls ${subDir} | wc -l
done

