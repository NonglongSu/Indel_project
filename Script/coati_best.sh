dir1=$1
dir2=$2

for file in ${dir1}/*
do
  seed=$(basename "${file}")
  #avoid repeated running

  if test -f "${dir2}/${seed}"; then
    continue
  else
    echo ${seed}
  fi
  ~/coati/builddir/src/coati-alignpair ${file} -o ${dir2}/${seed} -k 3
done
