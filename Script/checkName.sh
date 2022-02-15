
# add .fa to nameList.txt 
awk 'NF{print $0 ".fa"}' $1 > nameTmp.txt

while IFS= read -r file;do
	if ! [[ -e $2/$file ]]; then
        	printf '%s is missing in %s\n' "$file" "$2"
    	fi
done < "nameTmp.txt"

rm nameTmp.txt
