for dir in $2/*
do
	goodn=0;

	for file in $dir/*.png
	do
		if [ "$(./readsym $1 "$file")" = "$(basename "$dir")" ]
		then
			((++goodn))
		fi
	done
	
	echo "$(basename "$dir"): $goodn / $(ls "$dir" | wc -l)"
done
