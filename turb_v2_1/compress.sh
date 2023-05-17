for dir in *256*M_0.9*hydro; 
do tar -czvf "$dir.tar.gz" "$dir/"; # >> output_2.txt 2>&1; 
done
