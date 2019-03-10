for f in results_perth/*; do 
    div=$(echo $f | cut -f2 -d'.')
    p=$((100 - $div))
    res=$(cat "$f")
    echo $p","$res
done
