for i in `ls */do_*`
do
        echo $i
        sed -f change.cmd $i > tmp
        mv tmp $i
done
