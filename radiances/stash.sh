DIR=./backups/tmp
if [ ! -e $DIR ]; then mkdir $DIR; fi
mv job_* $DIR
mv Edir* $DIR
mv radiances* $DIR
mv input_params.txt $DIR
mv numus.txt $DIR
mv nphis.txt $DIR
mv wc3D.dat $DIR

