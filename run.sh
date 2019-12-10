echo $BASHPID; Rscript integrationMain.R 1 500 $1 $2&
echo $BASHPID;Rscript integrationMain.R 2 500 $1 $2 &
echo $BASHPID; Rscript integrationMain.R 3 500 $1 $2 &
echo $BASHPID; Rscript integrationMain.R 5 500 $1 $2 &
echo $BASHPID; Rscript integrationMain.R 10 500 $1 $2 &
echo $BASHPID; Rscript integrationMain.R 20 500 $1 $2 
# $1 gives the genz function - from 1 to 6
# $2 gives whether it is sequential or not - 0 or 1
