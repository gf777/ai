iQV=40
igsize=4657790
ictgn=1

awk -v iQV="$iQV" -v igsize="$igsize" -v ictgn="$ictgn"  '{cost=( ( $3 - iQV )^2 / iQV^2) + (( $6 - ictgn )^2 / ictgn^2) + (($4 - igsize )^2 / igsize^2); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"cost}' $1

