for((n=48; n<=47; ++n))
do
process num 65,65,65 range -0.2,1.2,-0.2,1,0,5.3 dump rawpoint${n}
$NEKBIN/FieldConvert -f -e -v -m interppoints:fromxml=wing.xml,wingc.xml:fromfld=../../wing_${n}.chk:topts=rawpoint${n}.csv data${n}.plt
done
rm .swap.center
for((n=52; n<=52; n+=1))
do
process num 129,129,129 range 0.2,1.1,-0.2,0.6,0,5 load data${n} process 0.03 output sdata${n}
echo "job finished"
done