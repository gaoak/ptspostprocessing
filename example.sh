process num 65,65,65 range -0.2,1.2,-0.2,1,0,5.3 dump rawpoint${n}
$NEKBIN/FieldConvert -f -e -v -m interppoints:fromxml=wing.xml,wingc.xml:fromfld=../../wing_${n}.chk:topts=rawpoint${n}.csv data${n}.plt

process num 65,65,65 range -0.2,1.2,-0.2,1,0,5.3 load data${n} process 0.03 output sdata${n}
echo "job finished"