process num 65,65,65 range 0.2,1.1,-0.2,0.6,0,5 dump small_rawpoint.csv
$NEKBIN/FieldConvert -f -e -v -m interppoints:fromxml=wing.xml,wingc.xml:fromfld=../../wing_32.chk:topts=small_rawpoint.csv small_data.csv
process num 65,65,65 range 0.2,1.1,-0.2,0.6,0,5 load small_data.csv output_def small_rawdata.plt
process num 65,65,65 range 0.2,1.1,-0.2,0.6,0,5 load small_data.csv dosmooth 0.03 output_def small_sdata.plt
echo "job finished"
