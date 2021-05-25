list_csv=""
list_name="gene"
for i in  ./*/*genes.results
do
	echo $i
	list_csv=${list_csv}" "${i}
	i_name1=${i%/*}
	i_name2=${i_name1##*/}
	list_name=${list_name}"       "${i_name2}
	echo ${list_name}
done
echo $list_name > ./gene-tpm-matrix.txt
python ./rsem-count-extract.py ${list_csv} >> ./gene-tpm-matrix.txt