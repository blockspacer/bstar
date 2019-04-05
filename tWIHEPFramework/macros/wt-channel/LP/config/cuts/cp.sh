for((k=1;k<11;k++)){
	j=$[${k}*40]
	echo bstar.cuts.MET${j}
	cp bstar.cuts bstar.cuts.MET${j}
	sed -i "s/MissingEt.Min:  0.0/MissingEt.Min:  ${j}/g" bstar.cuts.MET${j} 
}
