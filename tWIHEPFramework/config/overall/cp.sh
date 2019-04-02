for((k=0;k<1;k++)){
#	j=$[${k}*40]
#	echo bstar.config.BoostjetPt700.MET${j}
#	cp bstar.config bstar.config.MET${j}
#	sed -i "s/bstar.cuts/bstar.cuts.MET${j}/g" bstar.config.MET${j}
#	cp bstar.config.MET280.BoostjetPt700 bstar.config.BoostjetPt700.MET${j}
#	sed -i "s/bstar.cuts.MET280/bstar.cuts.MET${j}/g" bstar.config.BoostjetPt700.MET${j}
#	j=$[200+${k}*50]
#	echo bstar.config.MET160.BoostjetPt${j}
#	cp bstar.config.MET160 bstar.config.MET160.BoostjetPt${j}
#	sed -i "s/bstar.config/bstar.config.BoostjetPt${j}/g" bstar.config.MET160.BoostjetPt${j}
	j=$[${k}*10]
	echo bstar.config.MET280.BoostjetPt700.JetPt${j}
	cp bstar.config.MET280.BoostjetPt700 bstar.config.MET280.BoostjetPt700.JetPt${j}
	sed -i "s/bstar.config.BoostjetPt700/bstar.config.BoostjetPt700.JetPt${j}/g" bstar.config.MET280.BoostjetPt700.JetPt${j}
}
