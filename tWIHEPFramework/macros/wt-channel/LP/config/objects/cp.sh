for((k=0;k<1;k++)){
	#j=$[400+${k}*50]
	#echo bstar.config.BoostjetPt${j}
	#cp bstar.config bstar.config.BoostjetPt${j}
	#sed -i "s/BoostedJet.MinPt: 200.0/BoostedJet.MinPt: ${j}/g" bstar.config.BoostjetPt${j}
	j=$[${k}*10]
	echo bstar.config.BoostjetPt700.JetPt${j}
	cp bstar.config.BoostjetPt700 bstar.config.BoostjetPt700.JetPt${j}
	sed -i "s/Jet.MinPt: 30.0/Jet.MinPt: ${j}/g" bstar.config.BoostjetPt700.JetPt${j}
}
