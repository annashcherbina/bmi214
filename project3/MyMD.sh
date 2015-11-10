#python MyMD.py --iF test.rvc #1FW4_cold.rvc --kB 40000 --kN 400.0 --nbCutoff 0.50 --m 12.0 --dt 0.001 --n 1000

#python MyMD.py --iF 1FW4_med.rvc --kB 40000 --kN 400.0 --nbCutoff 0.50 --m 12.0 --dt 0.001 --n 1000
#python MyMD.py --iF 1FW4_cold.rvc # --kB 40000 --kN 400.0 --nbCutoff 0.50 --m 12.0 --dt 0.001 --n 1000
#python MyMD.py --iF 1FW4_hot.rvc 

#Weak interactions
#python MyMD.py --iF 1FW4_med.rvc --kB 1000.0 --kN 100.0 --out kB1000kN100
#python convertRVCtoCRD.py kB1000kN100_out.rvc kB1000kN100_out.crd 

#echo "Near Interactions"
#python MyMD.py --iF 1FW4_med.rvc --nbCutoff 0.25  --out nb25
#python convertRVCtoCRD.py nb25_out.rvc nb25_out.crd 

#echo "Far Interactions"
#python MyMD.py --iF 1FW4_med.rvc --nbCutoff 0.75  --out nb75
#python convertRVCtoCRD.py nb75_out.rvc nb75_out.crd 

#echo "Longer Time Steps"
#python MyMD.py --iF 1FW4_med.rvc --dt 0.01  --out t01
#python convertRVCtoCRD.py t01_out.rvc t01_out.crd 

#echo "Much Longer Time Steps"
python MyMD.py --iF 1FW4_med.rvc --dt 0.05  --out t05
python convertRVCtoCRD.py t05_out.rvc t05_out.crd 

#echo "Very Long Time Steps"
#python MyMD.py --iF 1FW4_med.rvc --dt .1  --out t1
#python convertRVCtoCRD.py t1_out.rvc t1_out.crd 

#echo "No Calcium"
#python MyMD.py --iF 1FW4_noCa_med.rvc
#python convertRVCtoCRD.py 1FW4_noCa_med_out.rvc 1FW4_noCa_med_out.crd
#python computeFEATUREScore1FW4.py 1FW4_noCa_med_out.rvc 1FW4_noCa_med_out.features

#echo "Mutated"
#python MyMD.py --iF 1FW4_MUT_med.rvc
#python convertRVCtoCRD.py 1FW4_MUT_med_out.rvc 1FW4_MUT_med_out.crd
#python computeFEATUREScore1FW4.py 1FW4_MUT_med_out.rvc 1FW4_MUT_med_out.features
