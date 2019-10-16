WD=$PWD

cd $WD/1_Static
wolf_sweep m1=0.001 m2=20 K1=0.01 K2=100 R2_rep_T=0.9 nSamples=1e9 > a1.txt

cd $WD/2_Expansion
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=1 > a1.txt
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=2 > a2.txt
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=3 > a3.txt
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=4 > a4.txt
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=5 > a5.txt
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=6 > a6.txt
wolf_expansion m1=0.001 m2=20 K1=0.01 K2=100 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 R2_rep_T=0.9 nSamples=1e9 firstDeme=7 > a7.txt

cd $WD/3_Bottleneck
wolf_bottleneck m1=0.001 m2=20 K1=0.01 K2=100 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 > a1.txt


cd $WD/4_Bottleneck_Expansion
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=1 > a1b.txt
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=2 > a2b.txt
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=3 > a3b.txt
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=4 > a4b.txt
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=5 > a5b.txt
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=6 > a6b.txt
wolf_expansion_bottleneck m1=0.001 m2=20 xBottle1=0 xBottle2=1 tStart1=5 tStart2=40 dT1=1e-3 dT2=1 t_ab=15 t_bc=40 R2_rep_T=0.9 nSamples=1e9 firstDeme=7 > a7b.txt

cd $WD
