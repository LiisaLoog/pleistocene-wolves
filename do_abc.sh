WD = $PWD

############################################## Static ##############################################

cd $WD/1_Static/ABC

rm -f ../ABC_RESULTS_STATIC.txt
touch ../ABC_RESULTS_STATIC.txt

# run ABCestimator
../../ABCestimator ../params.txt simName=abc.out_static_1.0.txt > result_static_1.0.txt

# extract marginal posterior probability density and number of accepted simulations
echo ABC1 `grep -r 'marginal density:' ./result_*.txt` `wc -l abc.out_*` >> ../ABC_RESULTS_STATIC.txt



############################################## Expansion ##############################################

cd $WD/2_Expansion

rm -f ../ABC_RESULTS_EXPANSION.txt
touch ../ABC_RESULTS_EXPANSION.txt

for i in {1..7}
do
	echo $i
	cd ABC$i

	# run ABCestimator
	../../ABCestimator ../params.txt simName=abc.out_$i'_expansion_1.0.txt' > result_$i'_expansion_1.0.txt'

	# extract marginal posterior probability density and number of accepted simulations
	echo ABC$i `grep -r 'marginal density:' ./result_*.txt` `wc -l abc.out_*` >> ../ABC_RESULTS_EXPANSION.txt

	cd ..
done

############################################## Bottleneck ##############################################
cd $WD/3_Bottleneck/ABC

rm -f ../ABC_RESULTS_BOTTLENECK.txt
touch ../ABC_RESULTS_BOTTLENECK.txt

# run ABCestimator
../../ABCestimator ../params.txt simName=abc.out_$i'_bottelneck_1.0.txt' > result_$i'_expansion_1.0.txt'

# extract marginal posterior probability density and number of accepted simulations
echo ABC1 `grep -r 'marginal density:' ./result_*.txt` `wc -l abc.out*` >> ../ABC_RESULTS_BOTTLENECK.txt


############################################## Bottleneck Expansion #####################################

cd $WD/4_Bottleneck_Expansion

rm -f ../ABC_RESULTS_BOTTLENECK_EXPANSION.txt
touch ../ABC_RESULTS_BOTTLENECK_EXPANSION.txt

for i in {1..7}
do
	echo $i
	cd ABC$i

	# run ABCestimator
	../../ABCestimator ../params.txt simName=abc.out_$i'_bottelneck_expansion_1.0.txt' > result_$i'_bottelneck_expansion_1.0.txt'

	# extract marginal posterior probability density and number of accepted simulations
	echo ABC$i `grep -r 'marginal density:' ./result_*.txt` `wc -l abc.out_*` >> ../ABC_RESULTS_BOTTLENECK_EXPANSION.txt

	cd ..
done
