./test_generator

for i in {1..7} 
do
	./parse 'data'$i'.txt'
	rm 'data'$i'.txt'

	./../bin/jellyfish count -m 31 -s 100M -t 10 -C -L 1 'data'$i'.fasta' -o 'data'$i'.jf'
	./../bin/jellyfish dump 'data'$i'.jf' > 'data'$i'.fa'
	rm 'data'$i'.jf'
	./create_s 'data'$i'.fa'
	rm 'data'$i'.fa'
	echo $i'. TEST RESULT: '
	./test 'data'$i'S.txt' 'data'$i'.fasta'
	rm 'contigs.txt'
	rm 'data'$i'S.txt'
	rm 'data'$i'.fasta'

done