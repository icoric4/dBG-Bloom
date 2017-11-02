./synthetic_raw_data

for i in {1..9} 
do
	./parse 'data'$i'.txt'
	rm 'data'$i'.txt'

	./jellyfish count -m 7 -s 100M -t 10 -C 'data'$i'.fasta' -o 'data'$i'.jf'
	./jellyfish dump 'data'$i'.jf' > 'data'$i'.fa'
	rm 'data'$i'.jf'
	./create_s 'data'$i'.fa'
	rm 'data'$i'.fa'

done