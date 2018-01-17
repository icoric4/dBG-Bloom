# dBG-Bloom

Space-efficient and exact de Bruijn graph representation based on a Bloom filter.

This program was created as a project assignment for Bioinformatics class at Faculty of Electrical Engineering and Computing, University of Zagreb (http://www.fer.unizg.hr/predmet/bio).

## Instructions

    # get a local copy of dBG-Bloom source code
    git clone --recursive https://github.com/icoric4/dBG-Bloom.git
    
    # compile the code on your computer
    cd dBG-Bloom
    sh INSTALL


# User manual	 

Type `./dbg --help` for usage instructions.

# Testing
   
    # to run small test execute following commands
    cd bin
    ./dbg -in ../test/test_small.fasta -abundance-min 1
    
    # to run tests on real data execute following commands
    cd test
    unzip SRR6472718.fasta.zip
    cd ..
    cd bin
    ./dbg -in ../test/SRR6472718.fasta
    
    # to run test on synthetic data run test.sh script
    # the script creates synthetic data, calculates k-mers using jellyfish program, and tries to reconstruct original data
    cd test
    ./test.sh
