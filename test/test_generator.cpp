#include <iostream>
#include <fstream>
#include <string> 
#include <random>

using namespace std;


int main(void) {

	mt19937 rng;
    rng.seed(random_device()());
    uniform_int_distribution<mt19937::result_type> dist(0,3);
    char alphabet[4] = {'A','C','T','G'};

	for(int i=1e2, b=1, j=1; i <= 1e6; i = (b) ? (i*5) : (i/5*10), b=!b, ++j) {
		
		ofstream f;
		if (j == 6) {j++;}
		f.open ("data"+to_string(j)+".txt");
		for (int k = 0; k != i; ++k) {
			char c = alphabet[dist(rng)];
			f << c;
		}
		f.close();
	}

	ofstream f;
		f.open ("data"+to_string(6)+".txt");
		for (int k = 0; k != 25000; ++k) {
			char c = alphabet[dist(rng)];
			f << c;
		}
		f.close();
}


