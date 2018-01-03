#include "dbg_bloom.hpp"





int main (int argc, char *argv[]) {

	dbg d(argv[1]);
	d.compute_cFP(100);
	d.traverse_graph();
	cout << d.contigs.size() << endl;
	for (auto s : d.contigs) {
		cout << s << " " << d.reverse_complement(s) << endl;
	}

}