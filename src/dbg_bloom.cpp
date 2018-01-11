#include "dbg_bloom.hpp"

/*
 * Constructor of dbg class.
 * Argument name represents the name of the file which contains solid k-mers.
 * Set S which holds all solid k-mers (true positives) is olso created.
 */
dbg::dbg(string name) {
    ifstream ifs(name);
    string s;
    while (ifs >> s) {
        S.insert(s);
    }
    vector<bool> bla(S.size(), false);
    k_size = (*S.begin()).size();
    create_BF();
}

/*
 * Method for creation of bloom filter and adding all solid k-mers in the filter.
 */
void dbg::create_BF() {
    bloom_parameters parameters;
    parameters.projected_element_count = S.size();
    parameters.false_positive_probability = 0.0001;
    parameters.random_seed = 0xA5A5A5A5;
    parameters.compute_optimal_parameters();
    filter = bloom_filter(parameters);
    for (set<string>::iterator it = S.begin(); it != S.end(); ++it)	{
        filter.insert(*it);
        filter.insert(reverse_complement(*it));
    }
}

/*
 * Method for calculation of set P which represents the set of all extensions of
 * set S for which the Bloom filter answers yes.
 */
void dbg::compute_P() {
    string h[] = {"A","C","G","T"};
    set<string> nucleotides(h, h+4);
    for (set<string>::iterator it = S.begin(); it != S.end(); ++it)	{
        string prefix = (*it).substr(0, (*it).size()-1);
        string sufix = (*it).substr(1);
        for (set<string>::iterator nucl = nucleotides.begin(); nucl != nucleotides.end(); ++nucl) {

            if (filter.contains(sufix + *nucl)) {
                P.insert(canonical(sufix + *nucl));
            }

            if (filter.contains(*nucl + prefix)) {
                P.insert(canonical(*nucl + prefix));
            }
        }
    }
    //cout << S.size() << " " << P.size() << endl;
}

/*
 * Method for computation of set cFP which is
 * formally defined as cFP = P\S
 */
void dbg::compute_cFP(unsigned int M) {
    compute_P();
    set<string> D(P);
    set<string>::iterator it = S.begin();
    //int i = 0;
    while(it != S.end()) {
        set<string> Pi;
        set<string> Dn;
        //cout << i << "/" << S.size() << endl;
        while (Pi.size() < M && it != S.end()) {
            Pi.insert(*it);
            it++;
            //	++i;
        }

        for (set<string>::iterator m = D.begin(); m != D.end(); ++m) {
            if (Pi.count(canonical(*m)) == 0) {
                Dn.insert(canonical(*m));
            } //else delete
        }

        D = set<string>(Dn);
    }

    cFP = set<string>(D);
    //cout << cFP.size() << endl;
    vector<string> tmp(S.begin(), S.end());
    nodes = tmp;
}

/*
 * Method for graph traversal. It starts from a solid k-mer
 * and, using Bloom filter and cFP structure, looks for its
 * neighbours using bounded-breadth, bounded-depth BFS. Method
 * takes the name of the file in which contigs will be saved as
 * an argument.
 */
void dbg::traverse_graph(string name) {

    int index = 0;

    while (index < nodes.size()) {
        int depth = 0;
        int reason = -1;
        for (int tt = 0; tt != 2; ++tt) { // tt=0 go only left from nodes[index], tt=1 only right 
            set<string> front;
            front.insert(nodes[index]);
            bool l = !tt, r = tt;

            do {
                set<string> new_front;

                for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
                    set<string> ret = branch(*it, l, r);
                    new_front.insert(ret.begin(), ret.end());
                    if (r) {
                        marked.insert((*it).substr((*it).size()-k_size));
                        marked.insert(reverse_complement((*it).substr((*it).size()-k_size)));
                    }
                    if (l) {
                        marked.insert((*it).substr(0,k_size));
                        marked.insert(reverse_complement((*it).substr(0,k_size)));
                    }
                }

                if (new_front.empty()) {
                    reason = 2;
                    break;
                }
                depth++;

                if (depth > 1e6) {
                    reason = 0;
                    break;
                }
                front = new_front;
                if (new_front.size() == 1) {
                    set<string> ret = branch(*(new_front.begin()), l, r);
                    if (ret.size() >= 1) {
                        continue;
                    }
                    reason = 1;
                    break;
                }
            } while (1);
            if (reason != -1) {
                for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
                    //cout << "contig: " << *it << " " <<reverse_complement(*it) <<  endl;

                    if ((*it).size() > 2 * k_size + 1) {
                        //  cout << "contig: " << *it << " " <<reverse_complement(*it) <<  endl;
                        contigs.insert(*it);
                    }
                }
            }
        }

        index += 1;

        for (int i = index; i < nodes.size(); ++i) {
            if(!marked.count(nodes[i])) {
                index = i;
                break;
            }
            if (i == nodes.size()-1) {
                index = nodes.size();
            }
        }
    }

    ofstream myfile;
    myfile.open (name);
    for (set<string>::iterator it = contigs.begin(); it != contigs.end(); ++it) {
        myfile << *it << endl;
    }
    myfile.close();
}

/*
 * Method that takes a contig s and tries to expand it
 * to the left and right if corresponding flags: left and right
 * are true
 */
set<string> dbg::branch(string s, bool left, bool right) {

    char bla[4] = {'A', 'C', 'G', 'T'};
    list<char> l(bla, bla + sizeof(bla) / sizeof(char));
    set<string> ret;
    if (left) {
        for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
            string h = *c + s.substr(0, k_size - 1);
            if (left && filter.contains(h) && !cFP.count(canonical(h)) && !marked.count(h)) {
                ret.insert(*c + s);
            }
        }
    }
    if (right) {
        for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
            string h = s.substr(s.size() - k_size + 1) + *c;
            if (right && filter.contains(h) && !cFP.count(canonical(h)) && !marked.count(h)) {
                ret.insert(s + *c);
            }
        }
    }

    return ret;
}

/*
 * Method for calcuating caonical k-mer.
 * K-mer is canonical if it is less than
 * its reverse complement in means of ASCII.
 */
string dbg::canonical (string s1) {
    string s2 = reverse_complement(s1);
    return s1 < s2 ? s1 : s2;
}

/*
 * Method for calculation reverse complement of a k-mer.
 */
string dbg::reverse_complement (string s) {
    string tmp = "";
    for(string::iterator it = s.begin(); it != s.end(); ++it) {
        tmp = complement(string(1,*it)) + tmp;
    }
    return tmp;
}

/*
 * Method for calculation of a complement of a nucleotid.
 */
string dbg::complement(string c) {
    if (c == "C")
        return "G";
    if (c == "T")
        return "A";
    if (c == "A")
        return "T";
    if (c == "G")
        return "C";
    return "";
}