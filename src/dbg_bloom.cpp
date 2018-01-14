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
        S.push_back(s);
    }
    k_size = (*S.begin()).size();
    create_BF();
}

/*
 * Method for creation of bloom filter and adding all solid k-mers in the filter.
 */
void dbg::create_BF() {
    const unsigned k = k_size;
    const unsigned numHashes = 10;
    const unsigned size = 1e9;
    filter = KmerBloomFilter(size, numHashes, k);
    for (int i = 0; i != S.size(); ++i)	{
        filter.insert(S[i].c_str());
        filter.insert(reverse_complement(S[i]).c_str());
    }
}

/*
 * Method for calculation of set P which represents the set of all extensions of
 * set S for which the Bloom filter answers yes.
 */
vector<string> dbg::compute_P() {
    string h[] = {"A","C","G","T"};
    char hc[] = {'A', 'C', 'G', 'T'};
    unordered_set<string> P;
    stringstream ss;
    char end; 
    string start;
    for (int i = 0; i != S.size(); ++i) { // a bit uglier code to avoid usage of string::operator+ to save time
        end = S[i][(S[i].size()-1)];
        ss << S[i][0];
        ss >> start;
        for (int j = 0; j != 4; ++j) {
            S[i].pop_back();
            S[i].insert(0,h[j]);
            if (filter.contains(S[i].c_str())) {
                P.insert(canonical(S[i]));
            }
            S[i].erase(S[i].begin());
            S[i].push_back(end);

            S[i].erase(S[i].begin());
            S[i].push_back(hc[j]);
            if (filter.contains(S[i].c_str())) {
                P.insert(canonical(S[i]));
            }
            S[i].pop_back();
            S[i].insert(0,start);
        }
    }
    vector<string> ret(P.begin(),P.end());
    return ret;
}

/*
 * Method for computation of set cFP which is
 * formally defined as cFP = P\S
 */
void dbg::compute_cFP(unsigned int M) {
    clock_t begin = clock();
    vector<string> D = compute_P();
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << elapsed_secs << endl;


    int k = 0;
    while(k < S.size()) {
        unordered_set<string> Pi;
        vector<string> Dn;
        while (Pi.size() < M && k != S.size()) {
            Pi.insert(S[k]);
            k++;
        }

        for (int i = 0; i != D.size(); ++i) {
            if (Pi.count(canonical(D[i])) == 0) {
                Dn.push_back(canonical(D[i]));
            } //else delete
        }

        D = Dn;
    }
    unordered_set<string> d(D.begin(), D.end());
    cFP = d;
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

    while (index < S.size()) {
        int depth = 0;
        int reason = -1;
        vector<string> tmp;
        for (int right = 0; right != 2; ++right) { // tt=0 go only left from S[index], tt=1 only right 
            vector<string> front;
            front.push_back(S[index]);

            do {
                vector<string> new_front;

                for (int i = 0; i != front.size(); ++i) {
                    vector<string> ret = branch(front[i], right);

                    new_front.insert(new_front.end(), ret.begin(), ret.end());
                    if (right) {
                        marked.insert(front[i].substr(front[i].size()-k_size));
                        marked.insert(reverse_complement(front[i].substr(front[i].size()-k_size)));
                    }
                    if (!right) {
                        marked.insert(front[i].substr(0,k_size));
                        marked.insert(reverse_complement(front[i].substr(0,k_size)));
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
                    vector<string> ret = branch(new_front[0], right);
                    if (ret.size() >= 1) {
                        continue;
                    }
                    reason = 1;
                    break;
                }
            } while (1);

            if (right == 0 && reason != -1) { // remember contigs found going left
                for (int i = 0; i != front.size(); ++i) {
                    tmp.push_back(front[i]);
                }
            }

            if (right == 1 && reason != -1) { // merge right with left contigs
                for (int i = 0; i != front.size(); ++i) {
                    for (int j = 0; j != tmp.size(); ++j) {
                        if (front[i].size() + tmp[j].size()-k_size > 2 * k_size + 1) {
                            contigs.push_back(tmp[j] + front[i].substr(k_size));
                        }
                    }   
                }
            }
        }

        index += 1;

        for (int i = index; i < S.size(); ++i) {
            if(!marked.count(S[i])) {
                index = i;
                break;
            }
            if (i == S.size()-1) {
                index = S.size();
            }
        }
    }

    ofstream myfile;
    myfile.open (name);
    for (int i = 0; i != contigs.size(); ++i) {
        myfile << contigs[i] << endl;
    }
    myfile.close();
}

/*
 * Method that takes a contig s and tries to expand it
 * to the left and right if corresponding flags: left and right
 * are true
 */
vector<string> dbg::branch(string s, bool right) {

    char c[4] = {'A', 'C', 'G', 'T'};
    vector<string> ret;
    if (!right) {
        for (int i = 0; i != 4; ++i) {
            string h = c[i] + s.substr(0, k_size - 1);
            if (!right && filter.contains(h.c_str()) && !cFP.count(canonical(h)) && !marked.count(h)) {
                ret.push_back(c[i] + s);
            }
        }
    }

    if (right) {
        for (int i = 0; i != 4; ++i) {
            string h = s.substr(s.size() - k_size + 1) + c[i];
            if (right && filter.contains(h.c_str()) && !cFP.count(canonical(h)) && !marked.count(h)) {
                ret.push_back(s + c[i]);
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

void dbg::test_size() {
    // size_t bits_per_kmer = cFP.size()*sizeof(*cFP.begin())*8/float(S.size());
    // cout << "Critical false positive size: " << bits_per_kmer << endl;
    // cout << "Bloom filter size: " << sizeof(filter) << endl;
}

