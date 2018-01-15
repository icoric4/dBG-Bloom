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
    const unsigned numHashes = 4;
    const unsigned size = 1e9;
    filter = KmerBloomFilter(size, numHashes, k);
    for (int i = 0; i != S.size(); ++i) {
        filter.insert(S[i].c_str());
        filter.insert(reverse_complement(S[i]).c_str());
    }
}

/*
 * Method for calculation of set P which represents the set of all extensions of
 * set S for which the Bloom filter answers yes.
 */
vector<string> dbg::compute_P() {
    string c[] = {"A","C","G","T"};
    unordered_set<string> P;
    for (int i = 0; i != S.size(); ++i) {
        for (int j = 0; j != 4; ++j) {
            string h = S[i].substr(1) + c[j];
            if (filter.contains(h.c_str())) {
                P.insert(canonical(h));
            }
            h = c[j]+S[i].substr(0,S[i].size()-1);
            if (filter.contains(h.c_str())) {
                P.insert(canonical(h));
            }
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
    
    vector<string> D = compute_P();
    
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

        string start;
        vector<string> vs;
        for (int right = 0; right != 2; ++right) { // right=0 go only left from S[index], tt=1 only right
            int depth = 0, first_not_single = -1;
            vector<string> front; front.push_back(S[index]); start = S[index];
            vector<set<char> > mem;
            do { // branch in one direction, store new kmers in new_front, remember nucleotides in mem
                vector<string> new_front;
                for (int i = 0; i != front.size(); ++i) {
                    vector<char> ret = branch(front[i], right);
                    if (i == 0) { // init mem
                        set<char> sc(ret.begin(),ret.end());
                        mem.push_back(sc);
                    }

                    for (int j = 0; j != ret.size(); ++j) { // construct new kmers
                        string h;
                        if (right) {
                            h = front[i] + ret[j];
                            h.erase(0, 1);
                        } else {
                            h = ret[j] + front[i];
                            h.pop_back();
                        }
                        new_front.push_back(h); // save them
                        if (j != 0) {
                            mem[depth].insert(ret[j]); // save nucleotide
                        }
                    }
                    marked.insert(front[i]); // mark as used 
                    marked.insert(reverse_complement(front[i])); // mark as used
                }

                if (new_front.size() > 1 && first_not_single == -1) { // if more than one option, remember to update later when correct one is chosen
                        first_not_single = depth;
                } else if (first_not_single != -1) {
                    if (mem[first_not_single].size() > new_front.size()) { // we got less options than before and make an update
                        for(int j = first_not_single; j < depth && mem[j].size() != 1; ++j) { // all of them must have more options than current new_front, so trim them
                            set<char> sc;
                            for(int k = 0; k != new_front.size(); ++k) {
                                sc.insert(new_front[k][new_front[k].size()-1-(depth-j)]);
                            }
                            mem[j] = sc;
                            if (sc.size() == 1) {
                                first_not_single++;
                            }
                        }
                        if (mem[first_not_single].size()==1) {
                            first_not_single = -1;
                        }
                    }
                }


                if (new_front.empty()) {
                    break;
                }
                depth++;
                int width = new_front.size();
                if (depth > 1e6 || width > 20) {
                    break;
                }
                front = new_front;
                if (new_front.size() == 1) {
                    vector<char> ret = branch(new_front[0], right);
                    if (ret.size() >= 1) {
                        continue;
                    }
                    break;
                }
            } while (1);

            if (right == 0) { // remember contigs found going left in vs
                vs.push_back(start);
                for (int i = 0; i != mem.size(); ++i) {
                    vector<string> tmp;
                    for (int j = 0; j != vs.size(); ++j) {
                        for (set<char>::iterator it = mem[i].begin(); it != mem[i].end(); ++it) {
                            tmp.push_back(*it+vs[j]);
                        }
                    }
                    vs = tmp;
                }
            }

            if (right == 1) { // merge right with left contigs
                vector<string> rs = vs;
                for (int i = 0; i != mem.size(); ++i) {
                    vector<string> tmp;
                    for (int j = 0; j != rs.size(); ++j) {
                        for (set<char>::iterator it = mem[i].begin(); it != mem[i].end(); ++it) {
                            tmp.push_back(rs[j] + *it);
                        }
                    }
                    rs = tmp;
                }


                for (int i = 0; i != rs.size(); ++i) {
                        if (rs[i].size() > 2 * k_size + 1) {
                            contigs.push_back(rs[i]);
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
vector<char> dbg::branch(string s, bool right) {

    char c[4] = {'A', 'C', 'G', 'T'};
    vector<char> ret;
    if (!right) {
        for (int i = 0; i != 4; ++i) {
            string h = c[i] + s.substr(0, k_size - 1);
            if (!right && filter.contains(h.c_str()) && !cFP.count(canonical(h)) && !marked.count(h)) {
                ret.push_back(c[i]);
            }
        }
    }

    if (right) {
        for (int i = 0; i != 4; ++i) {
            string h = s.substr(s.size() - k_size + 1) + c[i];
            if (right && filter.contains(h.c_str()) && !cFP.count(canonical(h)) && !marked.count(h)) {
                ret.push_back(c[i]);
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