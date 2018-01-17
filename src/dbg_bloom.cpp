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
    if (S.size() == 0) {
        cout << "S has no content!" << endl;
        exit(0);
    }
    k_size = (*S.begin()).size();
    create_BF();
}

/*
 * Method for creation of bloom filter and adding all solid k-mers in the filter.
 */
void dbg::create_BF() {
    const unsigned k = k_size;
    const float fp = 2.08/(16*k);
    bloom_size = ceil(1.44*S.size()*log2(1/fp));
    const unsigned numHashes = ceil(bloom_size/float(S.size())*log(2));
    while (bloom_size%numHashes != 0) {bloom_size++;}
    filter = KmerBloomFilter(bloom_size, numHashes, k);
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
void dbg::traverse_graph(string name, int depth_bound, int width_bound) {

    int index = 0;

    while (index < S.size()) {

        string start;
        vector<string> vs;
        //cout << "traverse GRAPH: " << index / float(S.size()) << endl;
        for (int right = 0; right != 2; ++right) { // right=0 go only left from S[index], tt=1 only right
            int depth = 0;
            vector<string> front;
            front.push_back(S[index]);
            start = S[index];
            vector<vector<vector<char> > > mem;
            if (index / float(S.size()) > 5.1635*1e-5) {

            }
            do { // branch in one direction, store new kmers in new_front, remember nucleotides in mem
                vector<string> new_front;
                mem.push_back(vector<vector<char> >());
                vector<int> to_delete;
                for (int i = 0; i != front.size(); ++i) {
                    vector<char> ret = branch(front[i], right);

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
                    }
                    mem[depth].push_back(ret); // save nucleotide
                    if (ret.empty()) {          // dead line, mark for deletion
                        to_delete.push_back(i);
                    }
                    marked.insert(front[i]); // mark as used
                    marked.insert(reverse_complement(front[i])); // mark as used
                }

                if (new_front.empty()) {
                    break;
                }

                int deleted = 0;
                for (int i = 0; i != to_delete.size(); ++i) { // kill deleted branches
                    int ind = to_delete[i]-deleted;
                    mem[depth].erase(mem[depth].begin() + ind);
                    bool done = false;
                    for (int j = depth-1; !done ; --j) {
                        int position_change = 0;
                        for (int k = 0, tmp = 0; k != mem[j].size();k++) {
                            if (tmp + mem[j][k].size() - 1 < ind) {
                                tmp += mem[j][k].size();
                                position_change+= mem[j][k].size()-1;
                            } else {
                                if (mem[j][k].size() == 1) {
                                    mem[j].erase(mem[j].begin() + k);
                                    break;
                                } else {
                                    mem[j][k].erase(mem[j][k].begin() + (ind - tmp));
                                    done = true;
                                    break;
                                }
                            }
                        }
                        ind -= position_change;
                    }
                    deleted++;
                }


                depth++;
                int width = new_front.size();
                if (depth > depth_bound || width > width_bound) {
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
                vs = construct_contigs(mem, depth-1, right);
            }

            if (right == 1) { // merge right with left contigs
                vector<string> rs = construct_contigs(mem, depth-1, right);;
                for (int i = 0; i != vs.size(); ++i) {
                    for (int j = 0; j != rs.size(); ++j) {
                        if (vs[i].size() + k_size+ rs[j].size() > 2*k_size+1) {
                            contigs.push_back(vs[i] + start + rs[j]);
                        }
                    }
                }
            }

        }

        index += 1;

        for (int i = index; i < S.size(); ++i) {
            if (!marked.count(S[i])) {
                index = i;
                break;
            }
            if (i == S.size() - 1) {
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
 * Method for constructing contigs left and right part from mem structure.
 * It connects chars on different depths consistently.
 */
vector<string> dbg::construct_contigs(vector<vector<vector<char> > > mem, int depth, bool right) {
    vector<string> ret;
    if (depth == -1) {
        ret.push_back("");
        return ret;
    }

    vector<string> r = construct_contigs(mem, depth-1,right);
    for (int i = 0; i != r.size(); ++i) {
        for (int j = 0; j != mem[depth][i].size(); ++j) {
            if (right) {
                ret.push_back(r[i]+mem[depth][i][j]);
            } else {
                ret.push_back(mem[depth][i][j]+r[i]);
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
    char a,b;
    for(int i = 0; i != s.size()/2; ++i) {
        a = complement(s[i]); b = complement(s[s.size()-1-i]);
        s[i]= b; s[s.size()-1-i]=a;
    }
    if (s.size() %2 == 1) {
        s[s.size()/2] = complement(s[s.size()/2]);
    }
    return s;
}

/*
 * Method for calculation of a complement of a nucleotid.
 */
char dbg::complement(char c) {
    if (c == 'C')
        return 'G';
    if (c == 'T')
        return 'A';
    if (c == 'A')
        return 'T';
    if (c == 'G')
        return 'C';
    return 'U';
}

void dbg::test_size() {
    double space = cFP.size()*k_size/(1024.*1024);
    cout << setw (30) << "cFP space : " << space << " MB." << endl;
    space = bloom_size/(8*1024.*1024.);
    cout << setw (30) << "Bloom filter space : " << space << " MB." << endl;
}

