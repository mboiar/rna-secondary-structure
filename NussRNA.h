#ifndef nussrna
#define nussrna

#include <iostream>
#include<vector>
#include<utility>
#include<stack>

class NussRNA{
    std::string seq, bra;
    size_t n;
    int score;
    std::stack<std::pair<size_t,size_t>> base_pairs;
    std::vector<std::vector<int> > nrg;

    public:
        static char get_pair(char&);
        NussRNA(std::string& s, char init);
        bool is_valid_pair(size_t i, size_t j){
            return (i<j) && (seq[j] == get_pair(seq.at(i)));
        }
        int get_nrg(size_t i, size_t j);
        void compute();
        void print_nrg();
        void traceback(size_t, size_t);
        void reset();
        void get_bracket_seq();
        bool is_correct_seq();

        friend std::ostream& operator<<(std::ostream& os, NussRNA& R);
};

#endif