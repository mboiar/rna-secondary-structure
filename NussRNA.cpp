/*
Implementation of the NussRNA class for predicting RNA secondary structure based on the Nussinov algorithm with custom restrictions:
  - Minimal hairpin loop size: 3
  - Energy of A-U base pair: -2
  - Energy of G-C base pair: -3

Algorithm has O(n^3) time complexity and O(n^2) space complexity. Pseudoknots are not allowed. 

Original paper: R Nussinov and A B Jacobson, 'Fast algorithm for predicting the secondary structure of single-stranded RNA.' (doi: 10.1073/pnas.77.11.6309)
*/

#include <string>
#include<fstream>
#include<stdexcept>
#include<iomanip>
#include<utility>
#include<algorithm>
#include"NussRNA.h"

using std::cout;
using std::endl;
using std::min;

// initialize class variables to default values
void NussRNA::reset(){
    n = seq.size();
    nrg.resize(n,{0});
	for (auto& v: nrg){
        v.resize(n,0);
	}
    bra.resize(n, '.');
}

NussRNA::NussRNA(std::string& s, char init){

    if (init=='f'){
        std::ifstream str ( s );
        if ( ! str.is_open() ){
            throw std::runtime_error ( "File not found.\n" );
        } else{
            str >> seq;
        }
    } else if (init=='s'){
        seq = s;
    } else { throw std::invalid_argument("Unknown option."); }
    reset();
}

std::ostream& operator<<(std::ostream& os, NussRNA& R){
    R.compute();
    os << endl;
    os << "Optimal energy of the configuration: " << R.score << endl;
    os << "Sequence in the dot-bracket notation: " + R.bra << endl;
    os << "Original sequence: " + R.seq << endl;
    //os << "Energy matrix:\n";
    //R.print_nrg();
    return os;
}

// output matrix with minimal energy values to a standard ostream
void NussRNA::print_nrg(){
    cout << endl;
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cout << std::setw((int)(80-n)/2/n) << nrg.at(i).at(j) << "\t";
        }
        cout <<endl;
    }
    cout << endl;
}

// compute binding energy for sequence elements at given indices
int NussRNA::get_nrg(size_t i, size_t j){
            if (j-i<4) return 0;
            if (!is_valid_pair(i, j)){
                return 0;
            }
            return (seq.at(i)=='A' || seq.at(i)=='U') ? -2: -3;
}

// dynamically retrace sequence based on the energy matrix to find base pairs
void NussRNA::traceback(size_t i, size_t j){
    if (i<j){
        if (nrg.at(i).at(j) == nrg.at(i+1).at(j)){
            traceback(i+1, j);
        } else if (nrg.at(i).at(j) == nrg.at(i).at(j-1)){
            traceback(i, j-1);
        } else if (nrg.at(i).at(j) == nrg.at(i+1).at(j-1)+get_nrg(i,j)){
            traceback(i, j-1);
            base_pairs.push(std::make_pair(i,j));
            traceback(i+1,j-1);
        } else {
            for (size_t l=i+1; l<j; l++){
                if (nrg.at(i).at(j)==nrg.at(i).at(l)+nrg.at(l+1).at(j)){
                    traceback(i, l);
                    traceback(l+1, j);
                    break;
                }
            }
        }
    }
}

// the main part of the algorithm
void NussRNA::compute(){
    size_t j;
    for (size_t k=1; k<n; k++){
        for (size_t i=0; i<n-k; i++){
            j = k+i;
            nrg.at(i).at(j) = nrg.at(i+1).at(j-1)+get_nrg(i,j);
            for (size_t t=i+1; t<j; t++){
                if (nrg.at(i).at(t)+nrg.at(t+1).at(j) < nrg.at(i).at(j)){
                    nrg.at(i).at(j) = nrg.at(i).at(t)+nrg.at(t+1).at(j);
                }
            }
            nrg.at(i).at(j) = min(min(nrg.at(i).at(j), nrg.at(i+1).at(j)), nrg.at(i).at(j-1));
        }
    }
    traceback(0, n-1);
    get_bracket_seq();
    if (!is_correct_seq()){
        throw std::runtime_error("Unable to determine optimal structure. Ensure that there are no forbidden structures like knots.");
    }
    score= nrg.at(0).at(n-1);
}

// construct sequence in dot-bracket notation based on base pairs
void NussRNA::get_bracket_seq(){
    while(!base_pairs.empty()){
    std::pair <size_t,size_t> p = base_pairs.top();
    base_pairs.pop();
    size_t i = std::get<0>(p);
    size_t j = std::get<1>(p);
    bra[i] = '(';
    bra[j] = ')';
    }

}

bool NussRNA::is_correct_seq(){
    return ( std::count(bra.begin(),bra.end(), '(') == std::count(bra.begin(),bra.end(), ')'));
}

char NussRNA::get_pair(char& base){
    switch(base){
        case 'A': return 'U';
        case 'U': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: throw std::invalid_argument("Given argument is not one of: A, U, G, C.");
    }
}