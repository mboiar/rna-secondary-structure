/*
A terminal client for the RNA secondary structure prediction algorithm (based on the Nussinov algorithm).

Author: Maksym Boiar
Year: 2022
*/

#include"NussRNA.h"

// using std::cerr;
// using std::exception;
// using std::string;
// using std::cout;
// using std::endl;
using namespace std;

int main(int argc, char* argv[]){
    string DFT = "GACGAAAGUC";                // default sequence to analyze
    if (argc==1){
        cout << "Using default sequence:\n"
             << DFT << endl;
        NussRNA r(DFT, 's');
        cout << r << endl;
    } else if (argc!=3 || (*argv[1]!='f' && *argv[1]!='s')) {
        cout << "Usage: .\\NussRNA[.exe] OPTION ARGUMENT\n"
             << "OPTIONS\nf\tread sequence from the file ARGUMENT\n"
             << "s\tuse ARGUMENT as an input sequence\n";
        return -1;

    } else {
        try{
            string s(argv[2]);
            NussRNA r(s, *argv[1]);
            cout << r << endl;
        } catch (const exception& e){
            cerr << e.what();
            return -1;
        }
    }
    return 0;
}

