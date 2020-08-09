#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// calcolo dell'agreement fra bitset
NumericVector compute_bitcor_CPP(NumericMatrix bitvects) {
    std::vector<std::vector<double>*>* COR = new std::vector<std::vector<double>*>();
    for(int i=0; i<bitvects.nrow(); i++){
        COR->push_back(new std::vector<double>());
        for(int j=0; j<=i; j++){
            COR->at(i)->push_back(0);
            for(int k=0;k<bitvects.ncol();k++){
                COR->at(i)->at(j) += (bitvects(i,k) == bitvects(j,k))?1:0;     
            }
            COR->at(i)->at(j) /= bitvects.ncol();
        }
    }
    NumericMatrix out(bitvects.nrow(),bitvects.nrow());
    for(int i=0; i<bitvects.nrow(); i++){
        for(int j=0; j<=i; j++){
            out(i,j) = COR->at(i)->at(j);
            out(j,i) = COR->at(i)->at(j);
        }
    }
    
    return out;
}

