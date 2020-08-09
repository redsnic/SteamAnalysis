#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_shared_userbase_similarity(NumericMatrix A, NumericMatrix E, int size, NumericVector selection) {
  
    // Inizializza la matrice dei pesi
    std::vector<std::vector<double>*>* W = new std::vector<std::vector<double>*>();
    for(int i=0; i<size; i++){
        W->push_back(new std::vector<double>());
        for(int j=0; j<size; j++){
            W->at(i)->push_back(0);
        }
    }
    
    // calcolo della matrice dei pesi
    for(int i=0; i<E.nrow(); i++){
      W->at((int)E(i,0)-1)->at((int)E(i,1)-1) = E(i,2)-1;
    }

    // calcolo della lista di adiacenza
    std::vector<std::vector<double>*>* AL = new std::vector<std::vector<double>*>();
    for(int i=0; i<size; i++){
      AL->push_back(new std::vector<double>());
      for(int j=0;j<size; j++){
        if(A(i,j)){
          AL->at(i)->push_back(j);
        }
      }
    }

    // algoritmo BFS-based
    NumericVector out(size);
    
    for(int s:selection){
      int iter = 0;
      NumericVector flow(size);
      NumericVector considered(size);
      NumericVector queued(size);
      std::vector<int> Q = std::vector<int>();
      Q.push_back(s);
      while(Q.size() > 0){
        int source = Q.back(); 
        Q.pop_back();
        for(int destination:*(AL->at(source))){
          if(!considered[destination]){
            flow[destination] += flow[source]/2 + W->at(source)->at(destination)/(pow(2,iter));
            iter++;
            if(!queued[destination]){
              Q.push_back(destination);
              queued[destination] = 1;
            }
          }
          considered[source] = 1;
        }
      }
      
      // salva il flusso 
      for(int i=0; i<size;i++){
        out[i] += flow[i];
      }
    }
    
    return out;
}

