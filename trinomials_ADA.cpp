
#include "is_prime.h"
#include <vector>
#include <iostream>
#include <set>
#include <fstream>
#include <algorithm>
#include <string>
#include <bitset>
#include "common_factor_rt.hpp"

#define debug(x) cerr << #x << " -> " << (x) << "\n"; 
using namespace std;
typedef long long lint;

//------------------------------------------------------------------------------
//utility functions

//compute exponents by repeated squaring
lint powm(lint b, lint r, lint p) {
  r = r % (p - 1);
  const bitset<64> pow2(r);
  lint result = 1;
  lint count = 0;
  for (int i = 0; count < pow2.count(); ++i) {
    if (pow2.test(i)) {
      ++count;
      lint a = b;
      for (int j = 0; j < i; ++j) {
        a = (a * a) % p;
      }
      result = (result * a) % p;
    }
  }
  return result;
}

lint next_prime(lint p) {
  p += 2; 
  while (!is_prime_mr(p)) {
    p += 2;
  }
  return p;
}

//------------------------------------------------------------------------------
//trinomial testing

//global data for all parallel threads
vector<lint> A;
vector<lint> B;
lint kiran_A;
lint degree;
lint mult_d;

lint trinomial_count(lint p) {
  lint c = 1;
  degree = -1;
  fill(A.begin(), A.end(), 0);  
  while (A.size() < p) {
    A.push_back(0);
    B.push_back(0);
  }
  
  for (lint n = 2; n < (p - 1) / 2 + 2; ++n) {
    lint cn = 1;
    #pragma omp parallel
    {

      #pragma omp for
      for (lint i = 0; i < p; ++i) {
        B[i] = 1;
      }

      #pragma omp for
      for (lint i = 2; i < p; ++i) {
        lint next = (A[i] - 1); 
        if (next == -1)
          next = p - 1;         
        next *= i;              
        A[i] = next % p;       
        
        #pragma omp atomic
          ++B[A[i]];            
      }
      

      #pragma omp for reduction(max : cn)
      for (lint i = 2; i < p; ++i) {
        if (B[i] > cn) {
          cn = B[i];
        }
      }
    } //end parallel section

    if (cn > c) {
      c = cn;
      degree = n;
      mult_d = 1;
      kiran_A = -1;
	  
      #pragma omp parallel for
      for (int i = 2; i < p; ++i) {
        if (B[i] == cn) {
          if (kiran_A == -1)
            kiran_A = i;
        }
      }	//end parallel section
    }
    else if (cn == c) {
      ++mult_d;
    }
  }

  return c;
}

void print_data(lint p, lint n, ostream& out = cout) {
  //p prime
  //n = maximal roots of a trinomial
  //kiran_A = kiran's "A" parameter
  //d = minimal degree of trinomial with n roots
  //mult_d = number of distinct degrees that admit n roots
  out << p << "\t\t " << n << "\t " << kiran_A << "\t\t " <<"\t " << degree << "\t\t " << mult_d * 2 << "\n"; 
}

//begin testing at Fp; ignore all results of N or less
void search(lint p, lint N) {
  set<int> n_found;
  for (lint n = 0; n <= N; ++n)
    n_found.insert(n);

  while (true) {
    lint n = trinomial_count(p);
    ofstream ofs2("trinomials.txt", fstream::app);
    print_data(p, n, ofs2);
    //print_data(p, n);
    ofs2.close();

    if (n_found.count(n) == 0) {
      n_found.insert(n);
      ofstream ofs("extremal_pn.txt", fstream::app);
      print_data(p, n, ofs);
      print_data(p, n);
    }

    p = next_prime(p);
  }
}

//------------------------------------------------------------------------------

int main(int argc, char** argv) {
  lint P, N;
  if (argc >= 2) P = std::stoi(argv[1]);
  else P = 5;
  if (argc == 3) N = std::stoi(argv[2]);
  else N = 1;
  
  ofstream ofs("threads used.txt", fstream::trunc);
  #pragma omp parallel
  {
    ofs << "x";
  }
  ofs.close();

  search(P, N);
}
