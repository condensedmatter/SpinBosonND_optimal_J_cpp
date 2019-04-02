/*
Author: Jian Wang 
Date: Jan 16th, 2019
Goal: To find the optimal J, the jumpings step of Monte Carlo updating,
such that the calculated average quantity gets the smallest standard error.
*/

#include "Ising.h"  
/*
rng.mt.seed(seed);

Ising e1(N0,N1,K0,K1,alpha);     
e1.updating(n);     
e1.Ntotal
e1.s
*/


#include "calculate_powerSpectrum.h"  
/*
powerSpectrum pw(N0,N1);   
pw.calculate(e1.s);
pw.m_spectrum
*/


#include "file_generator.h"  
/*
Summation s(e1.Ntotal,sumN);      
s.zero();
s.add(pw.m_spectrum);
s.ave();
s.c;

writeFiles w(seed);
w.save(s.c);
*/


#include <ctime>
/*
clock_t begin;
clock_t end;
elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
*/

#include <iostream>

#include <fstream>
/*
f.open("timeSpend.txt");
f << elapsed_secs << '\n';
f.close();
*/


#include <string>

int main(int argc, char const *argv[]) {

  int thremoLOOPS=3000000;//thermolization steps
  int LOOPS=10;//LOOPS with fixed parameters
  int sumN=20;//summation steps to generate one saved file

  //this is the variable of interest.
  //seed and J are set to be the same, to identify different names of file groups
  std::vector<int> Js={  21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584 };
  int J,seed;
    
  /* the total numbers of file generated by this program is
    LOOPS x  len(Js)  + 1
   where 1 is the "timeSpend.txt"  file, recording len(Js) times as normalization factor
  */
    
    
  //this is the physics parameters  
  int N0=128;
  int N1=64;
  //int N2=16;
  double K0=0.136;
  double K1=0.2;
  //double K2=0.1;
  double alpha=0.2;

  Ising e1(N0,N1,  K0, K1 , alpha);
  powerSpectrum pw(N0,N1 );
    
  /***************************  parameters are above this line  ****************************/

  Summation s(e1.Ntotal,sumN);
    
  e1.updating(thremoLOOPS);

  clock_t begin;
  clock_t end;
  double elapsed_secs;
  ofstream f;
  f.open("timeSpend.txt");

  for (auto v:Js) {
    begin = clock();
    J=v;
    seed=v;
    rng.mt.seed(seed);
    writeFiles w(seed);

    for (size_t i = 0; i < LOOPS; i++) {
      std::cout << i << '\n';
      s.zero();
      for (size_t j = 0; j < sumN; j++) {
        e1.updating(J);
        pw.calculate(e1.s);
        s.add(pw.m_spectrum);
      }
      s.ave();
      w.save(s.c);
    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    f << elapsed_secs << '\n';
  }
  f.close();
  return 0;
}
