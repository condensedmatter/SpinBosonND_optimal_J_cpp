#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <cmath> //exp()
#include <random>  // class my_random{};
#include <numeric> // accumulate();
using namespace std;
class my_random{
    public:
    std::random_device rd;
    std::mt19937 mt;
    std::uniform_int_distribution<int> dis_bool;
    std::uniform_real_distribution<double> dis_01;

    my_random():mt(rd()), dis_bool(0,1),dis_01(0.0,1.0) {}
    my_random(int seed):mt(seed), dis_bool(0,1),dis_01(0.0,1.0){}

    bool rd_bool(){
        return dis_bool(mt);
    }
    double rd_01(){
        return dis_01(mt);
    }

    void rd_thermo(bool *s,int Ntotal){
        for(int i=0;i<Ntotal;i++)
            s[i]=rd_bool();
    }

    int rd_getone(int Ntotal){
        return Ntotal*rd_01();
    }
};
my_random rng(2);//global object rng
// to reset new seed call this:
// rng.mt.seed(3);

class my_random_neighbour{
    public:
    int M; //numbers to neighbhour to be selected,  M>=1  M is integer
    double **Cum; //Cum[M][M+2]
    my_random_neighbour(vector<double> Couplings, int M_input){
        M=M_input;
        Cum=new double*[M];
        for(int i=0;i<M;i++){
            Cum[i]=new double[M+2];
        }
        for(int i=0;i<M;i++){
            Cum[i][i+1]=Couplings[i];
            Cum[i][M+1]=1.0;
            for(int j=0;j<i+1;j++){
                Cum[i][j]=0.0;
            }
        }
        for(int i=0;i<M-1;i++){
            for(int j=i+2;j<M+1;j++){
                Cum[i][j]=Cum[i][j-1]+Couplings[j-1];
            }
        }
        for(int i=0;i<M;i++){
            for(int j=i+1;j<M+1;j++){
                Cum[i][j]=1.0-exp(-2.0*Cum[i][j]);
            }
        }
    }
    ~my_random_neighbour(){
        for(int i=0;i<M;i++){
            delete [] Cum[i];
        }
        delete [] Cum;
    }

    vector<int> rd_neighbour(){
        vector<int> output;
        int index=0;
        for(int i=0;i<M; ){
            /*
            index =insert_rand( Cum[i],rng.rd_01(),M+2);
            if(index<M) {output.push_back(index);}//cout<<"push "<<index<<endl;}
            i=index+1;
           */
            index =insert_rand( Cum[i]+i,rng.rd_01(),M+2-i);
            if(index+i <M) {output.push_back(index+i);}
            i=index+1+i;
        }
        return output;   // { output= {x} | x in 0,1,...,M-1}
    }

    int insert_rand(double *a,double r,int M){
        int i=0;
        int j=M;
        int mid=0;
        while(i<j-1){
            mid=(i+j)/2;
            if(r< *(a+mid)){j=mid;}
            else{i=mid;}
        }
    return i;
    }
}; //Cumulative probability method

class My_random{
    public:
    my_random_neighbour a;
    My_random(int N0,double K0,double A):a(couples(N0,K0,A),N0-1){
    }
    My_random(int N0,int N1,double K0,double K1,double A):a(couples(N0,N1,K0,K1,A),N0+1){
    }
    My_random(int N0,int N1,int N2,double K0,double K1,double K2,double A):a(couples(N0,N1,N2,K0,K1,K2,A),N0+3){
    }


    std::vector<double> couples(int N0,double K0,double A){
        vector<double> output;
        double pi=atan(1.0)*4.0;//3.14159265
        for(int i=1;i<N0;i++){
            output.push_back(   0.5*A*pow(pi/N0,2)/pow(sin(pi*i/N0),2) );
        }
        output[0]+=K0;
        output[N0-2]+=K0;
        return output;
    }


    vector<double> couples(int N0,int N1,double K0,double K1,double A){
        vector<double> output;
        double pi=atan(1.0)*4.0;//3.14159265
        for(int i=1;i<N0;i++){
            output.push_back(   0.5*A*pow(pi/N0,2)/pow(sin(pi*i/N0),2) );
        }
        output[0]+=K0;
        output[N0-2]+=K0;
        output.push_back(K1);
        output.push_back(K1);
        return output;
    }

    vector<double> couples(int N0,int N1,int N2,double K0,double K1, double K2,double A){
        vector<double> output;
        double pi=atan(1.0)*4.0;//3.14159265
        for(int i=1;i<N0;i++){
            output.push_back(   0.5*A*pow(pi/N0,2)/pow(sin(pi*i/N0),2) );
        }
        output[0]+=K0;
        output[N0-2]+=K0;
        output.push_back(K1);
        output.push_back(K1);
        output.push_back(K2);
        output.push_back(K2);


        return output;
    }

    vector<int> get_rand(){
        return a.rd_neighbour();
    }
};


// h[i][j]  the j-th neighbour of site i has index h[i][j]
// the coding is in pdf on Overleaf.com
void hashing_1d(int **h,int N0){
    int s;
    for(int i=0;i<N0;i++){
            for(int k=0;k<N0-1;k++){
                h[i][k]=(i+k+1)%N0;
            }
    }
}
void hashing_2d(int **h,int N0,int N1){
    int s;
    for(int i=0;i<N0;i++){
        for(int j=0;j<N1;j++){
            s=i+N0*j;
            for(int k=0;k<N0-1;k++){    // N0-1
                h[s][k]=(i+k+1)%N0+N0*j;
            }
            h[s][N0-1]=i+N0*((j-1+N1)%N1);
            h[s][N0]=i+N0*((j+1)%N1);
        }
    }
}
void hashing_3d(int **h,int N0,int N1,int N2){
    int s;
    for(int i=0;i<N0;i++){
        for(int j=0;j<N1;j++){
            for (int n = 0; n < N2; n++) {
                s=i+N0*j+N0*N1*n;

                for(int k=0;k<N0-1;k++){    // N0-1
                    h[s][k]=(i+k+1)%N0+N0*j+N0*N1*n;
                }
                h[s][N0-1]=i+N0*((j-1+N1)%N1)+N0*N1*n;
                h[s][N0]=i+N0*((j+1)%N1)+N0*N1*n;
                h[s][N0+1]=i+N0*j+N0*N1*((n-1+N2)%N2);
                h[s][N0+2]=i+N0*j+N0*N1*((n+1+N2)%N2);
            }
        }
    }


}




#include <queue>   // class Ising.updating()  Q
#include <vector> // class Ising.Ns
class Ising{
    public:
    int dimension;
    int Ntotal;
    vector<int> Ns;
    vector<double> Ks;//next neighbour coupling strength
    double A;//longe range coupling strength
    bool *s;//pointer to the Ising field

    int **hashing;//pointer to a 2D array  hashing[i][j]  site i's j-th neighbour
                  //j-th neighbour

    My_random pp;
        // 0+1 D long range Ising problem
	Ising(int N0,double K0,double A0):pp(N0,K0,A0){
	    //neither work for N1=1 nor N1=2,  we should set N1>=3
	    //N0 is imaginary time, should be larger
	    dimension=1; Ntotal=N0;
	    Ns.push_back(N0);
	    Ks.push_back(K0);A=A0;
	    s=new bool[Ntotal];
	    rng.rd_thermo(s,Ntotal);// RANDOM 1
	    hashing=new int*[Ntotal];
	    for(int i=0;i<Ntotal;i++){
		hashing[i]=new int[(dimension-1)*2+N0-1];
	    }
	    hashing_1d(hashing,N0);
	}
    // 1+1 D long range Ising problem
    Ising(int N0,int N1,double K0,double K1,double A0):pp(N0,N1,K0,K1,A0){
        //neither work for N1=1 nor N1=2,  we should set N1>=3
        //N0 is imaginary time, should be larger
        dimension=2; Ntotal=N0*N1;
        Ns.push_back(N0);Ns.push_back(N1);
        Ks.push_back(K0);Ks.push_back(K1);A=A0;
        s=new bool[Ntotal];
        rng.rd_thermo(s,Ntotal);// RANDOM 1
        hashing=new int*[Ntotal];
        for(int i=0;i<Ntotal;i++){
            hashing[i]=new int[(dimension-1)*2+N0-1];
        }
        hashing_2d(hashing,N0,N1);
    }
    // 2+1 D long range Ising problem
    Ising(int N0,int N1,int N2,double K0,double K1,double K2,double A0):pp(N0,N1,N2,K0,K1,K2,A0){
        dimension=3; Ntotal=N0*N1*N2;
        Ns.push_back(N0);Ns.push_back(N1);Ns.push_back(N2);
        Ks.push_back(K0);Ks.push_back(K1);Ks.push_back(K2);A=A0;
        s=new bool[Ntotal];
        rng.rd_thermo(s,Ntotal);// RANDOM 1
        hashing=new int*[Ntotal];
        for(int i=0;i<Ntotal;i++){
            hashing[i]=new int[(dimension-1)*2+N0-1];
        }
        hashing_3d(hashing,N0,N1,N2);
    }
    ~Ising(){
        delete [] s;
        for(int i=0;i<Ntotal;i++){
            delete [] hashing[i];
        }
        delete [] hashing;
    }

    // Wolff cluster updating algorithm
    void updating(){
          queue<int> Q;
          vector<int> bb;
          int temp=rng.rd_getone(Ntotal);// RANDOM 2

          Q.push(temp);
          s[temp]=!s[temp];
          int site, neighbour;
          while(!Q.empty()){
              site=Q.front();
              Q.pop();
              bb=pp.get_rand(); // RANDOM 3
              for(int i=0;i<bb.size();i++){
                  neighbour=hashing[site][bb[i]];
                  if(s[neighbour]!=s[temp]){         //*hash[site][bb[i]]==color;
                      Q.push(neighbour);
                      s[neighbour]=!s[neighbour];
                  }
              }
          }
    }

    // updating "loops" times
    void updating(int loops){
      for (size_t i = 0; i < loops; i++) {
        updating();
      }
    }

};



#endif
