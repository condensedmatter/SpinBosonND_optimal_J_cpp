#ifndef H_FILE
#define  H_FILE


#include <iostream>
#include <fstream>

#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;


class writeFiles{
    public:
    int m_i;
    static const int m_i_DIGITS=5;
    int m_seed;
    static const int m_seed_DIGITS=3;

    ofstream m_file;
    string m_name;

    writeFiles(int seed){
        m_seed=seed;
        m_i=0;
    }

    // method 1, save vector
    void save(vector<double> vd){
          m_name=get_name(m_seed,m_i);
          std::cout<<"creating file: "<<m_name<<std::endl;

          m_file.open (m_name.c_str());
          write_file(m_file,vd);
          m_file.close();
          m_i=m_i+1;
    }
    // method 2, save from pointer and size
    void save(double *s, int spinsize) {
      m_name=get_name(m_seed,m_i);
      std::cout<<"creating file: "<<m_name<<std::endl;
      m_file.open (m_name.c_str());
      for(int i=0;i<spinsize;i++){
          m_file<<s[i]<<endl;
      }
      m_file.close();
      m_i=m_i+1;
    }

    void write_file(ofstream & f,vector<double> vd){
        for(int i=0;i<vd.size();i++){
            f<<vd[i]<<endl;
        }
    }

    string get_name(int seed, int index){
        std::stringstream ss;
        ss <<std::setw(m_seed_DIGITS) << std::setfill('0') <<seed <<"_" <<  std::setw(m_i_DIGITS) << std::setfill('0') << index <<".dat";
        std::string s = ss.str();
        return s;
    }

};





class Summation{
  public:
  int spinsize;
  int sumN;
  std::vector<double> c;
  Summation(int i_spinsize,int i_sumN){
    spinsize=i_spinsize;
    sumN=i_sumN;
    for (size_t i = 0;i < spinsize; i++) {
      c.push_back(0.0);
    }
  }

  void add(double *s) {
    for(int i=0;i<spinsize;i++){
        c[i]=c[i]+s[i];
    }
  }

  void zero() {
    for(int i=0;i<spinsize;i++){
        c[i]=0.0;
    }
  }

  void ave() {
    for(int i=0;i<spinsize;i++){
        c[i]=c[i]/sumN;
    }
  }

};





#endif
