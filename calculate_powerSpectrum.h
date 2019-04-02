#ifndef powerSpectrum_H
#define powerSpectrum_H

#include <fftw3.h>

class powerSpectrum{
  public:
    fftw_complex *in, *out;
    fftw_plan p;

    int m_N0,m_N1,m_N2;
    int m_spinsize;
    double * m_spectrum;


    //1D
    powerSpectrum(int N0 ){
       m_N0=N0; 
       m_spinsize=N0 ;
       m_spectrum=new double[m_spinsize];
       in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_spinsize);
       out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_spinsize);
       p=  fftw_plan_dft_1d( N0,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    }

    //2D
    powerSpectrum(int N0,int N1){
       m_N0=N0;
       m_N1=N1;
       m_spinsize=N0*N1;
       m_spectrum=new double[m_spinsize];
       in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_spinsize);
       out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_spinsize);
       p=  fftw_plan_dft_2d(N1,N0,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    }

    //3D
    powerSpectrum(int N0,int N1,int N2){
       m_N0=N0;
       m_N1=N1;
       m_N2=N2;
       m_spinsize=N0*N1*N2;
       m_spectrum=new double[m_spinsize];
       in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_spinsize);
       out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_spinsize);
       p=  fftw_plan_dft_3d(N2,N1,N0,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    }

    ~powerSpectrum(){
      delete [] m_spectrum;
    }

    void calculate(bool *s){
        for (size_t i = 0; i < m_spinsize; i++) {
           in[i][0]=(s[i])?1.0:-1.0;
           in[i][1]=0.0;
        }
        fftw_execute(p); /* FFT */
        for (size_t i = 0; i < m_spinsize; i++) {
            m_spectrum[i]=out[i][0]*out[i][0]+out[i][1]*out[i][1];
        }
    }

};
#endif
