#ifndef PFAPACK_FORTRAN_H
#define PFAPACK_FORTRAN_H

#include <complex>

extern "C" {
  //prototypes for the PFAPACK routines S,D,C,ZSKTRF

  void ssktrf_(const char *, const char *, const int *,
               float *, const int *, int *, float *,
               const int *, int *);

  void dsktrf_(const char *, const char *, const int *,
               double *, const int *, int *, double *,
               const int *, int *);

  void csktrf_(const char *, const char *, const int *,
               std::complex<float> *, const int *, int *,
               std::complex<float> *, const int *, int *);

  void zsktrf_(const char *, const char *, const int *,
               std::complex<double> *, const int *, int *,
               std::complex<double> *, const int *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKTF2

  void ssktf2_(const char *, const char *, const int *,
               float *, const int *, int *, int *);

  void dsktf2_(const char *, const char *, const int *,
               double *, const int *, int *, int *);

  void csktf2_(const char *, const char *, const int *,
               std::complex<float> *, const int *, int *,
               int *);

  void zsktf2_(const char *, const char *, const int *,
               std::complex<double> *, const int *, int *,
               int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKTRD

  void ssktrd_(const char *, const char *, const int *,
               float *, const int *, float *, float *,
               float *, const int *, int *);

  void dsktrd_(const char *, const char *, const int *,
               double *, const int *, double *, double *,
               double *, const int *, int *);

  void csktrd_(const char *, const char *, const int *,
               std::complex<float> *, const int *, float *,
               std::complex<float> *,
               std::complex<float> *, const int *, int *);

  void zsktrd_(const char *, const char *, const int *,
               std::complex<double> *, const int *,
               double *, std::complex<double> *,
               std::complex<double> *, const int *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKTD2

  void ssktd2_(const char *, const char *, const int *,
               float *, const int *, float *, float *, int *);

  void dsktd2_(const char *, const char *, const int *,
               double *, const int *, double *, double *,
               int *);

  void csktd2_(const char *, const char *, const int *,
               std::complex<float> *, const int *, float *,
               std::complex<float> *, int *);

  void zsktd2_(const char *, const char *, const int *,
               std::complex<double> *, const int *,
               double *, std::complex<double> *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKPFA

  void sskpfa_(const char *, const char *, const int *,
               float *, const int *, float *, int *,
               float *, const int *, int *);

  void dskpfa_(const char *, const char *, const int *,
               double *, const int *, double *, int *,
               double *, const int *, int *);

  void cskpfa_(const char *, const char *, const int *,
               std::complex<float> *, const int *,
               std::complex<float> *, int *,
               std::complex<float> *, const int *, float *,
               int *);

  void zskpfa_(const char *, const char *, const int *,
               std::complex<double> *, const int *,
               std::complex<double> *, int *,
               std::complex<double> *, const int *,
               double *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKPF10

  void sskpf10_(const char *, const char *, const int *,
                float *, const int *, float *, int *,
                float *, const int *, int *);

  void dskpf10_(const char *, const char *, const int *,
                double *, const int *, double *, int *,
                double *, const int *, int *);

  void cskpf10_(const char *, const char *, const int *,
                std::complex<float> *, const int *,
                std::complex<float> *, int *,
                std::complex<float> *, const int *,
                float *, int *);

  void zskpf10_(const char *, const char *, const int *,
                std::complex<double> *, const int *,
                std::complex<double> *, int *,
                std::complex<double> *, const int *,
                double *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKBTRD

  void sskbtrd_(const char *, const char *, const char *,
                const int *, const int *, float *,
                const int *, float *, float *, const int *,
                float *, int *);

  void dskbtrd_(const char *, const char *, const char *,
                const int *, const int *, double *,
                const int *, double *, double *,
                const int *, double *, int *);

  void cskbtrd_(const char *, const char *, const char *,
                const int *, const int *,
                std::complex<float> *, const int *,
                float *, std::complex<float> *,
                std::complex<float> *, const int *,
                std::complex<float> *, float *, int *);

  void zskbtrd_(const char *, const char *, const char *,
                const int *, const int *,
                std::complex<double> *, const int *,
                double *, std::complex<double> *,
                std::complex<double> *, const int *,
                std::complex<double> *, double *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKBPFA

  void sskbpfa_(const char *, const int *, const int *,
                float *, const int *, float *, float *,
                int *);

  void dskbpfa_(const char *, const int *, const int *,
                double *, const int *, double *, double *,
                int *);

  void cskbpfa_(const char *, const int *, const int *,
                std::complex<float> *, const int *,
                std::complex<float> *,
                std::complex<float> *, float *, int *);

  void zskbpfa_(const char *, const int *, const int *,
                std::complex<double> *, const int *,
                std::complex<double> *,
                std::complex<double> *, double *, int *);

  //prototypes for the PFAPACK routines S,D,C,ZSKBPF10

  void sskbpf10_(const char *, const int *, const int *,
                 float *, const int *, float *, float *,
                 int *);

  void dskbpf10_(const char *, const int *, const int *,
                 double *, const int *, double *, double *,
                 int *);

  void cskbpf10_(const char *, const int *, const int *,
                 std::complex<float> *, const int *,
                 std::complex<float> *,
                 std::complex<float> *, float *, int *);

  void zskbpf10_(const char *, const int *, const int *,
                 std::complex<double> *, const int *,
                 std::complex<double> *,
                 std::complex<double> *, double *, int *);
}

#endif
