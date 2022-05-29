%module mymodule

%{
    #define SWIG_FILE_WITH_INIT
    #include "qr.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double *IN_ARRAY2, int DIM1, int DIM2) {(double *a, int ma, int na)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* q, int nq)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* r, int mr, int nr)};

%include "qr.h"
