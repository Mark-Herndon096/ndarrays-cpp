#ifndef _NDARRAY_UTILS_H_
#define _NDARRAY_UTILS_H_

#include <cstring>
#include <fstream>

// -------- Function declarations ---------- //
template<class Type>
ndarray<Type> zeros(int ni);
template<class Type>
ndarray<Type> zeros(int ni, int nj);
template<class Type>
ndarray<Type> zeros(int ni, int nj, int nk);
template<class Type>
ndarray<Type> zeros(int ni, int nj, int nk, int nl);
template<class Type>
ndarray<Type> zeros(int ni, int nj, int nk, int nl, int nm);

template<class Type>
ndarray<Type> ones(int ni);
template<class Type>
ndarray<Type> ones(int ni, int nj);
template<class Type>
ndarray<Type> ones(int ni, int nj, int nk);
template<class Type>
ndarray<Type> ones(int ni, int nj, int nk, int nl);
template<class Type>
ndarray<Type> ones(int ni, int nj, int nk, int nl, int nm);

template<class Type>
ndarray<Type> linspace(Type xstart, Type xend, int nsteps);

ndarray<int> slice(int is, int ie);
ndarray<int> slice(int is, int ie, int js, int je);
ndarray<int> slice(int is, int ie, int js, int je, int ks, int ke);
ndarray<int> slice(int is, int ie, int js, int je, int ks, int ke, int ls, int le);
ndarray<int> slice(int is, int ie, int js, int je, int ks, int ke, int ls, int le, int ms, int me);

template<class Type>
ndarray<Type> squeeze(ndarray<Type> x);

// ----------- Function definitions --------------- //
template<class Type>
ndarray<Type> zeros(int ni){
    ndarray<Type> return_array(ni);
    return return_array;
};
template<class Type>
ndarray<Type> zeros(int ni, int nj){
    ndarray<Type> return_array(ni,nj);
    return return_array;
};

template<class Type>
ndarray<Type> zeros(int ni, int nj, int nk){
    ndarray<Type> return_array(ni,nj,nk);
    return return_array;
};

template<class Type>
ndarray<Type> zeros(int ni, int nj, int nk, int nl){
    ndarray<Type> return_array(ni,nj,nk,nl,nm);
    return return_array;
};

template<class Type>
ndarray<Type> zeros(int ni, int nj, int nk, int nl, int nm){
    ndarray<Type> return_array(ni,nj,nk,nl,nm);
    return return_array;
};

template<class Type>
ndarray<Type> ones(int ni){
    ndarray<Type> return_array(ni);
    return_array = (Type) 1;
    return return_array;
};
template<class Type>
ndarray<Type> ones(int ni, int nj){
    ndarray<Type> return_array(ni,nj);
    return_array = (Type) 1;
    return return_array;
};

template<class Type>
ndarray<Type> ones(int ni, int nj, int nk){
    ndarray<Type> return_array(ni,nj,nk);
    return_array = (Type) 1;
    return return_array;
};

template<class Type>
ndarray<Type> ones(int ni, int nj, int nk, int nl){
    ndarray<Type> return_array(ni,nj,nk,nl,nm);
    return_array = (Type) 1;
    return return_array;
};

template<class Type>
ndarray<Type> ones(int ni, int nj, int nk, int nl, int nm){
    ndarray<Type> return_array(ni,nj,nk,nl,nm);
    return_array = (Type) 1;
    return return_array;
};

template<class Type>
ndarray<Type> linspace(Type xstart, Type xend, int nsteps){
    ndarray<Type> return_array(nsteps);
    Type xspan = xend - xstart;
    Type dx = (Type) xspan/(nsteps-1);
    for (int i=0; i<nsteps; i++){
        return_array.elements[i] = xstart + (Type)i*dx;
    }; 
    return return_array;
};

ndarray<int> slice(int is, int ie){
    int js = 0; int je = 0;
    int ks = 0; int ke = 0;
    int ls = 0; int le = 0;
    int ms = 0; int me = 0;
    ndarray<int> return_array(10);
    return_array(0) = is;
    return_array(1) = ie;
    return_array(2) = js;
    return_array(3) = je;
    return_array(4) = ks;
    return_array(5) = ke;
    return_array(6) = ls;
    return_array(7) = le;
    return_array(8) = ms;
    return_array(9) = me;
    return return_array;
};
ndarray<int> slice(int is, int ie, int js, int je){
    int ks = 0; int ke = 0;
    int ls = 0; int le = 0;
    int ms = 0; int me = 0;
    ndarray<int> return_array(10);
    return_array(0) = is;
    return_array(1) = ie;
    return_array(2) = js;
    return_array(3) = je;
    return_array(4) = ks;
    return_array(5) = ke;
    return_array(6) = ls;
    return_array(7) = le;
    return_array(8) = ms;
    return_array(9) = me;
    return return_array;
};
ndarray<int> slice(int is, int ie, int js, int je, int ks, int ke){
    int ls = 0; int le = 0;
    int ms = 0; int me = 0;
    ndarray<int> return_array(10);
    return_array(0) = is;
    return_array(1) = ie;
    return_array(2) = js;
    return_array(3) = je;
    return_array(4) = ks;
    return_array(5) = ke;
    return_array(6) = ls;
    return_array(7) = le;
    return_array(8) = ms;
    return_array(9) = me;
    return return_array;
};
ndarray<int> slice(int is, int ie, int js, int je, int ks, int ke, int ls, int le){
    int ms = 0; int me = 0;
    ndarray<int> return_array(10);
    return_array(0) = is;
    return_array(1) = ie;
    return_array(2) = js;
    return_array(3) = je;
    return_array(4) = ks;
    return_array(5) = ke;
    return_array(6) = ls;
    return_array(7) = le;
    return_array(8) = ms;
    return_array(9) = me;
    return return_array;
};
ndarray<int> slice(int is, int ie, int js, int je, int ks, int ke, int ls, int le, int ms, int me){
    ndarray<int> return_array(10);
    return_array(0) = is;
    return_array(1) = ie;
    return_array(2) = js;
    return_array(3) = je;
    return_array(4) = ks;
    return_array(5) = ke;
    return_array(6) = ls;
    return_array(7) = le;
    return_array(8) = ms;
    return_array(9) = me;
    return return_array;
};

// Squeeze function to move singleton dimensions from array
template<class Type>
ndarray<Type> squeeze(ndarray<Type> x){
    int ni = 1;
    int nj = 1;
    int nk = 1;
    int nl = 1;
    int nm = 1;
    bool idim = true;
    bool jdim = true;
    bool kdim = true;
    bool ldim = true;
    bool mdim = true;
    int ndims;
    int datalength = sizeof(x.elements[0])*x.nelms;
    ndarray<Type> return_array;
    if ( x.ni > 1 ){
        ni = x.ni;
        ndims = 1;
    } else if ( x.nj > 1){
        ni = x.nj;
        ndims = 1;
        jdim = false;
    } else if (x.nk > 1){
        ni = x.nk;
        ndims = 1;
        jdim = false;
        kdim = false;
    } else if (x.nl > 1){
        ni = x.nl;
        ndims = 1;
        jdim = false;
        kdim = false;
        ldim = false;
    } else if (x.nm > 1){
        ni = x.nm;
        ndims = 1;
        jdim = false;
        kdim = false;
        ldim = false;
        mdim = false;
    };

    if ( x.nj > 1 && jdim == true ){
        nj = x.nj;
        ndims = 2;
    } else if (x.nk > 1 && kdim == true){
        nj = x.nk;
        ndims = 2;
        kdim = false;
    } else if (x.nl > 1 && ldim == true){
        nj = x.nl;
        ndims = 2;
        kdim = false;
        ldim = false;
    } else if (x.nm > 1 && mdim == true){
        nj = x.nm;
        ndims = 2;
        kdim = false;
        ldim = false;
        mdim = false;
    };

    if ( x.nk > 1 && kdim == true ){
        nk = x.nk;
        ndims = 3;
    } else if (x.nl > 1 && ldim == true){
        nk = x.nl;
        ndims = 3;
        ldim = false;
    } else if (x.nm > 1 && mdim == true){
        nk = x.nm;
        ndims = 3;
        ldim = false;
        mdim = false;
    };

    if ( x.nl > 1 && ldim == true ){
        nl = x.nl;
        ndims = 4;
    } else if (x.nm > 1 && mdim == true){
        nl = x.nm;
        ndims = 4;
        mdim = false;
    };

    if ( x.nm > 1 && mdim == true ){
        ndims = 5;
        nm = x.nm;
    };

    if ( ndims == 1 ){
        return_array.alloc(ni);
        memcpy(&return_array.elements[0],&x.elements[0],datalength);
    };
    if ( ndims == 2 ){
        return_array.alloc(ni,nj);
        memcpy(&return_array.elements[0],&x.elements[0],datalength);
    };
    if ( ndims == 3 ){
        return_array.alloc(ni,nj,nk);
        memcpy(&return_array.elements[0],&x.elements[0],datalength);
    };
    if ( ndims == 4 ){
        return_array.alloc(ni,nj,nk,nl);
        memcpy(&return_array.elements[0],&x.elements[0],datalength);
    };
    if ( ndims == 5 ){
        return_array.alloc(ni,nj,nk,nl,nm);
        memcpy(&return_array.elements[0],&x.elements[0],datalength);
    };
    return return_array;    
};

// ------ write array to file ------- //
template<class Type>
void write_array_to_file(const char* filename, const ndarray<Type> &x){
    char* headerBuffer;
    char* arrayBuffer;
    int nints = 7;
    headerBuffer = new char[nints*sizeof(int)]; 
    arrayBuffer = new char[x.nelms*sizeof(Type)]; 

    memcpy(&headerBuffer[0],  &x.ndims, sizeof(int));
    memcpy(&headerBuffer[4],  &x.nelms, sizeof(int));
    memcpy(&headerBuffer[8],  &x.ni,    sizeof(int));
    memcpy(&headerBuffer[12], &x.nj,    sizeof(int));
    memcpy(&headerBuffer[16], &x.nk,    sizeof(int));
    memcpy(&headerBuffer[20], &x.nl,    sizeof(int));
    memcpy(&headerBuffer[24], &x.nm,    sizeof(int));
    memcpy(&arrayBuffer[0], &x.elements[0], x.nelms*sizeof(Type));
    
    std::ofstream outFile (filename, std::ios::out | std::ios::binary);
    outFile.write(headerBuffer,nints*sizeof(int));
    outFile.write(arrayBuffer,x.nelms*sizeof(Type));
    outFile.close();
    
};

#endif
