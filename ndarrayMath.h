#ifndef _NDARRAY_MATH_H
#define _NDARRAY_MATH_H_
#include "ndarray.h"
#include <cmath>

// ---------- Declarations ---------- //

// Element wise functions
template<class Type>
ndarray<Type> vSin(const ndarray<Type> &x);

template<class Type>
ndarray<Type> vCos(const ndarray<Type> &x);

// These templates should be restricted to floating point types //
template<class Type>
void cubic_spline_interpolation(ndarray<Type> &x, ndarray<Type> &y, ndarray<Type> &xq, ndarray<Type> &yq);

template<class Type>
void cubic_spline_interpolation_2D(ndarray<Type> &x, ndarray<Type> &y, ndarray<Type> &Z, ndarray<Type> &xq, ndarray<Type> &yq, ndarray<Type> &ZQ);

template<class Type>
void cubic_spline_coefs(ndarray<Type> &x, ndarray<Type> &y, ndarray<Type> &coefs);

template<class Type>
void look_up_inds(ndarray<Type> &x, Type xq, int &ileft, int &iright, int &mOld);

template<class Type>
void initial_loop_up_guess(int &ileft,int &iright,int &mOld,const ndarray<Type> &x,const Type &xq);

// ------------ Definitions -------------- //
template<class Type>
ndarray<Type> vSin(const ndarray<Type> &x){
    ndarray<Type> y;
    // Copy x to y to initalize to same shape
    y = x;
    for (int i=0; i<x.nelms; i++){
        //y.elements[i] = std::sin(x.elements[i]);
        y.elements[i] = sin(x.elements[i]);
    };
    return y;
};

template<class Type>
ndarray<Type> vCos(const ndarray<Type> &x){
    ndarray<Type> y;
    // Copy x to y to initalize to same shape
    y = x;
    for (int i=0; i<x.nelms; i++){
        //y.elements[i] = std::cos(x.elements[i]);
        y.elements[i] = cos(x.elements[i]);
    };
    return y;
};

template<class Type>
void cubic_spline_interpolation(ndarray<Type> &x, ndarray<Type> &y, ndarray<Type> &xq, ndarray<Type> &yq){
    ndarray<Type> coefs;
    Type xqq, xl, xr, xt1, xt2, xt3, xt4;
    double ts;
    double te;
    int ileft, iright, ind;
    int mOld = 0;   // Initialize guess for lookup
    ileft = 0;
    iright = x.ni-1;
    
    {
    cubic_spline_coefs(x, y, coefs);
    for (int n=0; n<xq.ni; n++){
        look_up_inds(x,xq(n),ileft,iright,mOld);
        if ( ileft == iright ){
            yq(n) = y(ileft);
            continue;
        };
        ind = ileft;
        xqq = xq(n);
        xr  = x(iright);
        xl  = x(ileft);
        xt1 = (xqq-xl);
        xt2 = xt1*xt1*xt1;
        xt3 = (xr-xqq);
        xt4 = xt3*xt3*xt3;
        yq(n) = coefs(ind,0)*(xqq-xl) + 
                coefs(ind,1)*(xr-xqq) + 
                coefs(ind,2)*xt2      + 
                coefs(ind,3)*xt4;
        };
    };
};

template<class Type>
void cubic_spline_interpolation_2D(ndarray<Type> &x, ndarray<Type> &y, ndarray<Type> &Z, ndarray<Type> &xq, ndarray<Type> &yq, ndarray<Type> &ZQ){
    
    int ni  = Z.ni;
    int nj  = Z.nj;
    int niq = ZQ.ni;
    int njq = ZQ.nj;

    ndarray<Type> ZQTMP(niq,nj);
    #pragma omp parallel shared(ZQTMP)
    {
    ndarray<Type> xqtmp(niq);
    ndarray<Type> yqtmp(njq);
    ndarray<Type> zqtmp(niq);
    ndarray<Type> zqtmp2(njq);
    ndarray<Type> ztmp(ni);
    ndarray<Type> ztmp2(nj);

    #pragma omp for
        for (int j=0; j<nj; j++){
            memcpy(&ztmp.elements[0],&Z.elements[j*ni],sizeof(Type)*ni);
            cubic_spline_interpolation(x,ztmp,xq,zqtmp);
            memcpy(&ZQTMP.elements[j*niq],&zqtmp.elements[0],sizeof(Type)*niq);
        };

    #pragma omp barrier

    #pragma omp for
        for (int i=0; i<niq; i++){
            for (int j=0; j<nj; j++){
                ztmp2.elements[j] = ZQTMP.elements[j*niq+i]; 
            };

            cubic_spline_interpolation(y,ztmp2,yq,zqtmp2);

            for (int j=0; j<njq; j++){
                ZQ.elements[j*niq+i] = zqtmp2.elements[j];
            };
        };
    };
};

template<class Type>
void cubic_spline_coefs(ndarray<Type> &x, ndarray<Type> &y, ndarray<Type> &coefs){
    ndarray<Type> h(x.ni);
    ndarray<Type> mu(x.ni);
    ndarray<Type> lb(x.ni);
    ndarray<Type> di(x.ni);
    ndarray<Type> d(x.ni);
    ndarray<Type> M(x.ni);
    di = (Type) 2.0;
    coefs.alloc(x.ni-1,4);
    int ii;
    int np = x.ni;
    double ts, te;
    
    Type xdiff, xi, yi;
    
   // #pragma omp for
    for (int i=1; i<np; i++){
        h(i) = x(i) - x(i-1);
    };
    for (int i=1; i<np-1; i++){
        xdiff = x.elements[i+1]-x.elements[i-1];
        xi    = x.elements[i];
        yi    = y.elements[i];
	    mu.elements[i] = h.elements[i]/(h.elements[i]+h.elements[i+1]);
	    lb.elements[i] = 1.0 - mu.elements[i];
	    d.elements[i]  =  (y.elements[i+1] - yi) / ((x.elements[i+1]-xi)*xdiff) -
                          (yi-y.elements[i-1]) / ((xi-x.elements[i-1])*xdiff);
    };
    
    Type w;
    for (int i=1; i<np; i++){
	    w = mu(i)/di(i-1);
	    di(i) = di(i) - w*lb(i-1);
	    d(i)  = d(i)  - w*d(i-1);
    }; 

    M(np) = d(np)/di(np);
    
    for (int i=np-2; i>=0; i--){
        M.elements[i] = (d.elements[i] - mu.elements[i]*M.elements[i+1])/di.elements[i];
    };

    Type hh, hh2, Mi, Mi1, yi1;
    for (int i=0; i<np-1; i++){
        ii = i+1;
        hh = h.elements[ii];
        hh2 = hh*hh;
        Mi  = M.elements[ii];
        Mi1 = M.elements[ii-1];
        yi  = y.elements[ii];
        yi1 = y.elements[ii-1];
        coefs.elements[i] = (1/hh)*(yi-(Mi*hh2)/6.0); 
        coefs.elements[(np-1)+i] = (1/hh)*(yi1-(Mi1*hh2)/6.0);
        coefs.elements[2*(np-1)+i] = Mi/(6.0*hh);
        coefs.elements[3*(np-1)+i] = Mi1/(6.0*hh);
    };
};

template<class Type>
void look_up_inds(ndarray<Type> &x, Type xq, int &ileft, int &iright, int &mOld){
    ileft = mOld;
    iright = x.ni-1;
    int m;
    
    if ( mOld == -1 ){
        ileft = 0;
        mOld = ((ileft+iright)/2);
    } else {
        // Calculate initial guess for smaller search range
        initial_loop_up_guess(ileft,iright,mOld,x,xq);
    };
    if ( xq == x(ileft) ){
        iright = ileft;
        mOld = ileft;
        return;
    } else if ( xq == x(iright) ) {
        ileft = iright;
        mOld = iright;
        return;
    } else if ( xq == x(mOld) ) {
        ileft = mOld;
        iright = mOld;
        return;
    }; 
    int loop_counter = 0; 
    while ( iright-ileft != 1 && iright != ileft){
        if ( loop_counter == 0 ){
            m = mOld;
        } else {
            m = ((ileft+iright)/2);
        };
        if ( xq < x(m) ){
            iright = m;
            mOld = ileft;
        } else if ( xq > x(m) ){
            ileft = m;
            mOld = ileft;
        } else if ( xq == x(m) ){
            ileft = m;
            iright = m;
            mOld = ileft;
            return;
        };
        loop_counter = loop_counter + 1;
    };
};

template<class Type>
void initial_loop_up_guess(int &ileft,int &iright,int &mOld,ndarray<Type> &x,const Type &xq){
    if ( mOld < x.ni-3 ){  
        if ( xq <= x(mOld+1) && xq >= x(mOld)){
            ileft = mOld;
            iright = mOld + 1;
        } else if ( xq <= x(mOld+2) && xq >= x(mOld+1) ){
            ileft = mOld+1;
            iright = mOld+2;
            mOld = mOld+1;
        };
    };
    if ( mOld >= x.ni-3 ){
        ileft = x.ni-3;
        iright = x.ni-1;
        mOld = ileft;
    };
};
#endif
