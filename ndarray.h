#ifndef ND_ARRAY_H
#define ND_ARRAY_H
#include <iostream>
#include <cstring>

/*  Template for contiguous nd-arrays (up to 5D) 
    some examples of initializing:
    1D: ndarray<float> x(ni);       --> 1-dimensional array of floats with ni elements
    2D: ndarray<int> x(ni,nj);      --> 2-dimensional array of ints  with ni*nj elements
    3D: ndarray<float> x(ni,nj,nk); --> 3-dimensional array of doubles with ni*nj*nk elements
    ndarray<double> x --> null ndarray of type double (needs to be assigned from existing ndarray)
    where elements are accessed by x(i,j,k) (for 3D), etc.
    ndarray<ndarry<Type>> --> ndarray of ndarrays of ...
    
    Function taking ndarray as input require taking array by reference, for example,
        void some_function(ndarray<Type> &input_array){};
    which will allow direct modification of original array within this function.
    
    *******
    DATA IS STORED AS A 1D CONTIGUOUS ARRAY WITH SPECIAL ACCESS FUNCTIONS TO GET
    ELEMENTS AT INDEPENDENT INDICES. THE DATA IS STORED IN 'COLUMN MAJOR' FORMAT
    FOR EASILY TRANSITIONING BETWEEN Fortran AND MATLAB 
    
    EXAMPLE: 
    3D LOOPS SHOULD BE STRUCTURES AS
    for (int k=0; k<nk; k++){
        for (int k=0; k<nk; k++){
            for (int k=0; k<nk; k++){
                x(i,j,k) = expression;
            }
        }
    }
    WHICH IS OPPOSITE OF TYPICAL C/C++ CONVENTIONS
    
    THIS IS NOT IDEAL FOR VECTORIZATION -- IF YOU WANT VECTORIZED CODE
    -->
    for (int k=0; k<nk; k++){
        for (int k=0; k<nk; k++){
            for (int k=0; k<nk; k++){
                x.elements[(k*x.nj+j)*x.ni+i] = expression;
            }
        }
    }
    *******

    Written by: Mark A. Herndon
    January 12, 2023
*/

// ---------- ndarray declarations ------------- //

template<class Type>
class ndarray {
    
    private:
        void allocate(); // Allocates memory for array elements

    public:
        int ni; // Number of elements in 1st index position
        int nj; // Number of elements in 2nd index position
        int nk; // Number of elements in 3rd index position
        int nl; // Number of elements in 4th index position
        int nm; // Number of elements in 5th index position
        
        int ndims; // Number of non-zero dimensions        
        int nelms; // Total number of elements
        int* __restrict__ size; // Gives shape of array {ni, nj, nk, nl, nm}
                   // size array is 5D for all array sizes
    
        Type* __restrict__ elements; // Pointer to elements to be defined 
                        // by new (1D array of contiguous memory)
    
        typedef ndarray<Type> ndArrayType; // Defines generic type of ndarray

        // --------- Constructors ---------- //
        inline ndarray(); // null constructor
        ndarray(int ni);  // 1D array constructor
        ndarray(int ni,int nj); // 2D array constructor
        ndarray(int ni,int nj,int nk); // 3D array constructor
        ndarray(int ni,int nj,int nk, int nl); // 4D array constructor
        ndarray(int ni,int nj,int nk, int nl, int nm); // 5D array constructor
        // Can add 6D, 7D, ... , nD array constructors but that is overkill 
        ndarray(const ndArrayType&); // Copy constructor

        // --------- Local operators --------- //
        Type& operator()(int i);
        Type& operator()(int i,int j);
        Type& operator()(int i, int j, int k);
        Type& operator()(int i, int j, int k, int l);
        Type& operator()(int i, int j, int k, int l, int m);

        void operator=(const ndarray<Type>&); // Assignment operator
        void operator=(const Type&);          // Assignment to value 

        ndarray<Type> operator()(ndarray<int> inds);
        
        ~ndarray(); // Destructor

       // ------------ allocate null object explicitly ------- //
        void alloc(int nii);
        void alloc(int nii, int nji);
        void alloc(int nii, int nji, int nki);
        void alloc(int nii, int nji, int nki, int nli);
        void alloc(int nii, int nji, int nki, int nli, int nmi);
};

// ---------------- Declarations of Global Functions ------------------ //
// Element wise addition
template<class Type>
ndarray<Type> operator+(const ndarray<Type>&, const ndarray<Type>&);
// Element wise subtraction
template<class Type>
ndarray<Type> operator-(const ndarray<Type>&, const ndarray<Type>&);
// Element wise multiplication
template<class Type>
ndarray<Type> operator*(const ndarray<Type>&, const ndarray<Type>&);
// Element wise division
template<class Type>
ndarray<Type> operator/(const ndarray<Type>&, const ndarray<Type>&);

// -----------------  Store in ndarray.cpp ? ----------------------//
// ----------------- Constructor Definitions ---------------------//

// null constructor will call assignment operator when creating an
// ndarray of an ndarray . . . yikes!
template<class Type>
ndarray<Type>::ndarray(){
    elements = nullptr;
    size = nullptr;
};

// 1D array constructor
template<class Type>
ndarray<Type>::ndarray(int ni) : ni(ni){
    nelms = ni;
    allocate();
    nj=1; nk=1; nl=1; nm=1;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    ndims = 1;
};

// 2D array
template<class Type>
ndarray<Type>::ndarray(int ni,int nj) : ni(ni), nj(nj){
    nelms = ni*nj;
    allocate();
    nk=1; nl=1; nm=1;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    ndims = 2;
};

// 3D array
template<class Type>
ndarray<Type>::ndarray(int ni,int nj,int nk) : ni(ni), nj(nj), nk(nk){
    nelms = ni*nj*nk;
    allocate();
    nl=1; nm=1;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    ndims = 3;
};

// 4D array
template<class Type>
ndarray<Type>::ndarray(int ni,int nj,int nk, int nl) : ni(ni), nj(nj), nk(nk), nl(nl){
    nelms = ni*nj*nk*nl;
    allocate();
    nm=1;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    ndims = 4;
};

// 5D array
template<class Type>
ndarray<Type>::ndarray(int ni,int nj,int nk, int nl, int nm) : ni(ni), nj(nj), nk(nk), nl(nl), nm(nm){
    nelms = ni*nj*nk*nl*nm;
    allocate();
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    ndims = 5;
};

// Copy constructor
template<class Type>
ndarray<Type>::ndarray(const ndarray<Type>& v) : nelms(v.nelms), ni(v.ni), nj(v.nj), nk(v.nk), nl(v.nl), nm(v.nm), ndims(v.ndims){
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    if (v.elements)
    {
        allocate();
        for (int i=0; i<nelms; i++){
            elements[i] = v.elements[i];
        }
    };
};
template<class Type>
void ndarray<Type>::alloc(int nii){
    ni = nii;
    nj = 1;
    nk = 1;
    nl = 1;
    nm = 1;
    nelms = ni;
    ndims = 1;
    size = new int[5];
    size[0] = ni;
    size[1] = 1;
    size[2] = 1;
    size[3] = 1;
    size[4] = 1;
    allocate();
};
template<class Type>
void ndarray<Type>::alloc(int nii, int nji){
    ni = nii;
    nj = nji;
    nk = 1;
    nl = 1;
    nm = 1;
    nelms = ni*nj;
    ndims = 2;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = 1;
    size[3] = 1;
    size[4] = 1;
    allocate();
};
template<class Type>
void ndarray<Type>::alloc(int nii, int nji, int nki){
    ni = nii;
    nj = nji;
    nk = nki;
    nl = 1;
    nm = 1;
    nelms = ni*nj*nk;
    ndims = 3;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    allocate();
};
template<class Type>
void ndarray<Type>::alloc(int nii, int nji, int nki, int nli){
    ni = nii;
    nj = nji;
    nk = nki;
    nl = nli;
    nm = 1;
    nelms = ni*nj*nk*nl;
    ndims = 4;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    allocate();
};

template<class Type>
void ndarray<Type>::alloc(int nii, int nji, int nki, int nli, int nmi){
    ni = nii;
    nj = nji;
    nk = nki;
    nl = nli;
    nm = nmi;
    nelms = ni*nj*nk*nl*nm;
    ndims = 5;
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;
    allocate();
};

// ------------- Private member functions -------------- //
// Memory allocation for elements of ndarray
template<class Type>
void ndarray<Type>::allocate(){
    elements = new Type[nelms];
    const int n = nelms;
    for (int i=0; i<n; i++){
        elements[i] = (Type) 0;
    };
};

// Destructor (Free pointers to dynamically allocated memory)
template<class Type>
ndarray<Type>::~ndarray(){
    if (elements)
    {
        delete[] elements;
       // _mm_free(elements);
    };
    if (size)
    {
        delete[] size;
    }
};

// ------------------ Member operators for array access ---------- //
// WARNING: THIS DOES NOT CHECK FOR SIZE CONSISTENCY, ONLY NUMBER OF ELEMENTS //

// 1D array access array(i)
template<class Type>
Type& ndarray<Type>::operator()(int i){
    if ( i > nelms ){
        std::cout << "WARNING, INDEX " << i << " EXCEEDS (FLATTENED) ARRAY DIMENSIONS\n";
    }
    return elements[i];
};

// 2D array access array(i,j)
template<class Type>
Type& ndarray<Type>::operator()(int i,int j){
    int ind = j*ni+i;
    if ( ind > nelms ){
        std::cout << "WARNING, INDEX " << ind << " EXCEEDS (FLATTENED) ARRAY DIMENSIONS\n";
        std::cout << "CHECK ARRAY SIZE\n";
    }
    return elements[ind];
};

// 3D array access array(i,j,k)
template<class Type>
Type& ndarray<Type>::operator()(int i, int j, int k){
    int ind = (k*nj+j)*ni+i;
    if ( ind > nelms ){
        std::cout << "WARNING, INDEX " << ind << " EXCEEDS (FLATTENED) ARRAY DIMENSIONS\n";
        std::cout << "CHECK ARRAY SIZE\n";
    }
    return elements[ind];
};

// 4D array access array(i,j,k,l)
template<class Type>
Type& ndarray<Type>::operator()(int i, int j, int k, int l){
    int ind = ((l*nk+k)*nj+j)*ni+i;
    if ( ind > nelms ){
        std::cout << "WARNING, INDEX " << ind << " EXCEEDS (FLATTENED) ARRAY DIMENSIONS\n";
        std::cout << "CHECK ARRAY SIZE\n";
    }
    return elements[ind];
};

// 5D array access array(i,j,k,l,m)
template<class Type>
Type& ndarray<Type>::operator()(int i, int j, int k, int l, int m){
    int ind = (((m*nl+l)*nk+k)*nj+j)*ni+i;
    if ( ind > nelms ){
        std::cout << "WARNING (5D), INDEX " << ind << " EXCEEDS (FLATTENED) ARRAY DIMENSIONS\n";
        std::cout << "CHECK ARRAY SIZE\n";
    }
    return elements[ind];
};

// Assignment to value 
template<class Type>
void ndarray<Type>::operator=(const Type& val){
    if (elements){
        for (int i=0; i<nelms; i++){
            elements[i] = val;
        }
    }
};

/* inds strucuture:
    inds.elements[0] = istart
    inds.elements[1] = iend
    inds.elements[2] = jstart
    inds.elements[3] = jend
    inds.elements[4] = kstart
    inds.elements[5] = kend
    inds.elements[6] = lstart
    inds.elements[7] = lend
    inds.elements[8] = mstart
    inds.elements[9] = mend
*/
template<class Type>
ndarray<Type> ndarray<Type>::operator()(ndarray<int> inds){
    int is = inds.elements[0];
    int ie = inds.elements[1];
    int js = inds.elements[2];
    int je = inds.elements[3];
    int ks = inds.elements[4];
    int ke = inds.elements[5];
    int ls = inds.elements[6];
    int le = inds.elements[7];
    int ms = inds.elements[8];
    int me = inds.elements[9];
    int nitmp = ie-is+1;
    int njtmp = je-js+1;
    int nktmp = ke-ks+1;
    int nltmp = le-ls+1;
    int nmtmp = me-ms+1;
    int ind;
    
    if (ie == ni){
        ie=ie-1;
        nitmp = nitmp-1;
    };
    if (je == nj){
        je=je-1;
        njtmp = njtmp-1;
    };
    if (ke == nk){
        ke=ke-1;
        nktmp = nktmp-1;
    };
    if (le == nl){
        le=le-1;
        nltmp = nltmp-1;
    };
    if (me == nm){
        me=me-1;
        nmtmp = nmtmp-1;
    };
    
    ndarray<Type> tmpArray(nitmp, njtmp, nktmp, nltmp, nmtmp);
    ndarray<Type> return_array;
    
    int ii, jj, kk, ll, mm;
    for(int m=ms; m<=me; m++){
        for(int l=ls; l<=le; l++){
            for(int k=ks; k<=ke; k++){
                for(int j=js; j<=je; j++){
                    for(int i=is; i<=ie; i++){
                        ii = i - is;
                        jj = j - js;
                        kk = k - ks;
                        ll = l - ls;
                        mm = m - ms;
                        ind = (((m*nl+l)*nk+k)*nj+j)*ni+i;
                        tmpArray(ii,jj,kk,ll,mm) = elements[ind];
                    }
                }
            }
        }
    };
    return_array = squeeze(tmpArray);
    return return_array;
};

// ----------------- Global functions -----------------//

// Assignment operator 
/* Assignment will allocate memory for a null ndarray to 
   correspond to data structure of assigning array */
template<class Type>
void ndarray<Type>::operator=(const ndarray<Type>& v){
    if (elements && nelms == 0) { 
        delete[] elements;
        delete[] size;
        elements = nullptr;
        size = nullptr;
    };
    if (!elements){
    ni = v.ni;
    nj = v.nj;
    nk = v.nk;
    nl = v.nl;
    nm = v.nm;
    nelms = v.nelms;
    ndims = v.ndims;
    allocate();
    size = new int[5];
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = nl;
    size[4] = nm;

        memcpy(&elements[0],&v.elements[0],sizeof(Type)*nelms);
    } else if (elements && nelms > 0) { 
        // If elements exist, need to check dimension consistency  
        memcpy(&elements[0],&v.elements[0],sizeof(Type)*nelms);
    };
};

template<class Type>
ndarray<Type> operator+(const ndarray<Type>& A, const ndarray<Type>& B){
    ndarray<Type> Z(A.nelms);
    if ( A.size[0] != B.size[0] || A.size[1] != B.size[1] ||
         A.size[2] != B.size[2] || A.size[3] != B.size[3] ||
         A.size[4] != B.size[4] ) {
        std::cout << "WARNING, ARRAY SHAPES DO NOT MATCH\n" << "\n";
        abort();
    }
    for (int i=0; i<Z.nelms; i++){
        Z.elements[i] = A.elements[i] + B.elements[i];
    }
    return Z;
};

template<class Type>
ndarray<Type> operator-(const ndarray<Type>& A, const ndarray<Type>& B){
    ndarray<Type> Z(A.nelms);
    if ( A.size[0] != B.size[0] || A.size[1] != B.size[1] ||
         A.size[2] != B.size[2] || A.size[3] != B.size[3] ||
         A.size[4] != B.size[4] ) {
        std::cout << "WARNING, ARRAY SHAPES DO NOT MATCH\n" << "\n";
        abort();
    }
    for (int i=0; i<Z.nelms; i++){
        Z.elements[i] = A.elements[i] - B.elements[i];
    }
    return Z;
};

template<class Type>
ndarray<Type> operator*(const ndarray<Type>& A, const ndarray<Type>& B){
    ndarray<Type> Z(A.nelms);
    if ( A.size[0] != B.size[0] || A.size[1] != B.size[1] || 
         A.size[2] != B.size[2] || A.size[3] != B.size[3] ||
         A.size[4] != B.size[4] ) {
        std::cout << "WARNING, ARRAY SHAPES DO NOT MATCH\n" << "\n";
        abort();
    }
    for (int i=0; i<Z.nelms; i++){
        Z.elements[i] = A.elements[i] * B.elements[i];
    }
    return Z;
};

template<class Type>
ndarray<Type> operator/(const ndarray<Type>& A, const ndarray<Type>& B){
    ndarray<Type> Z(A.nelms);
    if ( A.size[0] != B.size[0] || A.size[1] != B.size[1] ||
         A.size[2] != B.size[2] || A.size[3] != B.size[3] ||
         A.size[4] != B.size[4] ) {
        std::cout << "WARNING, ARRAY SHAPES DO NOT MATCH\n" << "\n";
        abort();
    }
    for (int i=0; i<Z.nelms; i++){
        Z.elements[i] = A.elements[i] / B.elements[i];
    }
    return Z;
};

// -------------- ndarrayUtils.h ------------ //
#include "ndarrayUtils.h"

#endif
