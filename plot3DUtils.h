#ifndef _PLOT_3D_UTILS_H_
#define _PLOT_3D_UTILS_H_

#include "ndarray.h"
#include <fstream>

// ------------- Declarations ------------- //
void read_grid_file(const char* filename, ndarray<float> &x, ndarray<float> &y, ndarray<float> &z, int &ni, int &nj, int &nk, int &ng);
void write_grid_file(const char* filename, const ndarray<float> &x,const ndarray<float> &y, const ndarray<float> &z, const int &ni, const int &nj, const int &nk, const int &ng);
void read_solution_file(const char* filename, ndarray<float> &q, ndarray<float> &dat);

// ------------ Definitions -------------- //
void write_grid_file(const char* filename, const ndarray<float> &x,const ndarray<float> &y, const ndarray<float> &z, const int &ni, const int &nj, const int &nk, const int &ng){
    char* headerBuffer;
    char* gridBuffer;
    int nints=4;
    std::ofstream gridFile (filename, std::ios::out | std::ios::binary);
    headerBuffer = new char[sizeof(int)*4];
    memcpy(&headerBuffer[0],  &ng, sizeof(int));
    memcpy(&headerBuffer[4],  &ni, sizeof(int));
    memcpy(&headerBuffer[8],  &nj, sizeof(int));
    memcpy(&headerBuffer[12], &nk, sizeof(int));
    gridFile.write(headerBuffer,sizeof(int)*4);
    int nelms = ni*nj*nk*ng;
    int dataLength = 4*3*nelms;
    gridBuffer = new char[dataLength];
    memcpy(&gridBuffer[0],&x.elements[0],sizeof(float)*nelms);    
    int offset = sizeof(float)*nelms;
    memcpy(&gridBuffer[offset],&y.elements[0],sizeof(float)*nelms);    
    offset = 2*offset;
    memcpy(&gridBuffer[offset],&z.elements[0],sizeof(float)*nelms);    
    gridFile.write(gridBuffer,dataLength);
    gridFile.close();

    delete [] headerBuffer;
    delete [] gridBuffer;
    headerBuffer = nullptr;
    gridBuffer = nullptr;
};

void read_grid_file(const char* filename, ndarray<float> &x, ndarray<float> &y, ndarray<float> &z, int &ni, int &nj, int &nk, int &ng){
    char* headerBuffer;
    char* gridBuffer;
    int nints=4;
    std::ifstream gridFile (filename, std::ios::in | std::ios::binary);
    headerBuffer = new char[sizeof(int)*4];
    gridFile.read(headerBuffer,sizeof(int)*nints);
    memcpy(&ng,&headerBuffer[0], sizeof(int));
    memcpy(&ni,&headerBuffer[4], sizeof(int));
    memcpy(&nj,&headerBuffer[8], sizeof(int));
    memcpy(&nk,&headerBuffer[12],sizeof(int));
    int nelms = ni*nj*nk*ng;
    int dataLength = 4*3*nelms;
    gridBuffer = new char[dataLength];
    gridFile.read(gridBuffer,dataLength);
    gridFile.close();
    
    x.alloc(ni,nj,nk);
    y.alloc(ni,nj,nk);
    z.alloc(ni,nj,nk);

    memcpy(&x.elements[0],&gridBuffer[0],sizeof(float)*nelms);    
    int offset = sizeof(float)*nelms;
    memcpy(&y.elements[0],&gridBuffer[offset],sizeof(float)*nelms);    
    offset = 2*offset;
    memcpy(&z.elements[0],&gridBuffer[offset],sizeof(float)*nelms);    
    
    delete [] headerBuffer;
    delete [] gridBuffer;
    headerBuffer = nullptr;
    gridBuffer = nullptr;
};

void read_solution_file(const char* filename, ndarray<float> &q, ndarray<float> &dat){
    char* headerBuffer;
    char* datBuffer;
    char* solutionBuffer;
    int ng;
    int ni;
    int nj;
    int nk;
    int nvars = 5;
    int nints = 4;

    std::ifstream solutionFile (filename, std::ios::in | std::ios::binary);
    headerBuffer = new char[sizeof(int)*4];
    datBuffer = new char[sizeof(float)*4];

    solutionFile.read(headerBuffer,sizeof(int)*nints);
    memcpy(&ng,&headerBuffer[0], sizeof(int));
    memcpy(&ni,&headerBuffer[4], sizeof(int));
    memcpy(&nj,&headerBuffer[8], sizeof(int));
    memcpy(&nk,&headerBuffer[12],sizeof(int));
    int nelms = ni*nj*nk*ng;
    int dataLength = sizeof(float)*nvars*nelms;
    solutionBuffer = new char[dataLength];
    solutionFile.read(datBuffer, sizeof(float)*4);
    solutionFile.read(solutionBuffer,dataLength);
    solutionFile.close(); 
    
    dat.alloc(4); 
    q.alloc(ni,nj,nk,nvars);
    memcpy(&dat.elements[0],&datBuffer[0],4*sizeof(float));
    memcpy(&q.elements[0],&solutionBuffer[0],dataLength);
    
    delete [] headerBuffer;
    delete [] datBuffer;
    delete [] solutionBuffer;
    headerBuffer = nullptr;
    datBuffer = nullptr;
    solutionBuffer = nullptr;
};

#endif
