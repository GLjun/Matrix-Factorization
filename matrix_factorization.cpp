/*************************************************************************
	> File Name: matrix_factorization.cpp
	> Author: 
	> Mail: 
	> Created Time: Sun Jul  8 16:42:13 2018
 ************************************************************************/

#include <iostream>
#include <cstdio>
#include <random>
#include <functional>
#include <cstdlib>
#include <cmath>
#include "sparse_matrix_csr.h"
#include "util.h"

#define K 100
#define delta 0.15

using namespace std;

template <typename T>
void calculate_pq(SparseMatrix<T>& sm, T p[][K], T qt[][K], const T DELTA)
{
    int nrow = sm.nrow, ncol = sm.ncol;

    cout << nrow <<" " << ncol;

    std::random_device rd;
    std::minstd_rand gen(rd());
    std::uniform_int_distribution<> dis(1, 10000);
    auto dice = std::bind(dis, gen);

    //init p
    for(int i = 0;i < nrow; i++)
        for(int j = 0;j < K; j++)
            p[i][j] = dice()/10000.0;

    //init qt
    for(int i = 0;i < ncol; i++)
        for(int j = 0;j < K; j++)
            qt[i][j] = dice()/10000.0;

    const double alpha = 0.0002;
    const double beta = 0.01;


    T eij = T();
    int col = 0;
    double error = DELTA;
    int cnt = 0;
    //while(error >= DELTA)
    TIME_T start, end;
    MARK_TIME(start);
    while(error >= DELTA)
    {
        error = 0.0;
        for(int i = 0;i < nrow; i ++)
        {
            for(int j = sm.row_ptr[i];j < sm.row_ptr[i+1]; j++)
            {
                eij = sm.values[j];
                col = sm.col_idx[j];
                for(int k = 0; k < K; k ++)
                    eij -= p[i][k] * qt[col][k];
                error += eij * eij;

                for(int k = 0; k < K; k ++)
                {
                    p[i][k] = p[i][k] + alpha * 2.0*eij*qt[col][k] - alpha* beta*p[i][k];
                    qt[col][k] = qt[col][k] + alpha * 2.0*eij*p[i][k] - alpha * beta*qt[col][k];
                }
            }
        }
        
        cnt ++;
        printf("\riteration %d, diff error %.7f <-> %.7f", cnt, error, DELTA);
        fflush(stdout);
    }
    MARK_TIME(end);
    LOG("total time %.5f s, time per iter %.5f s", DIFF_TIME(end, start), DIFF_TIME(end, start)/cnt);
}

int main(int argc, char** argv)
{
    transfer_movieslen_to_csr(argv[1], argv[2]);


    SparseMatrix<double> sm;
    sm.read_csr(argv[2]);

    //matrix P
    double (*P)[K] = new double[sm.nrow][K];
    //permutation of Q
    double (*QT)[K] = new double[sm.ncol][K];

    const double DELTA = sm.nnz * delta;
    calculate_pq(sm, P, QT, DELTA);

    delete [] P;
    delete [] QT;

    return 0;
}

