/*this demo is the same as demo.cpp except for that this demo calls the C interface*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "nicslu.h"
#include <sys/time.h>
#include <stdbool.h>
#define MICRO_IN_SEC 1000000.00

double microtime()
{
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv,&tz);

	return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}

const char *const ORDERING_METHODS[] = { "AMD", "mAMD", "AMF","mAMF1","mAMF2","mAMF3" };

int main( int argc[], char *argv[])
{
    _double_t *ax = NULL, *b = NULL, *x = NULL;
    _uint_t *ai = NULL, *ap = NULL;
    _uint_t n, row, col, nz, nnz, i, j;
    _handle_t solver = NULL;
    _double_t res[4], cond, det1, det2, fflop, sflop;
    size_t mem;
    _double_t *cfg;
    const _double_t *stat;
    const char *last_err;
	_uint_t *row_perm, *row_perm_inv, *col_perm, *col_perm_inv, *ui, *li;
	_double_t *lx, *ux, *row_scale, *col_scale;
	_size_t *lp, *up;
	_bool_t sort = 1;

    /*read matrix A*/
    if (__FAIL(ReadMatrixMarketFile(argv[1], &row, &col, &nz, NULL, NULL, NULL, NULL, NULL, NULL)))
    {
        printf("Failed to read matrix A\n");
        goto EXIT;
    }
    n = row;
    nnz = nz;
    ax = (_double_t *)malloc(sizeof(_double_t)*nnz);
    ai = (_uint_t *)malloc(sizeof(_uint_t)*nnz);
    ap = (_uint_t *)malloc(sizeof(_uint_t)*(1 + n));
	row_perm = (_uint_t *)malloc(sizeof(_uint_t)*n);
	col_perm = (_uint_t *)malloc(sizeof(_uint_t)*n);
	col_perm_inv = (_uint_t *)malloc(sizeof(_uint_t)*n);
	row_perm_inv = (_uint_t *)malloc(sizeof(_uint_t)*n);
	row_scale = (_double_t *)malloc(sizeof(_double_t)*n);
	col_scale = (_double_t *)malloc(sizeof(_double_t)*n); 
	lp = (_size_t *)malloc(sizeof(_size_t)*(1+n));
	up = (_size_t *)malloc(sizeof(_size_t)*(1+n));
	
    ReadMatrixMarketFile(argv[1], &row, &col, &nz, ax, ai, ap, NULL, NULL, NULL);
    printf("***********%s: row %d, col %d, nnz %d\n", argv[1], n, n, nnz);

    /*read RHS B*/
    b = (_double_t *)malloc(sizeof(_double_t)*n);
	for ( i = 0; i < n; i++ ) b[i] = 1.0;
    x = (_double_t *)malloc(sizeof(_double_t)*n);
    memset(x, 0, sizeof(_double_t) * n);

    /*initialize solver*/
    if (__FAIL(NicsLU_Initialize(&solver, &cfg, &stat, &last_err)))
    {
        printf("Failed to initialize\n");
        goto EXIT;
    }
    cfg[0] = 1.; /*enable timer*/
//    cfg[3] = 4;
    /*pre-ordering (do only once)*/
    NicsLU_Analyze(solver, n, ax, ai, ap, MATRIX_COLUMN_REAL, NULL, NULL, NULL, NULL);

    /*create threads (do only once)*/
    NicsLU_CreateThreads(solver, 0); /*use all physical cores*/
    /*factor & solve (first-time)*/
    NicsLU_FactorizeMatrix(solver, ax, 0); /*use all created threads*/
	
	
 	double *l_gp, *l_gp_sn_v6, *l_gp_sn_v5_row_computing, *l_gp_sn_v5_multi_row_computing;
    int l_gp_error = 0, l_gp_sn_v6_error = 0, l_gp_sn_v5_row_computing_error = 0, l_gp_sn_v5_multi_row_computing_error = 0;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	double start_GP, end_GP;
	int thresold = atoi(argv[3]);
	
	 double time_fact_start, time_fact_end;
	   double sum_time_nic = 0;
         for ( i = 0; i < 10; i++ )
	 {
	 time_fact_start = microtime();
		 NicsLU_ReFactorize(solver, ax, 0); 
	 time_fact_end = microtime() - time_fact_start;
	 printf("Time of nicslu: %g\n", time_fact_end);
	 sum_time_nic += time_fact_end;
	 }
	 printf("Average Time of nicslu: %g\n", sum_time_nic/10);
	////printf("***********SuperNode************\n");
	 printf("sn of NICSLU is: %lf\n", stat[12]);
 
    /*finally, print some statistical information*/

EXIT:
    free(ax);
    free(ai);
    free(ap);
    free(b);
    free(x);
    NicsLU_Free(solver);
    return 0;
}
