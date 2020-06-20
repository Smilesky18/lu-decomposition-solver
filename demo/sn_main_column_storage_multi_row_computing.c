/*this demo is the same as demo.cpp except for that this demo calls the C interface*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "nicslu.h"
#include <sys/time.h>
#include <stdbool.h>
#define MICRO_IN_SEC 1000000.00

double Abs(double x)
{
	  return x < 0 ? -x : x;
} 
bool equal( double a, double b )
{
	  if ( Abs(a-b) < 0.1 )
	  {
		  return true;
	  }
	  else
	  {
		  return false;
	  }
}
int min ( int a, int b )
{
    if ( a > b ) return b;
    return a;
}
int detect ( int *asub, int *xa, int lower_col, int higher_col )
{
    int lower_col_ptr = xa[lower_col + 1] - 1;
    int higher_col_ptr = xa[higher_col + 1] - 1;
    int lower_col_count = xa[lower_col + 1] - xa[lower_col];
    int higher_col_count = xa[higher_col + 1] - xa[higher_col];
    int count = min(lower_col_count, higher_col_count) - 1;
    int i;

    for ( i = 0; i < count; i++ )
    {
        if ( asub[lower_col_ptr] == asub[higher_col_ptr] )
        {
            lower_col_ptr--;
            higher_col_ptr--;
            continue;
        }
        else
        {
          break;
        }
    }

    if ( i == count )
    {
        if ( asub[lower_col_ptr] == higher_col && asub[higher_col_ptr] == higher_col ) return 1;
    }
    return 0;
}
 

/* Time Stamp */
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

    /*pre-ordering (do only once)*/
    NicsLU_Analyze(solver, n, ax, ai, ap, MATRIX_COLUMN_REAL, NULL, NULL, NULL, NULL);

    /*create threads (do only once)*/
    NicsLU_CreateThreads(solver, 0); /*use all physical cores*/
    /*factor & solve (first-time)*/
    NicsLU_FactorizeMatrix(solver, ax, 0); /*use all created threads*/

    /* Get L/U structural information from NicSLU */
    lx = (_double_t *)malloc(sizeof(_double_t) * stat[9]);
    ux = (_double_t *)malloc(sizeof(_double_t) * stat[10]);
    li = (_uint_t *)malloc(sizeof(_uint_t) * stat[9]);
    ui = (_uint_t *)malloc(sizeof(_uint_t) * stat[10]);
    NicsLU_GetFactors(solver, lx, li, lp, ux, ui, up, sort, row_perm, col_perm, row_scale, col_scale);
    for (i = 0; i < n; i++) col_perm_inv[col_perm[i]] = i;
    for (i = 0; i < n; i++) row_perm_inv[row_perm[i]] = i;
    
    int lnz = (int)stat[10]; // Number of non-zeros in L 
    int unz = (int)stat[9];  // Numbe of non-zeros in U
    int *offset_U = (int *)malloc(sizeof(int) * (n+1)); // This array stores the column offset in U 
    int *offset_L = (int *)malloc(sizeof(int) * (n+1)); // This array stores the column offset in L
    for ( i = 0; i <= n; i++ ) offset_U[i] = lp[i];	// assignment: get lp from NicSLU 
    for ( i = 0; i <= n; i++ ) offset_L[i] = up[i];	// assignment: get up from NicSLU

    int *row_ptr_L = (int *)malloc(sizeof(int) * lnz);	// This array stores the row ptr in L
    int *row_ptr_U = (int *)malloc(sizeof(int) * unz);	// This array stores the row ptr in U
    for ( i = 0; i < lnz; i++ ) row_ptr_L[i] = ui[i];	// assignment: get ui from NicSLU
    for ( i = 0; i < unz; i++ ) row_ptr_U[i] = li[i];	// assignment: get li from NicSLU
    
    int *sn_record = (int *)malloc(sizeof(int) * n);	
    memset(sn_record, -1, sizeof(int) * n);
    int *sn_start = (int *)malloc( sizeof( int ) * n);
    int *sn_end = (int *)malloc( sizeof( int ) * n);
    int col_thresold = atoi(argv[2]);
    memset(sn_start, 0, sizeof(int)*n);
    memset(sn_end, 0, sizeof(int)*n);
    int sn_sum = 0;
    int sn_sum_final = 0;
    
    int lower = 0, higher;
   
    /* detect supernode */ 
    for ( i = 0; i < n-1; i++ )
    {
	higher = i + 1;
	if ( detect(row_ptr_L, offset_L, lower, higher) )
        {
		sn_start[sn_sum] = lower;
		sn_end[sn_sum] = higher;
        } 
	else
	{
		lower = higher;
		sn_sum++;
	}
    }
    
    //int sum = 0, sum1 = 0, sum2 = 0;
    int *sn_column_start = (int *)malloc(sizeof(int) * sn_sum);
    memset(sn_column_start, 0, sizeof(int) * sn_sum);
    int *sn_column_end = (int *)malloc(sizeof(int) * sn_sum);
    memset(sn_column_end, 0, sizeof(int) * sn_sum);
    int *sn_num_record = (int *)malloc(sizeof(int) * n);
    memset(sn_num_record, 0, sizeof(int) * n); 
    
    /* record supernode information */
    	// sn_sum_final: Number of SuperNodes.
	// sn_record[j] = the distance between column j and the end column in this supernode.
	// sn_num_record[j]: This array records which supernode column j belongs.
	// sn_column_start[sn]: The start column of supernode sn.
	// sn_column_end[sn]: The end column of supernode sn. 
    for ( i = 0; i <= sn_sum; i++ )
    {
	if ( sn_end[i] - sn_start[i] >= col_thresold ) 
	{
		for ( j = sn_start[i]; j <= sn_end[i]; j++ )
            	{
                	sn_record[j] = sn_end[i] - j;
			sn_num_record[j] = sn_sum_final;
            	}
		sn_column_start[sn_sum_final] = sn_start[i];
		sn_column_end[sn_sum_final] = sn_end[i];
            	sn_sum_final++;
	}
    }

    double *x_result, *x_real_result;
    int error_result = 0;
    double t_s, t_e;
    int thresold = atoi(argv[3]);
    
    double sum_time = 0;	
    for ( i = 0; i < 10; i++ )
    {
	    t_s = microtime();
	    x_result = lu_gp_sparse_supernode_dense_column_computing_v5_multi_row_computing(ax, ai, ap, n, lnz, unz, row_perm_inv, col_perm, row_ptr_L, offset_L, row_ptr_U, offset_U, sn_record, thresold, sn_num_record, sn_column_start, sn_column_end, sn_sum_final);
	    t_e = microtime() - t_s;
            printf("Time of LU_GP_ssn_column_storage_multi_row_computing: %g\n", t_e);
	    sum_time += t_e;
    }
    printf("Average Time of LU_GP_ssn_column_storage_multi_row_computing: %g\n", sum_time/10); 
    
    NicsLU_Solve(solver, b, x);
   // double *x_v5_multi_row_computing = (double *)malloc(sizeof(double) * n);
    for ( i = 0; i < n; i++ )
    {
	    x_real_result[row_perm_inv[i]] = x_result[i];
    } 
	
    for ( i = 0; i < n; i++ )
    {
	    if ( !equal(x_real_result[i], x[i]) )
	    {
		    //printf("nicslu[%d] = %.16f me[%d] = %.16f\n", i, x[i], i, x_v5_multi_row_computing[i]);
		    error_result++;
	    }
    } 

    printf("Errors are: %d\n", error_result);
    printf("***********SuperNode************\n");
    printf("Number of NicSLU's sn is: %lf\n", stat[12]);
    printf("Number of lu_gp's sn is: %d\n", sn_sum_final);
   // printf("nnz of L = %d nnz of U = %d\n", lnz, unz);
 
EXIT:
    free(ax);
    free(ai);
    free(ap);
    free(b);
    free(x);
    NicsLU_Free(solver);
    return 0;
}
