# include <stdio.h>
# include <stdlib.h>
# include <float.h>
# include <string.h>
# include <immintrin.h>
#include <sys/time.h>
#define MICRO_IN_SEC 1000000.00

  typedef __attribute__((aligned(64))) union
  {
    __m512d vec;
    double ptr_vec[8];
  }v2df_t;

  typedef union
  {
    __m256i vec;
    int ptr_vec[8];
  }v2if_t;

int mmm_multi_value( int a, int b )
{
	if ( a < b ) return a;
	return b;
}

/*double microtime()
{
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv,&tz);

	return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}*/


double* lu_gp_sparse_supernode_dense_column_computing_v5_multi_row_computing(double *a, int *row_ptr, int *offset, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *row_ptr_L, int *offset_L, int *row_ptr_U, int *offset_U, int *sn_record, int thresold, int *sn_num_record, int *sn_column_start, int *sn_column_end, int sn_sum)
{
  double *L, *U, *xx;
  double U_diag;
  int j, k, current_column;
  int i;
  L = ( double * )_mm_malloc(sizeof(double) * nzl, 64);
  U = ( double * )_mm_malloc(sizeof(double) * nzu, 64);
  xx = ( double *)_mm_malloc(sizeof(double) * n, 64);
  
  for ( i = 0; i < n; i++ )
  {
    xx[i] = 0;
    L[offset_L[i]] = 1.0;
  }

  int  val, columns, column_end, pack, column_start, column_number_sn, column_sn_end, row_num, row_count, column_divid, divid, row_sn_start, qqq, row_else1, row_else2, dense_vec_counter, sss, dividd, row_record, row_record_2;
  double temp, *dense_vec;
  dense_vec = ( double * )_mm_malloc(sizeof(double) * n, 64);
  memset(dense_vec, 0, sizeof(double) * n);

  v2df_t v_l, v_row, v_mul, v_sub, v, add_sum, zero, v_l_2, v_sub_2, v_2, add_sum_2;
  v2if_t vi, vi_2;
  zero.vec =  _mm512_setzero_pd();
  
  for ( k = 0; k < n; k++ )
  {
	  current_column = perm_c[k];
	  for ( j = offset[current_column]; j < offset[current_column+1]; j++ )
	  {
		  xx[perm_r[row_ptr[j]]] = a[j];
	  }
	  columns = offset_U[k+1] - offset_U[k] - 1;
	  column_end = row_ptr_U[offset_U[k+1] - 2];
	  
	  for ( j = 0; j < columns; j+=pack )
	  {
		  val = j+offset_U[k];
		  column_start = row_ptr_U[val];
		  temp = xx[column_start];
		  if ( column_end-column_start > sn_record[column_start] ) column_number_sn = sn_record[column_start]+1;
		  else column_number_sn = column_end-column_start+1;
		  column_sn_end = row_ptr_U[val+column_number_sn-1];
		  row_num = offset_L[column_sn_end+1] - offset_L[column_sn_end] - 1;
		  row_count = offset_L[column_start + 1] - offset_L[column_start] - 1;
		  column_divid = row_num % 16;
		  divid = row_num / 16;
		  if ( column_number_sn < 4 || row_num < thresold )
		  {
			  for ( i = offset_L[column_start]+1; i < offset_L[column_start+1]; i++ )
			  {
				  xx[row_ptr_L[i]] -=  temp*L[i];
			  }
			  pack = 1;
		  }
		  else
		  {
			  row_sn_start = offset_L[column_sn_end] + column_divid + 1;
			 /* for ( qqq = 0; qqq < column_number_sn; qqq++ )
			  {
				  row_else1 = row_ptr_U[val+qqq];
				  temp = xx[row_else1];
				  for ( sss = offset_L[row_else1]+1; sss < offset_L[row_else1+1]-divid*16; sss++ )
				  {
					  xx[row_ptr_L[sss]] -= temp * L[sss];
				  }
			  }*/
			  for ( qqq = 0; qqq < column_number_sn-1; qqq++ )
			  {
				  row_else1 = row_ptr_U[val+qqq];
				  row_else2 = row_ptr_U[val+qqq+1];
				  temp = xx[row_else1];
				  dense_vec_counter = qqq;
				  for ( sss = offset_L[row_else1]+1; sss < offset_L[row_else1+1]-divid*16; sss++ )
				  {
					  dense_vec[dense_vec_counter++] += temp*L[sss];
				  }
				  xx[row_else2] -= dense_vec[row_else2-column_start-1];
				  dense_vec[row_else2-column_start-1] = 0;
			  }
			  temp = xx[column_sn_end];
			  dense_vec_counter = column_number_sn - 1;
			  for ( i = offset_L[column_sn_end]+1; i < offset_L[column_sn_end+1]-divid*16; i++ )
			  {
				  dense_vec[dense_vec_counter++] += temp*L[i];
			  }
			  dense_vec_counter = column_number_sn - 1;
			  for ( sss = offset_L[column_start+1]-row_num; sss < offset_L[column_start+1]-divid*16; sss++ )
			  {
				  xx[row_ptr_L[sss]] -= dense_vec[dense_vec_counter];
				  dense_vec[dense_vec_counter] = 0;
				  dense_vec_counter++;
			  }
			  dividd = divid;
			  for ( qqq = 0; qqq < divid; qqq++ )
			  {
				  add_sum.vec = zero.vec; 
				  add_sum_2.vec = zero.vec;
				  
				  for ( sss = 0; sss < column_number_sn; sss++ )
				  {
					  row_else1 = row_ptr_U[val+sss];
					  row_record = offset_L[row_else1+1] - dividd*16;
					  row_record_2 = offset_L[row_else1+1] - dividd*16 + 8;
					  
					  v_row.vec = _mm512_set1_pd(xx[row_else1]);
					  v_l.vec = _mm512_load_pd(&L[row_record]);
					  v_l_2.vec = _mm512_load_pd(&L[row_record_2]);
					  add_sum.vec = _mm512_fmadd_pd(v_row.vec, v_l.vec, add_sum.vec);	
					  add_sum_2.vec = _mm512_fmadd_pd(v_row.vec, v_l_2.vec, add_sum_2.vec);	
				  }

				  vi.vec = _mm256_load_si256((__m256i const *)&row_ptr_L[row_sn_start+qqq*16]);
				  v.vec = _mm512_i32gather_pd(vi.vec, &xx[0], 8);
				  v_sub.vec = _mm512_sub_pd(v.vec, add_sum.vec);
				  _mm512_i32scatter_pd(&xx[0], vi.vec, v_sub.vec, 8);
				  
				  vi_2.vec = _mm256_load_si256((__m256i const *)&row_ptr_L[row_sn_start+qqq*16+8]);
				  v_2.vec = _mm512_i32gather_pd(vi_2.vec, &xx[0], 8);
				  v_sub_2.vec = _mm512_sub_pd(v_2.vec, add_sum_2.vec);
				  _mm512_i32scatter_pd(&xx[0], vi_2.vec, v_sub_2.vec, 8); 
				  
				  dividd--;
			  }
			  pack = column_number_sn;
		  }    
	  }
	  /* solve for U[:,k] and L[:, k] */
	  for ( i = offset_U[k]; i < offset_U[k+1]; i++ )
	  {
		  U[i] = xx[row_ptr_U[i]];
		  xx[row_ptr_U[i]] = 0;
	  }
	  
	  U_diag = U[i-1];
	  for ( i = offset_L[k]+1; i < offset_L[k+1]; i++ )
	  {
		  L[i] = xx[row_ptr_L[i]] / U_diag;
		  xx[row_ptr_L[i]] = 0;
	  }
  }

  /* solve for Ly = b and Ux = y */
  double *y, *x;
  y = ( double *)_mm_malloc( sizeof( double ) * n, 64 );
  x = ( double *)_mm_malloc( sizeof( double ) * n, 64 );

  for ( i = 0; i < n; i++ )
  {
	  y[i] = 1.0;
  }

  for ( i = 0; i < n; i++ )
  {
	  for ( j = offset_L[i]+1; j < offset_L[i+1]; j++ )
	  {
		  y[row_ptr_L[j]] -= y[i] * L[j];
	  }
  }

  //x[n-1] = y[n-1];
  for ( i = 0; i < n; i++ )
  {
	  x[i] = y[i];
  }

  x[n-1] = y[n-1]/U[nzu-1];
  for ( i = n-1; i > 0; i-- )
  {
	  for ( j = offset_U[i]; j < offset_U[i+1]-1; j++ )
	  {
		  x[row_ptr_U[j]] -= x[i] *U[j];
	  }
	  x[i-1] = x[i-1]/U[offset_U[i]-1];
  }
  return x;
}

