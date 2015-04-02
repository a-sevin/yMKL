plug_in,"libmkl";

extern potri_cpu;
/* DOCUMENT potri_cpu
   potri_cpu, h_mat
     
   SEE ALSO:
 */
extern gemm_cpu;
/* DOCUMENT gemm_cpu
   gemm_cpu, C, A, B 
   C = gemm_cpu(A, B)
      
   SEE ALSO:
 */
extern symm_cpu;
/* DOCUMENT symm_cpu
   symm_cpu, C, A, B [, side, alpha, beta]
   C = symm_cpu(A, B [, side, alpha])
   
   preform C = alpha*A*B + beta*C when it's called by procedure
   preform C = alpha*A*B when it's called by function

   where A is symmetric (or B if side is set to 'R')

   defauts :
     - side = 'L'
     - alpha = 1.
     - beta = 0.

   SEE ALSO:
 */
extern syevd_cpu;
/* DOCUMENT syevd_cpu
   evd_cpu, mat, ev 
      
   SEE ALSO:
 */
func check_mkl(n){
  tmp = random(n, n)-0.5;
  write, "doing mat = tmp(,+)*tmp(,+)...";
  mat = tmp(,+)*tmp(,+);
  write, "doing mat2= gemm_cpu(tmp, tmp, 'n', 't')...";
  mat2= gemm_cpu(tmp, tmp, 'n', 't');
  write, "max(abs(mat-mat2))";
  max(abs(mat-mat2));
  write, "doing sytmp = gemm_cpu(mat, mat, 'n', 'n')...";
  sytmp = gemm_cpu(mat, mat);
  write, "doing mat3= gemm_cpu(sytmp, mat)...";
  mat3= gemm_cpu(sytmp, mat);
  write, "doing mat4= symm_cpu(sytmp, mat)...";
  mat4= symm_cpu(sytmp, mat);
  write, "max(abs(mat3-mat4))";
  max(abs(mat3-mat4));
  write, "doing mat5= symm_cpu(mat, sytmp, 'R')...";
  mat5= symm_cpu(sytmp, mat, 'R');
  write, "max(abs(mat3-mat5))";
  max(abs(mat3-mat5));
  write, "potri_cpu, mat";
  potri_cpu, mat2;
  write, "max(abs( mat(,+)*inv_mat(+,) - unit(n) ));";
  max(abs( gemm_cpu(mat, mat2) - unit(n) ));
  write, "ev=SVdec(mat)";
  ev=SVdec(mat);
  write, "ev2=syevd_cpu(mat)";
  ev2=syevd_cpu( mat );
  write, "max(abs(ev-ev2(::-1))";
  max(abs(ev-ev2(::-1)));
  error;
}
func bench_potri(n, niter){
  tmp = random(n, n)-0.5;
  mat = gemm_cpu(tmp, tmp, 'n', 't');
  write, "potri_cpu, mat";
  tmp=0.;
  for(i=0; i<niter; i++){
    mat2= mat;
    tic; potri_cpu, mat2;tmp+=tac()/niter;
  }
  write, format="temps moyen : %fs\n",tmp;
  write, "max(abs( mat(,+)*inv_mat(+,) - unit(n) ));";
  max(abs( gemm_cpu(mat, mat2) - unit(n) ));
  error;
}
// 1.156205s

func bench_mkl(n){
  tmp = random(n, n)-0.5;
  write, "doing mat= gemm_cpu(tmp, tmp, 'n', 't')...";
  mat= gemm_cpu(tmp, tmp, 'n', 't');
  sytmp = 0*tmp;
  write, "doing gemm_cpu, sytmp, mat, tmp...";
  tic; gemm_cpu, sytmp, mat, tmp;
  write, format= "in %0.3fs\n", tac();
  write, "doing symm_cpu, sytmp, mat, tmp...";
  tic; symm_cpu, sytmp, mat, tmp;
  write, format= "in %0.3fs\n", tac();
  write, "doing symm_cpu, sytmp, tmp, mat, 'R'...";
  tic; symm_cpu, sytmp, tmp, mat, 'R';
  write, format= "in %0.3fs\n", tac();
  error;
}

/*
    if (yType == Y_FLOAT) {
      LAPACKE_ssyevd( n, 'V', 'L', n, (float*)h_mat, n, (float*)h_ev);
    } else if (yType == Y_DOUBLE) {
      LAPACKE_dsyevd( n, 'V', 'L', n, (double*)h_mat, n, (double*)h_ev);
    } else {
      y_error("carma_potri not implemented for this type");
    }
*/
