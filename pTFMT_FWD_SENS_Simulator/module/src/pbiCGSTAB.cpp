/******************************************************************************

	module:		pbicgstab.cpp

	function that applies a D-ILU preconditioned BiConjugateStab (BiCGSTAB) to 
	solve the linear system A*u = rhs iteratively.

******************************************************************************/

/*	
	a(k): = nonezero values of sparse matrix using a compress storage row (CSR)
	csrRowPtr(k) = row pointer
	csrColInd(k) = column index for k-th entry in A(k)

	n:= size of solution vector
	u:= solution vector of size n
	rhs:= right-hand size
	nnz:= number of nonezero entries in sparse matrix A
	iterMax: = number of maximum iterations

*/
void pbicgstab(int prob, int preconditioner, double *a, int *csrRowPtr,
           int *csrColInd, vector<double> &x, vector<double> &b,
           int nnz, int iterMax, double eps, int progressMonitor)

{
  unsigned int n = x.size();
 	unsigned int iter;
 	double rho,rho_old;
 	double omega;
 	double alpha;
 	double beta;
  double bnorm2;
  double resid;
	
  vector <double> cd(n,0.0);
  vector <double> q(n,0.0);
  vector <double> r(n,0.0);
  vector <double> p(n,0.0),phat(n,0.0);  
  vector <double> w(n,0.0); // matrix-vector product
  vector <double> rt(n,0.0);
  vector <double> s(n,0.0),shat(n,0.0);
  vector <double> t(n,0.0);
	
	//initialize
	fill(x.begin(),x.end(),0.0);
	bnorm2 = vecNorm(b);
	if(preconditioner == 1) constructDILU (n, nnz, a, csrRowPtr, csrColInd, cd);

	//r0=(b-A*x0)
  ax_cr(prob, a, csrRowPtr, csrColInd, x, w, n, nnz);
	r = b - w;
	p = r;
	rt = r;
	if(iter == 0 && progressMonitor == 1){
		cout << "\tinitial residual = " << vecNorm(r) << "\n";
	}

	/* main loop */
  for(iter = 0; iter < iterMax; iter++ ){
		// rho = rt*r
		rho = rt * r;    
		if(iter > 0){
			beta = (rho/rho_old)*(alpha/omega);
			p = r + beta*(p - omega*q);		
		}

		// M*phat = q
    if(preconditioner == 1) solveDILU(prob, n, nnz, a, csrRowPtr, csrColInd, cd, p, phat);
    ax_cr(prob, a, csrRowPtr, csrColInd, phat, q, n, nnz);

		// alpha = rho/(rt*q);
		alpha = rho / (rt*q);

		// s = r - alpha * q
    s = r - alpha *q;		
    resid = vecNorm(s);
		if(resid <= eps) break;
		
		// x = x + alpha*phat
    x = x + alpha *phat;
		
		// M*shat = s
    if(preconditioner == 1) solveDILU(prob, n, nnz, a, csrRowPtr, csrColInd, cd, s, shat);
    ax_cr(prob, a, csrRowPtr, csrColInd, shat, t, n, nnz);
		omega = (t*s)/(t*t);
		
		x = x + omega*shat;
		r = s - omega*t;		  			
		rho_old = rho;
  }
	
	if(progressMonitor == 1) {
		cout << "\tnumber of iterations = " << ( iter + 1 ) << "\n";
		cout << "\tfinal residual = " << resid << "\n\n";
	}

  vector<double> ().swap(cd);
  vector<double> ().swap(shat);
  vector<double> ().swap(r);
  vector<double> ().swap(phat);
  vector<double> ().swap(p);
  vector<double> ().swap(q);  
  vector<double> ().swap(w); // matrix-vector product
  vector<double> ().swap(rt);
  vector<double> ().swap(s);
  vector<double> ().swap(t);
  
  return;
}

