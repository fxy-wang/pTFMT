/******************************************************************************

	module:		pmgmres.cpp

	function that applies a D-ILU preconditioned modified restarted GMRES method to solve the linear
	system Ax = b iteratively.

******************************************************************************/

/*	function that applies the preconditioned modified version of a restarted GMRES (pmgmres) solver to
	a linear system Ax=b.

	a(k) = nonezero values of sparse matrix using a row-compress scheme.
	rowIndex(k) = row index for k-th entry in A(k)
	columnIndex(k) = column index for k-th entry in A(k)

	u(n) = complex-valued intensity in unit W/cm^2/sr.
	unknowns = number of unknowns for which Ax=b is solved
	nonezeroNum = number of nonezero entries in sparse matrix A
	iterMax = number of maximum outer iterations
	mr = number of maximum inner iterations(should be smaller than unknowns)
	eps  = tolerance
*/


using namespace std;

void pmgmres(int prob, int preconditioner, double *a, int *rowIndex,
           int *columnIndex, vector <double> &u, vector<double> &rhs, int n,
           int nonezeroNum, int iterMax, int mr, double eps, int forwardMonitor)

{
  int i;
  int itr;
  int j;
  int k;
	
  double av;
  double delta = 1.0e-3;
  double htmp;
	double mu;
  double rho;
  double rho_tol;
	
	vector<double> c(mr+1,0.0);
  vector<double> g(mr+1,0.0);
  vector<vector<double>> h(mr+1,vector<double> (mr,0.0));
  vector<double> r(n,0.0);
  vector<double> x(n,0.0);
  vector<double> s(mr+1,0.0);
  vector<vector<double>> v(mr+1,vector<double> (n,0.0));
  vector<double> y(mr+1,0.0);
  vector<double> cd(n,0.0);


	if(preconditioner == 1) constructDILU (n, nonezeroNum, a, rowIndex, columnIndex, cd);

	fill(x.begin(),x.end(),0.0);
  rho = vecNorm(r);
	rho_tol = eps * rho;
  if(forwardMonitor == 1) cout << "\tinitial residual = " << rho << "\n";

  //begin::outer iteration
  for ( itr = 0; itr < iterMax; itr++ )
  {
 		
    ax_cr(prob, a, rowIndex, columnIndex, x, r, n, nonezeroNum);

    r = rhs - r;
		
    if(preconditioner == 1) solveDILU(prob, n, nonezeroNum, a, rowIndex, columnIndex, cd, r, r);

    rho = residual(r, n);
		
    for ( i = 0; i < n; i++)
    {
      v[0][i] = r[i] / rho;
    }

    for ( i = 0; i <= mr; i++ )
    {
      g[i] = 0.0;
      for ( j = 0; j < mr; j++ )
      {
        h[i][j] = 0.0;
      }
    }

    g[0] = rho;
    k = 0;
	 //begin::inner iteration
    while ( ( rho_tol < rho ) && ( k < mr ) )
    {
      ax_cr (prob, a, rowIndex, columnIndex, v[k], v[k+1], n, nonezeroNum );
      if(preconditioner == 1) solveDILU(prob, n, nonezeroNum, a, rowIndex, columnIndex, cd, v[k+1], v[k+1]);
      av = residual ( v[k+1], n );

      for ( j = 0; j <= k; j++ )
      {
        h[j][k] = dot_product ( v[k+1], v[j], n );
        for ( i = 0; i < n; i++ )
        {
          v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
        }
      }

      h[k+1][k] = residual ( v[k+1], n );

      if ( ( av + delta * h[k+1][k] ) == av )
      {
         for ( j = 0; j <= k; j++ )
         {
           htmp = dot_product ( v[k+1], v[j], n );
           h[j][k] = h[j][k] + htmp;
           for ( i = 0; i < n; i++ )
           {
             v[k+1][i] = v[k+1][i] - htmp * v[j][i];
           }
         }
         h[k+1][k] = residual ( v[k+1], n );
      }

      if ( h[k+1][k] != 0.0 ){
        for ( i = 0; i < n; i++ ){
          v[k+1][i] = v[k+1][i] / h[k+1][k];
        }
      }

      if ( 0 < k )
      {
        for ( i = 0; i <= k+1; i++ )
        {
          y[i] = h[i][k];
        }
        for ( j = 0; j < k; j++ )
        {
          mult_givens ( c[j], s[j], j, y );
        }
        for ( i = 0; i <= k+1; i++ )
        {
          h[i][k] = y[i];
        }
      }

      mu = sqrt ( pow ( h[k][k], 2 ) + pow ( h[k+1][k], 2 ) );
      c[k] = h[k][k] / mu;
      s[k] = -h[k+1][k] / mu;
      h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
      h[k+1][k] = 0;
      mult_givens ( c[k], s[k], k, g );
      rho = fabs ( g[k+1] );
      k = k + 1;

    }

    //end::inner iteration

    k = k-1;
    y[k] = g[k] / h[k][k];

    for ( i = k-1; 0 <= i; i-- )
    {
      y[i] = g[i];
      for ( j = k; i < j; j-- )
      {
        y[i] = y[i] - h[i][j] * y[j];
      }
      y[i] = y[i] / h[i][i];
    }

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < k; j++ )
      {
        x[i] = x[i] + v[j][i] * y[j];
      }
    }

    if (rho <= rho_tol)
    {
      break;
    }
    //*
	 	if(forwardMonitor == 1)
	 	{
	 		//cout << "\touter iteration #" << itr <<"\n";
	 		//cout << "\tresidual: " << setprecision(6) << scientific << rho << endl;
	 	}
		//*/
    /// outer iteration (end)
  }

	if(forwardMonitor == 1) {
		cout << "\tnumber of iterations = " << ( itr + 1 ) * mr << "\n";
		cout << "\tfinal residual = " << rho << "\n\n";
	}

	for(int i=0; i<n; i++) {
		u[i] = x[i];		
	}	

  vector<double> ().swap(c);
  vector<double> ().swap(g);
  vector<vector<double>> ().swap(h);
  vector<double> ().swap(r);
  vector<double> ().swap(x);
  vector<double> ().swap(s);  
  vector<vector<double>> ().swap(v);
  vector<double> ().swap(y);
  vector<double> ().swap(cd);

  return;
}

