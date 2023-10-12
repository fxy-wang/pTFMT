/*******************************************************************/
bool sortbyfirst(const pair<double,int> &a,  
               const pair<double,int> &b) 
/*******************************************************************/
/*
 Driver function to sort the vector elements by  
 first element of pair in descending order
*/
{ 
       return (a.first > b.first); 
} 

/*******************************************************************/
string splitFilename(const string str)
/*******************************************************************/
{
		string filename;
		string filenameOnly;
		
		size_t folder_pos = str.find_last_of("/\\");
		filename = (folder_pos == string::npos) ? str : str.substr(folder_pos+1);

		size_t dot_pos = filename.find_last_of(".");

		filenameOnly = (dot_pos == string::npos) ? filename : filename.substr(0,dot_pos);
		
		//cout << filename << "\t" << filenameOnly << endl;	
		
		return filenameOnly;
}

/*******************************************************************/
double dot_product(vector<double> v0, vector<double> v1, int n)
/*******************************************************************/
{
  int i;
  double s;

  s = 0.0;
  for (i=0; i<n; i++)
  {
    s = s + v0[i] * v1[i];
  }

  return s;
}

/********************************************************************/
double residual(vector<double> r, int n)
/********************************************************************/
/* function that performs a vector norm
*/
{
  int i;
  double s = 0.0;

  for (i=0; i<n; i++){
    s += r[i] * r[i];
  }

  return sqrt(s);
}

/*******************************************************************/
double vecDot(vector <double> v1, vector <double> v2)
/*******************************************************************/
{
  double s;
  s = 0.0;
  for (int i=0; i<v1.size(); i++){
    s = s + v1[i] * v2[i];
  }
  return s;
}

/********************************************************************/
double vecNorm(vector <double> v)
/********************************************************************/
{
  double s;
  s = 0.0;
  for (int i=0; i<v.size(); i++){
    s = s + v[i] * v[i];
  }
  return sqrt(s);
}

/*******************************************************************/
double vecDist(vector <double> v1, vector <double> v2)
/*******************************************************************/
{
  double d;
  vector <double> tmp;
  tmp = v1-v2;
  if(v2[0] > 1.0e10){
		d = sqrt(tmp[1]*tmp[1] + tmp[2]*tmp[2]);  	
  }else if(v2[1] > 1.0e10){
		d = sqrt(tmp[0]*tmp[0] + tmp[2]*tmp[2]);  	
	}else if(v2[2] > 1.0e10){
		d = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);  	
	}else{
  	d = sqrt(tmp*tmp);
	}
  return d;
}

// ******************************************************************/
// vector<double> ompMXV(int N, vector <vector<double>> M, vector <double> V)
// ******************************************************************/
// {
//   vector<double> MV;
// 	int i, j, m=M.size(),n = M[0].size();
// 	assert(n==V.size());
// 	MV.resize(m,0.0);	
// 	
// 	# pragma omp parallel for schedule(static) shared(M,V,MV,m,n) private(i,j) num_threads(N) 
// 		for ( i = 0; i < m; i++ )
// 			for ( j = 0; j < n; j++ )
// 					MV[i] += M[i][j] * V[j];
// 					
//   return MV;
// }
// 
// ******************************************************************/
// double ompVXV(int N, vector <double> v1, vector <double> v2)
// ******************************************************************/
// {
//   double s;
// 	int i, n = v1.size();
// 	assert(n==v2.size());
// 	s = 0.0;
// 	# pragma omp parallel for reduction( + : s ) shared(v1,v2,n) private(i) num_threads(N) 
// 		for ( i = 0; i < n; i++ ) 
// 			s += v1[i] * v2[i];
// 					
//   return s;
// }


/********************************************************************/
double gaussrand1()
/********************************************************************/
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

/**********************************************************************/
double gaussrand2()
/**********************************************************************/
{
	static double U, V;
	static int phase = 0;
	double Z;
	double onepi = 2.0*asin(1.0);

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * onepi * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * onepi * V);

	phase = 1 - phase;

	return Z;
}


/*******************************************************************/
void mult_givens(double c, double s, int k, vector<double> &g)
/*******************************************************************/
/*	function from a GMRES module that applies a Given rotation to the
	two successive entries of a vector.
*/
{
	double g1;
	double g2;

	g1     = c * g[k] - s * g[k+1];
	g2     = s * g[k] + c * g[k+1];

	g[k]	 = g1;
	g[k+1] = g2;

  return;
}

////////////////////////////////////////////////////////////////
double sfn(vector <double> ska, vector <double> skb, double g)
////////////////////////////////////////////////////////////////
{
//nonlinear Henyey-Greenstein scattering phase function
	double	cosine 	= 0.0;
	double	phaseFunc = 0.0;

	cosine	 =	vecDot(ska, skb);
	phaseFunc =	(1.0 - g * g) / pow((1.0 + g * g - 2.0 * g * cosine),1.5);

	return phaseFunc;
}

/*******************************************************************/
double reflect(double dot, double nindex)
/*******************************************************************/
/*	function that computes reflectivity at boundary due to air-tissue
	refractive index mismatch problem.
*/
{
	double small=1.0e-25;
	double onepi = 2.0*asin(1.0);
	double angIn,angOut,angCut;
	double rho,rho1,rho2,c1,c2;
	
	angIn    =   acos(dot);
	angOut   =   asin(nindex*sin(angIn));
	angCut   =   asin(1.0 / nindex);

	if(angIn >= angCut && angIn < (onepi / 2.0)){
		rho  =  1.0;
	}
	else{
		rho1  =   (nindex * cos(angOut) - cos(angIn))
							/(nindex * cos(angOut) + cos(angIn)+small);
		rho2  =   (nindex * cos(angIn) - cos(angOut))
							/(nindex * cos(angIn) + cos(angOut)+small);
		rho   =   (1.0/2.0) * pow(rho1,2) + (1.0/2.0)*pow(rho2,2);
	}
	return rho;
}

/*******************************************************************/
void reflect_momentum(User *user)
/*******************************************************************/
/*	function that computes reflection moments for SP3 */
{
	int n=100;
	double onepi = 2.0*asin(1.0);
	double angle;
	double dot;
	double totalAngle=onepi/2.0;
	double dtheta=totalAngle/double(n);
	double rho;

	user->R[0]=user->R[1]=user->R[2]=user->R[3]=user->R[4]=user->R[5]=0.0;
 	for(int k=0;k<n;k++){
 	  angle=double(k)*dtheta;
 	  dot=cos(angle);
 	  rho = reflect(dot,user->nindex);
 	  //rho = 1.0;
 	  //cerr << angle*180/onepi << " " << dot << " " << rho << " " << endl;
		if(k == 0 || k == n-1){
		  user->R[0] += cos(angle)*sin(angle)*rho*(dtheta/2.0);
		  user->R[1] += pow(cos(angle),2.0)*sin(angle)*rho*(dtheta/2.0);
		  user->R[2] += pow(cos(angle),3.0)*sin(angle)*rho*(dtheta/2.0);
		  user->R[3] += pow(cos(angle),4.0)*sin(angle)*rho*(dtheta/2.0);
		  user->R[4] += pow(cos(angle),5.0)*sin(angle)*rho*(dtheta/2.0);
		  user->R[5] += pow(cos(angle),6.0)*sin(angle)*rho*(dtheta/2.0);
		  
		}else{
		  user->R[0] += cos(angle)*sin(angle)*rho*dtheta;
		  user->R[1] += pow(cos(angle),2.0)*sin(angle)*rho*dtheta;
		  user->R[2] += pow(cos(angle),3.0)*sin(angle)*rho*dtheta;
		  user->R[3] += pow(cos(angle),4.0)*sin(angle)*rho*dtheta;
		  user->R[4] += pow(cos(angle),5.0)*sin(angle)*rho*dtheta;
		  user->R[5] += pow(cos(angle),7.0)*sin(angle)*rho*dtheta;
	  }
	}
	//cerr << R[0] << " " << R[1] << endl;
	user->R_eff = (2.0*user->R[0] + 3.0*user->R[1])
	  					/ (2.0 - 2.0*user->R[0] + 3.0*user->R[1]);
	
	//cerr << user->R[0] << " " << user->R[1] << endl;

	return;
}

/*******************************************************************/
double dir_reflect(vector <double> sfnorm, vector <double> skb,
 double nindex)
/*******************************************************************/ 
/*	function that computes reflectivity at boundary due to air-tissue
	refractive index mismatch problem.
*/
{
	double onepi = 2.0*asin(1.0);
	double dot,angIn,angOut,angCut;
	double rho,rho1,rho2,c1,c2;
		
	dot      =  vecDot(sfnorm, skb) / vecNorm(skb);
	angIn    =   acos(dot);
	angCut   =   asin(1.0 / nindex);

	if(angIn >= angCut && angIn < (onepi / 2.0)){
		rho  =  1.0;
	}	
	else{
		angOut   =   asin(nindex*sin(angIn));		
		c1       =   cos(angIn);
		c2       =   cos(angOut);
		rho1     =   (nindex * c2 - c1)/(nindex * c2 + c1);
		rho2     =   (nindex * c1 - c2)/(nindex * c1 + c2);
		rho      =   0.5 * (rho1 * rho1 + rho2 * rho2);
	}
	return rho;
}	

/*******************************************************************/
double get_coeff(int order, int od, int ele, Param *param, User *user)
/*******************************************************************/
/*	function that computes coefficients appearing on the boundary condition for SP3 */
{
	double A1,B1,C1,D1;
  double A2,B2,C2,D2;
  double a1,b1,c1,d1;
  double a2,b2,c2,d2;  
  double coefficient;
  
  A1= -user->R[0];
  B1 = 3.*user->R[1];
  C1 = (-3./2.)*user->R[1]-(5./2.)*user->R[3];
  D1 = (3./2.)*user->R[1]-(5./2.)*user->R[3];
  
  A2= (-9./4.)*user->R[0]+(15./2.)*user->R[2]-(25./4.)*user->R[4];
  B2 = (63./4.)*user->R[1]-(105./2.)*user->R[3]+(175./4.)*user->R[5];
  C2 = (-3./2.)*user->R[0]+(5./2.)*user->R[2];
  D2 = (3./2.)*user->R[1]-(5./2.)*user->R[3];      			  

  a1=(1./2.)+A1;
  b1=(1.+B1)/(3.0*(param->xka[ele]+(1.-user->g)*param->xsig[ele]));
  c1=(1./8.)+C1;
  d1=D1/(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele]);
  
  a2=(7./24.)+A2;
  b2=(1.+B2)/(7.0*(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele]));
  c2=(1./8.)+C2;
  d2=D2/(param->xka[ele]+(1.-user->g)*param->xsig[ele]); 
  
  if(order == 0 && od == 0) coefficient=(a1*b2-d1*c2)/(b1*b2-d1*d2)/(param->xka[ele]+(1.-user->g)*param->xsig[ele])/3.;
  if(order == 0 && od == 1) coefficient=(d1*a2-b2*c1)/(b1*b2-d1*d2)/(param->xka[ele]+(1.-user->g)*param->xsig[ele])/3.;
  if(order == 1 && od == 0) coefficient=(a1*d2-b1*c2)/(b1*b2-d1*d2)/(param->xka[ele]+(1.-pow(user->g,3.))*param->xsig[ele])/7.;
  if(order == 1 && od == 1) coefficient=(a2*b1-c1*d2)/(b1*b2-d1*d2)/(param->xka[ele]+(1.-pow(user->g,3.))*param->xsig[ele])/7.;
             
	return coefficient;
}		

/*******************************************************************/
double get_source_coeff(int eqn, int ele, double s1, double s2, 
Param *param, User *user)
/*******************************************************************/
/*	function that computes coefficients for source term for SP3 */
{
	double A1,B1,C1,D1;
  double A2,B2,C2,D2;
  double a1,b1,c1,d1;
  double a2,b2,c2,d2;  
  double coefficient;
  
  A1= -user->R[0];
  B1 = 3.*user->R[1];
  C1 = (-3./2.)*user->R[1]-(5./2.)*user->R[3];
  D1 = (3./2.)*user->R[1]-(5./2.)*user->R[3];
  
  A2= (-9./4.)*user->R[0]+(15./2.)*user->R[2]-(25./4.)*user->R[4];
  B2 = (63./4.)*user->R[1]-(105./2.)*user->R[3]+(175./4.)*user->R[5];
  C2 = (-3./2.)*user->R[0]+(5./2.)*user->R[2];
  D2 = (3./2.)*user->R[1]-(5./2.)*user->R[3];      			  

  a1=(1./2.)+A1;
  b1=(1.+B1)/(3.0*(param->xka[ele]+(1.-user->g)*param->xsig[ele]));
  c1=(1./8.)+C1;
  d1=D1/(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele]);
  
  a2=(7./24.)+A2;
  b2=(1.+B2)/(7.0*(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele]));
  c2=(1./8.)+C2;
  d2=D2/(param->xka[ele]+(1.-user->g)*param->xsig[ele]); 
  
  if(eqn == 0 ) coefficient=(d1*s2+s1*b2)/(b1*b2-d1*d2)/(3.0*(param->xka[ele]+(1.-user->g)*param->xsig[ele]));
  if(eqn == 1 ) coefficient=(b1*s2+d2*s1)/(b1*b2-d1*d2)/(7.0*(param->xka[ele]+(1.-pow(user->g,3.))*param->xsig[ele]));
              
	return coefficient;
}	

/*******************************************************************/
void get_sp3_Q(int lw, int ele, Meas *meas, Param *param, User *user)
/*******************************************************************/
/*	function that computes coefficients for source term for SP3 */
{
	double A1,B1,C1,D1;
  double A2,B2,C2,D2;
	double J0,J1,J2,J3;

  double a1,b1,g1,d1; // alpha1,beta1,gamma1,delta1 based on SP3 paper(2016)
  double a2,b2,g2,d2; // alpha2,beta2,gamma2,delta2 based on SP3 paper(2016)
	double eta11,eta12,eta21,eta22;
  
  A1= -user->R[0];
  B1 = 3.*user->R[1];
  C1 = (-3./2.)*user->R[1]-(5./2.)*user->R[3];
  D1 = (3./2.)*user->R[1]-(5./2.)*user->R[3];
  
  A2= (-9./4.)*user->R[0]+(15./2.)*user->R[2]-(25./4.)*user->R[4];
  B2 = (63./4.)*user->R[1]-(105./2.)*user->R[3]+(175./4.)*user->R[5];
  C2 = (-3./2.)*user->R[0]+(5./2.)*user->R[2];
  D2 = (3./2.)*user->R[1]-(5./2.)*user->R[3];      			  

	J0 = -user->R[0]/2.;
	J1 = -3.*user->R[1]/2.;
	J2 = 5.*user->R[0]/4.-15.*user->R[2]/4.;
	J3 = 21.*user->R[1]/4.-35.*user->R[3]/4.;

  a1=(1./2.)+A1;
  a2=(7./24.)+A2;
  b1=(1.+B1);
  b2=(1.+B2);
  g1=(1./8.)+C1;
  g2=(1./8.)+C2;  
  
  d1=1./(3.0*(param->xka[ele]+(1.-user->g)*param->xsig[ele]));
  d2=1./(7.0*(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele]));
  
  eta11 = (a1*b2-7.*d1*g2)/(b1*b2-21.*d1*d2);
  eta21 = (3.*d2*eta11-g2)/b2;
  eta12 = (7.*d1*a2-g1*b2)/(b1*b2-21.*d1*d2);
  eta22 = (3.*d2*eta12+a2)/b2;
  
  user->nu1 = 1./4.+J0+(1./2.+J1)*eta11+J3*eta21;
  user->nu2 = -3./48.-2.*J0/3.+J2/3.+(1./2.+J1)*eta12+J3*eta22;

}	

//****************************************************************************80
void csrRearrange ( int n, int nz_num, int *ia, int *ja, double *a)
//****************************************************************************80
{
  double temp;
  int i;
  int is;
  int itemp;
  int j;
  int j1;
  int j2;
  int k;

  for ( i = 0; i < n; i++ ){
    j1 = ia[i];
    j2 = ia[i+1];
    is = j2 - j1;

    for ( k = 1; k < is; k++ ){
      for ( j = j1; j < j2 - k; j++ ){
        if ( ja[j+1] < ja[j] ){
          itemp = ja[j+1];
          ja[j+1] =  ja[j];
          ja[j] =  itemp;

          temp = a[j+1];
          a[j+1] =  a[j];
          a[j] = temp;
        }
      }
    }
  }
  return;
}

//****************************************************************************80

void constructDILU (int nx, int nz_num, double *a, int *ia, int *ja, vector<double> &cd)

//****************************************************************************80

{
/*
The incomplete preconditoner (ILU) is constructed as follows:
the coefficient matrix A into (L+D+U) where L and U are strict lower and upper
triangle matrices and D is a diagonal matrix containing the pivots based on:

C. POMMERELL, Solution of Large Unsymmetric Systems of Linear Equations,
vol. 17 of Series in Micro-electronics, volume 17, Hartung-Gorre Verlag, Konstanz, 1992

*/

  int irow2;
  int jcol2;
  int	nRow = nx;
	
  for (int irow = 0; irow < nRow; irow++ )
  {
    for (int nz = ia[irow]; nz <ia[irow+1]; nz++ )
    {
      if(irow == ja[nz]) cd[irow] = a[nz];
    }
	}
	
	for (int irow = 0; irow < nRow; irow++ )
	{
		cd[irow] = 1.0 / (cd[irow]);
		for (int nz = ia[irow]; nz <ia[irow+1]; nz++ )
		{
			irow2 = ja[nz];
			jcol2 = irow;
			for (int nz2 = ia[irow2]; nz2 <ia[irow2+1]; nz2++ )
			{
				if(ja[nz] > irow && ja[nz2] == jcol2)
				{
					cd[ja[nz]] = cd[ja[nz]] - a[nz2] * cd[ja[nz2]] * a[nz];

				}
			}
		}
	}

  return;
}


/*******************************************************************/
void  ax_cr(int matType, double *a, int *rowIndex,
			int *columnIndex, vector<double> &x, vector<double> &ax,
			int nx, int nonezeroNum)
/*******************************************************************/
/*
	subroutine that does the matrix-vector multiplication with the sparse matrix A
	and the vector x, and it returns the vector ax
*/
{
	int nzBegin,nzEnd;
	int nRow = nx;

	for(int row=0; row< nRow; row++){
		ax[row] = 1.0e-25;
	}

	// A * x
	if(matType == 1)
	{
		for(int row=0; row< nRow; row++)
		{
	   	nzBegin = rowIndex[row];
	    nzEnd = rowIndex[row+1];
	    for(int nz=nzBegin;nz<nzEnd;nz++)
	    {
	    	ax[row] = ax[row]+a[nz]*x[columnIndex[nz]];
	  	}
	  }
	}
	// At(transpose of A) * vector x
	else if(matType == 2){
		for(int row=0; row< nRow; row++){
	   	nzBegin = rowIndex[row];
	    nzEnd = rowIndex[row+1];
	    for(int nz=nzBegin;nz<nzEnd;nz++){
	    	ax[columnIndex[nz]] = ax[columnIndex[nz]]+a[nz]*x[row];
	  	}
	  }
	}

	return;

}

//****************************************************************************80

void solveDILU(int matType, int nx, int nz_num, double *a, int *ia, int *ja, 
			vector<double> &cd, vector<double> &r, vector<double> &z)

//****************************************************************************80
/*
	solve the ILU system given by M * z = r where M = (D + L) * D^(-1) * (D + U) as follows:

				(D + L) * y = r, and

				D^(-1) * (D + U) z = y

	and returns the vector z as a solution.
	The routine has been improved in calculation speed (HK,12/13/2019)
*/
//
{
  int irow2;
  int jcol2;

  int nRow = nx;
  double sum;
  double *tmp = new double [nRow];

  double *y = new double [nRow];
  double *zc = new double [nRow];

  for (int irow = 0; irow < nRow; irow++ )
  {
    y[irow] = r[irow];
  }
	// system associated with A
	if(matType == 1){
		// solve (D + L ) * y = r => y = D^(-1)*(r - L*y)
		for (int irow = 0; irow < nRow; irow++){
			for (int nz = ia[irow]; nz < ia[irow+1]; nz++ ){
				if(ja[nz] < irow){
						y[irow] = y[irow] - a[nz]* y[ja[nz]]; // strict L-matrix of A
				}
			}
			y[irow] = cd[irow] * y[irow];
		}

		// solve D^(-1)*(D+U) * z = y => z = y - D^(-1)*U*z
		for (int irow = nRow-1; irow >= 0; irow-- ){
			for (int nz = ia[irow]; nz < ia[irow+1]; nz++ ){
				if(ja[nz] > irow) {
					y[irow] = y[irow] - cd[irow]*a[nz]*y[ja[nz]]; // strict U-matrix] of A
				}
			}
		}
	}
	// system associated with At (transpose of A)
	if(matType == 2){
		// solve (D+Ut) * y = r => y = D^(-1)*(y-Ut*y)
		for (int irow = nRow-1; irow >= 0;irow-- ){
			for (int nz = ia[irow]; nz < ia[irow+1]; nz++ ){
				if(ja[nz] > irow) {
					//y[ja[nz]] = y[ja[nz]] - cd[irow]*a[nz]*y[irow]; // strict U-matrix] of A
					y[ja[nz]] = y[ja[nz]] - a[nz]* y[irow];
				}
			}
			y[irow] = cd[irow] * y[irow];
		}

		// solve D^(-1)*(D + Lt ) * z = y => y = (y - D^(-1)*Lt*y)
		for (int irow = 0; irow < nRow; irow++){
		
			for (int nz = ia[irow]; nz < ia[irow+1]; nz++ ){
				if(ja[nz] < irow){
						//y[ja[nz]] = y[ja[nz]] - a[nz]* y[irow]; // strict L-matrix of A
						y[ja[nz]] = y[ja[nz]] - cd[irow]*a[nz]*y[irow];
				}
			}
			//y[irow] = cd[irow] * y[irow];
		}

	}

	// copy w into z with complex value format
  for (int irow = 0; irow < nRow; irow++ )
  {
    z[irow] = y[irow];
  }

  delete [] y;
  delete [] zc;
	delete [] tmp;

  return;
}

