/******************************************************************************

	module:		sensingRHSEm


	function that incorporates the adjoint boundary condition into the rhs vector

******************************************************************************/

void sensingRHSEm(int lw, int lt, int ls, int ld, Mesh *mesh, Meas *meas,
Param *param, User *user, vector<double> &arhs)
{
	
	/*
	right-hand sides of the adjoint equation for sensitivity coefficient
	*/
	int size = param->fwd_size;
  int nmax = mesh->nmax;
  int nelem = mesh->nelem;

  double totArea;
  double dot;
  double dFdP;
  
  int detectNode;
  int vectorIndex;
  int lld;
  
	// initialize
 	fill(arhs.begin(),arhs.end(),0.0);
 	
	//lld = meas->detpairs[ls][ld];
	
	dFdP = 1.0;
	
	detectNode = meas->nodeldetector[ld];	
	
	/* diffusion model (begin) */
	if(user->lightModel == 1){
			vectorIndex = detectNode;
			arhs[vectorIndex] =  (1.0/2.0)*(1.0-2.0*user->R[0])/(1.0+3.0*user->R[1]);
	}//

	/* RTE model(begin) */
	else if(user->lightModel == 3){
		totArea = 0.0;
		for(int j=0;j<mesh->elem[detectNode].snum;j++){
			if(mesh->elem[detectNode].sbound[j] == 1){
				totArea = totArea + mesh->elem[detectNode].sarea[j];
			}
		}
		for(int k=0;k<nmax;k++){
			vectorIndex = nelem*k+detectNode;
			for(int j=0;j<mesh->elem[detectNode].snum;j++){
				if(mesh->elem[detectNode].sbound[j] == 1){
					dot = vecDot(mesh->s[k],mesh->elem[detectNode].snorm[j]);
					if(dot > 0.0){
						// specular reflection						
						arhs[vectorIndex] = arhs[vectorIndex] +  dot 
										* (1.0-dir_reflect(mesh->elem[detectNode].snorm[j],mesh->s[k],user->nindex))
										* mesh->w[k] * (mesh->elem[detectNode].sarea[j] / totArea);
									
						
						// diffuse reflection
						/*
						arhs[vectorIndex] = arhs[vectorIndex] 
															-  dFdP * (1.0-user->R_eff)
															* dot * mesh->w[k] * (mesh->elem[detectNode].sarea[j] / totArea);
									
					 */				
					}				
				}
			}
		}// RTE model (end)
	//cerr << ls << "\t" << lld << "\n";
	}					
	
	if(vecNorm(arhs) < 1.0e-25){
		cerr << "\tERROR in the vector b of Ax=b for the sensitivity problem\n";
		cerr << "\tMesh node " << meas->detpairs[ls][ld] << " is not detected:\n";
		cerr << "\tnode should be included in the measurement operator\n";
		exit(-1);
	}
	
} //

/******************************************************************************

	module:		sensingRHSEx.cpp

	release:	version 1.3
	author:		Hyun Keol Kim, Columbia University, NY 10027, USA
	e-mail:		hkk2107@columbia.edu
	date:			May 10, 2005
	remarks:

	function that incorporates the adjoint boundary condition into the rhs vector

******************************************************************************/

void sensingRHSEx(int lw, int lt, int ls, int ld, Mesh *mesh, Meas *meas,
Param *param, User *user, vector<double> &arhs)
{
	
	/*
	right-hand sides of the adjoint equation for sensitivity coefficient
	*/
	int nelem = mesh->nelem;
	int	nmax = mesh->nmax;
	int vectorIndex;
	int size = param->fwd_size;
	double sum;
	double onepi = 2.0*asin(1.0);

 	fill(arhs.begin(),arhs.end(),0.0);

	/* diffusion model (begin) */
	if(user->lightModel == 1){
		for(int ele=0;ele<nelem;ele++){
				arhs[ele] = 
										- param->axim[lw][lt][ls][ele] 
										* (param->qs[ele]*user->qYield/(1.0 - meas->Ls[lt] * user->lifeTime)) 
										* mesh->elem[ele].vol;
		}
	}//

	/* sp3 model (begin) */
	else if(user->lightModel == 2){
		for(int ele=0;ele<nelem;ele++){
			//arhs[ele] = arhs[ele + nelem] = 0.0;// not complete
		}//
	}//

	/* RTE model(begin) */
	else if(user->lightModel == 3){
		for(int ele=0;ele<nelem;ele++){
			for(int ka=0;ka<nmax;ka++){
				vectorIndex = nelem*ka+ele;
				sum = 0.0;
				for(int kb=0;kb<nmax;kb++){
					sum += param->axim[lw][lt][ls][nelem*kb+ele];
				}
				arhs[vectorIndex] = 
									- sum * mesh->w[ka] 
									* (param->qs[ele]*user->qYield/(1.0 - meas->Ls[lt] * user->lifeTime)) 
									* mesh->elem[ele].vol/4.0/onepi;
			}
		}
	}// RTE model(end)				
	
	if(vecNorm(arhs) < 1.0e-25){
		cerr << "\tERROR in the vector b of Ax=b for the sensitivity problem\n";
		cerr << "\tMesh node " << meas->detpairs[ls][ld] << " is not detected:\n";
		cerr << "\tnode should be included in the measurement operator\n";
		exit(-1);
	}
	
} //

