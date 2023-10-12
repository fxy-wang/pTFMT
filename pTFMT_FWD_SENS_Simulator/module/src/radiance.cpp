/******************************************************************************

	module:		radianceEx.cpp

	function that predicts radiative flux at detector positions on the surface.

******************************************************************************/

void radianceEx(int lw, int lt, int ls,
							Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int	nelem = mesh->nelem;
	int nmax = mesh->nmax;
	int vectorIndex,detectNode,detectIndex;
	double dot,a1,b1,c1,totArea;
	
	/*
	calculate the outgoing radiative flux at detectors
	*/
	
	for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
		meas->prediction[lw][lt][ls][meas->detpairs[ls][ld]] = 0.0;
	}

	/* 
	SP1: diffusion model 
	*/
	if(user->lightModel == 1){
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){	
			detectNode = meas->nodeldetector[meas->detpairs[ls][ld]];
			meas->prediction[lw][lt][ls][meas->detpairs[ls][ld]] = 
										(1.0/2.0)*((1.0-2.0*user->R[0])/(1.0+3.0*user->R[1]))
										*param->xis[lw][lt][ls][detectNode];
		}				
	}

	/* 
	SP3 
	*/
	if(user->lightModel == 2){
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){	
			detectNode = meas->nodeldetector[meas->detpairs[ls][ld]];
			get_sp3_Q(lw,detectNode,meas,param,user);
			
			meas->prediction[lw][lt][ls][meas->detpairs[ls][ld]]	
					= (user->nu1*param->xis[lw][lt][ls][detectNode] 
					+ user->nu2*param->xis[lw][lt][ls][detectNode+nelem]);
		}
	}

	/* 
	RTE model 
	*/
	if(user->lightModel == 3){
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
			detectNode = meas->nodeldetector[meas->detpairs[ls][ld]];
			totArea = 0.0;
			for(int j=0;j<mesh->elem[detectNode].snum;j++){
				if(mesh->elem[detectNode].sbound[j] == 1){
					totArea = totArea + mesh->elem[detectNode].sarea[j];
				}
			}
			for(int k=0;k<nmax;k++){
				vectorIndex=nelem*k+detectNode;			
				for(int j=0;j<mesh->elem[detectNode].snum;j++){
				if(mesh->elem[detectNode].sbound[j] == 1){
					dot=vecDot(mesh->s[k],mesh->elem[detectNode].snorm[j]);   				
					if(dot > 0.0){
						/* ideal specular reflection */
						/*
						meas->prediction[lw][lt][ls][meas->detpairs[ls][ld]] += param->xis[lw][lt][ls][vectorIndex]
															* dot * (1.0-dir_reflect(mesh->elem[detectNode].snorm[j],mesh->s[k],user->nindex))
															* mesh->w[k] * mesh->elem[detectNode].sarea[j]
															/ totArea;


						/*/
						
						/* more realistic diffuse transmission */
						meas->prediction[lw][lt][ls][meas->detpairs[ls][ld]] += param->xis[lw][lt][ls][vectorIndex]
															* dot * (1.0 - user->R_eff)
															* mesh->w[k] * mesh->elem[detectNode].sarea[j]
															/ totArea;	
															
						//*/													
					}
				}
			}

			//cerr << "\tprediction = " << ld << " " << meas->prediction[lw][lf][ls][ld] << endl;
		}
		}		
	}	

}

/******************************************************************************

	module:		computeExcitationField.cpp

	release:		version 1.3
	author:		Hyun Keol Kim, Columbia University, NY 10027, USA
	e-mail:		hkk2107@columbia.edu
	date:			May 10, 2005
	remarks:

	function that predicts radiative flux at detector positions on the surface.

******************************************************************************/
void computeExcitationField(int lw, int lt, int ls,
							Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int nelem	= mesh->nelem;
	int nmax 	= mesh->nmax;

	/* 
	SP1: diffusion model 
	*/
	if(user->lightModel == 1){
		for(int ele=0;ele<nelem;ele++){
			param->excitation[lw][lt][ls][ele] = param->xis[lw][lt][ls][ele];
		}
	}

	/* 
	SP3 
	*/
	if(user->lightModel == 2){
		for(int ele=0;ele<nelem;ele++){
			param->excitation[lw][lt][ls][ele]	
						= param->xis[lw][lt][ls][ele] - (2/3)*param->xis[lw][lt][ls][ele+nelem];
		}
	}

	/* 
	RTE 
	*/
	if(user->lightModel == 3){	
		int vectorIndex;	
		double sum;
		for(int ele=0;ele<nelem;ele++){
			sum = 0.0;
			for(int k=0;k<nmax;k++){
				vectorIndex = nelem*k+ele;
				sum += param->xis[lw][lt][ls][vectorIndex] * mesh->w[k];
			}
			param->excitation[lw][lt][ls][ele] = sum;
		}
	}

}//

/******************************************************************************

	module:		radianceEm.cpp

	release:		version 1.3
	author:		Hyun Keol Kim, Columbia University, NY 10027, USA
	e-mail:		hkk2107@columbia.edu
	date:			May 10, 2005
	remarks:

	function that predicts radiative flux at detector positions on the surface.

******************************************************************************/

void radianceEm(int lw, int lt, int ls,
							Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int	nelem = mesh->nelem;
	int nmax = mesh->nmax;
	int vectorIndex,detectNode,detectIndex;
	double dot,a1,b1,c1,totArea;
	
	/*
	calculate the outgoing radiative flux at detectors
	*/
	
	for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
		meas->prediction_m[lw][lt][ls][meas->detpairs[ls][ld]] = 0.0;
	}

	/* 
	SP1: diffusion model 
	*/
	if(user->lightModel == 1){
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){	
			detectNode = meas->nodeldetector[meas->detpairs[ls][ld]];
			meas->prediction_m[lw][lt][ls][meas->detpairs[ls][ld]] = 
							(1.0/2.0)*((1.0-2.0*user->R[0])/(1.0+3.0*user->R[1]))
							 * param->xism[lw][lt][ls][detectNode];
		}				
	}

	/* 
	SP3 
	*/
	if(user->lightModel == 2){
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){	
			detectNode = meas->nodeldetector[meas->detpairs[ls][ld]];
			get_sp3_Q(lw,detectNode,meas,param,user);
			
			meas->prediction_m[lw][lt][ls][meas->detpairs[ls][ld]] = 
							(user->nu1*param->xism[lw][lt][ls][detectNode] 
								+ user->nu2*param->xism[lw][lt][ls][detectNode+nelem]);
		}
	}

	/* 
	RTE model 
	*/
	if(user->lightModel == 3){
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
			detectNode = meas->nodeldetector[meas->detpairs[ls][ld]];
			totArea = 0.0;
			for(int j=0;j<mesh->elem[detectNode].snum;j++){
				if(mesh->elem[detectNode].sbound[j] == 1){
					totArea += mesh->elem[detectNode].sarea[j];
				}
			}
			for(int k=0;k<nmax;k++){
				vectorIndex=nelem*k+detectNode;			
				for(int j=0;j<mesh->elem[detectNode].snum;j++){
				if(mesh->elem[detectNode].sbound[j] == 1){
					dot=vecDot(mesh->s[k],mesh->elem[detectNode].snorm[j]);   				
					if(dot > 0.0){
						/* ideal specular reflection */
						/*
						meas->prediction[lw][lt][ls][meas->detpairs[ls][ld]] +=
							param->xis[lw][lt][ls][vectorIndex]
							* dot 
							* (1.0-dir_reflect(mesh->elem[detectNode].snorm[j],mesh->s[k],user->nindex))
							* mesh->w[k] * mesh->elem[detectNode].sarea[j]
							/ totArea;


						/*/
						
						/* more realistic diffuse transmission */
						meas->prediction_m[lw][lt][ls][meas->detpairs[ls][ld]] += 
							param->xism[lw][lt][ls][vectorIndex]
							* dot * (1.0 - user->R_eff)
							* mesh->w[k] * mesh->elem[detectNode].sarea[j]
							/ totArea;	
															
						//*/													
					}
				}
			}

			//cerr << "\tprediction = " << ld << " " << meas->prediction[lw][lf][ls][ld] << endl;
		}
		}		
	}	

}



