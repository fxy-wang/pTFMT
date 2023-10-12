/******************************************************************************

	module:		fwdrhsEx

	function that constructs the forward source vect b from the system Ax = b.

******************************************************************************/
void fwdRHSEx(int lw, int lt, int ls, Mesh *mesh, Meas *meas, Param *param, 
							User *user, vector<double> &rhs)
{
	int nelem	= mesh->nelem;
	int nmax 	= mesh->nmax;
	int size	= param->fwd_size;
	int sourceNode,vectorIndex;

	double dot;
	double onepi = 2.0*asin(1.0);
	double angularSourcePower;


	/*
	if meas->nsrcpersrc[ls] = 1 => point illumination
	else meas->nsrcpersrc[ls] > 1 => area illumination
	*/

	/*
	SP1(DE)
	*/
	if(user->lightModel == 1){
		fill(rhs.begin(),rhs.end(),0.0);		
		for(int src=0;src<meas->nsrcpersrc[ls];src++){
			sourceNode=meas->lsource[meas->srcpairs[ls][src]];
			//construct the forward rhs vector from AX=b (begin)
			for(int j=0;j<mesh->elem[sourceNode].snum;j++){
				if(mesh->elem[sourceNode].sbound[j] == 1){
					angularSourcePower = 0.0;//2.0 * meas->qin * onepi;
					//*
					for(int k=0;k<nmax;k++){
							dot=vecDot(mesh->s[k],mesh->elem[sourceNode].snorm[j]);
							if(dot < 0){
								angularSourcePower += meas->qin*fabs(dot)*2.0*mesh->w[k];
							}
					}
					//*/
					rhs[sourceNode] += (angularSourcePower/(1.0+3.0*user->R[1])) 
															* mesh->elem[sourceNode].sarea[j];
				}
			}
		}	
		
	}
	
	/*
	SP3
	*/
	int vectorIndex1,vectorIndex2;
	double s1,s2;
	if(user->lightModel == 2){
		fill(rhs.begin(),rhs.end(),0.0);
		for(int src=0;src<meas->nsrcpersrc[ls];src++){
			sourceNode=meas->lsource[meas->srcpairs[ls][src]];
			vectorIndex1 = sourceNode;
  		vectorIndex2 = sourceNode+nelem;
			//construct the forward rhs vector from AX=b (begin)
			for(int j=0;j<mesh->elem[sourceNode].snum;j++){
				if(mesh->elem[sourceNode].sbound[j] == 1){
      		s1=0.0;
      		s2=0.0;
      		for(int k=0;k<nmax;k++){     		  
     				dot=vecDot(mesh->s[k],mesh->elem[sourceNode].snorm[j]);
     				if(dot < 0){
   				  	s1 = s1+(meas->qin/onepi)*fabs(dot)*mesh->w[k];
   				  	s2 = s2+(meas->qin/onepi)*(5.0*pow(fabs(dot),3)-3.0*fabs(dot))*mesh->w[k];
  	   	  	}
  	 		 	}
   				rhs[vectorIndex1] = rhs[vectorIndex1]
   														+ get_source_coeff(0,sourceNode,s1,s2,param,user)
   														* mesh->elem[sourceNode].sarea[j];
   				rhs[vectorIndex2] = rhs[vectorIndex2]
   														+ get_source_coeff(1,sourceNode,s1,s2,param,user)
   														* mesh->elem[sourceNode].sarea[j];
				}
			}
		}		
	}

	/*
	RTE
	*/	
	if(user->lightModel == 3){
		fill(rhs.begin(),rhs.end(),0.0);
		//construct forward RHS vector from Ax=b (begin)	
		for(int src=0;src<meas->nsrcpersrc[ls];src++){		
			sourceNode=meas->lsource[meas->srcpairs[ls][src]];
			for(int k=0;k<nmax;k++){
				vectorIndex=nelem*k+sourceNode;
				for(int j=0;j<mesh->elem[sourceNode].snum;j++){
					if(mesh->elem[sourceNode].sbound[j] == 1){
						dot=vecDot(mesh->s[k],mesh->elem[sourceNode].snorm[j]);
						if(dot < 0.0){
							rhs[vectorIndex] = rhs[vectorIndex]
																- meas->qin * dot * mesh->elem[sourceNode].sarea[j];
						}
					}				
				}
			}	
		}
	}
	
	//check the sum of RHS to make sure it is not zero.
	if(vecNorm(rhs) < 1.0e-25){
		cerr << "\tError in the vector b of Ax=b for the forward problem\n";
		cerr << "\tCheck if the source location "<< ls << " is correct and try again!\n";
	}

}

/******************************************************************************

	module:		fwdrhsEm.cpp

	release:		version 1.3
	author:		Hyun Keol Kim, Columbia University, NY 10027, USA
	e-mail:		hkk2107@columbia.edu
	date:			May 10, 2005
	remarks:

	function that constructs the forward source vect b from the system Ax = b.

******************************************************************************/
void fwdRHSEm(int lw, int lt, int ls, Mesh *mesh, Meas *meas, Param *param, 
							User *user, vector<double> &rhs)
{
	int nelem	= mesh->nelem;
	int nmax 	= mesh->nmax;
	int size	= param->fwd_size;
	int sourceNode,vectorIndex;

	double dot;
	double onepi = 2.0*asin(1.0);
	double angularSourcePower;


	/*
	SP1(DE)
	*/
	if(user->lightModel == 1){
		fill(rhs.begin(),rhs.end(),0.0);
		for(int ele=0;ele<nelem;ele++){
			rhs[ele] = user->qYield * param->qs[ele] * param->excitation[lw][lt][ls][ele]
								* mesh->elem[ele].vol / (1.0 - meas->Ls[lt] * user->lifeTime);			
		}		
	}
	
	/*
	SP3
	*/
	int vectorIndex1,vectorIndex2;
	double s1,s2;
	if(user->lightModel == 2){
		fill(rhs.begin(),rhs.end(),0.0);
		for(int ele=0;ele<nelem;ele++){
   		//rhs[ele] = ;
   		//rhs[ele+nelem] = ;
		}		
	}

	/*
	RTE
	*/	
	if(user->lightModel == 3){
		fill(rhs.begin(),rhs.end(),0.0);
		for(int ele=0;ele<nelem;ele++){
			for(int k=0;k<nmax;k++){
				vectorIndex=nelem*k+ele;
				rhs[vectorIndex] = user->qYield * param->qs[ele] 
														* (param->excitation[lw][lt][ls][ele]/4.0/onepi)
														* mesh->elem[ele].vol / (1.0 - meas->Ls[lt] * user->lifeTime);
			}	
		}
	}
	
	//check the sum of RHS to make sure it is not zero.
	if(vecNorm(rhs) < 1.0e-25){
		cerr << "\tError in the vector b of Ax=b for the forward problem\n";
		cerr << "\tCheck if the source location "<< ls << " is correct and try again!\n";
	}

}


