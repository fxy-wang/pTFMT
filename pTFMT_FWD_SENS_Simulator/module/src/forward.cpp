/******************************************************************************

	module:		forward.cpp


	forward model that solves radiative transfer equation for angular intensity
	in discrete ordinates by using unstructured finite-volume discrete-ordinate
	method.

******************************************************************************/
void forward(Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int ntime	= meas->ntime;
	int nwave = meas->nwave;
	int nelem	= mesh->nelem;
	int nmax = mesh->nmax;
	int maxnei = mesh->maxnei;
	int ndetect	= meas->ndetect;
	int nsource	= meas->nsource;
	int itmax	= user->fwd_itmax;
	int	size = param->fwd_size;
	
	int fwdMatType,adjMatType;
	if(user->lightModel == 1){
		user->nz_max = (maxnei+1)*size;
		fwdMatType = 1;
		adjMatType = 2;
	}else if(user->lightModel == 2){
		user->nz_max = (maxnei+2)*size;
		fwdMatType = 1;
		adjMatType = 2;
	}else if(user->lightModel == 3){
		user->nz_max = (maxnei+nmax)*size;
		fwdMatType = 1;
		adjMatType = 2;
	}else{
		cerr << "\tERROR: light propagation model entry should be either dif or rte!!\n\n";
		exit(-1);
	}	

	int nz_max = user->nz_max;
	double *Aex = new double [nz_max];
	double *Aem = new double [nz_max];
	int *rowIndex  = new int [size+1];
	int *columnIndex = new int [nz_max];
	
	/*----------------------------------------------------------------------------------*/
	//	excitation field:
	//	sparse matrix Aex with row-compressed storage (RCS) scheme
	/*----------------------------------------------------------------------------------*/
	sparseMatrixEx(Aex, rowIndex, columnIndex, mesh, meas, param, user);
	nz_max = user->nz_max;
 	//csrRearrange(size,nz_max,rowIndex,columnIndex,Aex); 

	// declare variables
	double *AAex = new double [nz_max]; 
	// AA:= sparse matrix for Laplace transformed equation

	//csrRearrange(size,nz_max,rowIndex,columnIndex,AAex); 	
 	

	/*----------------------------------------------------------------------------------*/
	//	emission field: 
	//	sparse matrix Aem with row-compressed storage (RCS) scheme
	/*----------------------------------------------------------------------------------*/
	sparseMatrixEm(Aem, rowIndex, columnIndex, mesh, meas, param, user);
	nz_max = user->nz_max;
 	//csrRearrange(size,nz_max,rowIndex,columnIndex,Aem); 

	// declare variables
	double *AAem = new double [nz_max]; 
	//csrRearrange(size,nz_max,rowIndex,columnIndex,AAem);
 

	/*----------------------------------------------------------------------------------*/
	//	solve forward problems for sources and wavelengths with OpenMP
	/*----------------------------------------------------------------------------------*/
 	fill(param->xis.begin(),param->xis.end(),Vec3d(ntime,Vec2d(nsource,Vec1d(size,0.0))));
 	fill(param->xism.begin(),param->xism.end(),Vec3d(ntime,Vec2d(nsource,Vec1d(size,0.0))));

  int numDevs = user->numDevs;
	int lw,ls;
  
	double start_t,end_t,elapse_t;
	printf("\t...solving %d forward problems with %d CPUs\n\n",(nsource*ntime),numDevs);
	start_t = omp_get_wtime();
	
	#pragma omp parallel for collapse(2) num_threads(numDevs) private(lw,ls)
	for(lw=0;lw<nwave;lw++){
	for(ls=0;ls<nsource;ls++){
	for(int lt=0;lt<ntime;lt++){
	
		vector<double> rhs(size,0.0);		
		
		fwdRHSEx(lw, lt, ls, mesh, meas, param, user, rhs);
		
		for(int nz = 0; nz < user->nz_max; nz++ ){
		AAex[nz] = Aex[nz];
		}
		for(int irow = 0; irow < size; irow++ ){
			for(int nz = rowIndex[irow]; nz <rowIndex[irow+1]; nz++){
				if(irow == columnIndex[nz]){
					AAex[nz] = param->cdv[irow] * meas->Ls[lt] + AAex[nz];
				}
			}
		}		
		
		pbicgstab (fwdMatType, 1, AAex, rowIndex, columnIndex, param->xis[lw][lt][ls],
					rhs, user->nz_max, user->fwd_itmax, user->fwd_tol, user->fwd_monitor);
		
		computeExcitationField(lw, lt, ls, mesh, meas, param, user);
		
		fwdRHSEm(lw, lt, ls, mesh, meas, param, user, rhs);
		
		for(int nz = 0; nz < user->nz_max; nz++ ){
		AAem[nz] = Aem[nz];
		}
		for(int irow = 0; irow < size; irow++ ){
			for(int nz = rowIndex[irow]; nz <rowIndex[irow+1]; nz++){
				if(irow == columnIndex[nz]){
					AAem[nz] = param->cdv[irow] * meas->Ls[lt] + AAem[nz];
				}
			}
		}				
		
		pbicgstab (fwdMatType, 1, AAem, rowIndex, columnIndex, param->xism[lw][lt][ls],
					rhs, user->nz_max, user->fwd_itmax, user->fwd_tol, user->fwd_monitor);
					
		radianceEm(lw, lt, ls, mesh, meas, param, user);

		if(user->fwd_monitor == 1) saveFWD(lw, lt, ls, mesh, meas, param, user);
		vector<double> ().swap(rhs);
		
	}
	}
	}
			
	end_t = omp_get_wtime();
	elapse_t = end_t - start_t;
	printf("\tThe %d forward runs took %f seconds with %d CPUs\n\n",
					(nsource*ntime),elapse_t,numDevs);

	/*----------------------------------------------------------------------------------*/
	//	write predictions of detector readings for sources, frequencies, and wavelengths
	/*----------------------------------------------------------------------------------*/
	double amplitude_m;
  double std_dev_amp_m;
	size_t dot_pos = user->detector_readings.find_last_of(".");
	string filename, filenameOnly;
	filenameOnly = (dot_pos == string::npos) ? user->detector_readings : user->detector_readings.substr(0,dot_pos);
	filename = filenameOnly + "_mat.txt";		
	
	ofstream fwd_det;
	fwd_det.open(filename.c_str());

	for(int lw=0;lw<nwave;lw++){
	for(int ls=0;ls<nsource;ls++){
	for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
  for(int lt=0;lt<ntime;lt++){
		amplitude_m = meas->prediction_m[lw][lt][ls][meas->detpairs[ls][ld]];
		std_dev_amp_m = amplitude_m / pow(10.0,user->SNRAmp/10.0);
		amplitude_m  = amplitude_m + gaussrand2() * std_dev_amp_m;
		fwd_det << setprecision(6) << scientific 
						<< meas->prediction_m[lw][lt][ls][meas->detpairs[ls][ld]] << "\n";
	}
	}
  }
	}
	fwd_det.close();	
	
	delete [] Aex;
	delete [] AAex;
	delete [] Aem;
	delete [] AAem;
	delete [] rowIndex;
	delete [] columnIndex;	

}


