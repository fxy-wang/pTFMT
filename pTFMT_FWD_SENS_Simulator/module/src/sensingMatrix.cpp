/******************************************************************************

	module:		sensingMatrix.cpp


	module that computes Jacobian coefficients distribution, given a specific source
	-detector configuration.

******************************************************************************/

void sensingMatrix(Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int ntime	= meas->ntime;
	int nwave = meas->nwave;
	int nelem	= mesh->nelem;
	int nmax = mesh->nmax;
	int maxnei = mesh->maxnei;
	int ndetect	= meas->ndetect;
	int nsource	= meas->nsource;
	int	fwd_monitor = user->fwd_monitor;
	int itmax	= user->fwd_itmax;	
	int	fwd_size = param->fwd_size;
	int fwdMatType,adjMatType; // 1: sparse matrix A, 2: transpose of matrix A
	int numDevs = user->numDevs;

	double tol	= user->fwd_tol;
  
	if(user->lightModel == 1){
		user->nz_max = (maxnei+1)*fwd_size;
		fwdMatType = 1;
		adjMatType = 2;
	}else if(user->lightModel == 2){
		user->nz_max = (maxnei+2)*fwd_size;
		fwdMatType = 1;
		adjMatType = 2;
	}else if(user->lightModel == 3){
		user->nz_max = (maxnei+nmax)*fwd_size;
		fwdMatType = 1;
		adjMatType = 2;
	}else{
		cerr << "\tERROR: light propagation model entry should be either dif or rte!!\n\n";
		exit(-1);
	}	
	 
	int nz_max = user->nz_max;
	double *Aex = new double [nz_max];
	double *Aem = new double [nz_max];
	int *rowIndex  = new int [fwd_size+1];
	int *columnIndex = new int [nz_max];
	
	/*----------------------------------------------------------------------------------*/
	//	excitation field:
	//	sparse matrix Aex with row-compressed storage (RCS) scheme
	/*----------------------------------------------------------------------------------*/
	sparseMatrixEx(Aex, rowIndex, columnIndex, mesh, meas, param, user);
	nz_max = user->nz_max;
 	//csrRearrange(size,nz_max,rowIndex,columnIndex,Aex); 

	// AA:= sparse matrix for Laplace transformed equation
	double *AAex = new double [nz_max]; 

	//csrRearrange(fwd_size,nz_max,rowIndex,columnIndex,AAex); 	

	/*----------------------------------------------------------------------------------*/
	//	emission field: 
	//	sparse matrix Aem with row-compressed storage (RCS) scheme
	/*----------------------------------------------------------------------------------*/
	sparseMatrixEm(Aem, rowIndex, columnIndex, mesh, meas, param, user);
	nz_max = user->nz_max;
 	//csrRearrange(size,nz_max,rowIndex,columnIndex,Aem); 

	// declare variables
	double *AAem = new double [nz_max]; 
	// AA:= sparse matrix for Laplace transformed equation

	//csrRearrange(size,nz_max,rowIndex,columnIndex,AAem);
	
	/*
	initialization
	*/  
 	fill(param->xis.begin(),param->xis.end(),Vec3d(ntime,Vec2d(nsource,Vec1d(fwd_size,0.0))));
 	fill(param->xism.begin(),param->xism.end(),Vec3d(ntime,Vec2d(nsource,Vec1d(fwd_size,0.0))));

	/*
	solve for excitation field
	*/  
  int lw,ls,ld;	
  int detIndex;

	double start_t,end_t,elapse_t;
	printf("\t...solving %d forward problems with %d CPUs\n\n",(nsource*ntime),numDevs);
	start_t = omp_get_wtime();

	#pragma omp parallel for collapse(2) num_threads(numDevs) private(lw, ls)
	for(lw=0;lw<nwave;lw++){
	for(ls=0;ls<nsource;ls++){
	for(int lt=0;lt<ntime;lt++){
	
		vector<double> rhs(fwd_size,0.0);	

		fwdRHSEx(lw, lt, ls, mesh, meas, param, user, rhs);

		for(int nz = 0; nz < user->nz_max; nz++ ){
			AAex[nz] = Aex[nz];
		}
		for(int irow = 0; irow < fwd_size; irow++ ){
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
		for(int irow = 0; irow < fwd_size; irow++ ){
			for(int nz = rowIndex[irow]; nz <rowIndex[irow+1]; nz++){
				if(irow == columnIndex[nz]){
					AAem[nz] = param->cdv[irow] * meas->Ls[lt] + AAem[nz];
				}
			}
		}		
		
		pbicgstab (fwdMatType, 1, AAem, rowIndex, columnIndex, param->xism[lw][lt][ls],
					rhs, user->nz_max, user->fwd_itmax, user->fwd_tol, user->fwd_monitor);
					
		radianceEm(lw, lt, ls, mesh, meas, param, user);
 				
		vector<double> ().swap(rhs);
		
	}
	}
	}
	
	/* 
	solve for sensitivity coefficients
	*/
	int fwd_gmres_mr = 30;
 	fill(param->axi.begin(),param->axi.end(),Vec3d(ntime,Vec2d(nsource,Vec1d(fwd_size,0.0))));
 	fill(param->axim.begin(),param->axim.end(),Vec3d(ntime,Vec2d(nsource,Vec1d(fwd_size,0.0))));
 	fill(param->adjdet.begin(),param->adjdet.end(),Vec3d(ntime,Vec2d(ndetect,Vec1d(fwd_size,0.0))));
	
	printf("\t...solving %d sensitivity problems with %d CPUs\n\n",(ndetect*ntime),numDevs);
	#pragma omp parallel for collapse(2) num_threads(numDevs) private(lw,ld)
	for(lw=0;lw<nwave;lw++){
	for(ld=0;ld<ndetect;ld++){	
		for(int lt=ntime-1;lt>=0;lt--){

 			vector<double> arhs(fwd_size,0.0);

			sensingRHSEm(lw, lt, 0, ld, mesh, meas, param, user, arhs);
			
			pbicgstab(adjMatType, 1, AAem, rowIndex, columnIndex, param->adjdet[lw][lt][ld],
						arhs, user->nz_max, user->fwd_itmax, user->fwd_tol, user->fwd_monitor);

			//pmgmres (adjMatType, 1, AAem, rowIndex, columnIndex, param->adjdet[lw][lt][ld],
			//				arhs, fwd_size, user->nz_max, user->fwd_itmax, fwd_gmres_mr,
			//				user->fwd_tol, user->fwd_monitor);

			vector<double> ().swap(arhs);

		}
	}	
	}

	printf("\t...calculating sensitivity coefficients with %d CPUs\n\n",numDevs);
	#pragma omp parallel for collapse(2) num_threads(numDevs) private(lw,ls,ld)
	for(lw=0;lw<nwave;lw++) {
	for(ls=0;ls<nsource;ls++){
	for(ld=0;ld<meas->ndetpersrc[ls];ld++){	
		detIndex = meas->detpairs[ls][ld];
		for(int lt=ntime-1;lt>=0;lt--){
			sensitivity(lw, lt, ls, detIndex, mesh, meas, param, user);
		}
	}	
	}
	}

	end_t = omp_get_wtime();
	elapse_t = end_t - start_t;
	printf("\tThe %d sensitivity runs took %f seconds with %d CPUs\n\n",
					(ntime*ndetect),elapse_t,numDevs);

	/* 
	print sensing matrix in Matlab format
	*/
	size_t dot_pos = user->detector_readings.find_last_of(".");
	string filename, filenameOnly;	
	filenameOnly = (dot_pos == string::npos) ? user->detector_readings : user->detector_readings.substr(0,dot_pos);
	filename = filenameOnly + "_sensitivity_mat.txt";		
	ofstream fwd;
	fwd.open(filename.c_str());

	for(int lw=0;lw<meas->nwave;lw++) {
	for(int ls=0;ls<meas->nsource;ls++){
	for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
		detIndex  = meas->detpairs[ls][ld];
		for(int lt=0;lt<meas->ntime;lt++){	
			for(int ele=0;ele<mesh->nelem;ele++){
				if(mesh->elem[ele].center[2] < (param->minZ+1.0e-6)){
					fwd << setprecision(6) << scientific << param->S[lw][lt][ls][detIndex][ele] << "\t";
				}
			}
			fwd << "\n";
		}		
	}
	}
	}
	fwd.close();

	/* 
	release memory
	*/
 	delete		[] Aex;
 	delete		[] AAex;
 	delete		[] Aem;
 	delete		[] AAem; 	
 	delete		[] rowIndex;
	delete		[] columnIndex;

	//cerr << "\tcalculated sensitivity coefficients successfully!!\n\n";


}


