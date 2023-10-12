/******************************************************************************

	module:		allocate.cpp

	function that allocates memory

******************************************************************************/
void allocate(Mesh *mesh, Meas *meas, Param *param, User *user)
{
	//cerr<<"\t...allocating memory\n\n";

	int nwave = meas->nwave;
	int	ntime = meas->ntime;
	int nsource = meas->nsource;
	int ndetect = meas->ndetect;
	int nelem = mesh->nelem;
	int	nmax = mesh->nmax;

	if(user->lightModel == 1){
		param->fwd_size = nelem;
	}else if(user->lightModel == 2){
		param->fwd_size = 2 * nelem;
	}else if(user->lightModel == 3){
		param->fwd_size = nelem * nmax;
	}	

	// constant: dv / mspeed / dt	
	param->cdtdv.resize(param->fwd_size,0.0);
	param->cdv.resize(param->fwd_size,0.0);
	if(user->lightModel = 1){
		for(int ele=0;ele<nelem;ele++){
			param->cdv[ele] = mesh->elem[ele].vol/user->speed;
		}
	}else if(user->lightModel = 2){
		for(int ele=0;ele<nelem;ele++){
			param->cdv[ele] = mesh->elem[ele].vol/user->speed;
			param->cdv[ele+nelem] = mesh->elem[ele].vol/user->speed;
		}
	}else if(user->lightModel = 3){
		for(int dir=0;dir<nmax;dir++){
			for(int ele=0;ele<nelem;ele++){
				param->cdv[ele+dir*nelem] = mesh->elem[ele].vol/user->speed;
			}
		}
	}	
	
 	param->xis.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(param->fwd_size))));
 	param->xism.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(param->fwd_size))));
 	param->axi.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(param->fwd_size))));
 	param->axim.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(param->fwd_size))));
	param->excitation.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(param->fwd_size))));

	
	/* memory allocation for measurement */
	meas->prediction.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(ndetect))));
	meas->weight.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(ndetect))));
	meas->prediction_m.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(ndetect))));
	meas->weight_m.resize(nwave,Vec3d(ntime,Vec2d(nsource,Vec1d(ndetect))));
		
	/* memory allocation for parameters */	
	param->xka.resize(mesh->nelem,0.0);
	param->xsig.resize(mesh->nelem,0.0);
	param->xdif.resize(mesh->nelem,0.0);
	param->qs.resize(mesh->nelem,0.0);
	param->lifeTime.resize(mesh->nelem,0.0);
		
  /* sensing matrix */
  param->adjdet.resize(nwave,Vec3d(ntime,Vec2d(ndetect,Vec1d(param->fwd_size))));
  param->S.resize(nwave,Vec4d(ntime,Vec3d(nsource, Vec2d(ndetect, Vec1d(param->fwd_size)))));
  
  //cerr << "\tdone!\n";
}

