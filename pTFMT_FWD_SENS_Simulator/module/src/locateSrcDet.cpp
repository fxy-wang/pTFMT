/******************************************************************************

	module:		locatesrcdet.cpp

	function that loads all coordinate information about sources and detectors

******************************************************************************/

void locateSrcDet(Mesh *mesh, Meas *meas, Param *param, User *user)
{
	
	string comments;
	
	/* memory allocation for sources & detectors */
	meas->ldetector.resize(meas->ndetect);
	meas->nodeldetector.resize(meas->ndetect);
	meas->nsrcpersrc.resize(meas->nsource);
	meas->ndetpersrc.resize(meas->nsource);
	meas->detpairs.resize(meas->nsource);
	meas->srcpairs.resize(meas->nsource);

	meas->detector.resize(meas->ndetect, vector <double> (3));

	/*
	reading source-detector pairs
	*/
	ifstream src_det;
	src_det.open(user->src_det_cfg.c_str());
	if(!src_det.is_open()){
		cerr << "Error opening file:\n";
		cerr << user->src_det_cfg<<"\n";
		exit(-1);
	}

	// source coordinates
	int tmpInt;
	int srcPtsNum; // number of mesh nodes used for source illuminations
	src_det >> comments;
	src_det >> srcPtsNum;
	meas->lsource.resize(srcPtsNum);	
	meas->source.resize(srcPtsNum, vector <double> (3));	
	for(int ls=0;ls<srcPtsNum;ls++){
		src_det >> tmpInt;
		for(int k=0;k<3;k++){
		 src_det >> meas->source[ls][k];
		}	
	}


	// # of detectors
	int detPtsNum;
	src_det >> comments;
	meas->detector.resize(meas->ndetect, vector <double> (3));	
	for(int ld=0;ld<meas->ndetect;ld++){
		src_det >> tmpInt;
		for(int k=0;k<3;k++){
		 src_det >> meas->detector[ld][k];
		}	
	}
	
	// source-detector pair
	src_det >> comments;
	for(int ls=0;ls<meas->nsource;ls++){
		src_det >> meas->nsrcpersrc[ls]; // number of mesh nodes used for each illumination
		meas->srcpairs[ls].resize(meas->nsrcpersrc[ls]);
		for(int lls=0;lls<meas->nsrcpersrc[ls];lls++){
			src_det >> meas->srcpairs[ls][lls]; // source-node pairs with each illumination 
		}	
		src_det >> meas->ndetpersrc[ls]; // number of detector nodes associated with each illumination
		meas->detpairs[ls].resize(meas->ndetpersrc[ls]);
		for(int ld=0;ld<meas->ndetpersrc[ls];ld++){
			src_det >> meas->detpairs[ls][ld];
		}
	}
	while(src_det >> comments){
		if(src_det.eof()){ 
			src_det.close();
		}else{
			cerr << "\tERROR in reading the source detector configuration file!!\n";
			exit(-1);
		}	
	}


	/* find the mesh node closest to the user-specified source location */
		int nbound;
		double distance;
		vector <double> dist2Src(srcPtsNum,1.0e25);
		vector <double> dist2Det(meas->ndetect,1.0e25);
		for(int ele=0;ele<mesh->nelem;ele++){
			//sources
			nbound=0;
			for(int sf=0;sf<mesh->elem[ele].snum;sf++){
				if(mesh->elem[ele].sbound[sf] == 1) nbound=nbound+1;
			}
			for(int ls=0;ls<srcPtsNum;ls++){
				distance=vecDist(mesh->elem[ele].center,meas->source[ls]);		
				if(mesh->elem[ele].vol > 0.0 && distance < dist2Src[ls] && nbound > 0){
					meas->lsource[ls]=ele;
					dist2Src[ls] = distance;
				}
			}
			//detectors
			nbound=0;
			for(int sf=0;sf<mesh->elem[ele].snum;sf++){
				if(mesh->elem[ele].sbound[sf] == 1) nbound=nbound+1;
			}
			for(int ld=0;ld<meas->ndetect;ld++){
				distance=vecDist(mesh->elem[ele].center,meas->detector[ld]);		
				if(mesh->elem[ele].vol > 0.0 && distance < dist2Det[ld] && nbound > 0){
					meas->nodeldetector[ld]=ele;
					dist2Det[ld] = distance;
				}
			}
		}		
	// check the accuracy of found sources and detectors
	int el;
	for(int ls=0;ls<srcPtsNum;ls++){
		el = meas->lsource[ls];
		distance=vecDist(mesh->elem[el].center,meas->source[ls]);
		if(distance > 0.1){		
			cerr << "\t" << ls << "("<< el <<")" 
					<< " source found and error = " << distance << "\n";
		}
	}
	for(int ld=0;ld<meas->ndetect;ld++){
		el = meas->nodeldetector[ld];
		distance=vecDist(mesh->elem[el].center,meas->detector[ld]);	
		if(distance > 0.1){	
			cerr << "\t" << ld << "("<< el <<")" 
					<< " detector found and error = " << distance << "\n";
		}			
	}

}

