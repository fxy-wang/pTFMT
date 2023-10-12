/******************************************************************************

	module:		readmesh.cpp


	function that reads all mesh information including a finite volume and a GiD
	element.

******************************************************************************/

void readMesh(Mesh *mesh, Meas *meas, Param *param, User *user)
{
	//cerr <<"\t...reading finite-volume mesh\n\n";
	
	param->minX=param->minY=param->minZ=1.0e25;
	param->maxX=param->maxY=param->maxZ=-1.0e25;

	mesh->maxSurfNum = 0;

	string str;
	/*
	reading finite volume information(begin)
	*/
	ifstream fvm_file;
	fvm_file.open(user->finite_volume_mesh.c_str());
	if(!fvm_file.is_open()){
		cerr << setw(15) << "Error opening file: " << user->finite_volume_mesh 	<< endl;
		exit(-1);
	}
	fvm_file >> str;
	fvm_file >> str;
	if(str.compare("ver1.2+")==0){
		fvm_file >> str;
		//cerr << "\treading FV mesh in CVFEM-v1.1+ format ... \n";
	}else{
		cerr << "\tERROR: the input FV mesh format is not compatible with \n";
		cerr << "\t\tthe current version of Genecis code\n";
		cerr << "\tFV mesh MUST be created with CVFEM version 1.1 or higher\n";
		exit(-1); 
	}
	fvm_file >> mesh->nelem >> mesh->mshelem >> mesh->maxnei;
	mesh->elem.resize(mesh->nelem);
	mesh->center.resize(3);

	mesh->center[0] = mesh->center[1] = mesh->center[2] = 0.0;
	mesh->avg_vol = 0.0;
	mesh->bnode_num = 0;

	for (int ele=0;ele<mesh->nelem;ele++){

		mesh->elem[ele].center.resize(3);
		fvm_file  >>  mesh->elem[ele].center[0]  >>  mesh->elem[ele].center[1]  >>  mesh->elem[ele].center[2];

		mesh->center[0] = mesh->center[0] + mesh->elem[ele].center[0] / double(mesh->nelem);
		mesh->center[1] = mesh->center[1] + mesh->elem[ele].center[1] / double(mesh->nelem);
		mesh->center[2] = mesh->center[2] + mesh->elem[ele].center[2] / double(mesh->nelem);

		// information for conversion to structured data
		param->minX=min(mesh->elem[ele].center[0],param->minX);
		param->maxX=max(mesh->elem[ele].center[0],param->maxX);
		param->minY=min(mesh->elem[ele].center[1],param->minY);
		param->maxY=max(mesh->elem[ele].center[1],param->maxY);
		param->minZ=min(mesh->elem[ele].center[2],param->minZ);
		param->maxZ=max(mesh->elem[ele].center[2],param->maxZ);
		// end

		fvm_file  >>  mesh->elem[ele].vol;
		fvm_file  >>  mesh->elem[ele].snum;

		mesh->maxSurfNum = max(mesh->maxSurfNum,mesh->elem[ele].snum);
		mesh->avg_vol = mesh->avg_vol + mesh->elem[ele].vol/ double (mesh->nelem);

		int snum;
		snum=mesh->elem[ele].snum;

		mesh->elem[ele].sarea.resize(snum);
		mesh->elem[ele].sdist.resize(snum);
		mesh->elem[ele].sbound.resize(snum);
		mesh->elem[ele].snorm.resize(snum, vector <double> (3));
		mesh->elem[ele].snei.resize(snum, vector <int> (2));
		
		// determine whether a node is located on the boundary
		for(int sf=0;sf<snum;sf++){
			fvm_file  >>  mesh->elem[ele].sarea[sf];
			fvm_file  >>  mesh->elem[ele].sdist[sf];
			fvm_file  >>  mesh->elem[ele].sbound[sf];
			fvm_file  >>  mesh->elem[ele].snorm[sf][0]  >>  mesh->elem[ele].snorm[sf][1]  >>  mesh->elem[ele].snorm[sf][2];
			fvm_file  >>  mesh->elem[ele].snei[sf][0]  >>  mesh->elem[ele].snei[sf][1];
		}
		/*
		save indices of boundary nodes and compute DCT basis functions
		*/
		fvm_file  >>  mesh->elem[ele].layer;
		if(mesh->elem[ele].layer == 0) mesh->bnode_index.push_back(ele);
	}

	//read mesh element type and connectivity information
	int tmpType;
	fvm_file >> tmpType >> mesh->vertnum;
  if(tmpType == 2){
		mesh->elemtype = "FETriangle";
  }else if(tmpType == 3){
		mesh->elemtype = "FEQuadrilateral";
  }else if(tmpType == 4){
		mesh->elemtype = "FETetrahedron";
  }else if(tmpType == 5){
		mesh->elemtype = "FEBrick";
  }
	
	mesh->elm.resize(mesh->mshelem, vector <int> (mesh->vertnum));
	for(int msh=0;msh<mesh->mshelem;msh++){
		for(int k=0;k<mesh->vertnum;k++) {
			fvm_file >> mesh->elm[msh][k];
		}
	}
	fvm_file.close();

	param->Lx=fabs(param->minX-param->maxX);
	param->Ly=fabs(param->minY-param->maxY);
	param->Lz=fabs(param->minZ-param->maxZ);

}

