/******************************************************************************

	module:		readsetup.cpp


	function that reads all required information from the setup file

******************************************************************************/

void readSetup(Mesh *mesh, Meas *meas, Param *param, User *user)
{
	cout<<"\n\n";
	cout<<"\t...reading input information: " << user->input_ini << "\n\n";
	
	string str;
	ifstream fin;
	fin.open(user->input_ini.c_str());
	if(!fin.is_open()){
		cerr << setw(15) << "Error opening file: " << user->input_ini << endl;
		exit(-1);
	}

	/**************************************/
	// common settings
	/**************************************/
	fin >> str;
	fin >> str >> user->numCPUDevs;
 	fin >> str >> user->prob;
	fin >> str >> user->dim;
 	fin >> str >> user->lightModel;
	fin >> str >> mesh->ns;
	
	fin >> str >> meas->nsource; // number of source illuminations
	fin >> str >> meas->ndetect; // number of detector nodes 
	fin >> str >> user->nindex;
	fin >> str >> user->g;
	user->speed = 0.029987/user->nindex; // unit of [cm/picosecond]

	//calculate diffuse reflectivity momentum in case it is needed
	reflect_momentum(user); // calculate reflectivity momentum

	fin >> str >> user->lifeTime;
	fin >> str >> user->qYield;
	fin >> str >> meas->qin;

  // number of Laplace temporal filters
	fin >> str >> meas->ntime;
	// Laplace transform factor: exponential rise(early time) or decay(late time)
	fin >> str;
	meas->Ls.resize(meas->ntime);
	for(int i=0;i<meas->ntime;i++) {
		fin >> meas->Ls[i];
	}
	
	fin >> str >> meas->nwave;
	fin >> str;
	meas->wave.resize(meas->nwave);
	for(int lw=0;lw<meas->nwave;lw++){
		fin >> meas->wave[lw];
	}
	
	// background optical properties
	fin >> str >> user->xka_b;
  fin >> str >> user->xsig_b;
  fin >> str >> user->qs_b;
  user->qs_b = user->qs_b + 1.0e-6; // to avoid NaN error

	// gel layer properties
  fin >> str >> user->gel_layer_thickness;
  fin >> str >> user->gel_layer_mua;
  fin >> str >> user->gel_layer_mus;

  // min & max of absorption coefficient, xka
	fin >> str >> user->qs_min >> user->qs_max;
	user->qs_min = user->qs_min + 1.0e-6; // to avoid NaN error

  // scaling of variables using linear transformation: x=Dy+c
 	user->qs_d = (user->qs_max-user->qs_min)/2.0;
 	user->qs_c = (user->qs_max+user->qs_min)/2.0;

	// mesh and S-D config files
	fin >> str >> user->finite_volume_mesh;
	fin >> str >> user->src_det_cfg;

	/*************************************************/
	// forward & sensitivity settings
	/*************************************************/
  fin >> str;
	fin >> str >> user->fwd_monitor;
	fin >> str >> user->fwd_tol;
	fin >> str >> user->fwd_itmax;

	// embedded objects
	fin >> str >> user->nqs_o;
	if(user->nqs_o > 0){
		user->qs_o.resize(user->nqs_o,0);
		user->qs_o_location.resize(user->nqs_o,vector<double> (3));
		user->qs_o_radius.resize(user->nqs_o,0);
		for(int i=0;i<user->nqs_o;i++){
			fin >> user->qs_o[i];
			for(int j=0;j<3;j++) fin >> user->qs_o_location[i][j];
			fin >> user->qs_o_radius[i];
		}
	}

	
	// noise
	fin >> str >> user->SNRAmp;

  // file names
	fin >> str >> user->fwd_fluence;
	fin >> str >> user->fwd_optical_properties;
	fin >> str >> user->detector_readings;

	fin.close();

}



