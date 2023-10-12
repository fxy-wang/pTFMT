/******************************************************************************

	module:		saveforward.cpp


	function that prints the forward irradiance results on the file.

******************************************************************************/

void saveFWD(int lw, int lt, int ls,
                 Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int vectorIndex;
	int mshelem = mesh->mshelem;
	int nelem   = mesh->nelem;
	int nmax	  = mesh->nmax;

	string 	filename;
	double time = double (lt);
	double wave = meas->wave[lw];
	double small=1.0e-25;
	double excitation,emission;
	double onepi = 2.0*asin(1.0);

	filename = user->fwd_fluence + "_w" + to_string(int(wave))
							+ "_s" + to_string(ls+1) + "_t"+ to_string(int(time))+".tec";

	ofstream fwd;
	fwd.open(filename.c_str());
	fwd << "VARIABLES=\t" << "\"X[$cm$]\"" <<", \"Y[$cm$]\"" << ", \"Z[$cm$]\"" 
			<< ", \"excitation[$W{cm}^{-2}$]\"" << ", \"emission[$W{cm}^{-2}$]\"" << "\n";
	fwd << "ZONE N=\t" << nelem << ", E=\t" << mshelem 
			<< " DATAPACKING= POINT, ZONETYPE= " << mesh->elemtype << "\n";

	for(int ele=0;ele<mesh->nelem;ele++){
		
		if(user->lightModel = 1){
    	excitation = param->xis[lw][lt][ls][ele];
    	emission = param->xism[lw][lt][ls][ele];
    }else if(user->lightModel=2){
    	excitation = param->xis[lw][lt][ls][ele]-(2.0/3.0)*param->xis[lw][lt][ls][ele+nelem];
    	emission = param->xism[lw][lt][ls][ele]-(2.0/3.0)*param->xism[lw][lt][ls][ele+nelem];
    }else if(user->lightModel=3){
			excitation = 0.0;
			for(int dir=0;dir<nmax;dir++){
	   		vectorIndex	= nelem*dir+ele;
	   		excitation += param->xis[lw][lt][ls][vectorIndex]*mesh->w[dir];
	    }    
			emission = 0.0;
			for(int dir=0;dir<nmax;dir++){
	   		vectorIndex	= nelem*dir+ele;
	   		emission += param->xism[lw][lt][ls][vectorIndex]*mesh->w[dir];
	    }    
    }
    
    //fluence = log(fluence+1.0e-25);
		
		for(int k=0;k<3;k++) fwd	<<	mesh->elem[ele].center[k] << "\t";
		fwd	<<	excitation << "\t" << emission << "\n";
	}

	for(int ele=0;ele<mesh->mshelem;ele++){
		for(int vt=0;vt<mesh->vertnum;vt++){
			fwd	<<	mesh->elm[ele][vt]	<<	"\t";
		}
		fwd <<	"\n";
	}
	fwd.close();

}


