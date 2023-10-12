/******************************************************************************

	module:		setmedium.cpp

	function that allows perturbation in medium for numerical study

******************************************************************************/

void setMedium(int argc, char **filename, Mesh *mesh, Param *param, User *user)
{

		/*
		define the medium
		*/
		for(int	ele=0;ele<mesh->nelem;ele++){
			
			param->xka[ele] = user->xka_b;
			param->xsig[ele] = user->xsig_b;
			param->qs[ele] = user->qs_b;

			// gel layer
			if(mesh->elem[ele].center[2] > 0.05 && mesh->elem[ele].center[2] < (user->gel_layer_thickness+0.01)){
				param->xka[ele] = user->gel_layer_mua;
				param->xsig[ele] = user->gel_layer_mus;						
			}
			
			// embedded light sources
			if(user->nqs_o > 0){
				
				for(int i=0;i<user->nqs_o;i++){      	        		
					if(vecDist(mesh->elem[ele].center,user->qs_o_location[i]) 
						<= (user->qs_o_radius[i])){
						param->qs[ele]	=	user->qs_o[i];
					}
				}
			}//
										
		}

		
		/*
		write optical properties in tecplot format
		*/
		ofstream exact;
		exact.open(user->fwd_optical_properties.c_str());
		exact << "VARIABLES=\t" << "\"X[$cm$]\"" << ", \"Y[$cm$]\"" << ", \"Z[$cm$]\""
				 << ", \"$\\mu_{a}$[$cm^{-1}$]\"" << ", \"$\\mu_{s}$[$cm^{-1}$]\"" 
				 << ", \"$\\mu_{axf}$[$cm^{-1}$]\"" << "\n";
		exact << "ZONE N=\t" << mesh->nelem << ", E=\t"<< mesh->mshelem 
				<< ", DATAPACKING= POINT, ZONETYPE= "<<mesh->elemtype<<"\n";

		for(int ele=0;ele<mesh->nelem;ele++){
			for(int k=0;k<3;k++) exact<<mesh->elem[ele].center[k]<<"\t";
			exact << param->xka[ele] << "\t" << param->xsig[ele] << "\t" << param->qs[ele] << "\n";
		}

		for(int ele=0;ele<mesh->mshelem;ele++){
			for(int vt=0;vt<mesh->vertnum;vt++) exact << mesh->elm[ele][vt] << "\t";
			exact << "\n";
		}

		exact.close();

    //cout << "\tset the medium successfully!\n\n";

}

