/******************************************************************************

	module:		sensitivity.cpp

	function that obtains the gradients of the objective function with respect to
	the optical properties, such as the absorption and scattering coefficients.

******************************************************************************/

void sensitivity(int lw, int lt, int ls, int ld, Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int nelem	= mesh->nelem;
	int nmax = mesh->nmax;
	int vectorIndex;
	int nei;
		
	double onepi = 2.0*asin(1.0);	
	double sum;
	vector <double> sensitivity(nelem,0.0);

	/* SP1 model (begin) */
	if(user->lightModel == 1){		
		for(int ele=0;ele<nelem;ele++){
			sensitivity[ele] = param->adjdet[lw][lt][ld][ele]
										* user->qYield * param->excitation[lw][lt][ls][ele] 
										* mesh->elem[ele].vol / (1.0 - meas->Ls[lt] * user->lifeTime);
		}
	}
	/* SP1 model (end) */
	
	/* RTE model (begin) */
	else if(user->lightModel == 3){
		for(int ele=0;ele<nelem;ele++){
			sum = 0.0;
			for(int dir=0;dir<nmax;dir++){
				vectorIndex=nelem*dir+ele;
				sum += param->adjdet[lw][lt][ld][vectorIndex] 
								* param->excitation[lw][lt][ls][vectorIndex]
								/(1.0 - meas->Ls[lt] * user->lifeTime);
			}
			sensitivity[ele]=sum * mesh->elem[ele].vol*user->qYield/4.0/onepi;
		}			
	}
	/* RTE model (end) */
			
	/*
	copy sensitivity coefficients
	*/
	for(int ele=0;ele<nelem;ele++){
		param->S[lw][lt][ls][ld][ele]= sensitivity[ele];
	}

	// release memory
	vector<double> ().swap(sensitivity);

}
