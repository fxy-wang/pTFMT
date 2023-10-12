/******************************************************************************
                             pTFMT_FWD_SENS_Simulator
                    (Time-Domain Fluorescence Moleculary Tomography)

	module	:	main.cpp

	Main file that executes forward and sensitivty simulator 

	Forward model solver generates predicted measurements for a given set of optical
	properties based on the diffusion equation, SP3, or RTE model of light 
	propagation

	Sensing matrix solver calculates the sensitivites due to each source-detector pair of
	the model for a given medium using the adjoint theorem

******************************************************************************/
#define	ONEPI	3.1415926535898
#include <omp.h>
#include "./module/system.h"
using namespace std;
#include "./module/tofdot.h"

int	main(int argc, char	**argv)
{
	/* global	variables	*/
	Mesh		mesh;
	Meas		meas;
	Param		param;
	User		user;

	user.input_ini = argv[1];

	/* read	setup file */
	readSetup(&mesh, &meas,	&param,	&user);

	/* distribute N tasks on M devices (cpu or gpu) */
	int nsource = meas.nsource;
 	int maxNumCPUs = omp_get_max_threads();
 	int maxNumGPUs;
	cout << "\n\tNumber of max available CPUs = " << maxNumCPUs << "\n\n";
	int numDevs;
	if(user.numCPUDevs>0){
		numDevs = min(min(nsource,maxNumCPUs-1),user.numCPUDevs);
	}else{
		cerr << endl;
		cerr << "\t************************************************************\n"			
							"\tERROR: invalid entry - number of devices should be \n"
							"\tgreater than or equal to 1!\n";
		exit(-1);	
	}		
	user.numDevs = numDevs;
	
	/* read	control-volume finite-element mesh (CVFEM) */
	readMesh(&mesh, &meas, &param, &user);

	/* read angular meshes sn-quadrature */
	ordinates(&user, &mesh);

	/* allocate	memory for variables */
	allocate(&mesh, &meas, &param, &user);

	/* locate	sources	&	detectors	*/
	locateSrcDet(&mesh, &meas, &param, &user);

	/* set up the medium */
	setMedium(argc, argv, &mesh, &param, &user);
	
	/* main solver */
	if(user.prob == 1){
		forward(&mesh, &meas, &param, &user);		
	}else if(user.prob == 2){
		sensingMatrix(&mesh, &meas, &param, &user);		
	}
	return 0;

}


