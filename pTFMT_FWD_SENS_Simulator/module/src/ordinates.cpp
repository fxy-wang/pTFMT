/******************************************************************************

	module:		ordinates.cpp

	function that applies a level symetry angular discretization into a 4*onepi
	solid angle where ordinate direction vector is set to be (xmu,xxi,eta).

	discrete ordinate number for 2D problem is n(n+2)/2, but its weights are
	for 3D problem. So "2*w(k)" should be used instead of w(k) itself.

       s-2 : 4  ordinates
       s-4 : 12 ordinates
       s-6 : 24 ordinates
       s-8 : 40 ordinates
       s-10: 60 ordinates
       s-12: 84 ordinates
       s-14:112 ordinates
       s-16:144 ordinates

	mesh->nmax : total number of discrete ordinates for "2D" case.
	mesh->nmax(nmax+2)/2: total number of discrete ordinates for "3D" case

 	Also this module normalizes the scattering phase function and stores the
 	results onto a two-dimensional array. see below for details.

******************************************************************************/

void	ordinates(User *user, Mesh *mesh)
{
	double onepi = 2.0*asin(1.0);
	if(user->dim == 2) {
		mesh->nmax	= mesh->ns*(mesh->ns+2)/2;
	}
	else if(user->dim == 3) {
		mesh->nmax 	=	mesh->ns*(mesh->ns+2);
	}
 	/* memory allocation for mesh information */
	mesh->w.resize(mesh->nmax);
	mesh->s.resize(mesh->nmax, Vec1d(3));
	mesh->scatphase.resize(mesh->nmax, vector <double> (mesh->nmax));

//----------------------------------------------------------------------
//---- s-2 data
//----------------------------------------------------------------------
      if(mesh->ns == 2){
        mesh->w[0]=onepi/2.;
        mesh->s[0][0]=0.57735026; 					mesh->s[0][1]=0.57735026;
      }
//----------------------------------------------------------------------
//---- s-4 data
//----------------------------------------------------------------------
      if(mesh->ns == 4){
        mesh->w[0]=onepi/3./2.;         mesh->w[1]=onepi/3./2.;         mesh->w[2]=onepi/3./2.;
        mesh->s[0][0]=1./3.;         mesh->s[1][0]=0.88191710;         mesh->s[2][0]=1./3.;
        mesh->s[0][1]=0.88191710;         mesh->s[1][1]=1./3.;         mesh->s[2][1]=1./3.;
      }
//----------------------------------------------------------------------
//---- s-6 data
//----------------------------------------------------------------------
      if(mesh->ns == 6){
        mesh->w[0]=0.16086125*onepi/2.;					mesh->w[1]=0.17247209*onepi/2.;				mesh->w[2]=0.17247209*onepi/2.;
				mesh->w[3]=0.16086125*onepi/2.;         mesh->w[4]=0.17247209*onepi/2.;      mesh->w[5]=0.16086125*onepi/2.;

        mesh->s[0][0]=0.25819889;            mesh->s[1][0]=0.68313005;           mesh->s[2][0]=0.25819889;
				mesh->s[3][0]=0.93094934;							mesh->s[4][0]=0.68313005;						mesh->s[5][0]=0.25819889;

        mesh->s[0][1]=0.93094934;								mesh->s[1][1]=0.68313005;							mesh->s[2][1]=0.68313005;
				mesh->s[3][1]=0.25819889;								mesh->s[4][1]=0.25819889;							mesh->s[5][1]=0.25819889;
      }
//----------------------------------------------------------------------
//---- s-8 data
//----------------------------------------------------------------------
      if(mesh->ns == 8){
         mesh->w[0]=0.11678847*onepi/2.; 					 mesh->w[1]=0.09325523*onepi/2.;      		mesh->w[2]=0.09325523*onepi/2.;
				 mesh->w[3]=0.09325523*onepi/2.;          mesh->w[4]=0.09010320*onepi/2.;          mesh->w[5]=0.09325523*onepi/2.;
         mesh->w[6]=0.11678847*onepi/2.;          mesh->w[7]=0.09325523*onepi/2.;          mesh->w[8]=0.09325523*onepi/2.;
         mesh->w[9]=0.11678847*onepi/2.;

         mesh->s[0][0]=0.19232747;          mesh->s[1][0]=0.57735027;          mesh->s[2][0]=0.19232747;
         mesh->s[3][0]=0.79352178;          mesh->s[4][0]=0.57735027;         mesh->s[5][0]=0.19232747;
         mesh->s[6][0]=0.96229948;	         mesh->s[7][0]=0.79352178;	         mesh->s[8][0]=0.57735027;
         mesh->s[9][0]=0.19232747;

         mesh->s[0][1]=0.96229948;          mesh->s[1][1]=0.79352178;	         mesh->s[2][1]=0.79352178;
         mesh->s[3][1]=0.57735027;	         mesh->s[4][1]=0.57735027;	         mesh->s[5][1]=0.57735027;
         mesh->s[6][1]=0.19232747;	         mesh->s[7][1]=0.19232747;	         mesh->s[8][1]=0.19232747;
         mesh->s[9][1]=0.19232747;
      }
//----------------------------------------------------------------------
//---- s-10 data
//----------------------------------------------------------------------
      if(mesh->ns == 10){
			  mesh->w[0]=0.089842043*onepi/2.;         mesh->w[1]=0.067288705*onepi/2.;         mesh->w[2]=0.067288705*onepi/2.;
        mesh->w[3]=0.055780071*onepi/2.;         mesh->w[4]=0.053133809*onepi/2.;         mesh->w[5]=0.055780071*onepi/2.;
        mesh->w[6]=0.067288705*onepi/2.;         mesh->w[7]=0.053133809*onepi/2.;         mesh->w[8]=0.053133809*onepi/2.;
        mesh->w[9]=0.067288705*onepi/2.;         mesh->w[10]=0.089842043*onepi/2.;         mesh->w[11]=0.067288705*onepi/2.;
        mesh->w[12]=0.055780071*onepi/2.;         mesh->w[13]=0.067288705*onepi/2.;         mesh->w[14]=0.089842043*onepi/2.;

        mesh->s[0][0]=0.16962228;         mesh->s[1][0]=0.50714192;         mesh->s[2][0]=0.16962228;
        mesh->s[3][0]=0.69686020;         mesh->s[4][0]=0.50714192;         mesh->s[5][0]=0.16962228;
        mesh->s[6][0]=0.84500612;         mesh->s[7][0]=0.69686020;         mesh->s[8][0]=0.50714192;
        mesh->s[9][0]=0.16962228;         mesh->s[10][0]=0.97080202;         mesh->s[11][0]=0.84500612;
        mesh->s[12][0]=0.69686020;         mesh->s[13][0]=0.50714192;         mesh->s[14][0]=0.16962228;

        mesh->s[0][1]=0.97080202;         mesh->s[1][1]=0.84500612;         mesh->s[2][1]=0.84500612;
        mesh->s[3][1]=0.69686020;	        mesh->s[4][1]=0.69686020;         mesh->s[5][1]=0.69686020;
        mesh->s[6][1]=0.50714192;         mesh->s[7][1]=0.50714192;         mesh->s[8][1]=0.50714192;
        mesh->s[9][1]=0.50714192;         mesh->s[10][1]=0.16962228;         mesh->s[11][1]=0.16962228;
        mesh->s[12][1]=0.16962228;         mesh->s[13][1]=0.16962228;         mesh->s[14][1]=0.16962228;
      }
//----------------------------------------------------------------------
//  s-12 data
//----------------------------------------------------------------------
			if(mesh->ns == 12){
      mesh->w[0]=0.07332178*onepi/2.;       mesh->w[1]=0.05266740*onepi/2.;       mesh->w[2]=0.05266740*onepi/2.;
      mesh->w[3]=0.04161495*onepi/2.;      mesh->w[4]=0.03895667*onepi/2.;      mesh->w[5]=0.04161495*onepi/2.;
		mesh->w[6]=0.04161495*onepi/2.;       mesh->w[7]=0.03249018*onepi/2.;       mesh->w[8]=0.03249018*onepi/2.;
      mesh->w[9]=0.04161495*onepi/2.;       mesh->w[10]=0.05266740*onepi/2.;       mesh->w[11]=0.03895667*onepi/2.;
      mesh->w[12]=0.03249018*onepi/2.;      mesh->w[13]=0.03895667*onepi/2.;       mesh->w[14]=0.05266740*onepi/2.;
      mesh->w[15]=0.07332178*onepi/2.;       mesh->w[16]=0.05266740*onepi/2.;       mesh->w[17]=0.04161495*onepi/2.;
      mesh->w[18]=0.04161495*onepi/2.;       mesh->w[19]=0.05266740*onepi/2.;       mesh->w[20]=0.07332178*onepi/2.;

      mesh->s[0][0]=0.15395746;       mesh->s[1][0]=0.45769112;       mesh->s[2][0]=0.15395746;
      mesh->s[3][0]=0.62869660;       mesh->s[4][0]=0.45769112;       mesh->s[5][0]=0.15395746;
      mesh->s[6][0]=0.76225828;       mesh->s[7][0]=0.62869660;       mesh->s[8][0]=0.45769112;
      mesh->s[9][0]=0.15395746;       mesh->s[10][0]=0.87568027;       mesh->s[11][0]=0.76225828;
      mesh->s[12][0]=0.62869660;       mesh->s[13][0]=0.45769112;       mesh->s[14][0]=0.15395746;
      mesh->s[15][0]=0.97600932;       mesh->s[16][0]=0.87568027;       mesh->s[17][0]=0.76225828;
      mesh->s[18][0]=0.62869660;       mesh->s[19][0]=0.45769112;       mesh->s[20][0]=0.15395746;

      mesh->s[0][1]=0.97600932;       mesh->s[1][1]=0.87568027;       mesh->s[2][1]=0.87568027;
      mesh->s[3][1]=0.76225828;       mesh->s[4][1]=0.76225828;       mesh->s[5][1]=0.76225828;
      mesh->s[6][1]=0.62869660;       mesh->s[7][1]=0.62869660;       mesh->s[8][1]=0.62869660;
      mesh->s[9][1]=0.62869660;       mesh->s[10][1]=0.45769112;       mesh->s[11][1]=0.45769112;
      mesh->s[12][1]=0.45769112;       mesh->s[13][1]=0.45769112;       mesh->s[14][1]=0.45769112;
      mesh->s[15][1]=0.15395746;       mesh->s[16][1]=0.15395746;       mesh->s[17][1]=0.15395746;
      mesh->s[18][1]=0.15395746;       mesh->s[19][1]=0.15395746;       mesh->s[20][1]=0.15395746;
      }
//----------------------------------------------------------------------
//  s-14 data
//----------------------------------------------------------------------
			if(mesh->ns == 14){
      mesh->w[0]=0.062171628*onepi/2.;       mesh->w[1]=0.043325697*onepi/2.;       mesh->w[2]=mesh->w[1];
      mesh->w[3]=0.033217606*onepi/2.;       mesh->w[4]=0.030486324*onepi/2.;        mesh->w[5]=mesh->w[3];
      mesh->w[6]=0.031837060*onepi/2.;       mesh->w[7]=0.024545116*onepi/2.;       mesh->w[7]=mesh->w[7];
      mesh->w[9]=mesh->w[6];       mesh->w[10]=mesh->w[3];       mesh->w[11]=mesh->w[7];
      mesh->w[12]=0.019984453*onepi/2.;       mesh->w[13]=mesh->w[7];       mesh->w[14]=mesh->w[3];
      mesh->w[15]=mesh->w[1];       mesh->w[16]=mesh->w[4];       mesh->w[17]=mesh->w[7];
      mesh->w[18]=mesh->w[7];       mesh->w[19]=mesh->w[4];       mesh->w[20]=mesh->w[1];
      mesh->w[21]=mesh->w[0];       mesh->w[22]=mesh->w[1];       mesh->w[23]=mesh->w[3];
      mesh->w[24]=mesh->w[6];       mesh->w[25]=mesh->w[3];       mesh->w[26]=mesh->w[1];
      mesh->w[27]=mesh->w[0];

      mesh->s[0][0]=0.14238965;       mesh->s[1][0]=0.42048076;       mesh->s[2][0]=0.14238965;
      mesh->s[3][0]=0.57735027;       mesh->s[4][0]=0.42048076;       mesh->s[5][0]=0.14238965;
      mesh->s[6][0]=0.69990185;       mesh->s[7][0]=0.57735027;       mesh->s[8][0]=0.42048076;
      mesh->s[9][0]=0.14238965;       mesh->s[10][0]=0.80398498;       mesh->s[11][0]=0.69990185;
      mesh->s[12][0]=0.57735027;       mesh->s[13][0]=0.42048076;       mesh->s[14][0]=0.14238965;
      mesh->s[15][0]=0.89605866;       mesh->s[16][0]=0.80398498;       mesh->s[17][0]=0.69990185;
      mesh->s[18][0]=0.57735027;       mesh->s[19][0]=0.42048076;       mesh->s[20][0]=0.14238965;
      mesh->s[21][0]=0.97951538;       mesh->s[22][0]=0.89605866;       mesh->s[23][0]=0.80398498;
      mesh->s[24][0]=0.69990185;       mesh->s[25][0]=0.57735027;       mesh->s[26][0]=0.42048076;
      mesh->s[27][0]=0.14238965;

      mesh->s[0][1]=0.97951538;       mesh->s[1][1]=0.89605866;      mesh->s[2][1]=0.89605866;
      mesh->s[3][1]=0.80398498;       mesh->s[4][1]=0.80398498;       mesh->s[5][1]=0.80398498;
      mesh->s[6][1]=0.69990185;       mesh->s[7][1]=0.69990185;       mesh->s[8][1]=0.69990185;
      mesh->s[9][1]=0.69990185;       mesh->s[10][1]=0.57735027;       mesh->s[11][1]=0.57735027;
      mesh->s[12][1]=0.57735027;       mesh->s[13][1]=0.57735027;       mesh->s[14][1]=0.57735027;
      mesh->s[15][1]=0.42048076;       mesh->s[16][1]=0.42048076;       mesh->s[17][1]=0.42048076;
      mesh->s[18][1]=0.42048076;       mesh->s[19][1]=0.42048076;       mesh->s[20][1]=0.42048076;
      mesh->s[21][1]=0.14238965;       mesh->s[22][1]=0.14238965;       mesh->s[23][1]=0.14238965;
      mesh->s[24][1]=0.14238965;       mesh->s[25][1]=0.14238965;       mesh->s[26][1]=0.14238965;
      mesh->s[27][1]=0.14238965;
			}
//----------------------------------------------------------------------
//  s-16 data
//----------------------------------------------------------------------
      if(mesh->ns == 16){
      mesh->w[0]=0.05415425*onepi/2.;       mesh->w[1]=0.03679653*onepi/2.;       mesh->w[2]=mesh->w[1];
      mesh->w[3]=0.02777273*onepi/2.;       mesh->w[4]=0.02494275*onepi/2.;       mesh->w[5]=mesh->w[3];
      mesh->w[6]=0.02580284*onepi/2.;       mesh->w[7]=0.01962325*onepi/2.;       mesh->w[8]=mesh->w[7];
      mesh->w[9]=mesh->w[6];       mesh->w[10]=mesh->w[6];       mesh->w[11]=0.01879762*onepi/2.;
      mesh->w[12]=0.01544801*onepi/2.;       mesh->w[13]=mesh->w[11];       mesh->w[14]=mesh->w[6];
      mesh->w[15]=mesh->w[3];       mesh->w[16]=mesh->w[7];       mesh->w[17]=mesh->w[12];
      mesh->w[18]=mesh->w[12];       mesh->w[19]=mesh->w[8];       mesh->w[20]=mesh->w[3];
      mesh->w[21]=mesh->w[1];       mesh->w[22]=mesh->w[4];       mesh->w[23]=mesh->w[7];
      mesh->w[24]=mesh->w[11];       mesh->w[25]=mesh->w[7];       mesh->w[26]=mesh->w[4];
      mesh->w[27]=mesh->w[1];       mesh->w[28]=mesh->w[0];       mesh->w[29]=mesh->w[1];
      mesh->w[30]=mesh->w[3];       mesh->w[31]=mesh->w[6];       mesh->w[32]=mesh->w[6];
      mesh->w[33]=mesh->w[3];       mesh->w[34]=mesh->w[1];       mesh->w[35]=mesh->w[0];

      mesh->s[0][0]=0.13344572;       mesh->s[1][0]=0.39119433;       mesh->s[2][0]=0.13344572;
      mesh->s[3][0]=0.53689687;       mesh->s[4][0]=0.39119433;       mesh->s[5][0]=0.13344572;
      mesh->s[6][0]=0.65075610;       mesh->s[7][0]=0.53689687;       mesh->s[8][0]=0.39119433;
      mesh->s[9][0]=0.13344572;       mesh->s[10][0]=0.74746822;       mesh->s[11][0]=0.65075610;
      mesh->s[12][0]=0.53689687;       mesh->s[13][0]=0.39119433;       mesh->s[14][0]=0.13344572;
      mesh->s[15][0]=0.83302700;       mesh->s[16][0]=0.74746822;       mesh->s[17][0]=0.65075610;
      mesh->s[18][0]=0.53689687;       mesh->s[19][0]=0.39119433;       mesh->s[20][0]=0.13344572;
      mesh->s[21][0]=0.91058181;       mesh->s[22][0]=0.83302700;       mesh->s[23][0]=0.74746822;
      mesh->s[24][0]=0.65075610;       mesh->s[25][0]=0.53689687;       mesh->s[26][0]=0.39119433;
      mesh->s[27][0]=0.13344572;       mesh->s[28][0]=0.98203079;       mesh->s[29][0]=0.91058181;
      mesh->s[30][0]=0.83302700;       mesh->s[31][0]=0.74746822;       mesh->s[32][0]=0.65075610;
      mesh->s[33][0]=0.53689687;       mesh->s[34][0]=0.39119433;       mesh->s[35][0]=0.13344572;

      mesh->s[0][1]=0.98203079;       mesh->s[1][1]=0.91058181;       mesh->s[2][1]=0.91058181;
      mesh->s[3][1]=0.83302700;       mesh->s[4][1]=0.83302700;       mesh->s[5][1]=0.83302700;
      mesh->s[6][1]=0.74746822;       mesh->s[7][1]=0.74746822;       mesh->s[8][1]=0.74746822;
      mesh->s[9][1]=0.74746822;       mesh->s[10][1]=0.65075610;       mesh->s[11][1]=0.65075610;
      mesh->s[12][1]=0.65075610;       mesh->s[13][1]=0.65075610;       mesh->s[14][1]=0.65075610;
      mesh->s[15][1]=0.53689687;       mesh->s[16][1]=0.53689687;       mesh->s[17][1]=0.53689687;
      mesh->s[18][1]=0.53689687;       mesh->s[19][1]=0.53689687;       mesh->s[20][1]=0.53689687;
      mesh->s[21][1]=0.39119433;       mesh->s[22][1]=0.39119433;       mesh->s[23][1]=0.39119433;
      mesh->s[24][1]=0.39119433;       mesh->s[25][1]=0.39119433;       mesh->s[26][1]=0.39119433;
      mesh->s[27][1]=0.39119433;       mesh->s[28][1]=0.13344572;       mesh->s[29][1]=0.13344572;
      mesh->s[30][1]=0.13344572;       mesh->s[31][1]=0.13344572;       mesh->s[32][1]=0.13344572;
      mesh->s[33][1]=0.13344572;       mesh->s[34][1]=0.13344572;       mesh->s[35][1]=0.13344572;
      }
//----------------------------------------------------------------------
//  s-18 data
//----------------------------------------------------------------------
      if(mesh->ns == 18){
      mesh->w[0]=0.432081e-01*onepi/2.;       mesh->w[1]=0.363332e-01*onepi/2.;       mesh->w[2]=mesh->w[1];
      mesh->w[3]=0.146534e-01*onepi/2.;       mesh->w[4]=0.332555e-01*onepi/2.;       mesh->w[5]=mesh->w[3];
      mesh->w[6]=0.348583e-01*onepi/2.;       mesh->w[7]=0.685444e-02*onepi/2.;       mesh->w[8]=mesh->w[7];
      mesh->w[9]=mesh->w[6];       mesh->w[10]=0.137443e-03*onepi/2.;       mesh->w[11]=0.251561e-01*onepi/2.;
      mesh->w[12]=0.210216e-01*onepi/2.;       mesh->w[13]=mesh->w[11];       mesh->w[14]=mesh->w[10];
      mesh->w[15]=mesh->w[6];       mesh->w[16]=mesh->w[11];       mesh->w[17]=0.225641e-14*onepi/2.;
      mesh->w[18]=mesh->w[17];       mesh->w[19]=mesh->w[11];       mesh->w[20]=mesh->w[6];
      mesh->w[21]=mesh->w[3];       mesh->w[22]=mesh->w[7];       mesh->w[23]=mesh->w[12];
      mesh->w[24]=mesh->w[17];       mesh->w[25]=mesh->w[12]; 			mesh->w[26]=mesh->w[7];
      mesh->w[27]=mesh->w[3];       mesh->w[28]=mesh->w[1];       mesh->w[29]=mesh->w[4];
      mesh->w[30]=mesh->w[7];       mesh->w[31]=mesh->w[11];       mesh->w[32]=mesh->w[11];
      mesh->w[33]=mesh->w[7];       mesh->w[34]=mesh->w[4];       mesh->w[35]=mesh->w[1];
      mesh->w[36]=mesh->w[0];       mesh->w[37]=mesh->w[1];       mesh->w[38]=mesh->w[3];
      mesh->w[39]=mesh->w[6];       mesh->w[40]=mesh->w[10];       mesh->w[41]=mesh->w[6];
      mesh->w[42]=mesh->w[3];       mesh->w[43]=mesh->w[1];       mesh->w[44]=mesh->w[0];

      mesh->s[0][0]=0.13090100;       mesh->s[1][0]=0.36838759;       mesh->s[2][0]=0.13090100;
      mesh->s[3][0]=0.50426557;      mesh->s[4][0]=0.36838759;      mesh->s[5][0]=0.13090100;
      mesh->s[6][0]=0.61062109;      mesh->s[7][0]=0.50426557;      mesh->s[8][0]=0.36838759;
      mesh->s[9][0]=0.13090100;       mesh->s[10][0]=0.70102244;       mesh->s[11][0]=0.61062109;
      mesh->s[12][0]=0.50426557;      mesh->s[13][0]=0.36838759;      mesh->s[14][0]=0.13090100;
      mesh->s[15][0]=0.78102933;      mesh->s[16][0]=0.70102244;      mesh->s[17][0]=0.61062109;
      mesh->s[18][0]=0.50426557;      mesh->s[19][0]=0.36838759;      mesh->s[20][0]=0.13090100;
      mesh->s[21][0]=0.85356966;      mesh->s[22][0]=0.78102933;      mesh->s[23][0]=0.70102244;
      mesh->s[24][0]=0.61062109;      mesh->s[25][0]=0.50426557;      mesh->s[26][0]=0.36838759;
      mesh->s[27][0]=0.13090100;      mesh->s[28][0]=0.92041051;      mesh->s[29][0]=0.85356966;
      mesh->s[30][0]=0.78102933;      mesh->s[31][0]=0.70102244;      mesh->s[32][0]=0.61062109;
      mesh->s[33][0]=0.50426557;      mesh->s[34][0]=0.36838759;      mesh->s[35][0]=0.13090100;
      mesh->s[36][0]=0.98271555;      mesh->s[37][0]=0.92041051;      mesh->s[38][0]=0.85356966;
      mesh->s[39][0]=0.78102933;      mesh->s[40][0]=0.70102244;      mesh->s[41][0]=0.61062109;
      mesh->s[42][0]=0.50426557;      mesh->s[43][0]=0.36838759;      mesh->s[44][0]=0.13090100;

      mesh->s[0][1]=0.98271555;       mesh->s[1][1]=0.92041051;      mesh->s[2][1]=0.92041051;
      mesh->s[3][1]=0.85356966;      mesh->s[4][1]=0.85356966;      mesh->s[5][1]=0.85356966;
      mesh->s[6][1]=0.78102933;      mesh->s[7][1]=0.78102933;      mesh->s[8][1]=0.78102933;
      mesh->s[9][1]=0.78102933;      mesh->s[10][1]=0.70102244;      mesh->s[11][1]=0.70102244;
      mesh->s[12][1]=0.70102244;      mesh->s[13][1]=0.70102244;      mesh->s[14][1]=0.70102244;
      mesh->s[15][1]=0.61062109;      mesh->s[16][1]=0.61062109;      mesh->s[17][1]=0.61062109;
      mesh->s[18][1]=0.61062109;	     mesh->s[19][1]=0.61062109;       mesh->s[20][1]=0.61062109;
      mesh->s[21][1]=0.50426557;       mesh->s[22][1]=0.50426557;       mesh->s[23][1]=0.50426557;
      mesh->s[24][1]=0.50426557;      mesh->s[25][1]=0.50426557;      mesh->s[26][1]=0.50426557;
      mesh->s[27][1]=0.50426557;      mesh->s[28][1]=0.36838759;      mesh->s[29][1]=0.36838759;
      mesh->s[30][1]=0.36838759;      mesh->s[31][1]=0.36838759;      mesh->s[32][1]=0.36838759;
      mesh->s[33][1]=0.36838759;      mesh->s[34][1]=0.36838759;      mesh->s[35][1]=0.36838759;
      mesh->s[36][1]=0.13090100;      mesh->s[37][1]=0.13090100;      mesh->s[38][1]=0.13090100;
      mesh->s[39][1]=0.13090100;      mesh->s[40][1]=0.13090100;      mesh->s[41][1]=0.13090100;
      mesh->s[42][1]=0.13090100;      mesh->s[43][1]=0.13090100;      mesh->s[44][1]=0.13090100;
      }
//----------------------------------------------------------------------
//---- Reorder the ordinates and weights set with the given
//---- data above for each order of the S-N approximation.
//----------------------------------------------------------------------
/**************2D case**********************/
      if(user->dim == 2){
		//2nd quadrant
			int ll;
			for(int ia=mesh->nmax/4;ia<mesh->nmax/2;ia++){
      		ll=ia-mesh->nmax/4;
        	mesh->w[ia]=mesh->w[ll];
        	mesh->s[ia][0]=-mesh->s[ll][0];
         	mesh->s[ia][1]=mesh->s[ll][1];
			}

		//3rd quadrant
			for(int ia=mesh->nmax/2;ia<3*mesh->nmax/4;ia++){
        	ll=ia-mesh->nmax/2;
         	mesh->w[ia]=mesh->w[ll];
         	mesh->s[ia][0]=-mesh->s[ll][0];
         	mesh->s[ia][1]=-mesh->s[ll][1];
			}

		//4th quadrant
			for(int ia=3*mesh->nmax/4;ia<mesh->nmax;ia++){
        	 ll=ia-3*mesh->nmax/4;
         	 mesh->w[ia]=mesh->w[ll];
         	 mesh->s[ia][0]=mesh->s[ll][0];
          	 mesh->s[ia][1]=-mesh->s[ll][1];
			}

      	for(int k=0;k<mesh->nmax;k++){
      		mesh->w[k]=2.0*mesh->w[k];
      	}
		//compute eta from (xmu,xxi,eta) w.r.t 2D
      for(int ll=0;ll<mesh->nmax;ll++) mesh->s[ll][2]=sqrt(1. - mesh->s[ll][0]*mesh->s[ll][0]- mesh->s[ll][1]*mesh->s[ll][1]);

      }
/*****************end::2D***********************/


/***************3D case******************/
      if(user->dim == 3){
			int ll;
		//2nd quadrant
			for(int ia=mesh->nmax/8;ia<2*mesh->nmax/8;ia++){
      		ll=ia-mesh->nmax/8;
        	mesh->w[ia]=mesh->w[ll];
        	mesh->s[ia][0]=-mesh->s[ll][0];
         	mesh->s[ia][1]=mesh->s[ll][1];
			}

		//3rd quadrant
			for(int ia=2*mesh->nmax/8;ia<3*mesh->nmax/8;ia++){
        	ll=ia-2*mesh->nmax/8;
         	mesh->w[ia]=mesh->w[ll];
         	mesh->s[ia][0]=-mesh->s[ll][0];
         	mesh->s[ia][1]=-mesh->s[ll][1];
			}

		//4th quadrant
			for(int ia=3*mesh->nmax/8;ia<4*mesh->nmax/8;ia++){
        	 ll=ia-3*mesh->nmax/8;
         	 mesh->w[ia]=mesh->w[ll];
         	 mesh->s[ia][0]=mesh->s[ll][0];
          	 mesh->s[ia][1]=-mesh->s[ll][1];
			}

		//5nd quadrant
			for(int ia=4*mesh->nmax/8;ia<5*mesh->nmax/8;ia++){
      		ll=ia-4*mesh->nmax/8;
        	mesh->w[ia]=mesh->w[ll];
        	mesh->s[ia][0]=mesh->s[ll][0];
         	mesh->s[ia][1]=mesh->s[ll][1];
			}

		//6rd quadrant
			for(int ia=5*mesh->nmax/8;ia<6*mesh->nmax/8;ia++){
        	ll=ia-5*mesh->nmax/8;
         	mesh->w[ia]=mesh->w[ll];
         	mesh->s[ia][0]=-mesh->s[ll][0];
         	mesh->s[ia][1]=mesh->s[ll][1];
			}

		//7th quadrant
			for(int ia=6*mesh->nmax/8;ia<7*mesh->nmax/8;ia++){
        	 ll=ia-6*mesh->nmax/8;
         	 mesh->w[ia]=mesh->w[ll];
         	 mesh->s[ia][0]=-mesh->s[ll][0];
          	 mesh->s[ia][1]=-mesh->s[ll][1];
			}

		//7th quadrant
			for(int ia=7*mesh->nmax/8;ia<8*mesh->nmax/8;ia++){
        	 ll=ia-7*mesh->nmax/8;
         	 mesh->w[ia]=mesh->w[ll];
         	 mesh->s[ia][0]=mesh->s[ll][0];
          	 mesh->s[ia][1]=-mesh->s[ll][1];
			}

		//compute eta w.r.t 3D
			for(int ll=0;ll<4*mesh->nmax/8;ll++) mesh->s[ll][2]=sqrt(1. - mesh->s[ll][0]*mesh->s[ll][0]- mesh->s[ll][1]*mesh->s[ll][1]);
			for(int ll=4*mesh->nmax/8;ll<mesh->nmax;ll++) mesh->s[ll][2]=-sqrt(1. - mesh->s[ll][0]*mesh->s[ll][0]- mesh->s[ll][1]*mesh->s[ll][1]);

      }
/*********************end::3D******************************************/
/*
	from here we compute the scattering phase function that comes into
	the source term in radiative transfer equation. The Henyey-Greenstein
	phase function commonly used in tissue optics is applied here too.
*/	

	vector <double> normFactor(mesh->nmax,0.0);
		
	/* normalize phase function (begin) */
	for(int kb=0;kb<mesh->nmax;kb++) {
		normFactor[kb]=0.0;
		for(int ka=0;ka<mesh->nmax;ka++) {
			normFactor[kb] += sfn(mesh->s[kb],mesh->s[ka],user->g)*(mesh->w[kb])/4.0/onepi;	
		}
		normFactor[kb]=1.0/normFactor[kb];
		//normFactor[kb]=1.0;
	}
	
	for(int kb=0;kb<mesh->nmax;kb++) {
		for(int ka=0;ka<mesh->nmax;ka++) {
			mesh->scatphase[ka][kb]=normFactor[kb]*sfn(mesh->s[kb],mesh->s[ka],user->g);	
		}
	}
	/* normalize phase function (end) */	

}
