/*****************************************************************************

	module:	sparseMatrix.cpp

	function that constructs the sparse matrix A to setup the linear
	Ax=b which is solved iteratively by a GMRES solver. All nonezero entries
	in the sparse matrix A will be stored in the one-dimension row vector using
	a row-compress technique.

******************************************************************************/

void sparseMatrixEx(double *a, int *rowIndex, int *columnIndex,
									Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int nwave     = meas->nwave;
	int nelem		  = mesh->nelem;
	int nmax 		  = mesh->nmax;
	int maxnei		= mesh->maxnei;
	int ndetect		= meas->ndetect;
	int nsource		= meas->nsource;
  int nz_max		= user->nz_max;
	int	size 			= param->fwd_size;
	
	int nz,row,col,nei;
	int ll,kk;
	vector <int> used(mesh->maxSurfNum,0);
	vector <double> dist_center(mesh->maxSurfNum,0.0);
	vector <double> s_in(3,0.0); // exact incident direction at Fresnel interfact inside the medium
	double onepi = 2.0*asin(1.0);
	double dot,dot1,dot2;
	double small = 1.0e-25;
	double tmp;

	/* 
	sparse matrix A for SP1 model 
	*/
	///////////////////////////////////////////////
	if(user->lightModel == 1){
	///////////////////////////////////////////////
		for(int lw=0;lw<nwave;lw++){
			for(int ele=0;ele<nelem;ele++){
					param->xdif[ele] = 1.0/(3.0*(param->xka[ele]+(1.-user->g)*param->xsig[ele]));
			}

			// initialize
			nz = -1;

			for(int i=0;i<nz_max;i++){
					a[i]=0.0;
					columnIndex[i]=0;
			}
			for(int i=0;i<size+1;i++){
					rowIndex[i] = 0;
			}

			/* store only none-zero elements with row-compressed storage (RCS) format */
			for(int ele=0;ele<nelem;ele++){
					row = ele;
					col = ele;
					nz = nz + 1;
					columnIndex[nz] = col;

					// absorption by element
					a[nz] = (param->xka[ele]+param->qs[ele])*mesh->elem[ele].vol;

					// SP1 flux across the faces: SP1 strength from primary node (+)
					for(int j=0;j<mesh->elem[ele].snum;j++){
						if(mesh->elem[ele].sbound[j] == 0){
							nei = mesh->elem[ele].snei[j][0];
							a[nz] = a[nz] + 0.5 * (param->xdif[ele] + param->xdif[nei])*mesh->elem[ele].sdist[j];//
						} 
						else if(mesh->elem[ele].sbound[j] == 1){
							if(fabs(mesh->elem[ele].center[2] - param->maxZ) < 1.0e-6){
								a[nz] = a[nz] + mesh->elem[ele].sarea[j]
																		* (1.0/2.0-user->R[0])/(1.0+3.0*user->R[1]);
							}else if(mesh->elem[ele].center[2] < param->maxZ){
								a[nz] = a[nz] + mesh->elem[ele].sarea[j] * (1.0/2.0);
							}
						}
					}

					// SP1 flux across the faces: SP1 strength from neighbor nodes (-)
 					for(int j=0; j<mesh->elem[ele].snum; j++){
						if(mesh->elem[ele].sbound[j] == 0){
							nz = nz + 1;
							nei = mesh->elem[ele].snei[j][0];
							columnIndex[nz] = nei;
							a[nz] = -0.5 * (param->xdif[ele] + param->xdif[nei])*mesh->elem[ele].sdist[j];//
						}
					}

					rowIndex[row+1] = nz + 1;

			}

			//cerr << nz << " " << user->nz_max << endl;
			
			/*********************************************************/
			user->nz_max = nz + 1;
			/*********************************************************/
							
		} // end of wavelength loop

	}// end of SP1 sparse matrix A
	

	/* 
	sparse matrix A for SP3 model 
	*/
	int mua3_cent;
	int mua3_neib;
	int SPN = 2;
	double mua2;
	///////////////////////////////////////////////
	if(user->lightModel == 2){
	///////////////////////////////////////////////
		for(int lw=0;lw<nwave;lw++){
			for(int ele=0;ele<nelem;ele++){
					param->xdif[ele] = 1.0/(3.0*(param->xka[ele]+(1.0-user->g)*param->xsig[ele]));
			}

			// initialize
			nz = -1;

			for(int i=0;i<nz_max;i++){
					a[i]=0.0;
					columnIndex[i]=0;
			}
			for(int i=0;i<size+1;i++){
					rowIndex[i] = 0;
			}

			/* store only none-zero elements with row-compressed storage (RCS) format */
			for(int order=0;order<SPN;order++){
			
				for(int ele=0;ele<nelem;ele++){
				
				/* p1 */
 					if(order == 0){
						row = ele+order*nelem; 		  
						for(int od=0;od<SPN;od++){
							if(od == 0){
								col = ele+od*nelem;
								nz = nz + 1;
								columnIndex[nz] = col;

								// absorption by element
								a[nz] = param->xka[ele]*mesh->elem[ele].vol;

								// SP1 flux across the faces: SP1 strength from primary node (+)
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 0){
										nei = mesh->elem[ele].snei[j][0];
										a[nz] = a[nz] + 0.5 * (param->xdif[ele] + param->xdif[nei])*mesh->elem[ele].sdist[j];

									}else if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] + mesh->elem[ele].sarea[j]
																					* get_coeff(order,od,ele,param,user);
									}
								}

								// SP1 flux across the faces: SP1 strength from neighbor nodes (-)
								for(int j=0; j<mesh->elem[ele].snum; j++){
									if(mesh->elem[ele].sbound[j] == 0){

										nz = nz + 1;
										col = mesh->elem[ele].snei[j][0];
										columnIndex[nz] = col;

										a[nz] = -0.5 * (param->xdif[ele] + param->xdif[col])*mesh->elem[ele].sdist[j];//
									}
								}
							}else if(od == 1){
								col = ele+od*nelem;         	    
								nz=nz+1;
								columnIndex[nz]=col;	
								a[nz] = -(2.0/3.0)*param->xka[ele]*mesh->elem[ele].vol;  	    
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] + get_coeff(order,od,ele,param,user) * mesh->elem[ele].sarea[j];
									}
								}						
							}
						}				
						rowIndex[row+1] = nz + 1;
					}
				
					if(order == 1){
						row = ele+order*nelem; 	
						for(int od=0;od<SPN;od++){
							if(od == 1){
								col = ele+od*nelem;
								nz = nz + 1;
								columnIndex[nz] = col;	        
								mua2=param->xka[ele]+(1.-pow(user->g,2))*param->xsig[ele];
								a[nz] = ((4./9.)*param->xka[ele]+(5./9.)*mua2)*mesh->elem[ele].vol;        
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 0){
										nei = mesh->elem[ele].snei[j][0];
										mua3_cent=1./(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele])/7.;
										mua3_neib=1./(param->xka[nei]+(1.-pow(user->g,3))*param->xsig[nei])/7.;
										a[nz] = a[nz] + 0.5 * (mua3_cent+mua3_neib) * mesh->elem[ele].sdist[j];
									} else if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] + get_coeff(order,od,ele,param,user)*mesh->elem[ele].sarea[j];
									}
								}
					
								// neigbor nodes

								for(int j=0; j<mesh->elem[ele].snum; j++){
									if(mesh->elem[ele].sbound[j] == 0){
			
										nz = nz + 1;
										nei = mesh->elem[ele].snei[j][0];
										col = nei+od*nelem;
										columnIndex[nz] = col;
										mua3_cent=1./(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele])/7.;
										mua3_neib=1./(param->xka[nei]+(1.-pow(user->g,3))*param->xsig[nei])/7.;  			
										a[nz] = -0.5 * (mua3_cent + mua3_neib) * mesh->elem[ele].sdist[j];
									}
								}    		        
							}else if(od == 0){
								col = ele+od*nelem;
								nz = nz + 1;
								columnIndex[nz] = col;	        
								a[nz] = -(2./3.)*param->xka[ele]*mesh->elem[ele].vol; 
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] +get_coeff(order,od,ele,param,user)*mesh->elem[ele].sarea[j];
									}
								}      		
							}
						}				
						rowIndex[row+1] = nz + 1;
					}
				
				}
			}
			//cerr << nz << " " << user->nz_max << endl;
			
			/*********************************************************/
			user->nz_max = nz + 1;
			/*********************************************************/
							
		} // end of wavelength loop

	}// end of SP3 sparse matrix A

	/* 
	sparse matrix A for radiative transfer model (RTE) 
	*/
	///////////////////////////////////////////////
	if(user->lightModel == 3){
	///////////////////////////////////////////////
		for(int lw=0;lw<nwave;lw++){

			// initialize
			nz = -1;

			for (int i=0;i<nz_max;i++){
				a[i]=0.0;
				columnIndex[i]=0;
			}
			for(int i=0;i<size+1;i++){
				rowIndex[i] = 0;
			}

			/* store only none-zero elements with row-compressed storage (RCS) format */
			for(int ka=0;ka<nmax;ka++){
				for(int ele=0;ele<nelem;ele++){
				
					row = nelem * ka + ele;
				
					for(int kb=0;kb<nmax;kb++){
					
						col = nelem * kb + ele;
						nz = nz + 1;
						columnIndex[nz] = col;	
					
						if(kb == ka){
					
							a[nz] = -param->xsig[ele] * mesh->scatphase[ka][ka] 
													* mesh->elem[ele].vol * mesh->w[ka]/4.0/onepi
													+ (param->xka[ele] + param->qs[ele] + param->xsig[ele]) 
													* mesh->elem[ele].vol; 
					
							for(int j=0;j<mesh->elem[ele].snum;j++){
								if(mesh->elem[ele].sbound[j] == 0){
					
									dot = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
									a[nz] = a[nz] + max(0.0,dot) 
																							* mesh->elem[ele].sarea[j];
					
								}else if(mesh->elem[ele].sbound[j] == 1){
					
									dot = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
									a[nz] = a[nz] + max(0.0,dot)
																				      * mesh->elem[ele].sarea[j];
								}
							}
						}	
						else if(kb != ka){
					
							a[nz] = -param->xsig[ele] * mesh->scatphase[ka][kb]
													* mesh->elem[ele].vol * mesh->w[kb]/4.0/onepi;
							 
							//diffuse reflection
						
							/* 
							for(int j=0;j<mesh->elem[ele].snum;j++){
								if(mesh->elem[ele].sbound[j] == 1){
										dot1 = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
										dot2 = vecDot(mesh->s[kb],mesh->elem[ele].snorm[j]);
										if(dot2>0.0 && dot1<0.0){
											if(mesh->elem[ele].center[2] > param->minZ){
												a[nz] = a[nz] 
															+ (user->R_eff/onepi) 
															* dot2 * mesh->w[kb] * dot1 * mesh->elem[ele].sarea[j];
											}else{
												a[nz] = a[nz] 
															+ (1.0/onepi) 
															* dot2 * mesh->w[kb] * dot1 * mesh->elem[ele].sarea[j];	
											}								
										}				
								}
							}
 							*/
								
							
							
							//specular reflection
							//*
							for(int j=0;j<mesh->elem[ele].snum;j++){
								if(mesh->elem[ele].sbound[j] == 1){
									if(fabs(mesh->elem[ele].center[2] - param->maxZ) < 1.0e-6){
										dot1=vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
										if(dot1 < 0.0) {
			
												s_in = mesh->s[ka]-2.0 * mesh->elem[ele].snorm[j] * dot1;
												s_in = s_in / vecNorm(s_in);
												/* 
												check the correlation between 
												exact incident direction(s_in) and discrete ordinates (s_kb)
												as implemented below: vecDot(s_in,s_kb)). 
												When the correlation is greater than 0.9 
												(i.e., s_kb is close enough to s_in), the s_kb is taken as
												incidient direction.
												*/
												dot2 = vecDot(mesh->s[kb],mesh->elem[ele].snorm[j]);
												if(dot2 > 0.0 && fabs(vecDot(s_in,mesh->s[kb]) - 1.0) < 1.0e-6){
													/* 
													testing partial reflection at Fresnel interface
													*/
													/*
													if(ele == 0 && j == 0){
														cout << "kb = " << kb << endl;
														cout << "ele = " << ele << endl;
														cout << "surf = " << j << endl;
														cout << "s_norm: " << mesh->elem[ele].snorm[j][0] << " " << mesh->elem[ele].snorm[j][1] << " " << mesh->elem[ele].snorm[j][2] << endl; 
														cout << "s_ka: " << mesh->s[ka][0] << " " << mesh->s[ka][1] << " " << mesh->s[ka][2] << endl;
														cout << "s_in: " << s_in[0] << " " << s_in[1] << " " << s_in[2] << endl;
														cout << "s_kb: " << mesh->s[kb][0] << " " << mesh->s[kb][1] << " " << mesh->s[kb][2] << endl;
														cout << "error: " << fabs(vecDot(s_in,mesh->s[kb]) - 1.0) << endl;
														cout << endl;
												  }
												  */									
													a[nz] = a[nz] 
														+ (1.0-max(0.0,dot1/fabs(dot1+small)))
														* max(0.0,dot2/fabs(dot2+small))
												  	* dir_reflect(mesh->elem[ele].snorm[j],mesh->s[kb],user->nindex)
														* dot1 * mesh->elem[ele].sarea[j];
												}
										}	
									}
								}
							}
							//*/
							
						}
					}
					//neigbor nodes
					for(int j=0;j<mesh->elem[ele].snum;j++){
						if(mesh->elem[ele].sbound[j] == 0){
							nz = nz + 1;
							col = nelem * ka + mesh->elem[ele].snei[j][0];
							columnIndex[nz] = col;
							dot = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);				
							a[nz] = dot * max(0.0,-dot/fabs(dot+small)) * mesh->elem[ele].sarea[j];
						}
					}

					rowIndex[row+1] = nz + 1;
					
				}
			/*********************************************************/
			// set up sparse matrix with a row compress scheme (end)
			}		

			//cerr << nz << " " << user->nz_max << endl;
			
			/*********************************************************/
			//user->nz_max = nz + 1;
			/*********************************************************/
		} // end of wavelength loop
	
	}// end of RTE sparse matrix A

	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
}


void sparseMatrixEm(double *a, int *rowIndex, int *columnIndex,
									Mesh *mesh, Meas *meas, Param *param, User *user)
{
	int nwave     = meas->nwave;
	int nelem		  = mesh->nelem;
	int nmax 		  = mesh->nmax;
	int maxnei		= mesh->maxnei;
	int ndetect		= meas->ndetect;
	int nsource		= meas->nsource;
  int nz_max		= user->nz_max;
	int	size 			= param->fwd_size;
	
	int nz,row,col,nei;
	int ll,kk;
	vector <int> used(mesh->maxSurfNum,0);
	vector <double> dist_center(mesh->maxSurfNum,0.0);
	vector <double> s_in(3,0.0); // exact incident direction at Fresnel interfact inside the medium
	double onepi = 2.0*asin(1.0);
	double dot,dot1,dot2;
	double small = 1.0e-25;
	double tmp;

	/* 
	sparse matrix A for SP1 model 
	*/
	///////////////////////////////////////////////
	if(user->lightModel == 1){
	///////////////////////////////////////////////
		for(int lw=0;lw<nwave;lw++){
			for(int ele=0;ele<nelem;ele++){
					param->xdif[ele] = 1.0/(3.0*(param->xka[ele]+(1.-user->g)*param->xsig[ele]));
			}

			// initialize
			nz = -1;

			for(int i=0;i<nz_max;i++){
					a[i]=0.0;
					columnIndex[i]=0;
			}
			for(int i=0;i<size+1;i++){
					rowIndex[i] = 0;
			}

			/* store only none-zero elements with row-compressed storage (RCS) format */
			for(int ele=0;ele<nelem;ele++){
					row = ele;
					col = ele;
					nz = nz + 1;
					columnIndex[nz] = col;

					// absorption by element
					a[nz] = param->xka[ele]*mesh->elem[ele].vol;

					// SP1 flux across the faces: SP1 strength from primary node (+)
					for(int j=0;j<mesh->elem[ele].snum;j++){
						if(mesh->elem[ele].sbound[j] == 0){
							nei = mesh->elem[ele].snei[j][0];
							a[nz] = a[nz] + 0.5 * (param->xdif[ele] + param->xdif[nei])*mesh->elem[ele].sdist[j];//
						}else if(mesh->elem[ele].sbound[j] == 1){
							if(mesh->elem[ele].center[2] > param->minZ){
								a[nz] = a[nz] + mesh->elem[ele].sarea[j]
																		* (1.0/2.0-user->R[0])/(1.0+3.0*user->R[1]);
							}else if(fabs(mesh->elem[ele].center[2] - param->minZ) < 1.0e-6){
								a[nz] = a[nz] + mesh->elem[ele].sarea[j] * (1.0/(2.0*1.0e25));
							}
						}
						//*/
					}

					// SP1 flux across the faces: SP1 strength from neighbor nodes (-)
 					for(int j=0; j<mesh->elem[ele].snum; j++){
						if(mesh->elem[ele].sbound[j] == 0){
							nz = nz + 1;
							nei = mesh->elem[ele].snei[j][0];
							columnIndex[nz] = nei;
							a[nz] = -0.5 * (param->xdif[ele] + param->xdif[nei])*mesh->elem[ele].sdist[j];//
						}
					}

					rowIndex[row+1] = nz + 1;

			}

			//cerr << nz << " " << user->nz_max << endl;
			
			/*********************************************************/
			user->nz_max = nz + 1;
			/*********************************************************/
							
		} // end of wavelength loop

	}// end of SP1 sparse matrix A
	

	/* 
	sparse matrix A for SP3 model 
	*/
	int mua3_cent;
	int mua3_neib;
	int SPN = 2;
	double mua2;
	///////////////////////////////////////////////
	if(user->lightModel == 2){
	///////////////////////////////////////////////
		for(int lw=0;lw<nwave;lw++){
			for(int ele=0;ele<nelem;ele++){
					param->xdif[ele] = 1.0/(3.0*(param->xka[ele]+(1.0-user->g)*param->xsig[ele]));
			}

			// initialize
			nz = -1;

			for(int i=0;i<nz_max;i++){
					a[i]=0.0;
					columnIndex[i]=0;
			}
			for(int i=0;i<size+1;i++){
					rowIndex[i] = 0;
			}

			/* store only none-zero elements with row-compressed storage (RCS) format */
			for(int order=0;order<SPN;order++){
			
				for(int ele=0;ele<nelem;ele++){
				
				/* p1 */
 					if(order == 0){
						row = ele+order*nelem; 		  
						for(int od=0;od<SPN;od++){
							if(od == 0){
								col = ele+od*nelem;
								nz = nz + 1;
								columnIndex[nz] = col;

								// absorption by element
								a[nz] = param->xka[ele]*mesh->elem[ele].vol;

								// SP1 flux across the faces: SP1 strength from primary node (+)
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 0){
										nei = mesh->elem[ele].snei[j][0];
										a[nz] = a[nz] + 0.5 * (param->xdif[ele] + param->xdif[nei])*mesh->elem[ele].sdist[j];

									}else if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] + mesh->elem[ele].sarea[j]
																					* get_coeff(order,od,ele,param,user);
									}
								}

								// SP1 flux across the faces: SP1 strength from neighbor nodes (-)
								for(int j=0; j<mesh->elem[ele].snum; j++){
									if(mesh->elem[ele].sbound[j] == 0){

										nz = nz + 1;
										col = mesh->elem[ele].snei[j][0];
										columnIndex[nz] = col;

										a[nz] = -0.5 * (param->xdif[ele] + param->xdif[col])*mesh->elem[ele].sdist[j];//
									}
								}
							}else if(od == 1){
								col = ele+od*nelem;         	    
								nz=nz+1;
								columnIndex[nz]=col;	
								a[nz] = -(2.0/3.0)*param->xka[ele]*mesh->elem[ele].vol;  	    
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] + get_coeff(order,od,ele,param,user) * mesh->elem[ele].sarea[j];
									}
								}						
							}
						}				
						rowIndex[row+1] = nz + 1;
					}
				
					if(order == 1){
						row = ele+order*nelem; 	
						for(int od=0;od<SPN;od++){
							if(od == 1){
								col = ele+od*nelem;
								nz = nz + 1;
								columnIndex[nz] = col;	        
								mua2=param->xka[ele]+(1.-pow(user->g,2))*param->xsig[ele];
								a[nz] = ((4./9.)*param->xka[ele]+(5./9.)*mua2)*mesh->elem[ele].vol;        
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 0){
										nei = mesh->elem[ele].snei[j][0];
										mua3_cent=1./(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele])/7.;
										mua3_neib=1./(param->xka[nei]+(1.-pow(user->g,3))*param->xsig[nei])/7.;
										a[nz] = a[nz] + 0.5 * (mua3_cent+mua3_neib) * mesh->elem[ele].sdist[j];
														
									} else if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] + get_coeff(order,od,ele,param,user)*mesh->elem[ele].sarea[j];
									}
								}
					
								// neigbor nodes

								for(int j=0; j<mesh->elem[ele].snum; j++){
									if(mesh->elem[ele].sbound[j] == 0){
			
										nz = nz + 1;
										nei = mesh->elem[ele].snei[j][0];
										col = nei+od*nelem;
										columnIndex[nz] = col;
										mua3_cent=1./(param->xka[ele]+(1.-pow(user->g,3))*param->xsig[ele])/7.;
										mua3_neib=1./(param->xka[nei]+(1.-pow(user->g,3))*param->xsig[nei])/7.;  			
										a[nz] = -0.5 * (mua3_cent + mua3_neib) * mesh->elem[ele].sdist[j];
									}
								}    		        
							}else if(od == 0){
								col = ele+od*nelem;
								nz = nz + 1;
								columnIndex[nz] = col;	        
								a[nz] = -(2./3.)*param->xka[ele]*mesh->elem[ele].vol; 
								for(int j=0;j<mesh->elem[ele].snum;j++){
									if(mesh->elem[ele].sbound[j] == 1){
										a[nz] = a[nz] +get_coeff(order,od,ele,param,user)*mesh->elem[ele].sarea[j];
									}
								}      		
							}
						}				
						rowIndex[row+1] = nz + 1;
					}
				
				}
			}
			//cerr << nz << " " << user->nz_max << endl;
			
			/*********************************************************/
			user->nz_max = nz + 1;
			/*********************************************************/
							
		} // end of wavelength loop

	}// end of SP3 sparse matrix A

	/* 
	sparse matrix A for radiative transfer model (RTE) 
	*/
	///////////////////////////////////////////////
	if(user->lightModel == 3){
	///////////////////////////////////////////////
		for(int lw=0;lw<nwave;lw++){

			// initialize
			nz = -1;

			for (int i=0;i<nz_max;i++){
				a[i]=0.0;
				columnIndex[i]=0;
			}
			for(int i=0;i<size+1;i++){
				rowIndex[i] = 0;
			}

			/* store only none-zero elements with row-compressed storage (RCS) format */
			for(int ka=0;ka<nmax;ka++){
				for(int ele=0;ele<nelem;ele++){
				
					row = nelem * ka + ele;
				
					for(int kb=0;kb<nmax;kb++){
					
						col = nelem * kb + ele;
						nz = nz + 1;
						columnIndex[nz] = col;	
					
						if(kb == ka){
					
							a[nz] = -param->xsig[ele] * mesh->scatphase[ka][ka] 
													* mesh->elem[ele].vol * mesh->w[ka]/4.0/onepi
													+ (param->xka[ele] + param->qs[ele] + param->xsig[ele]) 
													* mesh->elem[ele].vol; 
					
							for(int j=0;j<mesh->elem[ele].snum;j++){
								if(mesh->elem[ele].sbound[j] == 0){
					
									dot = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
									a[nz] = a[nz] + max(0.0,dot) 
																							* mesh->elem[ele].sarea[j];
					
								}else if(mesh->elem[ele].sbound[j] == 1){
					
									dot = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
									a[nz] = a[nz] + max(0.0,dot)
																				      * mesh->elem[ele].sarea[j];
								}
							}
						}	
						else if(kb != ka){
					
							a[nz] = -param->xsig[ele] * mesh->scatphase[ka][kb]
													* mesh->elem[ele].vol * mesh->w[kb]/4.0/onepi;
							 
							//diffuse reflection

/* 
							for(int j=0;j<mesh->elem[ele].snum;j++){
								if(mesh->elem[ele].sbound[j] == 1){
										dot1 = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
										dot2 = vecDot(mesh->s[kb],mesh->elem[ele].snorm[j]);
										if(dot2>0.0 && dot1<0.0){
											if(mesh->elem[ele].center[2] > param->minZ){
												a[nz] = a[nz] 
															+ (user->R_eff/onepi) 
															* dot2 * mesh->w[kb] * dot1 * mesh->elem[ele].sarea[j];
											}else{
												a[nz] = a[nz] 
															+ (1.0/onepi) 
															* dot2 * mesh->w[kb] * dot1 * mesh->elem[ele].sarea[j];	
											}								
										}				
								}
							}
 */

							
							
							//specular reflection
							
							for(int j=0;j<mesh->elem[ele].snum;j++){
								if(mesh->elem[ele].sbound[j] == 1){
									if(fabs(mesh->elem[ele].center[2] - param->maxZ)<1.0e-6){
										dot1=vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);
										if(dot1 < 0.0) {
			
												s_in = mesh->s[ka]-2.0 * mesh->elem[ele].snorm[j] * dot1;
												s_in = s_in / vecNorm(s_in);
												/* 
												check the correlation between 
												exact incident direction(s_in) and discrete ordinates (s_kb)
												as implemented below: vecDot(s_in,s_kb)). 
												When the correlation is greater than 0.9 
												(i.e., s_kb is close enough to s_in), the s_kb is taken as
												incidient direction.
												*/
												dot2 = vecDot(mesh->s[kb],mesh->elem[ele].snorm[j]);
												if(dot2 > 0.0 && fabs(vecDot(s_in,mesh->s[kb]) - 1.0) < 1.0e-6){
													/* 
													testing partial reflection at Fresnel interface
													*/
													/*
													if(ele == 0 && j == 0){
														cout << "kb = " << kb << endl;
														cout << "ele = " << ele << endl;
														cout << "surf = " << j << endl;
														cout << "s_norm: " << mesh->elem[ele].snorm[j][0] << " " << mesh->elem[ele].snorm[j][1] << " " << mesh->elem[ele].snorm[j][2] << endl; 
														cout << "s_ka: " << mesh->s[ka][0] << " " << mesh->s[ka][1] << " " << mesh->s[ka][2] << endl;
														cout << "s_in: " << s_in[0] << " " << s_in[1] << " " << s_in[2] << endl;
														cout << "s_kb: " << mesh->s[kb][0] << " " << mesh->s[kb][1] << " " << mesh->s[kb][2] << endl;
														cout << "error: " << fabs(vecDot(s_in,mesh->s[kb]) - 1.0) << endl;
														cout << endl;
												  }
												  */									
													a[nz] = a[nz] 
														+ (1.0-max(0.0,dot1/fabs(dot1+small)))
														* max(0.0,dot2/fabs(dot2+small))
												  	* dir_reflect(mesh->elem[ele].snorm[j],mesh->s[kb],user->nindex)
														* dot1 * mesh->elem[ele].sarea[j];
												}
										}											
									}
								}
							}
							
							
						}
					}
					//neigbor nodes
					for(int j=0;j<mesh->elem[ele].snum;j++){
						if(mesh->elem[ele].sbound[j] == 0){
							nz = nz + 1;
							col = nelem * ka + mesh->elem[ele].snei[j][0];
							columnIndex[nz] = col;
							dot = vecDot(mesh->s[ka],mesh->elem[ele].snorm[j]);				
							a[nz] = dot * max(0.0,-dot/fabs(dot+small)) * mesh->elem[ele].sarea[j];
						}
					}

					rowIndex[row+1] = nz + 1;
					
				}
			/*********************************************************/
			// set up sparse matrix with a row compress scheme (end)
			}		

			//cerr << nz << " " << user->nz_max << endl;
			
			/*********************************************************/
			//user->nz_max = nz + 1;
			/*********************************************************/
		} // end of wavelength loop
	
	}// end of RTE sparse matrix A

	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
}