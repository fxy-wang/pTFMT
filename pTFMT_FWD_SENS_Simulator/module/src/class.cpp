/* Element class */
class Element
{
	public:
  	vector <double> center;
  	double 	vol;
  	double	dist2srcdet;
  	int    	snum;
  	vector <double> sarea;
  	vector <double> sdist;
  	vector <vector <double> > snorm;
    vector <double> weight;
    int			neibnum;
  	vector <int> neib;
  	vector <int> sbound;
  	vector <vector <int> > snei;
  	int	 		layer;
  	int			bnode;
  	Element();
  	~Element();
};
Element::Element(){
  	center.clear();
  	sarea.clear();
  	sdist.clear();
  	snorm.clear();
    weight.clear();
  	neib.clear();
  	sbound.clear();
  	snei.clear();
}
Element::~Element(){}

/* Mesh class */
class Mesh
{
	public:
  	int		mshelem; // number of finite elements
  	int		vertnum; // number of vertices for each mesh element
  	int		nelem;
  	int		nodeNumROI;
  	int		maxSurfNum;
  	int		maxnei;
  	int		maxsmnei; //max. number of nodal points used in the regularization
  	int 	ns;
  	int		nmax;
  	int		bnode_num;
  	vector <int>		bnode_index;
  	vector <vector <int>>		elm;
  	string		elemtype;
  	vector <Element> elem;
  	double  avg_vol;
  	vector <double> w;
  	vector <vector <double> > s; // discrete ordinates
  	vector <vector <double> >	scatphase;
  	vector <double> center;
  	Mesh();
  	~Mesh();
};

Mesh::Mesh(){
	  bnode_index.clear();
  	elm.clear();
  	elem.clear();
  	w.clear();
  	s.clear();
  	scatphase.clear();
  	center.clear();

}
Mesh::~Mesh(){}

/* Meas class */
class Meas
{
	public:
		int			nfreq;
		int     nwave;
		int 		nsource;
		int 		ndetect;
		int			NumOfUsedMeasurements;
		int			ntime;

		vector<int> filter;		
		vector<double>	Ls; // Laplace variable s
		
		vector <int>	lsource;
		vector <int>	ldetector;
		vector <int>	nodeldetector;
		vector <int>	nsrcpersrc; //number of source points used for each source illumination
		vector <int>	ndetpersrc; //number of detector points for each source illumination

		vector <vector <int> >	detpairs; // index of source-detector pair 
		vector <vector <int> >	srcpairs; // index of source points per each source illumination
		
 		Vec4d amplitude;
		Vec4d phase;
		Vec4d weight;
		Vec4d weight_m;
			
		Vec4d prediction;
		Vec4d prediction_m;
					
		vector <double> freq;
		vector <double> wave;
		vector <vector <double> > source;
		vector <vector <double> > detector;

		double qin; // heat flux input to source
	
		Meas();
		~Meas();
};

Meas::Meas(){
		lsource.clear();
		ldetector.clear();
		nsrcpersrc.clear(); //number of source points for each source illumination
		ndetpersrc.clear(); //number of measurements for each source illumination
		detpairs.clear(); // detectors associated with each source illumination 
		srcpairs.clear(); // source points associated with each source illumination
		amplitude.clear();
		filter.clear();
		phase.clear();
		prediction.clear();
		weight.clear();
		source.clear();
		detector.clear();
		freq.clear();
		wave.clear();
}
Meas::~Meas(){}

/* parameters */
class Param
{
  public:
		int		imax,jmax,kmax;
		int		fwd_size;
		int		inv_size;	
		double  dx,dy,dz;
		double  Lx,Ly,Lz;
		double  minX,minY,minZ;
		double  maxX,maxY,maxZ;
		vector <double> xka,xsig,xdif; // absorption,scattering,diffusion coefficients

		double  prefactor;
		vector <double> cdtdv;
		vector <double> cdv;

		vector <double> qs;
		vector <double> lifeTime;
		
		Vec4d xis; // forward solution with excitation
		Vec4d xism; // forward solution with emission
		Vec4d axi; // forward solution with excitation
		Vec4d axim; // forward solution with emission
		Vec4d excitation;
		Vec4d adjdet; // adjoint solution with source at detector location
		Vec5d S; // coefficients of sensing matrix
		
		Param();
		~Param();	
};
Param::Param(){
		xka.clear();
		xsig.clear();
		xdif.clear(); 	
		cdtdv.clear();
		cdv.clear();				
		qs.clear();
		adjdet.clear();
		xis.clear();
		xism.clear();
		axi.clear();
		axim.clear();
		excitation.clear();
		S.clear();
}
Param::~Param(){}

/* user defined variables and functions */
class User
{
	public:
  int 	prob; //type of problem to be solved
  int 	lightModel; // light propagation model
  int		dim; //dimensionality of the problem
  int		nz_max; //maximum number of nonzero elements in sparse matrix
  int		fwd_monitor;
  int		preconditioner;
  int		fwd_itmax; // maximum number of forward iteration	
	int		numDevs,numGPUDevs,numCPUDevs;
		
	double	lifeTime;
	double	qYield;
	double 	qs_b;
  double	xka_b,xsig_b;
  
  double	gel_layer_thickness;
  double	gel_layer_mua;
  double	gel_layer_mus;
  
  double	qs_min,qs_max;
  
  int		nqs_o;
  vector<double> qs_o;
  vector<double> qs_o_radius;
	vector<double> qs_center_max;
  vector<vector<double>> qs_o_location;

  double	qs_d,qs_c;
  
  double  SNRAmp;
  double  SNRArg;

  double	speed; // speed of light in the medium
  double	g; // anistropic factor
  double 	nindex;
  double	f; // f(0:explicit,1:implicit)
  double	R[6], R_eff; // medium reflectivity
  double	A1,B1,C1,D1,A2,B2,C2,D2; // coefficients for SP3
  double	nu1,nu2; // measurement operator coefficients for sp3 partial current
  double	alpha; //step length
  double	fwd_tol; // SOR solver tolerance
  double	noise; // noise level
	double	max_run_time; // max time limit for inverse problems

  string	input_ini;
  string	src_det_cfg; // file for source-detector information
  string	fwd_optical_properties; // file for exact distribution of xka and xsig

	string	filenameOnly;
  string	detector_readings;// file for prediction
  string	fwd_fluence; // file for predicted detector signal
  string	sens_matrix;
  string	progress_monitor; // file for residual monitor
  string	finite_volume_mesh; // fvm mesh
  time_t 	begintime,endtime;
  
  User();
  ~User();
};

User::User(){
  qs_o.clear();
  qs_o_location.clear();
  qs_o_radius.clear();
}

User::~User(){}

