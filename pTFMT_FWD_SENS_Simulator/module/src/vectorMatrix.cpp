using namespace std;

// define vector types
using Vec1d = vector <double>;
using Vec2d = vector <Vec1d>;
using Vec3d = vector <Vec2d>;
using Vec4d = vector <Vec3d>;
using Vec5d = vector <Vec4d>;

using comVec1d = vector <complex<double>>;
using comVec2d = vector <comVec1d>;
using comVec3d = vector <comVec2d>;
using comVec4d = vector <comVec3d>;
using comVec5d = vector <comVec4d>;

/*******************************************************************/
/* vector operation */
/*******************************************************************/
template <typename T>
vector <T> operator- (const vector <T> &v1){
	vector <T> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = -v1[i];
	}
	return v;
		
};

template <typename T>
vector <T> operator+ (const vector <T> &v1, const vector <T> &v2){
	if(v1.size() != v2.size())
		throw "can't add two vectors of different size!";
	vector <T> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] + v2[i];
	}
	return v;
		
};

template <typename T>
vector <T> operator- (const vector <T> &v1, const vector <T> &v2){
	if(v1.size() != v2.size())
		throw "can't subtract two vectors of different size!";
	vector <T> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] - v2[i];
	}
	return v;
		
};

template <typename T>
inline T operator* (const vector <T> &v1, const vector <T> &v2){
	if(v1.size() != v2.size())
		throw "can't multiply two vectors of different size!";
	T s;
	s=0;
	for(int i=0; i < v1.size(); i++){
		s = s + v1[i] * v2[i];
	}
	return s;
		
};

template <typename T1, typename T2>
vector <T1> operator* (const vector <T1> &v1, const T2 alpha){
	vector <T1> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] * ((T1) alpha);
	}
	return v;
		
};

template <typename T1, typename T2>
vector <T1> operator* (const T2 alpha, const vector <T1> &v1){
	vector <T1> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] * ((T1) alpha);
	}
	return v;
		
};

template <typename T1, typename T2>
vector <T1> operator/ (const vector <T1> &v1, const T2 alpha){
	vector <T1> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] / ((T1) alpha);
	}
	return v;
		
};

template <typename T1, typename T2>
vector <T1> operator/ (const T2 alpha, const vector <T1> &v1){
	vector <T1> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] / ((T1) alpha);
	}
	return v;
		
};

template <typename T1, typename T2>
vector <T1> operator/ (const vector <T1> &v1, const vector <T2> &v2){
	vector <T1> v(v1.size());
	for(int i=0; i < v1.size(); i++){
		v[i] = v1[i] / v2[i];
	}
	return v;
		
};

/* matrix*const, matrix-vector operaton */
template <typename T>
inline vector <vector <T> > transpose(const vector<vector <T> > &m1){
	int row = m1.size();
	int col = m1[0].size();
	vector <vector <T> > m(col,vector <T> (row));
	for(int i=0; i < row; i++){
		for(int j=0; j < col; j++){
			m[j][i] = m1[i][j];
		}
	}
	return m;		
};

template <typename T>
inline vector <vector <T> > operator- (const vector<vector <T> > &m1){
	vector <vector <T> > m(m1.size());
	for(int i=0; i < m1.size(); i++){
		m[i] = -m1[i];
	}
	return m;
		
};

template <typename T>
inline vector <T> operator* (const vector<vector <T> > &m, const vector <T> &v){
	if(m[0].size() != v.size())
		throw "can't multiply Matrix and vector of different size!";
	vector <T> mv(m.size());
	double init = 0.0;
	for(int i=0; i < m.size(); i++){
		mv[i] = inner_product(m[i].begin(),m[i].end(),v.begin(),init);
	}
	return mv;
		
};

template <typename T1, typename T2>
	vector <vector <T1> > operator/ (const vector<vector <T1> > &m, const T2 alpha){
	vector <vector <T1> > mv(m.size());
	for(int i=0; i < m.size(); i++){
			mv[i] = m[i] / ((T1) alpha);
	}
	return mv;
		
};

template <typename T1, typename T2>
inline	vector <vector <T1> > operator/ (const T2 alpha, const vector<vector <T1> > &m){
	vector <vector <T1> > mv(m.size());
	for(int i=0; i < m.size(); i++){
			mv[i] = m[i] / ((T1) alpha);
	}
	return mv;
		
};

template <typename T1, typename T2>
inline	vector <vector <T1> > operator* (const vector<vector <T1> > &m, const T2 alpha){
	vector <vector <T1> > mv(m.size());
	for(int i=0; i < m.size(); i++){
			mv[i] = m[i] * ((T1) alpha);
	}
	return mv;
		
};

template <typename T1, typename T2>
inline	vector <vector <T1> > operator* (const T2 alpha, const vector<vector <T1> > &m){
	vector <vector <T1> > mv(m.size());
	for(int i=0; i < m.size(); i++){
			mv[i] = m[i] * ((T1) alpha);
	}
	return mv;
			
};

/* matrix*const, matrix-vector operaton */
template <typename T>
inline vector <vector <T> > ompTranspose(const vector<vector <T> > &m1){
	int i,j;
	int row = m1.size();
	int col = m1[0].size();
	vector <vector <T> > m(col,vector <T> (row));
	# pragma omp parallel for schedule(static) shared(m1,m) private(i,j) 
	for(int i=0; i < row; i++){
		for(int j=0; j < col; j++){
			m[j][i] = m1[i][j];
		}
	}
	return m;		
};

template <typename T>
/*******************************************************************/
inline vector<T> ompMXV(int N, vector <vector<T>> M, vector <T> V)
/*******************************************************************/
{
  vector<T> MV;
	int i, j, m=M.size(),n = M[0].size();
	assert(n==V.size());
	MV.resize(m,0.0);	
	
	# pragma omp parallel for schedule(static) shared(M,V,MV,m,n) private(i,j) num_threads(N) 
		for ( i = 0; i < m; i++ )
			for ( j = 0; j < n; j++ )
					MV[i] += M[i][j] * V[j];
					
  return MV;
};

template <typename T>
/*******************************************************************/
inline double ompVXV(int N, vector <T> v1, vector <T> v2)
/*******************************************************************/
{
  T s;
	int i, n = v1.size();
	assert(n==v2.size());
	s = 0.0;
	# pragma omp parallel for reduction( + : s ) schedule(static) shared(v1,v2,n) private(i) num_threads(N) 
		for ( i = 0; i < n; i++ ) 
			s += v1[i] * v2[i];
					
  return s;
};



