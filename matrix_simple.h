// very simple matrix class based on an extension of the Vec template

#include "vec.h"

class matrix{
 public:
  int xdim,ydim;
  dmatrix mat;

  // this function assumes that xdim and ydim have reasonable values
  void init_mem(void){
    int i;
    mat.resize(ydim);
    for(i=0;i<ydim;i++) mat[i].resize(xdim);
  }
  void copy(const matrix& a){
    xdim=a.xdim; ydim=a.ydim;
    mat=a.mat; // note that "operator=" exists for dmatrix=Vec<dvec>
  }
  void test(const matrix& a){
    if((xdim!=a.xdim) || (ydim!=a.ydim)) 
      error("bad matrix dimensions");
  }
  int x(void){ return xdim;}
  int y(void){ return ydim;}
  int x(void) const{ return xdim;}
  int y(void) const{ return ydim;}
  // a very simple default constructor
  matrix(void){
    ydim=mat.dim;
    if(ydim>0){
      xdim=mat[0].dim;
    }else{
      xdim=0;
    }
  }
  // the main constructor
  matrix(int y,int x,double diag=0){
    int i,j;
    xdim=x; ydim=y;
    init_mem();
    for(i=0;i<ydim;i++) for(j=0;j<xdim;j++) mat[i][j]=0;
    if(diag!=0){
      int n=(xdim<ydim?xdim:ydim);
      for(i=0;i<n;i++) mat[i][i]=diag;
    }
  }
  ~matrix(){}
  // copy constructor (same as default)
  matrix(const matrix &a){
    copy(a);
  }
  // here only a trivial destructor is necessary
  matrix& operator=(const matrix& a){
    copy(a);
  }

  // access to matrix elements by a(i,j)
  matrix& operator=(double);
  double& operator()(int i,int j){
#ifdef RANGE_CHECK
    if ((i<0) || (i>=ydim) || (j<0) || (j>=xdim)){
      printf("i = %5d   j = %5d\n",i,j);
      error("Range error!");
    }
#endif
    return mat[i][j];
  }
  double& operator()(int i,int j) const{
#ifdef RANGE_CHECK
    if ((i<0) || (i>=ydim) || (j<0) || (j>=xdim)){
      printf("i = %5d   j = %5d\n",i,j);
      error("Range error!");
    }
#endif
    return mat[i][j];
  }

  // access to matrix lines by a[i] => access to elements by a[i][j]
  dvec& operator[](int i){
    return mat[i];
  }
  dvec& operator[](int i) const{
    return mat[i];
  }

  matrix& operator+=(const matrix&);
  matrix& operator+=(double);
  matrix& operator-=(const matrix&);
  matrix& operator-=(double);

  void print_size(const char *mess){
    printf("%s = %d x %d matrix\n",mess,ydim,xdim);
    fflush(stdout);
  }
};


// computing operations

matrix& matrix::operator+=(const matrix& a){
  test(a);
  int i,j;
  for(i=0;i<ydim;i++) for(j=0;j<xdim;j++) mat[i][j]+=a.mat[i][j];
  return *this;
}

matrix& matrix::operator+=(double lam){
  int i;

  if(xdim!=ydim) error("not a square matrix");
  for(i=0;i<xdim;i++) mat[i][i]+=lam;
  return *this;
}

matrix& matrix::operator-=(const matrix& a){
  test(a);
  int i,j;
  for(i=0;i<ydim;i++) for(j=0;j<xdim;j++) mat[i][j]-=a.mat[i][j];
  return *this;
}

matrix& matrix::operator-=(double lam){
  int i;

  if(xdim!=ydim) error("not a square matrix");
  for(i=0;i<xdim;i++) mat[i][i]-=lam;
  return *this;
}

// from here on non-class functions

matrix operator*(const matrix& a,const matrix& b){
  int i,j,k;
  double sum;
  if(a.x()!=b.y()) 
    error("non fitting matrix dimensions for matrix multiplication");
  matrix p(a.y(),b.x());
  for(i=0;i<a.y();i++) for(j=0;j<b.x();j++){
    sum=0;
    for(k=0;k<a.x();k++) sum+=(a(i,k)*b(k,j));
    p(i,j)=sum;
  }
  return p;
}

matrix operator*(double a,matrix p){
  int i,j;

  for(i=0;i<p.y();i++) for(j=0;j<p.x();j++) p(i,j)=(a*p(i,j));
  return p;
}

void operator*=(matrix& a,double lam){
  int i,j;
  
  for(i=0;i<a.y();i++) for(j=0;j<a.x();j++) a(i,j)*=lam;
}

void operator/=(matrix& a,double lam){
  int i,j;
  
  for(i=0;i<a.y();i++) for(j=0;j<a.x();j++) a(i,j)/=lam;
}

// computes the solution of the linear system a*z=d 
// or equivalently:  z=a^{-1}*d  corresponding to:  z=d/a

matrix operator/(matrix d,matrix a){
  int n,m,i,j,k,temp;
  double c;
  double eps_gleich=1e-14;

  if((a.x()!=a.y()) || (a.x()!=d.y())) 
    error("Wrong matrix dimension in linear system!");
  n=a.x(); m=d.x();
  ivec p(n);
  for(i=0;i<n;i++) p[i]=i;

  for(j=0;j<n;j++){
    k=j;
    for(i=j+1;i<n;i++)
      if(fabs(a(p[i],j))>fabs(a(p[k],j))) k=i;
    if(k>j){ temp=p[k]; p[k]=p[j]; p[j]=temp; }
    if(fabs(a(p[j],j))<=0.0) error("singular matrix!");
    for(i=j+1;i<n;i++){
      c=a(p[i],j)/a(p[j],j);
      if(c!=0){
	for(k=j+1;k<n;k++) a(p[i],k)-=(c*a(p[j],k));
	for(k=0;k<m;k++) d(p[i],k)-=(c*d(p[j],k));
      }
    }
  }
  // matrix for the solution 
  matrix z(n,m);
  for(j=0;j<m;j++){
    for(i=n-1;i>=0;i--){
      for(c=0,k=i+1;k<n;k++)
	c+=(a(p[i],k)*z(k,j));
      z(i,j)=(d(p[i],j)-c)/a(p[i],i);
    }
  }
  return z;
}

matrix adjoint(const matrix& a){
  int i,j,n,m;

  n=a.y(); m=a.x();
  matrix z(m,n);
  for(i=0;i<m;i++) for(j=0;j<n;j++) z(i,j)=a(j,i);
  return z;
}

double trace(const matrix& a){
  int i,n;
  double sum;

  if(a.xdim!=a.ydim) error("not a square matrix");
  for(sum=0,i=0;i<a.xdim;i++) sum+=a(i,i);
  return sum;
}

// now some inline functions to complete the set of usual operators

inline matrix operator+(const matrix& a,const matrix& b){
  matrix c(a);
  c+=b;
  return c;
}

inline matrix operator+(const matrix& a,double b){
  matrix c(a);
  c+=b;
  return c;
}

inline matrix operator+(double a,const matrix& b){
  matrix c(b);
  c+=a;
  return c;
}

inline matrix operator-(const matrix& a,const matrix& b){
  matrix c(a);
  c-=b;
  return c;
}

inline matrix operator-(const matrix& a,double b){
  matrix c(a);
  c-=b;
  return c;
}

inline matrix operator-(double a,const matrix& b){
  matrix c(b.x(),b.x(),a);
  c-=b;
  return c;
}

inline void operator*=(matrix& a,const matrix& b){
  a=a*b;
}

inline matrix operator*(const matrix& a,double b){
  matrix c(a);
  c*=b;
  return c;
}

// corresponds to: a=b^(-1)*a
inline void operator/=(matrix& a,const matrix& b){
  a=a/b;
}

// corresponds to: a/b with b = double number
inline matrix operator/(const matrix& a,double b){
  matrix c(a);
  c/=b;
  return c;
}

// corresponds to: b*(a^(-1)) with b = double number
inline matrix operator/(double b,const matrix& a){
  matrix c(a.x(),a.x(),b);
  return c/a;
}
