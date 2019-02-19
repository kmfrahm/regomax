#ifndef MY_VEC
#define MY_VEC

// a small class for template vectors with adaptable default constructor 
// useful in omp loops as private variables
// using vector of vector one also obtains a very simple matrix class

#include <stdio.h>
#include <stdlib.h>

inline void error(const char *s){
  puts(s);
  exit(1);
}


template <class T>
class Vec{
  static int default_size;
 public:
  void test(const Vec &a) const{
    if(dim!=a.dim) error("dimension error");
  }
  int dim;
  T *c;

  // default constructor
  Vec(void){
    dim=default_size;
    if(dim>0){
      c=new T[dim];
    }else{
      c=NULL;
    }
  }

  Vec(int n){
    if(n<0) n=0;
    dim=n;
    if(dim>0){
      c=new T[dim];
    }else{
      c=NULL;
    }
    default_size=n;
  }
  // constructor with initial value
  Vec(int n,T val){
    if(n==0) n=default_size;
    if(n<0) n=0; // allow to initialize default_size with 0
    dim=n;
    if(dim>0){
      c=new T[dim];
    }else{
      c=NULL;
    }
    default_size=n;
    for(int i=0;i<dim;i++) c[i]=val;
  }
  Vec(const Vec &a){
    dim=a.dim;
    if(dim>0){
      c=new T[dim];
    }else{
      c=NULL;
    }
    for(int i=0;i<dim;i++) c[i]=a.c[i];
  }
  ~Vec(){
    // note delete[] on NULL works also (with no effect)
    delete[] c;
  }
  int& size(void){ return dim;};
  void resize(int n){
    if(n!=dim){
      dim=n;
      delete[] c;
      if(dim>0){
	c=new T[dim];
      }else{
	c=NULL;
      }
    }
  }
  void put_value(T val){
    for(int i=0;i<dim;i++) c[i]=val;
  }
  T& operator[](int i){
    return c[i];
  }
  T& operator[](int i) const{
    return c[i];
  }
  Vec& operator=(const Vec &a){
    resize(a.dim);
    for(int i=0;i<dim;i++) c[i]=a.c[i];
    return *this;
  }
  Vec& operator=(T val){
    for(int i=0;i<dim;i++) c[i]=val;
    return *this;
  }
  // simple copy function
  void copy(const Vec &a){
    resize(a.dim);
    for(int i=0;i<dim;i++) c[i]=a.c[i];
  }
  void operator+=(const Vec& a){
    test(a);
    int i;
    for(i=0;i<dim;i++) c[i]+=a.c[i];
  }
  void operator-=(const Vec& a){
    test(a);
    int i;
    for(i=0;i<dim;i++) c[i]-=a.c[i];
  }
  void operator*=(T x){
    int i;
    for(i=0;i<dim;i++) c[i]*=x;
  }
  void operator/=(T x){
    int i;
    for(i=0;i<dim;i++) c[i]/=x;
  }
  static void set_size(int n);

  void lam_add(T sp,const Vec &a){
    test(a);
    int i;
    for(i=0;i<dim;i++) c[i]+=sp*a.c[i];
  }

  void lam_diff(T sp,const Vec &a){
    test(a);
    int i;
    for(i=0;i<dim;i++) c[i]-=sp*a.c[i];
  }


  // ATTENTION: bin_write and bin_read only work correctly if T=simple
  // variable type such as int, long, double etc.
  // but not a vec<...> itself !!
  void bin_write(char const *filename){
    int i;
    FILE *fp;

    fp=fopen(filename,"w");
    fwrite(&dim,sizeof(dim),(size_t)1,fp);
    fwrite(c,sizeof(T),(size_t)dim,fp);
    fclose(fp);
  }

  // returns 1 on success and 0 on non-success
  int bin_read(char const *filename){
    int n;
    FILE *fp;

    fp=fopen(filename,"r");
    if(fp==NULL) return 0;
    if(fread(&n,sizeof(n),(size_t)1,fp)<1){
      fclose(fp);
      return 0;
    }
    resize(n);
    if(fread(c,sizeof(T),(size_t)dim,fp)<(size_t)dim){
      fclose(fp);
      return 0;
    }
    fclose(fp);
    return 1;
  }

};

template<class T>
inline T abs(T x){
  return x>=0?x:(-x);
}

// some non member functions

// fast vector swap with pointer exchange
template<class T>
inline void swap(Vec<T> &a,Vec<T> &b){
  int t;
  T *t2;
  t=a.dim; a.dim=b.dim; b.dim=t;
  t2=a.c; a.c=b.c; b.c=t2;
}

template<class T>
inline Vec<T> operator+(const Vec<T> &a,const Vec<T> &b){
  Vec<T> res(a);
  res+=b;
  return res;
}

template<class T>
inline Vec<T> operator-(const Vec<T> &a,const Vec<T> &b){
  Vec<T> res(a);
  res-=b;
  return res;
}

// scalar product
//template<typename T>
template<class T>
inline T operator*(const Vec<T> &a,const Vec<T> &b){
  a.test(b);
  T sum;
  int i;
  sum=0;
  for(i=0;i<a.dim;i++) sum+=a[i]*b[i];
  return sum;
}

// scalar product
//template<typename T>
template<class T>
inline T scalar_product(const Vec<T> &a,const Vec<T> &b){
  a.test(b);
  T sum;
  int i;
  sum=0;
  for(i=0;i<a.dim;i++) sum+=a[i]*b[i];
  return sum;
}

// == operator
template<class T>
inline int operator==(const Vec<T> &a,const Vec<T> &b){
  a.test(b);
  int i;
  for(i=0;i<a.dim;i++) if(a[i]!=b[i]) return 0;
  return 1;
}

// != operator
template<class T>
inline int operator!=(const Vec<T> &a,const Vec<T> &b){
  return !(a==b);
}

template<class T>
inline Vec<T> operator*(const Vec<T> &a,T x){
  Vec<T> res(a);
  res*=x;
  return res;
}

template<class T>
inline Vec<T> operator*(T x,const Vec<T> &a){
  return a*x;
}

template<class T>
inline Vec<T> operator/(const Vec<T> &a,T x){
  Vec<T> res(a);
  res/=x;
  return res;
}

template<class T>
int Vec<T>::default_size=0;

template<class T>
void Vec<T>::set_size(int n){
  default_size=n;
}

template<class T>
inline T sum_vector(const Vec<T> &a){
  int i,n;
  T sum;

  n=a.dim;
  sum=0;
  for(i=0;i<n;i++) sum+=a[i];
  return sum;  
}

// norm definitons -------------------

template<class T>
inline T norm1(const Vec<T> &a){
  T sum;

  sum=0;
  for(int i=0;i<a.dim;i++) sum+=abs(a[i]);
  return sum;
}

template<class T>
inline T norm2(const Vec<T> &a){
  T sum;

  sum=0;
  for(int i=0;i<a.dim;i++) sum+=absquad(a[i]);
  return sqrt(sum);
}

template<class T>
inline T norm_inf(const Vec<T> &a){
  T sum,aa;

  sum=0;
  for(int i=0;i<a.dim;i++){
    aa=abs(a[i]);
    if(aa>sum){
      sum=aa;
    }
  }
  return sum;
}

// diff_norm definitons ----------------------

template<class T>
inline T diff_norm1(const Vec<T> &a,const Vec<T> &b){
  T sum;

  a.test(b);
  sum=0;
  for(int i=0;i<a.dim;i++) sum+=abs(a[i]-b[i]);
  return sum;
}

template<class T>
inline T diff_norm2(const Vec<T> &a,const Vec<T> &b){
  T sum;

  a.test(b);
  sum=0;
  for(int i=0;i<a.dim;i++) sum+=absquad(a[i]-b[i]);
  return sqrt(sum);
}

template<class T>
inline T diff_norm_inf(const Vec<T> &a,const Vec<T> &b){
  T sum,aa;

  a.test(b);
  sum=0;
  for(int i=0;i<a.dim;i++){
    aa=abs(a[i]-b[i]);
    if(aa>sum){
      sum=aa;
    }
  }
  return sum;
}

typedef Vec<char> cvec;
typedef Vec<int> ivec;
typedef Vec<long> lvec;
typedef Vec<double> dvec;
typedef Vec<ivec> imatrix;
typedef Vec<lvec> lmatrix;
typedef Vec<dvec> dmatrix;

#endif
