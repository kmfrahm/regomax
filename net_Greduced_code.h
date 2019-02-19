/* including file for different pagerank subroutines */

#ifndef __GREDUCED__
#define __GREDUCED__

#include <stdio.h>
#include <math.h>
#include "filesize.h"
#include "matrix_simple.h"
#include "network_class.h"
#include "read_ascii.h"
#include "utf_string.h"
#include "quicksort_page.h"

double eps_pagerank=1e-13;

/********************************/

/*
  copies the first part of in until the first '.' to out, 
  "in" and "out" maybe the same in which case "in" will be 
  effectively shortened
*/

void get_first_dot(char *out,const char *in){
  int i;
  for(i=0;in[i]!='\0' && in[i]!='.';i++){
    out[i]=in[i];
  }
  out[i]='\0';
}

/********************************/

/*
  copies the first part of in until the last '.' to out, 
  "in" and "out" maybe the same in which case "in" will be 
  effectively shortened
*/

void get_last_dot(char *out,const char *in){
  int i,l;
  for(i=0;in[i]!='\0';i++){
    out[i]=in[i];
  }
  out[i]='\0';
  for(l=i-1;l>=0;l--){
    if(in[l]=='.'){
      out[l]='\0';
      break;
    }
  }
}

/*****************************************/
/* computes projection on reduced space */

inline void project_sub(dvec &small,dvec &big,ivec &node){
  int len=node.dim;
  small.resize(len);
  int i;
  
  for(i=0;i<len;i++) small[i]=big[node[i]];
}

/*************************************/

void print_mat(matrix& a,const char *filename=NULL,
	       const char *node_file_names=NULL){
  int i,j,dimx,dimy,len=0,nlen,l;
  char **node_names=NULL;
  FILE *fp;

  if(node_file_names!=NULL){
    node_names=full_read_ascii_file(node_file_names,len);

    // determine max node_name length
    nlen=0;
    for(i=0;i<len;i++){
      l=strlen_utf(node_names[i]);
      if(l>nlen) nlen=l;
    }
    nlen++;
  }

  if(filename==NULL){
    fp=stdout;
  }else{
    fp=fopen(filename,"w");
  }
  
  dimx=a.x();
  dimy=a.y();
  for(i=0;i<dimy;i++){
    for(j=0;j<dimx;j++){
      fprintf(fp,"%5d\t  %5d\t  %24.16lg",i,j,a(i,j));
      //      if(node_file_names!=NULL && i<len && j<len){
      if(i<len && j<len){
	//	fprintf(fp," %20s %20s",node_names[i],node_names[j]);
	fprintf(fp,"\t");
	fprintf_utf(fp,node_names[i],nlen);
	fprintf(fp,"\t");
	fprintf_utf(fp,node_names[j],nlen);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  if(filename!=NULL) fclose(fp);
  clear_file_buff(node_names);
}


/*************************************/

void print_mat_dens(matrix& a,const char *filename,int reverse=0){
  int i,j,ii,dimx,dimy;
  char buff[200];
  FILE *fp;

  get_last_dot(buff,filename);
  strcat(buff,".mat");

  fp=fopen(buff,"w");
  dimx=a.x();
  dimy=a.y();
  fprintf(fp,"%d\n%d\n",dimx*dimy,dimy);
  for(i=0;i<dimy;i++){
    if(reverse) ii=dimy-1-i; else ii=i;
    for(j=0;j<dimx;j++){
      fprintf(fp,"%24.16lg\n",a(ii,j));
    }
  }
  fclose(fp);
}

/*************************************/

void read_mat_dens(matrix& a,const char *filename,int reverse=0){
  int ii,i,j,n,dimx,dimy;
  FILE *fp;

  fp=fopen(filename,"r");
  if(fp==NULL) error("Input file in read_mat_dens not found!");
  fscanf(fp,"%d%d",&n,&dimy);
  dimx=n/dimy;
  dvec::set_size(0);
  matrix b(dimy,dimx);
  fprintf(fp,"%d\n%d\n",dimx*dimy,dimy);
  for(i=0;i<dimy;i++){
    if(reverse) ii=dimy-1-i; else ii=i;
    for(j=0;j<dimx;j++){
      fscanf(fp,"%lf",&b(ii,j));
    }
  }
  fclose(fp);
  a=b;
}

/*************************************/

#define BUFF_READ_LEN 1001
void read_mat(matrix& a,const char *filename){
  FILE *fp;
  printf("Matrix input file = \"%s\"\n",filename);
  fflush(stdout);
  fp=fopen(filename,"r");
  if(fp==NULL) error("Input file in read_mat not found!");
  dvec::set_size(0);

  // use binary read mode if filename contains ".bin"
  if(strstr(filename,".bin")==NULL){
    printf("Using ascii mode for matrix read\n");
    fflush(stdout);
    // ascii file mode
    char buff[BUFF_READ_LEN];
    int ii,i,j,n,dimx,dimy;
    double val;

    // first reading to determine matrix size
    dimx=0; dimy=0;
    while(fgets(buff,BUFF_READ_LEN,fp)!=0){
      if(sscanf(buff,"%d%d%lf",&i,&j,&val)>=3){
	if(i>dimy) dimy=i;
	if(j>dimx) dimx=j;      
      }
    }
    fclose(fp);
    dimx++; dimy++;
    printf("matrix size: %d x %d\n",dimy,dimx);
    fflush(stdout);
    
    matrix b(dimy,dimx);
    // second reading 
    fp=fopen(filename,"r");
    if(fp==NULL) error("Input file in read_mat not found!");
    n=0;
    while(fgets(buff,BUFF_READ_LEN,fp)!=0){
      if(sscanf(buff,"%d%d%lf",&i,&j,&val)>=3){
	b(i,j)=val;
	n++;
      }
    }
    fclose(fp);
    printf("read n = %d  matrix elements: %d x %d - %d = %d\n",
	   n,dimy,dimx,n,dimy*dimx-n);
    fflush(stdout);
    a=b;
  }else{
    printf("Using binary mode for matrix read\n");
    fflush(stdout);
    // binary file mode
    int l,n,i;
    l=fread(&n,sizeof(int),1,fp);
    if(l!=1) error("binary file read error for matrix size");
    printf("Binary file: matrix size = %d x %d\n",n,n);
    fflush(stdout);
    if(n<=0) error("invalid matrix size");
    matrix b(n,n);
    for(i=0;i<n;i++){
      l=fread(&b[i][0],sizeof(double),n,fp);
      if(l!=n) error("binary file read error for matrix elements");
    }
    fclose(fp);    
    a=b;
  }
  printf("\n");
}

/*****************************************/
/* calculation of \sum_i abs(a(i)-b(i))/(abs(a(i)+b(i)) */

inline double diff_norm_rel(dvec &a,dvec &b){
  double sum,ss;
  int i,n;
  
  n=a.size(); sum=0.0;
  //#pragma omp parallel for reduction(+:sum)
  for(i=0;i<n;i++){
    ss=abs(a[i])+abs(b[i]);
    if(ss==0) continue;
    sum+=(abs(a[i]-b[i])/ss);
  }
  return sum;
}

/*****************************************/
/* calculation of \sum_i a(i) = e^T a */

inline double sum_vector(dvec &a){
  double sum;
  int i,n;
  
  n=a.size(); sum=0.0;
  //  #pragma omp parallel for reduction(+:sum)
  for(i=0;i<n;i++) sum+=a[i];
  return sum;
}

/*****************************************/
/* normalization of pagerank \sum_i a(i)=1 assuming a(i)>=0 
   return value = used 1-norm before normalization
*/

inline double pagerank_normalize(dvec &a){
  double sum;

  sum=sum_vector(a);
  a/=sum;
  return sum;
}

/*****************************************/
/* small matrix vector multiplication */

void mat_vec_mult(dvec &out,matrix& a,dvec &in){
  int n=in.size();
  if(n!=out.size() || n!=a.x() || n!=a.y()) 
    error("wrong mat-vec-size in mat_vec_mult");

  int i,j;
  for(i=0;i<n;i++){
    out[i]=0;
    for(j=0;j<n;j++) out[i]+=a(i,j)*in[j];
  }
}


/*****************************************/
/* calculation of pagerank with the power method 
   for small reduced google matrices 
*/

void calc_pagerank_small(dvec &pagerank,matrix& GR,int iprint){
  double quality,quality_rel,q1,qfak,pnorm;
  int i,max_iter;

  if(iprint<=0) iprint=1;
  max_iter=400;

  printf("max_iter = %d\n",max_iter);
  fflush(stdout);
  qfak=1.0+0.15/2.0;
  pnorm=pagerank_normalize(pagerank);
  dvec a(pagerank);
  quality_rel=1e40;
  for(i=0;i<=max_iter;i++){
    swap(a,pagerank);
    mat_vec_mult(pagerank,GR,a);
    if(i%iprint==0 || i==max_iter){
      quality=diff_norm1(pagerank,a);
      q1=quality_rel;
      quality_rel=diff_norm_rel(pagerank,a);
      //      pnorm=pagerank_normalize(pagerank);
      pnorm=sum_vector(pagerank);
      printf("%5d  %18.10lg  %18.10lg  %25.16lg\n",
	     i,quality,quality_rel,pnorm);
      fflush(stdout);
      if(quality_rel<eps_pagerank) break;
      if(quality_rel<1e-4){
	if(quality_rel*qfak>q1) break;
      }
    }
  }
  printf("Convergence at i = %d.\n",i);
  fflush(stdout);
}

/*****************************************/
/* calculation of pagerank with the power method */

void calc_pagerank_power(dvec &pagerank,network &net,
			 double delta_alpha,int iprint,int trans_flag=0){
  double quality,quality_rel,q1,qfak,pnorm;
  int i,max_iter;

  if(iprint<=0) iprint=1;
  max_iter=(int)(-log(eps_pagerank)/(delta_alpha+3E-7));
  max_iter*=2;

  printf("max_iter = %d\n",max_iter);
  fflush(stdout);
  qfak=1.0+delta_alpha/2.0;
  pnorm=pagerank_normalize(pagerank);
  dvec a(pagerank);
  quality_rel=1e40;
  for(i=0;i<=max_iter;i++){
    swap(a,pagerank);
    if(trans_flag){
      net.GTmult(delta_alpha,pagerank,a);
    }else{
      net.GGmult(delta_alpha,pagerank,a);
    }
    if(i%iprint==0 || i==max_iter){
      quality=diff_norm1(pagerank,a);
      q1=quality_rel;
      quality_rel=diff_norm_rel(pagerank,a);
      //      pnorm=pagerank_normalize(pagerank);
      pnorm=sum_vector(pagerank);
      printf("%5d  %18.10lg  %18.10lg  %25.16lg\n",
	     i,quality,quality_rel,pnorm);
      fflush(stdout);
      if(quality_rel<eps_pagerank) break;
      if(quality_rel<1e-4){
	if(quality_rel*qfak>q1) break;
      }
    }
  }
  printf("Convergence at i = %d.\n",i);
  fflush(stdout);
}

/*****************************************/
/* calculation of pagerank for "PG" with P=projector and 
   G=Google matrix, return value = corresponding eigenvalue,
   with the power method 

   node = array of length len, P=projector on nodes different from node[i]
*/

double calc_pagerank_project(dvec &pagerank,network &net,
			     double delta_alpha,int iprint,
			     ivec &node,int trans_flag=0){
  double quality,quality_rel,q1,qfak,pnorm,dlambda,dlambda_old;
  int i,max_iter,l;

  if(iprint<=0) iprint=1;
  max_iter=(int)(-log(eps_pagerank)/(delta_alpha+3E-7));
  max_iter*=2;

  printf("max_iter = %d\n",max_iter);
  fflush(stdout);
  qfak=1.0+delta_alpha/2.0;
  pnorm=pagerank_normalize(pagerank);
  dvec a(pagerank);
  quality_rel=1e40;
  dlambda=0;
  for(l=0;l<node.dim;l++){
    dlambda+=pagerank[node[l]];
    pagerank[node[l]]=0;
  }
  dlambda_old=dlambda;
  pnorm=pagerank_normalize(pagerank);
  if(trans_flag) dlambda=1.0-pnorm;
  for(i=0;i<=max_iter;i++){
    swap(a,pagerank);
    if(trans_flag){
      net.GTmult(delta_alpha,pagerank,a);
    }else{
      net.GGmult(delta_alpha,pagerank,a);
    }
    //    pnorm=pagerank_normalize(pagerank);
    //    printf("--> %5d  %25.16lg\n",i,pnorm);
    //    fflush(stdout);
    dlambda=0;
    for(l=0;l<node.dim;l++){
      dlambda+=pagerank[node[l]];
      pagerank[node[l]]=0;
    }
    pnorm=pagerank_normalize(pagerank);
    if(trans_flag) dlambda=1.0-pnorm;

    if(i%iprint==0 || i==max_iter){
      quality=diff_norm1(pagerank,a);
      q1=quality_rel;
      quality_rel=diff_norm_rel(pagerank,a);
      //      pnorm=pagerank_normalize(pagerank);
      //      pnorm=sum_vector(pagerank);
#pragma omp critical(print)
      {
	printf("%5d  %18.10lg  %18.10lg  %25.16lg  %18.10lg  %25.16lg\n",
	       i,quality,quality_rel,dlambda,abs(dlambda-dlambda_old),pnorm);
	fflush(stdout);
      }
      dlambda_old=dlambda;
      if(quality_rel<eps_pagerank) break;
      if(quality_rel<1e-3){
	if(quality_rel*qfak>q1) break;
      }
    }
  }
#pragma omp critical(print)
  {
    printf("Convergence at i = %d  with lambda = %25.16lg.\n",i,1.0-dlambda);
    fflush(stdout);
  }
  return dlambda;
}

/**************************************/
// print the pagerank vector to a file

void print_pagerank(dvec &pagerank,network &net,double delta_alpha,
		    const char *extra,int print_number=0,int ten_number=100){
  double sum;
  int j,n;
  char buff[300];
  FILE *fp;

  n=pagerank.size();
  ivec permut(n);

  init_permutation(permut.c,n);
  quicksort_down(pagerank.c,0,n-1,permut.c);

  // assure all elements of pagerank are positif
  if(pagerank[permut[0]]<0) pagerank*=(-1.0);

  // create file name
  sprintf(buff,"pagerank%s_%s_%lg.dat",extra,net.base_name.c,delta_alpha);

  // writing of pagerank file
  fp=fopen(buff,"w");
  fprintf(fp,"# size = %10d\n",n);

  if(print_number<=0){
    for(j=0;j<n;j++){
      fprintf(fp,"%6d %24.15lg %6d\n",j,pagerank[permut[j]],permut[j]);
    }
  }else{
    double ten_fak=exp(log(10.0)/(double)ten_number);
    int jlast=0;
    for(j=1;j<=n;j++){
      if(j<=print_number || (double)j>ten_fak*(double)jlast || j==n){
	jlast=j;
	fprintf(fp,"%6d %24.15lg %6d\n",j-1,pagerank[permut[j-1]],permut[j-1]);
      }
    }
  }
  fclose(fp);
}

/**************************************/
/* print the subpagerank vector to a file
*/

void print_subpagerank(dvec &pagerank,network &net,double delta_alpha,
		       const char *extra,ivec &node,
		       const char *node_file_names=NULL){
  double sum;
  int i,n,len2=0,nlen,l;
  char buff[300];
  char **node_names=NULL;
  FILE *fp;

  if(node_file_names!=NULL){
    node_names=full_read_ascii_file(node_file_names,len2);

    // determine max node_name length
    nlen=0;
    for(i=0;i<len2;i++){
      l=strlen_utf(node_names[i]);
      if(l>nlen) nlen=l;
    }
    nlen++;
  }

  n=pagerank.size();
  ivec permut(n),pinv(n);

  init_permutation(permut.c,n);
  quicksort_down(pagerank.c,0,n-1,permut.c);
  for(i=0;i<n;i++) pinv[permut[i]]=i;

  // assure all elements of pagerank are positif
  if(pagerank[permut[0]]<0) pagerank*=(-1.0);

  sum=0.0;
  for(i=0;i<node.dim;i++) sum+=pagerank[node[i]];
  
  // projection of subpagerank
  //  dvec pg(node.dim);
  //  project_sub(pg,pagerank,node);

  // create file name
  sprintf(buff,"subpagerank%s_%s_%lg.dat",extra,net.base_name.c,delta_alpha);

  // writing of pagerank file
  fp=fopen(buff,"w");
  fprintf(fp,"## size = %10d  norm = %25.15lg\n",node.dim,sum);
  fprintf(fp,"# reduced index, value, normalized value, original index, K index");
  if(len2>0) fprintf(fp,", name");
  fprintf(fp,"\n");
  sum=1.0/sum;

  for(i=0;i<node.dim;i++){
    fprintf(fp,"%6d\t %24.15lg\t %24.15lg\t %8d\t %8d",i,pagerank[node[i]],
	    sum*pagerank[node[i]],node[i],pinv[node[i]]);
    if(i<len2){
      fprintf(fp,"\t");
      fprintf_utf(fp,node_names[i],nlen);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  clear_file_buff(node_names);
}

/*****************************************/
/* also computes the usual PageRank since only 3 threads are used */

inline double compute_project(dvec &right,dvec &left,dvec &pg,
			      network &net,
			      double delta_alpha,ivec &node){
  int iprint=10;
  double sp,dlambda1,dlambda2,dlambda3;
  ivec node0(0);

  right.put_value(1.0);
  left.put_value(1.0);
  pg.put_value(1.0);

#pragma omp parallel sections
  {
#pragma omp section
    dlambda2=calc_pagerank_project(left,net,delta_alpha,iprint,node,1);
#pragma omp section
    dlambda1=calc_pagerank_project(right,net,delta_alpha,iprint,node);
#pragma omp section
    dlambda3=calc_pagerank_project(pg,net,delta_alpha,iprint,node0);
  }

  sp=1.0/scalar_product(left,right);
  left*=sp;
  sp=scalar_product(left,right);
#pragma omp critical(print)
  {
    printf("dlambda = %24.16lg   diff = %lg\n",
	   dlambda1,abs(dlambda1-dlambda2));
    printf("TEST: psi_left^T * psi_right = %26.16lg\n",sp);
    fflush(stdout);
  }

  return dlambda1;
}

/*****************************************/
/* computes: v = (1/f) P * v = (1/f) right * (left^T v) 
   with f being some inverse factor */

inline void projectP(dvec &right,dvec &left,dvec &v,double f=1){
  int i,n=v.size();
  double sp;

  sp=scalar_product(left,v)/f;
  v.test(right);
  for(i=0;i<n;i++) v[i]=sp*right[i];
}

/*****************************************/
/* computes: v = Q * v = (1 - P) * v = v - right * (left^T v) */

inline void projectQ(dvec &right,dvec &left,dvec &v){
  double sp;

  sp=scalar_product(left,v);
  v.lam_diff(sp,right);
}

/*****************************************/

inline void get_cnode(int n,ivec &cnode,ivec &node){
  int i,j;

  if(n<node.dim) error("n<len in get_cnode");
  ivec f(n);

  for(i=0;i<n;i++) f[i]=0;
  for(i=0;i<node.dim;i++){
    if(node[i]>=n || node[i]<0) error("node[i] out of range in get_cnode");
    f[node[i]]=1;
  }
  j=0;
  for(i=0;i<n;i++){
    if(!f[i]){
      cnode[j]=i;
      j++;
    }
  }
  //  printf("j = %d   cnode.dim = %d\n",j,cnode.dim);
  //  fflush(stdout);
  if(j!=cnode.dim) error("wrong cnode size in get_cnode");
}

/*****************************************/

#define MAX_GR_SIMPLE 3000 

void compute_GR_simple(matrix &G_R,network &net,double delta_alpha,
		      ivec &node){
  int n=net.size;
  if(n>MAX_GR_SIMPLE) error("Too large network for compute_GR_simple");
  int nr=node.dim;
  int ns=n-nr;
  ivec cnode(ns);
  get_cnode(n,cnode,node);
  dvec::set_size(0);
  matrix G_rr(nr,nr),G_rs(nr,ns),G_sr(ns,nr),G_ss(ns,ns);
  dvec in(n,0.0),out(n);

  int i,j;
  // filling of G_rr and G_sr
  for(i=0;i<nr;i++){
    in[node[i]]=1;
    net.GGmult(delta_alpha,out,in);
    in[node[i]]=0;

    for(j=0;j<nr;j++) G_rr(j,i)=out[node[j]];
    for(j=0;j<ns;j++) G_sr(j,i)=out[cnode[j]];
  }

  // filling of G_rs and G_ss
  for(i=0;i<ns;i++){
    in[cnode[i]]=1;
    net.GGmult(delta_alpha,out,in);
    in[cnode[i]]=0;

    for(j=0;j<nr;j++) G_rs(j,i)=out[node[j]];
    for(j=0;j<ns;j++) G_ss(j,i)=-out[cnode[j]];
  }
  for(j=0;j<ns;j++) G_ss(j,j)+=1.0;
  G_sr/=G_ss;
  G_R=G_rs*G_sr;
  G_R+=G_rr;

}

/*****************************************/
#ifndef ITER_MODE
// implementation with power series mode:
// s=(\sum_{n=0}^\infty g^n)]v
// using: s_0=v, f_0=v: and f_{n+1}=g f_n,  s_{n+1}=s_n+f_{n+1}
// with g= \bar G_{ss}


void compute_GR(matrix &G_R,matrix& G_rr,matrix& G_pr,matrix& G_qr,
		   matrix& G_I,dvec& psiL,dvec& psiR,dvec& pg,
		network &net,double delta_alpha,ivec node){
  int n=net.size;
  int nr=node.dim;
  int ns=n-nr;
  if(G_R.x()!=nr || G_R.y()!=nr) 
    error("Wrong matrix size of G_R  in comput_GR");
  if(G_rr.x()!=nr || G_rr.y()!=nr) 
    error("Wrong matrix size of G_rr in comput_GR");
  if(G_pr.x()!=nr || G_pr.y()!=nr) 
    error("Wrong matrix size of G_pr in comput_GR");
  if(G_qr.x()!=nr || G_qr.y()!=nr) 
    error("Wrong matrix size of G_qr in comput_GR");
  if(G_I.x()!=nr || G_I.y()!=nr) 
    error("Wrong matrix size of G_I  in comput_GR");
  double dlambda;

  int j,l;
  double quality;

  int i,max_iter;

  max_iter=(int)(-log(eps_pagerank)/(delta_alpha+3E-7));
  max_iter*=2;

  printf("Computation of left and right eigenvectors of G_ss\n");
  fflush(stdout);
  dlambda=compute_project(psiR,psiL,pg,net,delta_alpha,node);

  dvec in(n),out(n),s(n),t(n),f(n),f2(n);
  // note that the last line also fixes the default size of dvec to n 
  // which is important in the private declaration below which implicitely 
  // calls the default constructor of dvec for each thread

#pragma omp parallel for schedule(dynamic) private(in,out,s,t,f,f2,j,l,quality)
  for(i=0;i<nr;i++){
    in.put_value(0.0);
    in[node[i]]=1;
    net.GGmult(delta_alpha,out,in);
    in[node[i]]=0;
    for(j=0;j<nr;j++){
      G_R(j,i)=out[node[j]];
      G_rr(j,i)=out[node[j]];
      out[node[j]]=0;
    }
    s=out;
    projectP(psiR,psiL,out,dlambda);
    projectQ(psiR,psiL,s);
    f=s;

    for(l=0;l<max_iter;l++){
      t=s;
      net.GGmult(delta_alpha,f2,f,0);
      swap(f,f2);
      for(j=0;j<nr;j++) f[node[j]]=0;
      projectQ(psiR,psiL,f);
      s+=f;
      quality=diff_norm1(t,s);
#pragma omp critical(print)
      {
	if(l%10==0){
	  printf("%5d  %5d  %18.10lg  %18.10lg\n",i,l,quality,norm1(f));
	  fflush(stdout);
	}
      }
      //      if(quality<eps_pagerank) break;
      if(quality<=0) break;
    }
#pragma omp critical(print)
    {
      printf("%5d  ",i);
      printf("Convergence: %5d  %5d  %18.10lg  %18.10lg\n",
	     i,l,quality,norm1(f));
      fflush(stdout);
    }
    net.GGmult(delta_alpha,f,out,0);
    for(j=0;j<nr;j++){
      G_pr(j,i)=f[node[j]];
    }
    net.GGmult(delta_alpha,f,s,0);
    for(j=0;j<nr;j++){
      G_qr(j,i)=f[node[j]];
    }
    out+=s;
    net.GGmult(delta_alpha,f,out,0);
    for(j=0;j<nr;j++){
      G_I(j,i)=f[node[j]];
      G_R(j,i)+=f[node[j]];
    }
  }
}

/*****************************************/

#else 

// implementation with a modified iteration mode:
// s=(\sum_{n=0}^\infty g^n)]v
// using: s_0=v, s_{n+1}=v+g*s_n
// with g= \bar G_{ss}

void compute_GR(matrix &G_R,matrix& G_rr,matrix& G_pr,matrix& G_qr,
		   matrix& G_I,dvec& psiL,dvec& psiR,dvec& pg,
		network &net,double delta_alpha,ivec node){
  int n=net.size;
  int nr=node.dim;
  int ns=n-nr;
  if(G_R.x()!=nr || G_R.y()!=nr) 
    error("Wrong matrix size of G_R  in comput_GR");
  if(G_rr.x()!=nr || G_rr.y()!=nr) 
    error("Wrong matrix size of G_rr in comput_GR");
  if(G_pr.x()!=nr || G_pr.y()!=nr) 
    error("Wrong matrix size of G_pr in comput_GR");
  if(G_qr.x()!=nr || G_qr.y()!=nr) 
    error("Wrong matrix size of G_qr in comput_GR");
  if(G_I.x()!=nr || G_I.y()!=nr) 
    error("Wrong matrix size of G_I  in comput_GR");
  double dlambda;

  int j,l;
  double quality,old_quality,fak;
  //  fak=1-delta_alpha/2;
  //  fak=0.99;
  fak=1-delta_alpha/10;

  printf("### Convergence: FAK = %lg\n",fak);
  fflush(stdout);

  int i,max_iter;

  //  max_iter=(int)(-log(eps_pagerank)/(delta_alpha+3E-7));
  //  max_iter*=2;
  max_iter=(int)(-18*log(10.0)/(log(1-delta_alpha)+1e-7));
  printf("Using modified iteration mode with max_iter = %d\n",max_iter);

  printf("Computation of left and right eigenvectors of G_ss\n");
  fflush(stdout);
  dlambda=compute_project(psiR,psiL,pg,net,delta_alpha,node);

  dvec in(n),out(n),s(n),t(n),v(n),s2(n);
  // note that the last line also fixes the default size of dvec to n 
  // which is important in the private declaration below which implicitely 
  // calls the default constructor of dvec for each thread

#pragma omp parallel for schedule(dynamic) private(in,out,s,t,v,s2,j,l,quality,old_quality)
  for(i=0;i<nr;i++){
    in.put_value(0.0);
    in[node[i]]=1;
    net.GGmult(delta_alpha,out,in);
    in[node[i]]=0;
    for(j=0;j<nr;j++){
      G_R(j,i)=out[node[j]];
      G_rr(j,i)=out[node[j]];
      out[node[j]]=0;
    }
    s=out;
    projectP(psiR,psiL,out,dlambda);
    projectQ(psiR,psiL,s);
    v=s;
    quality=1e100;
    for(l=0;l<max_iter;l++){
      t=s;
      net.GGmult(delta_alpha,s2,s,0);
      for(j=0;j<nr;j++) s2[node[j]]=0;
      projectQ(psiR,psiL,s2);
      swap(s,s2);
      s+=v;
      old_quality=quality;
      quality=diff_norm1(t,s);
#pragma omp critical(print)
      {
	if(l%10==0){
	  printf("%5d  %5d  %18.10lg  %18.10lg\n",i,l,quality,old_quality);
	  fflush(stdout);
	}
      }
      if(quality<eps_pagerank && quality>old_quality*fak) break;
      //      if(quality<=0) break;
    }
#pragma omp critical(print)
    {
      printf("%5d  ",i);
      printf("Convergence: %5d  %5d  %18.10lg  %18.10lg\n",
	     i,l,quality,old_quality);
      fflush(stdout);
    }
    net.GGmult(delta_alpha,s2,out,0);
    for(j=0;j<nr;j++){
      G_pr(j,i)=s2[node[j]];
    }
    net.GGmult(delta_alpha,s2,s,0);
    for(j=0;j<nr;j++){
      G_qr(j,i)=s2[node[j]];
    }
    out+=s;
    net.GGmult(delta_alpha,s2,out,0);
    for(j=0;j<nr;j++){
      G_I(j,i)=s2[node[j]];
      G_R(j,i)+=s2[node[j]];
    }
  }
}

#endif 

/*****************************************/

void renorm_google(matrix &a){
  int n=a.x();
  if(n!=a.y()){
    printf("Warning: matrix not square in renorm_google.\n");
    fflush(stdout);
  }
  int i,j;
  double sum;
  for(i=0;i<n;i++){
    sum=0;
    for(j=0;j<n;j++){
      if(a(j,i)<0){
	//	printf("Warning: negative element a(%d,%d) in renorm_google.\n",j,i);
	//	printf("Taking abs-value.\n");
	a(j,i)=abs(a(j,i));
	fflush(stdout);
      }
      sum+=a(j,i);
    }
    sum=1.0/sum;
    for(j=0;j<n;j++) a(j,i)*=sum;    
  }
}

/*****************************************/
/* a = a + b for non-diagonal elements of b */

void add_non_diag(matrix &a,matrix &b){
  int n=a.x();
  int i,j;

  for(i=0;i<n;i++) for(j=0;j<n;j++){
    if(i!=j) a(i,j)+=b(i,j);
  }
}


#endif
