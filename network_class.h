/* simple C++ class for google-matrix networks */

#ifndef __NETWORK__
#define __NETWORK__

#include "vec.h" 
#include <string.h>
#include <time.h>
#include "quicksort_nonpermut.h"


#define BASE_NAME_LEN 201

class network{
  // the following function assumes reasonable values for size and link_len
  // and resizes the different vectors accordingly
  void init_mem(void){
    base_name.resize(BASE_NAME_LEN);
    from.resize(link_len);
    to.resize(link_len);
#ifdef USE_PROBS
    prob.resize(link_len);
#endif
    link_num.resize(size);
    firstpos.resize(size+1);
    // note that dangling vector will be resized later
    // when link_num and firstpos are constructed
  }

  // raw copy function
  void copy(const network &a){
    if(this==&a) return; // do nothing if copy of identical objects
    size=a.size;
    link_len=a.link_len;
    // note all complications of vector copy etc. are taken care of 
    // in Vec<> template defined in "vec.h"
    base_name=a.base_name;
    from=a.from;
    to=a.to;
#ifdef USE_PROBS
    prob=a.prob;
#endif
    link_num=a.link_num;
    firstpos=a.firstpos;
    dangling=a.dangling;
  }
 public:
  // number of nodes
  int size;
  // number of links
  int link_len;
  // base_name
  cvec base_name;
  // in- and out-going links (dim = link_len)
  ivec from,to;
#ifdef USE_PROBS
  // probability of links (dim = link_len)
  dvec prob;
#endif
  // link_num = vector of number of links per node (dim = size)
  // firstpos = vector of positions of new link positions (dim = size+1)
  // in principle prob[i]=1/link_num[from[i]];
  // and link_num[i]=firstpos[i+1]-firstpos[i]
  ivec link_num,firstpos;
  // vector of dangling nodes 
  ivec dangling;

  // ---- functions -----
  // default constructor => empty network 
  network(){}

  // calculate first-position array, dangling nodes 
  // eventually normalize prob-array 
  void complete(void);

  // constructor from file
  // this will be the main input function for networks
  network(const char *link_filename){
    //    read_netfile(link_filename);
    read_binfile(link_filename);
  }

  // very simple destructor, cleaning of memory is done in Vec<> destructor
  ~network(){}

  // "=" operator
  network& operator=(network &a){
    copy(a);
    return *this;
  };

  // copy constructor
  network(const network &a){
    copy(a);
  };

#ifdef USE_PROBS
  // the following function is only availabe if USE_PROBS is 
  // activated

  // renormalize the column sums to 1
  void renormalize(void){
    int i,j,a,b;
    double sum;

    b=firstpos[0];
    for(i=0;i<size;i++){
      sum=0;
      a=b; b=firstpos[i+1];
      for(j=a;j<b;j++){
	sum+=prob[j];
      }
      if(sum!=0){
	sum=1.0/sum;
	for(j=a;j<b;j++){
	  prob[j]*=sum;
	}
      }
    }
  }
#endif  

  // resorting of "from" and "to" links 
  // and reconstruction of firstpos and dangling nodes
  // the macro PROB_ARG is defined in "quicksort_nonpermut.h"
  void resort(){
    quicksort_up(from.c,to.c,0,link_len-1 PROB_ARG);
    complete();
  }

  void swap_in_out(void){
    swap(from,to);
  }

  // transformation to che rank
  void che_trans(void){
    int i,t;

    printf("starting of che transformation\n");
    fflush(stdout);
    swap_in_out();
    printf("link reverse done\n");
    fflush(stdout);
    // create new base name
    strcat(base_name.c,"che");
    resort();
#ifdef USE_PROBS
    //    renormalize();
#endif
    printf("resorting and completion done\n");
    fflush(stdout);
  }

  // multiplication with G
  void GGmult(double delta_alpha,dvec &out,const dvec &in,int norm_flag);
  // multiplication with G^T, transpose google matrix
  void GTmult(double delta_alpha,dvec &out,const dvec &in,int norm_flag);


  void write_netfile(const char *outfile);
  void read_netfile(const char *outfile);
  void write_binfile(const char *outfile);
  void read_binfile(const char *outfile);
  void write_fullnetfile(const char *outfile);
};


/********************************************/
/* diverse functions */

/****************************************/
// calculate first-position array, dangling nodes and prob-array
// from links
// the from array is assumed to be sorted 

void network::complete(void){
  int i,jj,a,b;
  int dangling_len;

  // building of 1stpos-structure
  jj=0;
  firstpos[jj]=0;
  for(i=0;i<link_len;i++){
    while(jj<from[i]){
      jj++;
      firstpos[jj]=i;
    }
  }
  while(jj<size){
    jj++;
    firstpos[jj]=i;
  }

  // calculation of dangling nodes
  dangling_len=0;
  for(i=0;i<size;i++){
    if(firstpos[i]==firstpos[i+1]){
      dangling_len++;
    }
  }
  dangling.resize(dangling_len);

  if(dangling_len>0){
    dangling_len=0;
    for(i=0;i<size;i++){
      if(firstpos[i]==firstpos[i+1]){
	dangling[dangling_len++]=i;
      }
    }
  }
  
  // compute link_num vector
  for(jj=0;jj<size;jj++){
    link_num[jj]=firstpos[jj+1]-firstpos[jj];
  }
}



/****************************************/

void network::read_netfile(const char *link_filename){
  int start,end,i,l;
  FILE *fp;

  printf("\n\n****** => Reading of data file \n");
  fflush(stdout);

  fp=fopen(link_filename,"r");
  if(fp==NULL) error("datafile file not found");
  if(fscanf(fp,"%d",&size)<1){
    error("unable to read size");
  }
  if(fscanf(fp,"%d",&link_len)<1){
    printf("--> %d   %d\n",size,link_len);
    fflush(stdout);
    error("unable to read link_len");
  }

  // create or adapt arrays
  init_mem();

  //----------------------------------------
  // create base name for output files
  //  for(i=0;link_filename[i]!='\0' && link_filename[i]!='.';i++){
  for(i=0;link_filename[i]!='\0';i++){
    base_name[i]=link_filename[i];
  }
  base_name[i]='\0';
  for(l=i-1;l>=0;l--){
    if(link_filename[l]=='.'){
      base_name[l]='\0';
      break;
    }
  }

  printf("****** => Reading of integer connection matrix\n");
  fflush(stdout);

  //----------------------------------------

  for(i=0;i<link_len;i++){
    if(fscanf(fp,"%d%d",&start,&end)<2) error("netfile error");
    // take into account that matrix indicies start at 0 and not 1
    from[i]=start-1;
    to[i]=end-1;
#ifdef USE_PROBS
    if(fscanf(fp,"%lf",&prob[i])<1) error("netfile error");
#endif
  }
  fclose(fp);

  complete();

  printf("size = %d   link_len = %d   dangling_len = %d\n",
	 size,link_len,dangling.dim);
  fflush(stdout);

  printf("****** => Reading of data file finished\n");
  fflush(stdout);
}

/****************************************/

void network::write_netfile(const char *outfile){
  int i;
  FILE *fp;

  fp=fopen(outfile,"w");
  fprintf(fp,"%d\n%d\n",size,link_len);
  for(i=0;i<link_len;i++){
    fprintf(fp,"%d %d",from[i]+1,to[i]+1);
#ifdef USE_PROBS
    fprintf(fp," %.17lg",prob[i]);
#endif
    fprintf(fp,"\n");
  }
  fclose(fp);
}

/****************************************/

void network::read_binfile(const char *link_filename){
  char buff[10];
  int l,i;
  FILE *fp;

  printf("\n\n****** => First try to read as binary data file \n");
  fflush(stdout);

  fp=fopen(link_filename,"r");
  if(fp==NULL) error("datafile file not found");
  
  l=fread(buff,sizeof(char),1,fp);
  if(l!=1) error("binary file read error with test charactor");
  if(buff[0]!='B'){
    fclose(fp);
    printf("****** => Data file is not a valid binary file\n");
    printf("****** => Try to read as ascii file\n");
    read_netfile(link_filename);
    return;
  }

  l=fread(&size,sizeof(int),1,fp);
  if(l!=1) error("binary file read error with size");
  l=fread(&link_len,sizeof(int),1,fp);
  if(l!=1) error("binary file read error with link_len");

  // create arrays
  init_mem();

  //----------------------------------------
  // create base name for output files
  //  for(i=0;link_filename[i]!='\0' && link_filename[i]!='.';i++){
  for(i=0;link_filename[i]!='\0';i++){
    base_name[i]=link_filename[i];
  }
  base_name[i]='\0';
  for(l=i-1;l>=0;l--){
    if(link_filename[l]=='.'){
      base_name[l]='\0';
      break;
    }
  }

  printf("****** => Reading of integer connection matrix\n");
  fflush(stdout);

  //----------------------------------------

  l=fread(from.c,sizeof(int),link_len,fp);
  if(l!=link_len) error("binary file read error with from array");

  l=fread(to.c,sizeof(int),link_len,fp);
  if(l!=link_len) error("binary file read error with to array");

#ifdef USE_PROBS
  l=fread(prob.c,sizeof(double),link_len,fp);
  if(l!=link_len) error("binary file read error with prob array");
#endif

  fclose(fp);
  complete();

  printf("size = %d   link_len = %d   dangling_len = %d\n",
	 size,link_len,dangling.dim);
  fflush(stdout);


  printf("****** => Reading of binary data file finished\n");
  fflush(stdout);
}

/****************************************/

void network::write_binfile(const char *outfile){
  int i;
  FILE *fp;

  fp=fopen(outfile,"w");
  // first write a marker for binary file
  fwrite("B",sizeof(char),1,fp);
  fwrite(&size,sizeof(int),1,fp);
  fwrite(&link_len,sizeof(int),1,fp);
  fwrite(from.c,sizeof(int),link_len,fp);
  fwrite(to.c,sizeof(int),link_len,fp);
#ifdef USE_PROBS
  fwrite(prob.c,sizeof(double),link_len,fp);
#endif
  fclose(fp);
}

/****************************************/

void network::write_fullnetfile(const char *outfile){
  int i;
  FILE *fp;

  fp=fopen(outfile,"w");
  fprintf(fp,"%d\n%d\n",size,link_len);
  for(i=0;i<link_len;i++){
#ifdef USE_PROBS
    fprintf(fp,"%d %d %.17lg\n",from[i]+1,to[i]+1,prob[i]);
#else
    fprintf(fp,"%d %d %.17lg\n",from[i]+1,to[i]+1,1.0/link_num[from[i]]);
#endif
  }
  fclose(fp);
}


/**************************************/
// note delta_alpha = 1-alpha !!
// if norm_flag is set the vector "in" is assumed to be sum normalized
// which increases the quality of convergence

void network::GGmult(double delta_alpha,dvec &out,const dvec &in,
		     int norm_flag=1){
  double sum,val;
  int i,a,b,jj;

  // contribution from dangling modes
  // ==> 1/N e d^T
  sum=0.0;
  if(dangling.dim>0){
    for(i=0;i<dangling.dim;i++){
      sum+=in[dangling[i]];
    }
    sum/=size;
  }
  for(i=0;i<size;i++) out[i]=sum;


  //  Computation of out=S*in
  for(jj=0;jj<size;jj++){
    a=firstpos[jj];
    b=firstpos[jj+1];
    if(a>=b) continue;
#ifndef USE_PROBS
    //    val=in[from[a]]/(b-a);
    //    val=in[jj]/(b-a);
    val=in[jj]/link_num[jj];
#else
    //    val=in[from[a]];
    val=in[jj];
#endif
    //    val=in[from[a]]*prob[a];
    for(i=a;i<b;i++){
#ifndef USE_PROBS
      out[to[i]]+=val;
#else
      out[to[i]]+=prob[i]*val;
#endif
    }
  }

  // computation of out=G*in, i.e. damping factor contributions
  // avoid comlications and rounding errors if delta_alpha==0
  if(delta_alpha==0) return;
  val=1.0-delta_alpha;
  for(i=0;i<size;i++) out[i]*=val;
  if(norm_flag){
    sum=1;
  }else{
    sum=0;
    for(i=0;i<size;i++) sum+=in[i];
  }
  sum*=(delta_alpha/(double)size);
  for(i=0;i<size;i++) out[i]+=sum;
}

/**************************************/
// multiplication with transposed google matrix
// note delta_alpha = 1-alpha !!
// if norm_flag is set the vector in is assumed to be sum normalized
// which increases the quality of convergence
// NOTE: this routine can be nicely parallized but it is not 
// done here because GTmult is used in the computation of the 
// reduced google matrix were parallization is done outside

void network::GTmult(double delta_alpha,dvec &out,const dvec &in,
		     int norm_flag=1){
  double sum,val;
  int i,a,b,jj;

  // contribution from dangling modes
  // note the modification with respect to GGmult
  // ==> 1/N d e^T
  if(dangling.dim>0){
    sum=0.0;
    for(i=0;i<size;i++){
      out[i]=0;
      sum+=in[i];
    }
    sum/=size;
    for(i=0;i<dangling.dim;i++) out[dangling[i]]+=sum;
  }else{
    for(i=0;i<size;i++){
      out[i]=0;
    }
  }

  //  Computation of out=S^T*in
  for(jj=0;jj<size;jj++){
    a=firstpos[jj];
    b=firstpos[jj+1];
    if(a>=b) continue;
    // note that from[a]=from[i]=jj for a<=i<b
#ifndef USE_PROBS
    sum=0;
    for(i=a;i<b;i++) sum+=in[to[i]];
    out[jj]+=sum/link_num[jj];
#else
    //    for(i=a;i<b;i++) sum+=prob[i]*in[to[i]];
    //    out[jj]+=sum;
    for(i=a;i<b;i++) out[jj]+=prob[i]*in[to[i]];
#endif
  }

  // computation of out=G^T*in, i.e. damping factor contributions
  // avoid comlications and rounding errors if delta_alpha==0
  if(delta_alpha==0) return;
  val=1.0-delta_alpha;
  for(i=0;i<size;i++) out[i]*=val;
  if(norm_flag){
    sum=1;
  }else{
    sum=0;
    for(i=0;i<size;i++) sum+=in[i];
  }
  sum*=(delta_alpha/(double)size);
  for(i=0;i<size;i++) out[i]+=sum;
}

#endif
