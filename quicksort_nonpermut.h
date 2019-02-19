/* two quicksort versions not using a permutation vector 
   for use with the network class
*/

/******************************************/
// simple quicksort version, 


void quicksort_up(int *a,int l,int r){
  int i,j;
  int m,t;
  
  m=a[l+(r-l)/2];
  i=l; j=r;
  do{
    while(a[i]<m) i++;
    while(a[j]>m) j--;
    if(i<=j){
      t=a[i]; a[i]=a[j]; a[j]=t;
      i++; j--;
    }
  }while(i<=j);
  if(j-l>r-i){
    if(i<r) quicksort_up(a,i,r);
    if(l<j) quicksort_up(a,l,j);
  }else{
    if(l<j) quicksort_up(a,l,j);
    if(i<r) quicksort_up(a,i,r);
  }
}

/******************************************/

/* quicksort with first key in a-vector and
   second key in b-vector
   ==> will be used to restruct the network class for in and out links
   version with simultaneous permutation for prob array
*/

#ifdef USE_PROBS
#define PROB_DECLARE ,double *prob
#define PROB_ARG ,prob
#else
#define PROB_DECLARE
#define PROB_ARG
#endif


void quicksort_up(int *a,int *b,int l,int r PROB_DECLARE){
#ifdef USE_PROBS
  double temp;
#endif
  int i,j,t;
  int ma,mb,mi;

  while(l<r){
    mi=l+(r-l)/2;
    ma=a[mi]; mb=b[mi];
    i=l; j=r;
    do{
      while((a[i]<ma) || ((a[i]==ma) && (b[i]<mb))) i++;
      while((a[j]>ma) || ((a[j]==ma) && (b[j]>mb))) j--;
      if(i<=j){
	if(i<j){
	  t=a[i]; a[i]=a[j]; a[j]=t;
	  t=b[i]; b[i]=b[j]; b[j]=t;
#ifdef USE_PROBS
	  temp=prob[i]; prob[i]=prob[j]; prob[j]=temp;
#endif
	}
	i++; j--;
      }
    }while(i<=j);
    if(j-l>r-i){
      if(i<r) quicksort_up(a,b,i,r PROB_ARG);
      r=j;
    }else{
      if(l<j) quicksort_up(a,b,l,j PROB_ARG);
      l=i;
    }
  }
}

/******************************************/
