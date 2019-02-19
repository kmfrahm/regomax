/* library of quick sort routines, for double and int and 
   in up and down directions, 
   versions with permutation vectors

   special version for pagerank sort choosing index as second 
   criterium
*/

/*****************************************/
/* small function to initialize a permutation vector, needs to be called 
   in main program before calling a quick-sort function
*/

inline void init_permutation(int *p,int n){
  int i;

  for(i=0;i<n;i++) p[i]=i;
}

/* returns -1 if a1<a2 and +1 if a2>a1 
   if a1==a2 returns -1 if i1>i2 and +1 if i1<i2
*/

inline int index_compare(double a1,int i1,double a2,int i2){
  if(a1<a2) return -1;
  if(a1>a2) return 1;
  if(i1>i2) return -1;
  if(i1<i2) return 1;
  return 0;
}

/*****************************************/
/* a = vector to sorted, l = left boundary, 
   r = right boundary, p = permutation vector 
*/


void quicksort_down(double *a,int l,int r,int *p){
  int i,j,t,mi;
  double m;
  
  //  m=a[p[(l+r)/2]];
  mi=p[l+(r-l)/2];
  m=a[mi];
  i=l; j=r;
  do{
    //    while(a[p[i]]>m) i++;
    //    while(a[p[j]]<m) j--;
    while(index_compare(a[p[i]],p[i],m,mi)>0) i++;
    while(index_compare(a[p[j]],p[j],m,mi)<0) j--;
    if(i<=j){
      t=p[i]; p[i]=p[j]; p[j]=t;
      i++; j--;
    }
  }while(i<=j);
  if(j-l>r-i){
    if(i<r) quicksort_down(a,i,r,p);
    if(l<j) quicksort_down(a,l,j,p);
  }else{
    if(l<j) quicksort_down(a,l,j,p);
    if(i<r) quicksort_down(a,i,r,p);
  }
}

void quicksort_up(double *a,int l,int r,int *p){
  int i,j,t,mi;
  double m;
  
  //  m=a[p[(l+r)/2]];
  mi=p[l+(r-l)/2];
  m=a[mi];
  i=l; j=r;
  do{
    //    while(a[p[i]]<m && i<=r) i++;
    //    while(a[p[j]]>m && j>=l) j--;
    while(index_compare(a[p[i]],p[i],m,mi)<0) i++;
    while(index_compare(a[p[j]],p[j],m,mi)>0) j--;
    if(i<=j){
      t=p[i]; p[i]=p[j]; p[j]=t;
      i++; j--;
    }
  }while(i<=j);
  if(j-l>r-i){
    if(i<r) quicksort_up(a,i,r,p);
    if(l<j) quicksort_up(a,l,j,p);
  }else{
    if(l<j) quicksort_up(a,l,j,p);
    if(i<r) quicksort_up(a,i,r,p);
  }
}

