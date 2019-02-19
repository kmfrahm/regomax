/* some test routines for simple Google matrices */

/*************************************/

double diff_mat(matrix &a,matrix &b){
  if(a.x()!=b.x() || a.y()!=b.y()) return 9999.9999;
  //  if((a.x()!=b.x()) || (a.y()!=b.y())) return 9999.9999;
  int i,j;
  double sum,t;

  sum=0;
  for(i=0;i<a.y();i++) for(j=0;j<a.x();j++){
    t=abs(a(i,j)-b(i,j));
    if(t>sum) sum=t;
  }
  return sum;
}

double sub_diff_mat(matrix &a,matrix &b){
  int dimx,dimy;
  dimx=(a.x()<b.x()?a.x():b.x());
  dimy=(a.y()<b.y()?a.y():b.y());
  int i,j;
  double sum,t;

  sum=0;
  for(i=0;i<dimy;i++) for(j=0;j<dimx;j++){
    t=abs(a(i,j)-b(i,j));
    if(t>sum) sum=t;
  }
  return sum;
}

double rel_diff_mat(matrix &a,matrix &b){
  int dimx,dimy;
  dimx=(a.x()<b.x()?a.x():b.x());
  dimy=(a.y()<b.y()?a.y():b.y());
  int i,j;
  double sum,t,s;

  sum=0;
  for(i=0;i<dimy;i++) for(j=0;j<dimx;j++){
    s=abs(a(i,j)+b(i,j));
    if(s==0) continue;
    t=2*abs(a(i,j)-b(i,j))/s;
    if(t>sum) sum=t;
  }
  return sum;
}

int test_google_mat(matrix &a){
  if(a.x()!=a.y()){
    printf("Test for Google matrix: not a square matrix\n");
    fflush(stdout);
    return 0;
  }

  int i,j,flag=1;
  double sum,diffmax;

  diffmax=0;
  for(i=0;i<a.x();i++){
    sum=0;
    for(j=0;j<a.y();j++){
      sum+=a(j,i);
    }
    sum-=1.0;
    //    printf("Test for Google matrix: column sum [%d] - 1 = %lg\n",i,sum);
    //    fflush(stdout);
    //    if(abs(sum)>1e-12) flag=0;
    sum=abs(sum)/a.y();
    if(diffmax<sum) diffmax=sum;
  }

  for(i=0;i<a.x();i++){
    sum=0;
    for(j=0;j<a.y();j++){
      if(a(j,i)<0){
	sum=abs(a(j,i));
	if(diffmax<sum) diffmax=sum;
	//	printf("Test for Google matrix: negative element: a(%d,%d) = %lg\n",
	//	       j,i,a(j,i));
      //	fflush(stdout);
      //	flag=0;
      }
    }
  }

  printf("diffmax = %24.16lg\n",diffmax);
  fflush(stdout);
  return diffmax<1e-10;
}
