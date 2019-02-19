// simple test of google matrix property of a (sum of) googlematrix(ces)

#include "net_Greduced_code.h"
#include "net_mattest.h"

int main(int argc,char **argv){
  matrix sum;
  int i;

  for(i=1;i<argc;i++){
    matrix g;
    printf("reading file: %s\n",argv[i]);
    fflush(stdout);
    read_mat(g,argv[i]);
    if(i==1) sum=g; else sum+=g;
  }


  if(test_google_mat(sum)){
       printf("sum is a Google matrix\n");
  }else{
       printf("sum is NOT a Google matrix\n");
  }
}
