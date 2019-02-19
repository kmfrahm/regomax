/* test programm for the reduced Google matrix with comparison 
   using full matrix inverse etc. for small networks */

//#define USE_PROBS
#include "net_Greduced_code.h"
#include "net_mattest.h"

int main(){
  char netfile[200],nodefile[200],buff[200],nodefilenames[200],*na;
  double delta_alpha,dlambda;
  int iprint,print_number,ten_number,len,i;
  FILE *fp;

  printf("input-file, 1-alpha, iprint, print_number, ten_number  = ");
  scanf("%s%lf%d%d%d",netfile,&delta_alpha,&iprint,&print_number,&ten_number);
  printf("%s  %lg  %d  %d  %d\n",netfile,delta_alpha,iprint,print_number,ten_number);
  printf("nodefile = "); scanf("%s",nodefile);
  printf("%s\n",nodefile);
  nodefilenames[0]='\0';
  printf("file of node names = "); scanf("%s",nodefilenames);
  printf("%s\n",nodefilenames);
  if(nodefilenames[0]=='\0') na=NULL; else na=nodefilenames;

  fp=fopen(nodefile,"r");
  if(fp==NULL) error("nodefile not found");
  fscanf(fp,"%d",&len);
  ivec node(len);
  for(i=0;i<len;i++){
    fscanf(fp,"%d",&node[i]);
  }
  printf("reading of nodefile finished: len = %d\n",len);
  fflush(stdout);

  network net(netfile);

  dvec::set_size(0);
  matrix GR1(len,len),GR2(len,len);
  matrix Grr(len,len),Gpr(len,len),Gqr(len,len),GI(len,len);
  // this also works: 
  /*
  dvec::set_size(len); dmatrix::set_size(len);
  matrix GR1,GR2,Grr,Gpr,Gqr,GI;
  GR1.print_size("GR1");
  GR2.print_size("GR2");
  Grr.print_size("Grr");
  Gpr.print_size("Gpr");
  Gqr.print_size("Gqr");
  GI.print_size(" GI");
  */
  int n=net.size;
  dvec psiL(n),psiR(n),pg(n);
  get_last_dot(nodefile,nodefile);

  if(net.size<=MAX_GR_SIMPLE){
    compute_GR_simple(GR1,net,delta_alpha,node);
    sprintf(buff,"GR1_%s_%s_%d.dat",net.base_name.c,nodefile,len);
    print_mat(GR1,buff,na);
    print_mat_dens(GR1,buff);
  }
  compute_GR(GR2,Grr,Gpr,Gqr,GI,psiL,psiR,pg,net,delta_alpha,node);

  sprintf(buff,"GR2_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(GR2,buff,na);
  print_mat_dens(GR2,buff);

  if(net.size<=MAX_GR_SIMPLE){
    double diff=diff_mat(GR1,GR2);
    printf("matrix diff = %24.16lg\n",diff);

    if(test_google_mat(GR1)){
      printf("GR1 is a Google matrix\n");
    }else{
      printf("GR1 is NOT a Google matrix\n");
    }
  }

  if(test_google_mat(GR2)){
       printf("GR2 is a Google matrix\n");
  }else{
       printf("GR2 is NOT a Google matrix\n");
  }

}
