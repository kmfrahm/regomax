/* computation of reduced Google matrix */

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

  matrix GR(len,len),Grr(len,len),Gpr(len,len),Gqr(len,len),GI(len,len);
  int n=net.size;
  dvec psiL(n),psiR(n),pg(n),a(n),small_pg(len,1.0),b(len);
  get_first_dot(nodefile,nodefile);


  compute_GR(GR,Grr,Gpr,Gqr,GI,psiL,psiR,pg,net,delta_alpha,node);

  ivec p1(len);
  init_permutation(p1.c,len);

  calc_pagerank_small(small_pg,GR,iprint);
  sprintf(buff,"_small_%s",nodefile);
  print_pagerank(small_pg,net,delta_alpha,buff,print_number,ten_number);
  print_subpagerank(small_pg,net,delta_alpha,buff,p1,na);

  // ------ extra part for Grr+Gqr with diagonal elements of Gqr
  matrix Grr_qr(Grr);
  add_non_diag(Grr_qr,Gqr);
  renorm_google(Grr_qr);
  dvec special_pg(len,1.0);

  calc_pagerank_small(special_pg,Grr_qr,iprint);
  sprintf(buff,"_GrrplusGqr_%s",nodefile);
  print_pagerank(special_pg,net,delta_alpha,buff,print_number,ten_number);
  print_subpagerank(special_pg,net,delta_alpha,buff,p1,na);
  // ------ end of extra part

  net.GTmult(delta_alpha,a,psiL,0);
  project_sub(b,a,node);
  sprintf(buff,"_left_%s",nodefile);
  print_pagerank(psiL,net,delta_alpha,buff,print_number,ten_number);
  print_subpagerank(b,net,delta_alpha,buff,p1,na);

  net.GGmult(delta_alpha,a,psiR,0);
  project_sub(b,a,node);
  sprintf(buff,"_right_%s",nodefile);
  print_pagerank(psiR,net,delta_alpha,buff,print_number,ten_number);
  print_subpagerank(b,net,delta_alpha,buff,p1,na);

  sprintf(buff,"_full_%s",nodefile);
  print_pagerank(pg,net,delta_alpha,buff,print_number,ten_number);
  print_subpagerank(pg,net,delta_alpha,buff,node,na);

  sprintf(buff,"GR_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(GR,buff,na);
  print_mat_dens(GR,buff);
  sprintf(buff,"Grr_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(Grr,buff,na);
  print_mat_dens(Grr,buff);
  sprintf(buff,"Gpr_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(Gpr,buff,na);
  print_mat_dens(Gpr,buff);
  sprintf(buff,"Gqr_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(Gqr,buff,na);
  print_mat_dens(Gqr,buff);
  sprintf(buff,"GI_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(GI,buff,na);
  print_mat_dens(GI,buff);
  sprintf(buff,"GrrplusGqr_%s_%s_%d.dat",net.base_name.c,nodefile,len);
  print_mat(Grr_qr,buff,na);
  print_mat_dens(Grr_qr,buff);

  if(test_google_mat(GR)){
       printf("GR is a Google matrix\n");
  }else{
       printf("GR is NOT a Google matrix\n");
  }

}
