// transformation to bin format

#include "network_class.h"

int main(int argc,char **argv){
  int i;
  char bin_name[200],*name;

  for(i=1;i<argc;i++){
    network net(argv[i]);

    sprintf(bin_name,"%s.bin",net.base_name.c);

    net.write_binfile(bin_name);
  }
}
