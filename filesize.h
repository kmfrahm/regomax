#include <sys/stat.h>

// determines the filesize
// return value = -1 if file not found

int file_size(const char *name){
  struct stat BUF;
  int res;
  
  res=stat(name,&BUF);
  if(res<0) return -1;
  return (int)BUF.st_size;
}
