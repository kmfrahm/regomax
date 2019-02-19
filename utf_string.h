#include <stdio.h>
#include <string.h>

/****************************************/

inline int strlen_utf(char *s){
  int i,j;

  i=0; j=0;
  while (s[i]) {
    if ((s[i] & 0xC0) != 0x80) j++;
    i++;
  }
  return j;
}

/****************************************/

inline void fprintf_utf(FILE *fp,char *s,int field=0){
  int i;
  int l=strlen_utf(s);
  //  printf("---> l = %d  field = %d\n",l,field);
  
  for(i=l;i<field;i++) fprintf(fp," ");
  fprintf(fp,"%s",s);
}

/****************************************/

inline void printf_utf(char *s,int field=0){
  fprintf_utf(stdout,s,field);
}

/****************************************/

