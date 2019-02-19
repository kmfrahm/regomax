// limit file size to 1 GB for full reading
#define MAX_FILE_LEN (1<<30)
#define BUFF_LEN (1<<15)


/****************************************/
/* Reading of ascii file in buffer and creating array 
   structure for lines, remove also the ^M characters which are \r in C/C++
 */

char** full_read_ascii_file(const char *name,int &line_len){
  int i,j,len;
  char **a,*buff;
  FILE *fp;

  len=file_size(name);
  if(len<0) error("Input file not found");
  if(len>MAX_FILE_LEN) error("Too large input file");

  //  printf("len = %d\n",len);
  fflush(stdout);

  buff=new char [len+10];
  fp=fopen(name,"r");
  fread(buff,sizeof(char),len,fp);
  fclose(fp);

  if(len>0){
    if(buff[len-1]!='\n'){
      buff[len]='\n';
      len++;
    }
  }else{
    buff[0]='\n';
    len=1;
  }

  line_len=0;
  for(i=0;i<len;i++){
    if(buff[i]=='\n') line_len++;
  }
  a=new char* [line_len+10];
  a[0]=buff; j=0;
  for(i=0;i<len;i++){
    if(buff[i]=='\n'){
      j++;
      buff[i]='\0';
      a[j]=buff+i+1;
    }
    if(buff[i]=='\r'){
      buff[i]='\0';
      if(a[j]==buff+i) buff++;
    }
  }
  //  printf("line_len = %d    j = %d\n",line_len,j);
  buff[len]='\0';

  return a;
}

/****************************************/

void clear_file_buff(char **a){
  delete[] a[0];
  delete[] a;
}

