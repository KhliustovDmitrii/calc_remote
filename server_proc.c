/*********************************************\
 * SERVER part of calc_remote program
 * 1. Passive socket creation
 * 2. General request processing
 * ---3. CLEAN
 * ------ source
 * ------ out
 * ------ configs
 * ------ progs
 * ---4. RECEIVE
 * ------ source
 * ------ configs
 * ------ progs
 * ---5. COMPILE
 * ---6. RUN
 * ---7. SEND STATUS
 * ---8. SEND
 * ------ out
\*********************************************/

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/signal.h>
#include <sys/wait.h>
#include <netinet/in.h>

#include <netdb.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <signal.h>

#define PORT 5556
#define PATHLN 1024
#define BUF_SIZE 3096

unsigned short portbase = 0;
char buf[BUF_SIZE];
pid_t calc_pid;
extern int errno;

int passiveTCP(int);
int clean_folders();
int compile_program();
int run_program();
int send_status(int);
int send_data(int);
int receive_data(int);

//Passive socket creation routine
int passiveTCP(int qlen)
{
   struct protoent *ppe;
   struct sockaddr_in sin;
   int s, type;
   
   memset(&sin, 0, sizeof(sin));
   sin.sin_family = AF_INET; //Set address family
   sin.sin_addr.s_addr = INADDR_ANY;
   sin.sin_port = htons(PORT);
   
   //Protocol name to number
   if((ppe = getprotobyname("tcp"))==0)
   {
      fprintf(stderr, "Could not get protocol name\n");
      exit(1);
   }
   
   type = SOCK_STREAM;
   
   //Create socket
   s = socket(PF_INET, type, ppe->p_proto);
   if(s<0)
   {
      fprintf(stderr, "Could not create socket\n");
      exit(1);
   }
   
   //Bind socket
   if(bind(s, (struct sockaddr *)&sin, sizeof(sin)) < 0)
   {
      fprintf(stderr, "Could not bind socket\n");
      exit(1);
   }
   
   if(listen(s, qlen) < 0)
   {
      fprintf(stderr, "Unable to listen on port %d\n", PORT);
      exit(1);
   }
   
   printf("Passive socket create - complete\n\n\n");
   return s;
}

//General processing procedure
void TCP_process_request(int fd)
{
   char c;
   int cc, ret_code;
   
   cc = read(fd, &c, 1); //Read first char of request
   
   if(cc < 0)
   {
      fprintf(stderr, "Unable to read request\n");
      exit(1);
   }
   
   switch(c)
   {
   case 'C': //CLEAN
      printf("Trying to clean folders...\n");
      ret_code = clean_folders();
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to clean folders\n");
         exit(1);
      }
      else
         printf("CLEAN folders - complete\n\n\n");
      break;
   case 'c': //COMPILE
      printf("Trying to compile program...\n");
      ret_code = compile_program();
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to compile program\n");
         exit(1);
      }
      else
         printf("COMPILE - complete\n\n\n");
      break;
   case 'r': //RUN
      printf("Trying to run program...\n");
      ret_code = run_program();
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to run program\n");
         exit(1);
      }
      else
         printf("RUN - complete\n\n\n");
      break;
   case 'u': //RECEIVE
      printf("Trying to receive data...\n");
      ret_code = receive_data(fd);
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to receive data\n");
         exit(1);
      }
      else
         printf("DATA receive - complete\n\n\n");
      break;
   case 'd': //SEND
      printf("Trying to send data...\n");
      ret_code = send_data(fd);
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to send data\n");
         exit(1);
      }
      else
         printf("DATA send - complete\n\n\n");
      break;
   case 's': //SEND STATUS
      printf("Trying to send status...\n");
      ret_code = send_status(fd);
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to send status\n");
         exit(1);
      }
      else
         printf("STATUS send - complete\n\n\n");
      break;
   }
}

//Folder cleaning routine
int clean_folders()
{
   DIR *dirp;
   struct dirent *dp;
   char rem_path[PATHLN];
   
   dirp = opendir("sources");
   if (dirp == NULL)
   {
     fprintf(stderr, "Opendir failed on sources \n");
     return -1;
   }
   for (;;)
   {
     dp = readdir(dirp);
     if (dp == NULL)
        break;
     if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
        continue;

     snprintf(rem_path, PATHLN, "sources/%s", dp->d_name);
     remove(rem_path);
   }
   
   dirp = opendir("out");
   if (dirp == NULL)
   {
     fprintf(stderr, "Opendir failed on out \n");
     return -1;
   }
   for (;;)
   {
     dp = readdir(dirp);
     if (dp == NULL)
        break;
     if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
        continue;

     snprintf(rem_path, PATHLN, "out/%s", dp->d_name);
     remove(rem_path);
   }
   
   dirp = opendir("configs");
   if (dirp == NULL)
   {
     fprintf(stderr, "Opendir failed on configs \n");
     return -1;
   }
   for (;;)
   {
     dp = readdir(dirp);
     if (dp == NULL)
        break;
     if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
        continue;

     snprintf(rem_path, PATHLN, "configs/%s", dp->d_name);
     remove(rem_path);
   }
   
   dirp = opendir("program");
   if (dirp == NULL)
   {
     fprintf(stderr, "Opendir failed on program\n");
     return -1;
   }
   for (;;)
   {
     dp = readdir(dirp);
     if (dp == NULL)
        break;
     if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
        continue;

     snprintf(rem_path, PATHLN, "program/%s", dp->d_name);
     remove(rem_path);
   }
   
   dirp = opendir("logs");
   if (dirp == NULL)
   {
     fprintf(stderr, "Opendir failed on logs\n");
     return -1;
   }
   for (;;)
   {
     dp = readdir(dirp);
     if (dp == NULL)
        break;
     if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
        continue;

     snprintf(rem_path, PATHLN, "logs/%s", dp->d_name);
     remove(rem_path);
   }
   return 0;
}

//Program compiling routine
int compile_program()
{
   int ret_code;
   
   ret_code = chdir("program");
   if(ret_code < 0)
   {
      fprintf(stderr, "Unable to change dir to program\n");
      return -1;
   }
   
   ret_code = system("make");
   if(ret_code < 0)
   {
      fprintf(stderr, "Unable to build executable\n");
      return -1;
   }
   
   ret_code = system("mv exec ../");
   if(ret_code < 0)
   {
      fprintf(stderr, "Unable to move executable\n");
      return -1;
   }
   
   ret_code = chdir("../");
   if(ret_code < 0)
   {
      fprintf(stderr, "Unable to change dir to home\n");
      return -1;
   }
   
   return 0;
}

//Program starting routine
int run_program()
{
   int ret_code;
   pid_t child_pid;

   switch(child_pid = fork())
   {
   case -1:
      fprintf(stderr, "Unable to fork() for calc execution\n");
      return -1;
   case 0:
      ret_code = system("./exec sources/data.XYZ out/result.XYZ 20 > logs/log.XYZ");
      if(ret_code < 0)
      {
         fprintf(stderr, "Unable to run program\n");
         return -1;
      }
      exit(1);
      break;
   default:
      calc_pid = child_pid;
   }
  
   printf("Program running with pid %d\n\n\n", calc_pid);
  
   return 0;
}

//Status sending routine
int send_status(int fd)
{
   int ld, flength, pos, i, ret_func;
   off_t fend, fstart;
   
   printf("Checking pid %d\n\n\n", calc_pid);
   
   errno = 0;
   ret_func = kill(calc_pid, 0); //Send zero signal to calc process
   if(ret_func != 0 && errno == ESRCH) //Calcuations finished
   { 
       printf("Pid %d does not exist, calculations finished\n\n\n", calc_pid);
       close(fd);
       return 0;
   }
   
   printf("Pid %d exists, calculations continue\n\n\n", calc_pid);
   
   ld = open("logs/log.XYZ", O_RDONLY); //Try to open log file
   if(ld < 0)
   {
      fprintf(stderr, "Unable to open log file\n");
      return -1;
   }
   
   fend = lseek(ld, 0, SEEK_END);
   fstart = lseek(ld, 0, SEEK_SET);
   
   flength = (int)fend - (int)fstart; //Get file size
   
   if(flength > BUF_SIZE)
      pos = -BUF_SIZE;
   else
      pos = -flength;
      
   lseek(ld, pos, SEEK_END); //Find end of file to send as message
   
   pos*=-1;
   
   if(pos != read(ld, buf, pos))
   {
      fprintf(stderr, "Error reading file\n");
      return -1;
   }
   
   buf[pos-1] = '\0';
   
   (void)write(fd, buf, pos);
   close(fd);
   close(ld);
   
   return 0;
}

//Data sending routine
int send_data(int fd)
{
   FILE *f;
   int str_len;
   
   f = fopen("out/result.XYZ", "r"); 
   
   while(fgets(buf, BUF_SIZE, f)) //Read next line of file
   {
      str_len = strlen(buf);
      (void)write(fd, buf, str_len); //Send to socket
   }
   
   close(fd);
   fclose(f);
   return 0;
}

//Data receiving routine
int receive_data(int fd)
{
   FILE *f;
   int n, len1, len2;
   char *ptr;
   char fname_buf[256];
   char progname[512];
   
   memset(fname_buf, 0, sizeof(fname_buf));
   
   while((n=read(fd, buf, sizeof(buf))) > 0)
   {
      buf[n] = '\0';
      ptr = strstr(buf, "$");
      if(ptr!=NULL) //Separator found
      {
         len1 = strlen(ptr); //Which part of buf is before separator?
         len2 = strlen(buf);
         
         
         strncat(fname_buf, buf, len2-len1);
         
         if(strcmp(fname_buf, "data.XYZ")==0)
         {
            f = fopen("sources/data.XYZ", "wt");
            break;
         }
         if(strcmp(fname_buf, "exec.conf")==0)
         {
            f = fopen("configs/exec.conf", "wt");
            break;
         }
         if(strcmp(fname_buf, "waveform.XYZ")==0)
         {
            f = fopen("configs/waveform.XYZ", "wt");
            break;
         }
         
         sprintf(progname, "program/%s", fname_buf);
         f = fopen(progname, "wt");
         break;
      }
      else //Keep reading file name
      {
         strcat(fname_buf, buf);
      }
   }
   
   ptr = &ptr[1];
   fprintf(f, "%s", ptr);//Write part of data in buf 
   
   while((n=read(fd, buf, BUF_SIZE)) > 0)
   {
      buf[n] = '\0';
      fprintf(f, "%s", buf);
   }
   
   close(fd);
   fclose(f);
   return 0;
}

void reaper(int sig)
{
   int status;
   while(wait3(&status, WNOHANG, (struct rusage *)0)>0)
   {
      
   }
}


//------------------------------- BOTTOM OF THE FILE -------------------------------
