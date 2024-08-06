/*********************************************\
 * CLIENT procedures
 * 1. Connection create
 * 2. CLEAN
 * --- source
 * --- out
 * --- configs
 * --- progs
 * 3. UPLOAD
 * --- source
 * --- configs
 * --- progs
 * 4. COMPILE
 * 5. RUN
 * 6. RECEIVE STATUS
 * 7. DOWNLOAD
 * --- out
\*********************************************/

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <netdb.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>

#ifndef INADDR_NONE
#define INADDR_NONE 0xffffffff
#endif
#define PORT 5556
#define NAME_LENGTH 1024

extern int errno;

//Socket connection routine
int connect_sock(const char *host)
{
   struct hostent *phe;
   struct protoent *ppe;
   struct sockaddr_in sin;
   int s, type;
   
   memset(&sin, 0, sizeof(sin));
   sin.sin_family = AF_INET; //Set address family
   
   sin.sin_port = htons(PORT); //Set port
   
   //Host name to IP
   if(phe=gethostbyname(host))
      memcpy(&sin.sin_addr, phe->h_addr, phe->h_length);
   else 
   {
      if((sin.sin_addr.s_addr = inet_addr(host)) == INADDR_NONE)
      {
         fprintf(stderr, "Could not get host IP\n");
         exit(1);  
      }
   }
   
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
   
   //Connect socket
   errno = 0;
   if(connect(s, (struct sockaddr *)&sin, sizeof(sin))<0)
   {
      fprintf(stderr, "Could not connect socket\n");
      fprintf(stderr, "%s\n", strerror(errno));
      exit(1);
   }
   
   return s;
}

/*Routine to do on server 
 *Controlled by single character c
 *Includes: 
 * CLEAN - c = "C"
 * COMPILE - c = "c"
 * RUN - c = "r"
 */
int TCP_simple(const char *host, char *c)
{
   int s, outchars, chsend;
   char buf[BUF_SIZE];
   
   s = connect_sock(host);
   sprintf(buf, "%s", c);
   outchars = strlen(buf);
   chsend = (int)write(s, buf, outchars);
   
   if(chsend <=0) //Error: no bytes send
      return -1;
   
   close(s);
   return 0;
}

//Status updating routine
int TCP_status_update(const char *host, char *msg)
{
   int s, n, tot_chars, chsend;
   
   s = connect_sock(host);
   chsend = (int)write(s, "s\0", 2); //Send command info to server
   if(chsend <=0) //Error: no bytes send
      return -1;
   
   tot_chars = 0;
   while((n=read(s, msg, MSG_SIZE)) > 0)
      tot_chars+=n;
      
   msg[tot_chars] = '\0';
   
   close(s);
   return tot_chars;
}

//Out downloading routine
int TCP_download(const char *host)
{
   int s, n, chsend, tot_chars;
   FILE *f;
   char buf[BUF_SIZE];
   
   s = connect_sock(host);
   chsend = (int)write(s, "d\0", 2); //Send command info to server
   if(chsend <=0) //Error: no bytes send
      return -1;
      
   f = fopen("calc_result.XYZ", "wt");
   
   tot_chars = 0;
   while((n=read(s, buf, BUF_SIZE)) > 0)
   {
      buf[n] = '\0';
      fprintf(f, "%s", buf);
   }
   
   close(s);
   fclose(f);
   return 0;
}

//Upload routine
int TCP_upload(const char *host)
{
   int s, n, str_len, chsend;
   FILE *f;
   char buf[BUF_SIZE];
   char fname[NAME_LENGTH];
   DIR *dirp;
   struct dirent *dp;
      
   dirp = opendir("upl"); //Open directory with files to upload
   if (dirp == NULL)
     return -1;
     
   for (;;) //Find all files in directory
   {
     dp = readdir(dirp);
     if (dp == NULL) //Have read all files
         break;
     if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0) //Current or parent directory
         continue;
        
     sprintf(fname, "upl/%s", dp->d_name); //Open next file
     f = fopen(fname, "r");
     
     s = connect_sock(host);
     chsend = (int)write(s, "u", 1); //Send command info to server
     if(chsend <=0) //Error: no bytes send
       return -1;
       
     str_len = strlen(dp->d_name); //Send file name
     chsend = (int)write(s, dp->d_name, str_len);
     if(chsend <=0) //Error: no bytes send
       return -1;
       
     printf("%s\n", dp->d_name);
     chsend = (int)write(s, "$", 1); //Send separator
     if(chsend <=0) //Error: no bytes send
       return -1;
     
     while(fgets(buf, BUF_SIZE, f)) //Read next line of file
     {
        str_len = strlen(buf);
        chsend = (int)write(s, buf, str_len); //Send to socket
        if(chsend <=0) //Error: no bytes send
          return -1;
     }
     
     fclose(f);
     close(s);
   }
   
   closedir(dirp);
   return 0;
}

//------------------------------- BOTTOM OF THE FILE -------------------------------
