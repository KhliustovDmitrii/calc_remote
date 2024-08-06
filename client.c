/*********************************************\
 * CLIENT part of calc_remote program
 * 1. CLEAN
 * --- source
 * --- out
 * --- configs
 * --- progs
 * 2. UPLOAD
 * --- source
 * --- configs
 * --- progs
 * 3. COMPILE
 * 4. RUN
 * 5. SHOW STATUS
 * 6. DOWNLOAD
 * --- out
\*********************************************/

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MSG_SIZE 1024
#define BUF_SIZE 1024

#include "client_proc.c"

int connect_sock(const char *);
int TCP_simple(const char *, char *);
int TCP_status_update(const char *, char *);
int TCP_download(const char *);
int TCP_upload(const char *);


int main(int argc, char **argv)
{
   printf("------------------------------- CLIENT PROGRAM STARTED -------------------------------\n");
   char *host = "localhost";
   char msg[MSG_SIZE];
   int ret_code;
   //+++++++ 1
   ret_code = TCP_simple(host, "C");
   if(ret_code < 0)
   {
      fprintf(stderr, "Could not clean remote folders\n");
      exit(1);
   }
   else
      printf("CLEAN remote folders - complete\n\n\n");
   
   sleep(5);
   
   //+++++++ 2
   ret_code = TCP_upload(host);
   if(ret_code < 0)
   {
      fprintf(stderr, "Could not upload files to server\n");
      exit(1);
   }
   else
      printf("UPLOAD files - complete\n\n\n");
      
   sleep(5);
   
   //+++++++ 3
   ret_code = TCP_simple(host, "c");
   if(ret_code < 0)
   {
      fprintf(stderr, "Could not compile remote program\n");
      exit(1);
   }
   else
      printf("COMPILE program - complete\n\n\n");
      
   sleep(5);
      
   //+++++++ 4
   ret_code = TCP_simple(host, "r");
   if(ret_code < 0)
   {
      fprintf(stderr, "Could not run remote program\n");
      exit(1);
   }
   else
      printf("RUN program - complete\n\n\n");
      
   sleep(3);
      
   //+++++++ 5
   do
   {
      ret_code = TCP_status_update(host, msg);
      if(ret_code < 0)
      {
         fprintf(stderr, "Could not update status\n");
         exit(1);
      }
      else
         printf("%s\n", msg);
         
      sleep(5);
      
   } while(ret_code>0);
   
   
   //+++++++ 6
   ret_code = TCP_download(host);
   if(ret_code < 0)
   {
      fprintf(stderr, "Could not download calc output\n");
      exit(1);
   }
   else
      printf("DOWNLOAD output - complete\n");
      
   sleep(3);
   
   printf("------------------------------- BOTTOM OF THE OUTPUT -------------------------------\n");
   
   exit(0);
}

//------------------------------- BOTTOM OF THE FILE -------------------------------
