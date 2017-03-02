#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>

#define PORT 9000 // the port client will be connecting to 

#define MAXDATASIZE 100 // max number of bytes we can get at once 

int cclient(void)//char *marker)
{
  int sockfd, numbytes;  
  char host_name[40];
  char *marker_send="go";//marker, 
  char marker_recv[MAXDATASIZE];
  struct hostent *he;
  struct sockaddr_in their_addr; // connector's address information 
  int len=strlen(marker_send);
  FILE *sockInfo;

  // read in the host name to connet to in sockets.optional
  if ((sockInfo=fopen("sockets.optional","r")) == NULL){
    perror("can't find sockets.optional");
    exit(1);}
  fscanf(sockInfo,"%s",host_name);
  fclose(sockInfo);
  printf("server is at: %s\n",host_name);

  if ((he=gethostbyname(host_name)) == NULL) {  // get the host info 
    perror("gethostbyname");
    exit(1);
  }

  their_addr.sin_family = AF_INET;    // host byte order 
  their_addr.sin_port = htons(PORT);  // short, network byte order 
  their_addr.sin_addr = *((struct in_addr *)he->h_addr);
  memset(&(their_addr.sin_zero), '\0', 8);  // zero the rest of the struct 

  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
    perror("socket");
    exit(1);
  }

  //connect to vacancyServer
  if (connect(sockfd, (struct sockaddr *)&their_addr,
	      sizeof(struct sockaddr)) == -1) {
    perror("connect");
    exit(1);
  }

  printf("connected to server\n");

  //transmit number or array
  if ((numbytes=send(sockfd, marker_send, len, 0)) == -1) {
    perror("send");
    exit(1);
  }
  
  printf("sent signal\n");

  if ((numbytes=recv(sockfd, marker_recv, MAXDATASIZE-1, 0)) == -1) {
    perror("recv");
    exit(1);
  }

  marker_recv[numbytes] = '\0';

  //printf("Received: %s\n",marker_recv);

  close(sockfd);

  return 0;
} 
