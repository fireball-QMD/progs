/*
  copyright info: 
  
                              @Copyright 2003
                             Fireball Committee
  Brigham Young University - James P. Lewis, Chair
  Arizona State University - Otto F. Sankey 
  Universidad de Madrid - Jose Ortega 
 
  Other contributors, past and present:
  Auburn University - Jian Jun Dong 
  Arizona State University - Gary B. Adams
  Arizona State University - Kevin Schmidt
  Arizona State University - John Tomfohr
  Lawrence Livermore National Laboratory - Kurt Glaesemann 
  Motorola, Physical Sciences Research Labs - Alex Demkov 
  Motorola, Physical Sciences Research Labs - Jun Wang 
  Ohio University - Dave Drabold
  University of Regensburg - Juergen Fritsch
 
 
  RESTRICTED RIGHTS LEGEND 
  Use, duplication, or disclosure of this software and its documentation
  by the Government is subject to restrictions as set forth in subdivision
  { (b) (3) (ii) } of the Rights in Technical Data and Computer Software 
  clause at 52.227-7013.
*/

#include "mysocks.h"

/*
 * -------------------------------------------------------------
 * soc_INIT: child processes call this to link to master socket
 * -------------------------------------------------------------
 */
void SOC_INIT(FORTINT *port)
{
int on=1;
struct sockaddr_in netsoc, remsoc;
struct hostent *hp;
char my_host[256];
getsockcast len ;

(void) strncpy(my_host, "localhost", 256);
hp = gethostbyname( my_host );         /* we are on same machine */ 

/*
 * Construct remote socket
 */
remsoc.sin_family = hp->h_addrtype;
(void) bcopy(hp->h_addr_list[0], (char *) &remsoc.sin_addr.s_addr, hp->h_length);
remsoc.sin_port = htons((ushort) *port);
/*
 *           create local socket 
 */
sock = socket( AF_INET, SOCK_STREAM, 0 );
(void) setsockopt( sock, IPPROTO_TCP, TCP_NODELAY, (char *) &on, sizeof(on));
/*
 *           construct the socket name
 */
netsoc.sin_family = AF_INET;
/*
 *           fill in sockaddr structure
 *           `hton' stands for host-to-network
 *           `s' indicates (unsigned) short integer data type
 *           htons converts host byte ordering to the 
 *           network byte ordering for ports
 */
netsoc.sin_addr.s_addr = htonl(INADDR_ANY);
netsoc.sin_port = htons((ushort) 0);
(void) bcopy( (char *) hp->h_addr, (char *) &netsoc.sin_addr, hp->h_length );
/*
 * Bind local socket
 */
len = (getsockcast) sizeof(netsoc) ;
(void) bind(sock, (struct sockaddr *) &netsoc, len );
(void) getsockname( sock, (struct sockaddr *) &netsoc, &len );
/*
 *           connect on the named sockets
 */
againcon1:
 if (connect( sock, (struct sockaddr *) &remsoc, sizeof(remsoc) )< 0)
   {
     if (errno == EINTR)
       goto againcon1;
     else
     {
         (void) printf("Connect to socket has failed, error %d\n",errno);
         abort();
     }
   }
}

