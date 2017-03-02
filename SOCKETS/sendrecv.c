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
 * -----------------------------------------------------------------------
 *           synchronous SEND subroutine with FORTRAN interface 
 * -----------------------------------------------------------------------
 */
void SOC_SEND(char *buff,FORTINT *msglen)
{
int nbytes;
char ack;
nbytes = *msglen;
(void) send( sock, buff, nbytes, 0 );
(void) recv( sock, &ack, 1, 0 );  /* synch acknowledgement */
}
/*
 * -----------------------------------------------------------------------
 *           synchronous RECEIVE subroutine with FORTRAN interface
 * -----------------------------------------------------------------------
 */
void SOC_RECV(char *buff,FORTINT *msglen)
{
int nsent, nbytes;
char ack;
nbytes = *msglen;
while( nbytes > 0 ) {
    nsent = recv( sock, buff, nbytes, 0 );
    buff   += nsent;
    nbytes -= nsent;
  }
(void) send( sock, &ack, 1, 0 );  /* synch acknowledgement */
}
