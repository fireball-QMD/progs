#include<sys/types.h>
#include<sys/socket.h>
#include<netinet/in.h>
#include<netinet/tcp.h>
#include<stdio.h>
#include<netdb.h>
#include<errno.h>
#include<signal.h>
#include<stdlib.h>
#include<unistd.h>
#include<inttypes.h>
#include<string.h>

/*   ----- begin machine dependency section ------ */

#ifdef COMPAQ
#define SOC_INIT      soc_init_
#define SOC_SEND      soc_send_
#define SOC_RECV      soc_recv_
#define SOC_FORK      soc_fork_
#define KILL_KID      kill_kid_
#define FORTINT long
#define getsockcast int
#endif

#ifdef HP9000
#define SOC_INIT      soc_init
#define SOC_SEND      soc_send
#define SOC_RECV      soc_recv
#define SOC_FORK      soc_fork
#define KILL_KID      kill_kid
#define FORTINT int
#define getsockcast int
#endif

#if (defined IBM32) || (defined IBM64)
#define SOC_INIT      soc_init
#define SOC_SEND      soc_send
#define SOC_RECV      soc_recv
#define SOC_FORK      soc_fork
#define KILL_KID      kill_kid
#endif
#ifdef IBM32
#define FORTINT int
#define getsockcast size_t
#endif
#ifdef IBM64
#define FORTINT long
#define getsockcast unsigned int
#endif

#ifdef LINUX
#define SOC_INIT      soc_init__
#define SOC_SEND      soc_send__
#define SOC_RECV      soc_recv__
#define SOC_FORK      soc_fork__
#define KILL_KID      kill_kid__
#define FORTINT int
#define getsockcast int
#endif

#ifdef ABSOFT
#define FORTINT int
#define getsockcast int
#endif

#ifdef LINUXIA64
#define SOC_INIT      soc_init_
#define SOC_SEND      soc_send_
#define SOC_RECV      soc_recv_
#define SOC_FORK      soc_fork_
#define KILL_KID      kill_kid_
#define FORTINT intptr_t
#define getsockcast socklen_t
#endif

#if (defined SGI32) || (defined SGI64)
#define SOC_INIT      soc_init_
#define SOC_SEND      soc_send_
#define SOC_RECV      soc_recv_
#define SOC_FORK      soc_fork_
#define KILL_KID      kill_kid_
#define getsockcast int
#endif
#ifdef SGI32
#define FORTINT int
#endif
#ifdef SGI64
#define FORTINT long
#endif

#ifdef SUN
#define SOC_INIT      soc_init_
#define SOC_SEND      soc_send_
#define SOC_RECV      soc_recv_
#define SOC_FORK      soc_fork_
#define KILL_KID      kill_kid_
#define FORTINT int
#define getsockcast int
#endif

#ifdef CRAY
#define SOC_INIT      SOC_INIT
#define SOC_SEND      SOC_SEND
#define SOC_RECV      SOC_RECV
#define SOC_FORK      SOC_FORK
#define KILL_KID      KILL_KID
#define FORTINT long
#define getsockcast int
#endif

/*
 * External stuff
 */

extern int errno;
 
/*
 * Globally accessible values
 */
 
int sock;

/*
 * Prototypes
 */
void SOC_SEND(char *buff,FORTINT *msglen) ;
void SOC_RECV(char *buff,FORTINT *msglen) ;
void SOC_INIT(FORTINT *port) ;
void KILL_KID(int sig) ;
void SOC_FORK(void) ;
