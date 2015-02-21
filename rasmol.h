/* rasmol.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */


/*===========================*/
/*  User Definable Options!  */
/*===========================*/

/* #define UNIX      */
/* #define APPLEMAC  */
/* #define IBMPC     */

#define MSWIN

#define DIALBOX
#define SOCKETS
#define TERMIOS
#define PROFILE
#define MITSHM
#define DDEIPC

/* #define HAVEZLIB     */
/* #define HAVELIBJPEG  */
/* #define HAVELIBPNG   */

/* #define MMIOLIB      */
/* #define CEXIOLIB     */

#if !defined(EIGHTBIT) && !defined(THIRTYTWOBIT) && !defined(SIXTEENBIT)
/* #define THIRTYTWOBIT */
/* #define SIXTEENBIT   */
/* #define EIGHTBIT     */
#endif


/*========================*/
/*  Default User Options! */
/*========================*/

#if defined(__MINGW32__) && !defined(_WIN32)
#define _WIN32
#endif

#if defined(__unix) || defined(__CYGWIN__)
#ifndef UNIX
#define UNIX
#endif
#undef APPLEMAC
#undef IBMPC
#endif

#ifdef VMS
#undef APPLEMAC
#undef IBMPC
#undef UNIX
#endif

#ifdef macintosh
#ifndef APPLEMAC
#define APPLEMAC
#endif
#undef  IBMPC
#undef  UNIX
#endif

#if defined(_WIN32) && !defined(__CYGWIN__)
#ifndef IBMPC
#define IBMPC
#endif
#undef APPLEMAC
#undef UNIX
#endif

#ifdef __MINGW32__
#ifndef IBMPC
#define IBMPC
#endif
#undef APPLEMAC
#undef DDEIPC
#undef UNIX
#endif

#ifndef IBMPC
#undef MSWIN
#endif

#ifndef RASMOLDIR
#ifdef IBMPC
#define RASMOLDIR  "C:\\RASWIN\\"
#endif

#if !defined(IBMPC) && !defined(APPLEMAC) && !defined(VMS)
#define RASMOLDIR "/usr/local/lib/rasmol/"
#endif
#endif

#if !defined(EIGHTBIT) && !defined(THIRTYTWOBIT) && !defined(SIXTEENBIT)
#ifdef __linux
#define SIXTEENBIT
#else
#define EIGHTBIT
#endif
#endif

#if !defined(INVERT) && defined(IBMPC)
#define INVERT
#endif

#if !defined(TIME) && defined(VMS)
#define TIME
#endif


/*==============================*/
/*  Application-wide Constants  */
/*==============================*/

#ifndef True
#define True  1
#define False 0
#endif

#ifndef PI   /* Avoid Linux Warnings! */
#define PI   3.14159265358979323846
#endif


typedef double Real;

#ifdef __STDC__
typedef signed char Char;
#else
typedef char Char;
#endif

#if defined(_LONGLONG) || defined(LP64) || defined(_LP64)
typedef unsigned int Card;
typedef int Long;
#else
typedef unsigned long Card;
typedef long Long;
#endif

#ifdef EIGHTBIT
typedef unsigned char Pixel;
#else
#ifdef THIRTYTWOBIT
typedef Long Pixel;
#else
typedef short Pixel;
#endif
#endif


#if defined(__sgi) || defined(__GNUC__)
#define UnusedArgument(x)  ((x)=(x))
#else
#define UnusedArgument(x)
#endif


#define Rad2Deg      (180.0/PI)
#define Deg2Rad      (PI/180.0)
#define AbsFun(a)    (((a)<0)? -(a) : (a))
#define MinFun(a,b)  (((a)<(b))? (a) : (b) )
#define MaxFun(a,b)  (((a)>(b))? (a) : (b) )

#if defined(__STDC__) || defined(IBMPC)  || defined(APPLEMAC) || defined(__sgi)
#define ToUpper(x)   (toupper((x)))
#else
#define ToUpper(x)   (islower((x))?toupper((x)):(x))
#endif


#if !defined(IBMPC) || defined(_WIN32)
#ifdef APPLEMAC
#include <MacMemory.h>
#define _fmalloc   NewPtrSys
#define _ffree(x)  DisposePtr((Ptr)(x))
#else
#define _fmalloc   malloc
#define _ffree     free
#endif
#define _fstrnicmp strnicmp
#define _fstrcmp   strcmp
#define _fmemset   memset
#define __huge
#define __far
#endif

#define AdvItemCount    8

#define AdvPickAtom     0
#define AdvPickNumber   1
#define AdvSelectCount  2
#define AdvName         3
#define AdvIdent        4
#define AdvClass        5
#define AdvImage        6
#define AdvPickCoord    7


#ifdef APPLEMAC
#include <MacTypes.h>
void SetFileInfo( char *ptr, OSType appl, OSType type, short icon );
#endif

void WriteChar( int );
void WriteString( char* );
void RasMolFatalExit( char* );
void AdviseUpdate( int );
void RefreshScreen( void );
void RasMolExit( void );

