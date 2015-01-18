/* rasmol.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */


/*===========================*/
/*  User Definable Options!  */
/*===========================*/

/* #define IBMPC        */
/* #define MSWIN        */
/* #define APPLEMAC     */
#define X11WIN
#define UNIX

#define DIALBOX
#define SOCKETS
#define TERMIOS
#define PROFILE
#define MITSHM

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

#ifdef IBMPC
#undef THIRTYTWOBIT
#undef SIXTEENBIT
#endif

#if !defined(EIGHTBIT) && !defined(THIRTYTWOBIT) && !defined(SIXTEENBIT)
#define EIGHTBIT
#endif

#ifndef RASMOLDIR
#ifdef IBMPC
#define RASMOLDIR  "C:\\RASWIN\\"
#endif

#if !defined(IBMPC) && !defined(APPLEMAC) && !defined(VMS)
#define RASMOLDIR "/usr/local/lib/rasmol/"
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
#ifndef APPLEMAC
typedef unsigned char Byte;
#endif

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


#if defined(__sgi)
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


#define ItemCount       8
#define AdvPickAtom     0
#define AdvPickNumber   1
#define AdvSelectCount  2
#define AdvName         3
#define AdvIdent        4
#define AdvClass        5
#define AdvImage        6
#define AdvPickCoord    7


void WriteChar( int );
void WriteString( char* );
void RasMolFatalExit( char* );
void AdviseUpdate( int );
void RefreshScreen( void );
void RasMolExit( void );

