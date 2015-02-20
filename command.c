/* command.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#include <Types.h>
#include <Errors.h>
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdio.h>

#if !defined(IBMPC) && !defined(VMS) && !defined(APPLEMAC)
#include <pwd.h>
#endif

#define COMMAND
#include "command.h"
#include "tokens.h"
#include "molecule.h"
#include "infile.h"
#include "abstree.h"
#include "transfor.h"
#include "cmndline.h"
#include "render.h"
#include "repres.h"
#include "graphics.h"
#include "pixutils.h"
#include "outfile.h"
#include "scripts.h"
#include "tmesh.h"


#if defined(__sun) && !defined(_XOPEN_SOURCE)
/* SUN doesn't always define popen in stdio! */
extern FILE *popen( const char*, const char* );
#endif

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(ptr=group->alist;ptr;ptr=ptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)


#define IsIdentChar(x)  ((isalnum((int)(x)))||((x)=='_')||((x)=='$'))


#ifndef VMS
#ifdef IBMPC
#define DirChar  '\\'
#else
#define DirChar  '/'
#endif
#endif


#define ErrSyntax        0
#define ErrBigNum        1
#define ErrBadOpt        2
#define ErrParam         3
#define ErrFilNam        4
#define ErrBadLoad       5
#define ErrNotNum        6
#define ErrNotSep        7
#define ErrNotBrac       8
#define ErrNoCol         9
#define ErrColour       10
#define ErrBadArg       11
#define ErrBadExpr      12
#define ErrParen        13
#define ErrScript       14
#define ErrFunc         15
#define ErrSetName      16
#define ErrBadSet       17
#define ErrInScrpt      18
#define ErrOutScrpt     19
#define ErrBadMolDB     20


static char *ErrorMsg[] = {
        "Invalid command syntax",            /* ErrSyntax   */
        "Parameter value too large",         /* ErrBigNum   */
        "Invalid parameter setting",         /* ErrBadOpt   */
        "Invalid parameter name",            /* ErrParam    */
        "Filename string expected",          /* ErrFilNam   */
        "Molecule database loaded",          /* ErrBadLoad  */
        "Integer value expected",            /* ErrNotNum   */
        "Comma separator missing",           /* ErrNotSep   */
        "Close bracket ']' expected",        /* ErrNotBrac  */
        "No colour specified",               /* ErrNoCol    */
        "Unknown or incorrect colour",       /* ErrColour   */
        "Invalid command argument",          /* ErrBadArg   */
        "Syntax error in expression",        /* ErrBadExpr  */
        "Close parenthesis ')' expected",    /* ErrParen    */
        "Script command stack too deep",     /* ErrScript   */
        "Open parenthesis '(' expected",     /* ErrFunc     */
        "Invalid or missing atom set name",  /* ErrSetName  */
        "Not enough memory to define set",   /* ErrBadSet   */
        "Command disabled in script file",   /* ErrInScrpt  */
        "Command invalid outside script",    /* ErrOutScrpt */
        "Molecule database not loaded"       /* ErrBadMolDB */
    };


typedef struct {
        unsigned char red;
        unsigned char grn;
        unsigned char blu;
    } RGBStruct;

static RGBStruct ColourTable[34] = {
    {   0,   0,   0 },  /* Black      */
    {   0,   0, 255 },  /* Blue       */
    { 175, 214, 255 },  /* BlueTint   */
    { 175, 117,  89 },  /* Brown      */
    {   0, 255, 255 },  /* Cyan       */
    { 255, 156,   0 },  /* Gold       */
    { 125, 125, 125 },  /* Gray       */
    {   0, 255,   0 },  /* Green      */
    {  46, 139,  87 },  /* GreenBlue  */
    { 152, 255, 179 },  /* GreenTint  */
    { 255,   0, 101 },  /* HotPink    */
    { 255,   0, 255 },  /* Magenta    */
    { 255, 165,   0 },  /* Orange     */
    { 255, 101, 117 },  /* Pink       */
    { 255, 171, 187 },  /* PinkTint   */
    { 160,  32, 240 },  /* Purple     */
    { 255,   0,   0 },  /* Red        */
    { 255,  69,   0 },  /* RedOrange  */
    {   0, 250, 109 },  /* SeaGreen   */
    {  58, 144, 255 },  /* SkyBlue    */
    { 238, 130, 238 },  /* Violet     */
    { 255, 255, 255 },  /* White      */
    { 255, 255,   0 },  /* Yellow     */
    { 246, 246, 117 }   /* YellowTint */
        };


typedef struct _HlpEntry {
                struct _HlpEntry __far *next;
                struct _HlpEntry __far *info;
                char __far *keyword;
                Long fpos;
                } HlpEntry;

#define HelpPool   16
static char *HelpFileName;
static char HelpFileBuf[80];
static HlpEntry __far *FreeInfo;
static HlpEntry __far *HelpInfo;


static char ResidueChar[29] = {
        'A', 'G', 'L', 'S', 'V', 'T', 'K', 'D', 'I', 'N',
        'E', 'P', 'R', 'F', 'Q', 'Y', 'H', 'C', 'M', 'W',
        'B', 'Z', '*', 'P',
        'A', 'C', 'G', 'T',
        'U'
    };


#define STACKSIZE  10
static char *NameStack[STACKSIZE];
static FILE *FileStack[STACKSIZE];
static int LineStack[STACKSIZE];

static int TokenLength;
static Long TokenValue;
static char TokenIdent[128];
static char *TokenStart;
static char *TokenPtr;
static int CurToken;


static int RVal, GVal, BVal;
static int SeqFormat;


/*=======================*/
/*  Function Prototypes  */
/*=======================*/

int ExecuteCommand( void );
int ExecuteIPCCommand( char __huge* );
static int PrefixString( char __far*, char __far* );



static void CommandError( char *error )
{
    register char *ptr;
    char buffer[40];

    if( TokenPtr )
    {   if( FileDepth > -1 )
        {   InvalidateCmndLine();
            WriteString(CurLine);
            WriteChar('\n');
        } else WriteString("        ");

        for( ptr=CurLine; ptr<TokenStart; ptr++ )
            WriteChar(' ');
        WriteString("^\n");
    }

    if( FileDepth > -1 )
    {   if( LineStack[FileDepth] )
        {   if( NameStack[FileDepth] )
            {   WriteChar('"');
                WriteString(NameStack[FileDepth]);
                WriteString("\",");
            }
            sprintf(buffer,"line %d: ",LineStack[FileDepth]);
            WriteString(buffer);
        } else
        {   WriteString(NameStack[FileDepth]);
            WriteString(": ");
        }
    }

    if( error )
    {   WriteString(error);
        WriteString("!\n");
    }
    CommandActive = False;
    CurToken = 0;
}



/*==========================*/
/*  File Handling Services  */
/*==========================*/

#ifdef IBMPC
static char *ProcessFileName( char *name )
{
    register char *ptr;

    while( *name==' ' )
        name++;

    ptr = DataFileName;
    while( *name )
        *ptr++ = *name++;

    /* Strip trailing spaces! */
    while( (ptr!=DataFileName) && (ptr[-1]==' ') )
        ptr--;
    *ptr = '\0';
    return ptr;
}
#endif

#ifdef APPLEMAC
static char *ProcessFileName( char *name )
{
    register char *ptr;

    while( *name==' ' )
        name++;

    ptr = DataFileName;
    while( *name )
        *ptr++ = *name++;

    /* Strip trailing spaces! */
    while( (ptr!=DataFileName) && (ptr[-1]==' ') )
        ptr--;
    *ptr = '\0';
    return ptr;
}
#endif

#ifdef VMS 
static char *ProcessFileName( char *name )
{
    register char *ptr;

    while( *name==' ' )
        name++;

    ptr = DataFileName;
    while( *name && (*name!=' ') )
    {   *ptr++ = ToUpper(*name);
        name++;
    }
    *ptr = '\0';
    return ptr;
}
#endif


#if !defined(IBMPC) && !defined(APPLEMAC) && !defined(VMS)
static int IsSecure( int ch )
{
    switch( ch )
    {   /* Dangerous characters in UNIX "popen"!  */
        case('<'):  case('>'):  case('('):  case(')'):
        case('{'):  case('}'):  case('['):  case(']'):
        case('\''): case(';'):  case('|'):  case('&'):
        case('"'):  case(' '):  case('\t'): case('\n'):
            return False;
    }
    return True;
}


static char *ProcessFileName( char *name )
{
    register struct passwd *entry;
    register char *temp;
    register char *ptr;
    char username[64];

    while( *name==' ' )
        name++;

    /* Perform filename globbing */
    if( *name=='~' )
    {   ptr = username;  name++;
        while( *name && (*name!=' ') && (*name!='/') )
            *ptr++ = *name++;
        *ptr = '\0';

        ptr = DataFileName;
        if( *username )
        {   if( (entry=getpwnam(username)) )
            {   temp = entry->pw_dir;
                endpwent();
            } else /* Unknown user! */
            {   temp = username;
                *ptr++ = '~';
            }

        } else if( !(temp=(char*)getenv("HOME")) )
            temp = ".";

        while( *temp )
            *ptr++ = *temp++;
    } else ptr = DataFileName;

    /* Strip dubious characters! */
    while( *name && (*name!=' ') )
        if( IsSecure(*name) )
        {   *ptr++ = *name++;
        } else name++;
    *ptr = '\0';
    return ptr;
}
#endif


#ifdef UNIX

#define MaxFileExt  4
/* UNIX Compressed Filename extensions! */
static char *FileExt[MaxFileExt] = { "", ".Z", ".gz", ".z" };

static FILE *OpenDataFile( char *begin, char *end )
{
    register char *src, *dst;
    register FILE *fp;
    register int i;
    
    for( i=0; i<MaxFileExt; i++ )
    {   dst = end; src = FileExt[i];
        while( (*dst++ = *src++) );
        if( (fp=fopen(begin,"rb")) )
            break;
    }
    fp = fopen(begin,"rb");
    *end = '\0';
    return fp;
}
#else /* !defined(UNIX) */

static FILE *OpenDataFile( char *begin, char *end )
{
    register FILE *fp;

    /* Avoid compiler warning! */
    UnusedArgument(end);

    fp = fopen(begin,"rb");
    return fp;
}
#endif


int ProcessFile( int format, int info, FILE *fp )
{
    register int done;

    switch( format )
    {   case(FormatPDB):      done = LoadPDBMolecule(fp,False);  break;
        case(FormatNMRPDB):   done = LoadPDBMolecule(fp,True);   break;
        case(FormatMacroMod): done = LoadMacroModelMolecule(fp); break;
        case(FormatAlchemy):  done = LoadAlchemyMolecule(fp);    break;
        case(FormatCharmm):   done = LoadCharmmMolecule(fp);     break;
        case(FormatBiosym):   done = LoadBiosymMolecule(fp);     break;
        case(FormatMOPAC):    done = LoadMOPACMolecule(fp);      break;
        case(FormatSHELX):    done = LoadSHELXMolecule(fp);      break;
        case(FormatMol2):     done = LoadMol2Molecule(fp);       break;
        case(FormatFDAT):     done = LoadFDATMolecule(fp);       break;
        case(FormatMDL):      done = LoadMDLMolecule(fp);        break;
        case(FormatXYZ):      done = LoadXYZMolecule(fp);        break;
#ifdef CEXIOLIB
        case(FormatCEX):      done = LoadCEXMolecule(fp);        break;
#endif
        default:              done = False;
    }

    if( !done )
    {   return( False );
    } else if( !Database )
        return( True );

    if( info )
        DescribeMolecule();
    DataFileFormat = format;
    AdviseUpdate(AdvName);
    AdviseUpdate(AdvClass);
    AdviseUpdate(AdvIdent);

    if( Interactive )
        FetchEvent(False);

    ReDrawFlag |= RFInitial;
    if( CalcBondsFlag )
    {   if( Info.bondcount < (MainAtomCount+HetaAtomCount)-Info.chaincount )
        {   if( MainAtomCount+HetaAtomCount > 255 )
            {   CreateMoleculeBonds(info,False);
            } else CreateMoleculeBonds(info,True);
        }
    }

    /* Explicit Hydrogen Bonds! */
    if( Info.hbondcount > 0 )
        SetHBondStatus(True,True,0);

    InitialTransform();

    VoxelsClean = False;
    ApplyTransform();
    return True;
}


int FetchFile( int format, int info, char *name )
{
#ifndef APPLEMAC
#ifdef UNIX
    register int comp;
#endif /* UNIX */
    register char *src,*dst;
    register char *tmp;
    char buffer[128];
#endif /* APPLEMAC */

    register int done;
    register FILE *fp;

    DataFileFormat = 0;
    name = ProcessFileName(name);
    fp = OpenDataFile(DataFileName,name);

#ifndef APPLEMAC
    /* Search for directory specification! */
#ifndef VMS
    if( !fp )
    {   src = DataFileName;
        while( *src && (*src!=DirChar) )
            src++;
        done = !(*src);
    } else done = False;
#else
    done = True;
#endif

    /* Try using a default file path! */
    if( !fp && done )
    {   switch( format )
        {   case(FormatNMRPDB):
            case(FormatPDB):     src = (char*)getenv("RASMOLPDBPATH");  break;
            case(FormatMol2):    src = (char*)getenv("RASMOLMOL2PATH"); break;
            case(FormatMMDB):    src = (char*)getenv("RASMOLMMDBPATH"); break;
            case(FormatAlchemy): src = (char*)getenv("RASMOLMOLPATH");  break;
            case(FormatMDL):     src = (char*)getenv("RASMOLMDLPATH");  break;
            case(FormatXYZ):     src = (char*)getenv("RASMOLXYZPATH");  break;
            default:             src = (char*)0;
        }

        if( src && *src )
        {   
#ifdef VMS
            dst = buffer;
            while( *src ) *dst++ = *src++;

            tmp = DataFileName;
            while( *dst = *tmp++ ) dst++;
            if( fp = OpenDataFile(buffer,dst) )
                strcpy(DataFileName,buffer);
#else
            while( *src )
            {   dst = buffer;
                while( *src && (*src!=':') ) 
                    *dst++ = *src++;
                if( *src == ':' ) 
                    src++;

                if( dst != buffer )
                {   if( *(dst-1) != DirChar )
                        *dst++ = DirChar;
                    tmp = DataFileName;
                    while( (*dst = *tmp++) ) dst++;
                    if( (fp = OpenDataFile(buffer,dst)) )
                    {   strcpy(DataFileName,buffer);
                        break;
                    }
                }
            }
#endif
        }
    }
#endif /* APPLEMAC */


#ifdef CEXIOLIB
    if( !fp && (format==FormatCEX) )
    {   if( ProcessFile(format,info,fp) )
            return True;
    }
#endif

    if( !fp )
    {   *name = '\0';
        InvalidateCmndLine();
        WriteString("Error: File '");
        WriteString(DataFileName);
        WriteString("' not found!\n\n");
        return( False );
    }

#ifdef UNIX
    done = getc(fp);
    if( done == 0x1f )
    {   done = getc(fp);
        fclose(fp);

        if( done == 0x9d )
        {   /* Should #include <signal.h> and trap "" SIGPIPE */
            sprintf(buffer,"trap \"\" 13; uncompress -c %s 2> /dev/null\n",
                                                              DataFileName);
        } else if( done == 0x8b )
        {   /* Should #include <signal.h> and trap "" SIGPIPE */
            sprintf(buffer,"trap \"\" 13; gzip -cdq %s 2> /dev/null\n",
                                                          DataFileName);
        } else /* bad magic number! */
        {   InvalidateCmndLine();
            WriteString("Error: Unrecognised compression format!\n\n");
            return( False );
        }
   
        comp = True;
        fp = popen(buffer,"r");
        if( !fp )  
        {   InvalidateCmndLine();
            WriteString("Error: Unable to decompress file!\n\n");
            return( False );
        }
    } else /* Uncompressed! */
    {   ungetc(done,fp);
        comp = False;
    }
#endif

    done = ProcessFile(format,info,fp);

#ifdef UNIX
    if( comp )
    {   if( pclose(fp) )
        {   InvalidateCmndLine();
            WriteString("Error: Unable to decompress file!\n\n");
            return(False);
        }
    } else fclose(fp);
#else /* !defined(UNIX) */
    fclose(fp);
#endif
    return done;
}


void LoadScriptFile( FILE *fp, char *name )
{
    register char *ptr;
    register int ch,len;
    register int stat;

    if( fp )
    {   FileDepth++;
        len = strlen(name)+1;
        ptr = (char*)malloc(len);
        memcpy(ptr,name,len);
        NameStack[FileDepth] = ptr;
        FileStack[FileDepth] = fp;
        LineStack[FileDepth] = 0;

        do {
            len = 0;
            ch = getc(fp);
            while( (ch!='\n') && (ch!='\r') && (ch!=EOF) )
            {   if( len<MAXBUFFLEN )
                    CurLine[len++] = ch;
                ch = getc(fp);
            }

            if( ch == '\r' )
            {   ch = getc(fp);
                if( ch != '\n' )
                    ungetc(ch,fp);
            }

            LineStack[FileDepth]++;
            if( len<MAXBUFFLEN )
            {   CurLine[len] = '\0';
                stat = ExecuteCommand();
                if( stat )
                {   if( stat == QuitTok )
                    {   while( FileDepth >= 0 )
                        {   fclose(FileStack[FileDepth]);
                            free(NameStack[FileDepth]);
                            FileDepth--;
                        }
                        RasMolExit();
                    } else /* ExitTok */
                        break;
                } else if( IsPaused )
                    return;

            } else CommandError("Script command line too long");
        } while( ch != EOF );
        free(NameStack[FileDepth]);
        fclose( fp );
        FileDepth--;
    } else 
    {   CommandError( (char*)0 );
        WriteString("Cannot open script file '");
        WriteString(name);  WriteString("'\n");
    }
}



/*====================================*/
/*  Command Line On-Line Help System  */
/*====================================*/

static int PrefixString( char __far *str1, char __far *str2 )
{
    while( *str1 )
        if( *str1++ != *str2++ )
            return False;
    return True;
}


static HlpEntry __far *EnterHelpInfo( char *text )
{
    register HlpEntry __far * __far *tmp;
    register HlpEntry __far *ptr;
    register int res,len,i;
    register char ch;
    char keyword[32];

    ptr = (HlpEntry __far*)0;
    while( *text && (*text!='\n') )
    {   while( *text && (*text!='\n') && (*text==' ') )
            text++;

        len = 0;
        while( *text && (*text!='\n') && (*text!=' ') )
            if( len<31 )
            {   ch = *text++;
                keyword[len++] = ToUpper(ch);
            } else text++;
        keyword[len]='\0';

        if( ptr )
        {   tmp = &ptr->info;
            ptr = (HlpEntry __far*)0;
        } else tmp = &HelpInfo;

        while( *tmp )
        {   res = _fstrcmp(keyword,(*tmp)->keyword);
            if( res==0 ) /* Exact Match */
            {   ptr = *tmp;
                break;
            } else if( res<0 )
                break;
            tmp = &(*tmp)->next;
        }

        if( !ptr )
        {   if( !FreeInfo )
            {   ptr = (HlpEntry __far*)_fmalloc(HelpPool*sizeof(HlpEntry));
                if( !ptr ) 
                    RasMolFatalExit("Command Error: Insufficient memory!");
                for( i=1; i<HelpPool; i++ )
                {   ptr->next = FreeInfo;
                    FreeInfo = ptr++;
                }
            } else
            {   ptr = FreeInfo;
                FreeInfo = ptr->next;
            }

            ptr->keyword = (char __far*)_fmalloc(len+1);
            for( i=0; i<=len; i++ )
                ptr->keyword[i] = keyword[i];

            ptr->info = (HlpEntry __far*)0;
            ptr->next = *tmp;
            ptr->fpos = 0;
            *tmp = ptr;
        }
    }
    return ptr;
}


static void InitHelpFile( void )
{
    register char *src,*dst;
    register HlpEntry __far *fix;
    register HlpEntry __far *ptr;
    register FILE *fp;
    register Long pos;
    char buffer[82];

    HelpFileName = "rasmol.hlp";
    fp=fopen(HelpFileName,"rb");

    if( !fp )
    {   src = (char*)getenv("RASMOLPATH");
        if( src )
        {   HelpFileName = dst = HelpFileBuf; 
            while( *src )
                *dst++ = *src++;
#ifndef VMS
            if( (dst!=HelpFileBuf) && (*(dst-1)!=DirChar) )
                *dst++ = DirChar;
#endif

            strcpy(dst,"rasmol.hlp");
            fp = fopen(HelpFileName,"rb");
        }
    }

#ifdef RASMOLDIR
    if( !fp )
    {   src = RASMOLDIR;
        HelpFileName = dst = HelpFileBuf;
        while( *src )
            *dst++ = *src++;
#ifndef VMS
        if( (dst!=HelpFileBuf) && (*(dst-1)!=DirChar) )
            *dst++ = DirChar;
#endif

        src = "rasmol.hlp"; 
        while( (*dst++ = *src++) );
        fp = fopen(HelpFileName,"rb");
    }
#endif

    if( !fp )
    {   InvalidateCmndLine();
        WriteString("Unable to find RasMol help file!\n");
        HelpFileName = (char*)0;
        return;
    }

    pos = 0;
    fgets(buffer,80,fp);
    while( !feof(fp) )
    {    fix = (HlpEntry __far*)0;
         while( *buffer=='?' )
         {   ptr = EnterHelpInfo(buffer+1);
             if( ptr )
             {   ptr->info = fix;
                 fix = ptr;
             }

             pos = ftell(fp);
             if( !fgets(buffer,80,fp) )
                 break;
         }

         while( fix )
         {   ptr = fix->info;
             fix->info = (HlpEntry __far*)0;
             fix->fpos = pos;
             fix = ptr;
         }

         while( fgets(buffer,80,fp) )
             if( *buffer=='?' )
                 break;
    }
    fclose(fp);
}


static void FindHelpInfo( void )
{
    register HlpEntry __far * __far *tmp;
    register HlpEntry __far *ptr;
    register int res,len;
    register Long pos;
    register FILE *fp;
    register char ch;
    char keyword[32];
    char buffer[82];

    while( *TokenPtr && (*TokenPtr==' ') )
        TokenPtr++;

    if( *TokenPtr )
    {   ptr = (HlpEntry __far*)0;
        do {
            len = 0;
            while( *TokenPtr && (*TokenPtr!=' ') )
                if( len<31 )
                {   ch = *TokenPtr++;
                    keyword[len++] = ToUpper(ch);
                } else TokenPtr++;
            keyword[len]='\0';

            if( ptr )
            {   tmp = &ptr->info;
                ptr = (HlpEntry __far*)0;
            } else tmp = &HelpInfo;

            while( *tmp )
            {   res = _fstrcmp(keyword,(*tmp)->keyword);
                if( res<0 )
                {   if( PrefixString(keyword,(*tmp)->keyword) )
                    {   ptr = *tmp;
                        if( ptr->next && 
                            PrefixString(keyword,ptr->next->keyword) )
                        {   InvalidateCmndLine();
                            WriteString("Ambiguous help topic requested!\n");
                            return;
                        } else break;
                    } else break;
                } else if( res==0 ) 
                {   ptr = *tmp;
                    break;
                }
                tmp = &(*tmp)->next;
            }

            while( *TokenPtr && (*TokenPtr==' ') )
                TokenPtr++;
        } while( *TokenPtr && ptr );

        if( !ptr || !ptr->fpos )
        {   InvalidateCmndLine();
            WriteString("No available help on requested topic!\n");
            return;
        } else pos=ptr->fpos;
    } else pos=0;


    if( !(fp=fopen(HelpFileName,"rb")) )
        RasMolFatalExit("Command Error: Unable to reopen help file!");

    InvalidateCmndLine();

    fseek(fp,pos,0);
    while( fgets(buffer,80,fp) )
        if( *buffer!='?' )
        {   WriteString(buffer);
        } else break;
    fclose(fp);
}



/*=================================*/
/*  Command Line Lexical Analysis  */
/*=================================*/

static int FetchToken( void )
{
    register char ch;

    CurToken = 0;
    ch = *TokenPtr++;
    while( ch && (ch!='#') )
    {   if( isspace((int)ch) )
        {   ch = *TokenPtr++;
            continue;
        }

        TokenStart = TokenPtr-1;
        if( isalpha((int)ch) )
        {   TokenLength = 1;
            *TokenIdent = ToUpper(ch);
            while( IsIdentChar(*TokenPtr) && (TokenLength<32) )
            {   ch = *TokenPtr++;
                TokenIdent[TokenLength++] = ToUpper(ch);
            }
            if( TokenLength==32 )
            {   CommandError("Identifier too long");
                return(0);
            } else TokenIdent[TokenLength] = '\0';
            return( CurToken = LookUpKeyword(TokenIdent) );

        } else if( isdigit((int)ch) )
        {   TokenValue = ch-'0';
            while( isdigit((int)*TokenPtr) )
                TokenValue = 10*TokenValue + (*TokenPtr++)-'0';
            return( CurToken = NumberTok );

        } else if( (ch=='\'') || (ch=='\"') || (ch=='`') )
        {   TokenLength = 0;
            while( *TokenPtr && (TokenLength<128) && (*TokenPtr!=ch) )
                TokenIdent[TokenLength++] = *TokenPtr++;

            if( ch != *TokenPtr )
            {   if( *TokenPtr )
                {   CommandError("String constant unterminated");
                } else CommandError("String constant too long");
                return( 0 );
            } else TokenPtr++;

            TokenIdent[TokenLength]='\0';
            return( CurToken = StringTok );
        } else if( ispunct((int)ch) )
            return( CurToken = ch );

        ch = *TokenPtr++;
    }
    TokenPtr--;
    return 0;
}


static int NextIf( int tok, int err )
{
    if( FetchToken() != tok )
    {   CommandError(ErrorMsg[err]);
        return True;
    } else return False;
}


static void FetchFloat( Long value, int scale )
{
    register int count;
    register int mant;

    if( !value && !isdigit((int)*TokenPtr) )
    {   CommandError("Invalid floating point number");
        TokenValue = 0;
        return;
    }

    mant = 0;
    count = 1;
    while( isdigit((int)*TokenPtr) )
    {   if( count < scale )
        {   mant = 10*mant + (*TokenPtr-'0');
            count *= 10;
        }
        TokenPtr++;
    }

    mant = (scale*mant)/count;
    TokenValue = value*scale + mant;
}


static int ParseColour( void )
{
    register RGBStruct *rgb;

    if( IsColourToken(CurToken) )
    {   rgb = ColourTable + Token2Colour(CurToken);
        RVal = rgb->red;
        GVal = rgb->grn;
        BVal = rgb->blu;
        return True;

    } else if( CurToken == '[' )
    {   RVal = GVal = BVal = 0;

        if( NextIf(NumberTok,ErrNotNum) ) 
        {   return False;
        } else if( TokenValue>255 )
        {   CommandError(ErrorMsg[ErrBigNum]); 
            return False;
        } else RVal = (int)TokenValue;

        if( NextIf(',',ErrNotSep) ) 
            return False;

        if( NextIf(NumberTok,ErrNotNum) ) 
        {   return False;
        } else if( TokenValue>255 )
        {   CommandError(ErrorMsg[ErrBigNum]); 
            return False;
        } else GVal = (int)TokenValue;

        if( NextIf(',',ErrNotSep) ) 
            return False;

        if( NextIf(NumberTok,ErrNotNum) ) 
        {   return False;
        } else if( TokenValue>255 )
        {   CommandError(ErrorMsg[ErrBigNum]);
            return False;
        } else BVal = (int)TokenValue;

        return !NextIf(']',ErrNotBrac);

    } else if( !CurToken && (*TokenPtr=='#') )
    {   RVal = 0;
        GVal = 0;
        BVal = 0;
    
    } else if( CurToken == IdentTok )
        if( Interactive )
            return LookUpColour(TokenIdent,&RVal,&GVal,&BVal);

    return False;
}


static Expr *ParseRange( int neg )
{
    register Expr *tmp1,*tmp2;
    register char ch;

    tmp1 = AllocateNode();
    tmp1->type = OpLftProp|OpRgtVal;
    tmp1->rgt.val = neg? -(int)TokenValue : (int)TokenValue;
    tmp1->lft.val = PropResId;

    if( *TokenPtr == '-' )
    {   TokenPtr++;
        neg = (*TokenPtr=='-');
        if( neg ) TokenPtr++;
        FetchToken();

        if( CurToken != NumberTok )
        {   CommandError(ErrorMsg[ErrNotNum]);
            DeAllocateExpr( tmp1 );
            return( (Expr*)0 );
        }

        tmp1->type |= OpMoreEq;
        tmp2 = AllocateNode();
        tmp2->rgt.ptr = tmp1;
        tmp2->type = OpAnd;

        tmp1 = AllocateNode();
        tmp1->type = OpLftProp|OpRgtVal|OpLessEq;
        tmp1->rgt.val = neg? -(int)TokenValue : (int)TokenValue;
        tmp1->lft.val = PropResId;
        tmp2->lft.ptr = tmp1;
        tmp1 = tmp2;
    } else tmp1->type |= OpEqual;

    if( *TokenPtr == ':' )
        TokenPtr++;

    ch = *TokenPtr;
    if( isalnum((int)ch) )
    {   ch = ToUpper(ch);
        TokenPtr++;

        tmp2 = AllocateNode();
        tmp2->type = OpAnd;
        tmp2->rgt.ptr = tmp1;

        tmp1 = AllocateNode();
        tmp1->type = OpEqual | OpLftProp | OpRgtVal;
        tmp1->lft.val = PropChain;               
        tmp1->rgt.val = ch;

        tmp2->lft.ptr = tmp1;
        tmp1 = tmp2;
    } else if( (ch=='?') || (ch=='%') || (ch=='*') )
        TokenPtr++;

    FetchToken();
    return tmp1;
}


static Expr *ParseExpression( int level )
{
    register Expr *tmp1,*tmp2;
    register int done, pred;
    register int neg;

    switch( level )
    {    case(0): /* Disjunctions */
                  tmp1 = ParseExpression(1);
                  while( (CurToken==OrTok) || (CurToken=='|') ||
                         (CurToken==',') )
                  {   if( CurToken=='|' )
                      {   if( FetchToken()=='|' )
                              FetchToken();
                      } else FetchToken();

                      tmp2 = AllocateNode();
                      tmp2->type = OpOr;
                      tmp2->lft.ptr = tmp1;
                      tmp2->rgt.ptr = (Expr*)0;
                      if( !(tmp1=ParseExpression(1)) )
                      {   DeAllocateExpr(tmp2);
                          return( tmp1 );
                      }
                      tmp2->rgt.ptr = tmp1;
                      tmp1 = tmp2;
                  }
                  return( tmp1 );

         case(1): /* Conjunctions */
                  tmp1 = ParseExpression(2);
                  while( (CurToken==AndTok) || (CurToken=='&') )
                  {   if( CurToken=='&' )
                      {   if( FetchToken()=='&' )
                              FetchToken();
                      } else FetchToken();

                      tmp2 = AllocateNode();
                      tmp2->type = OpAnd;
                      tmp2->lft.ptr = tmp1;
                      tmp2->rgt.ptr = (Expr*)0;
                      if( !(tmp1=ParseExpression(2)) )
                      {   DeAllocateExpr(tmp2);
                          return( tmp1 );
                      }
                      tmp2->rgt.ptr = tmp1;
                      tmp1 = tmp2;
                  }
                  return( tmp1 );

         case(2): /* Primitives */
                  if( IsPredTok(CurToken) || (CurToken==BackboneTok) )
                  {   /* Avoid compiler warning! */
                      pred = 0;
                      switch( CurToken )
                      {   case(HelixTok):    if( Info.helixcount < 0 )
                                                 DetermineStructure(False);
                                             pred = PredHelix;
                                             break;
                          case(SheetTok):    if( Info.laddercount < 0 )
                                                 DetermineStructure(False);
                                             pred = PredSheet;
                                             break;
                          case(TurnTok):     if( Info.turncount < 0 )
                                                 DetermineStructure(False);
                                             pred = PredTurn;
                                             break;
                          case(CystineTok):  if( Info.ssbondcount < 0 )
                                                 FindDisulphideBridges();
                                             pred = PredCystine;     
                                             break;
                          case(BackboneTok): pred = PredMainChain;   break;
                          case(SelectedTok): pred = PropSelect;      break;
                          default:  pred = PredAbsChr(PredTokOrd(CurToken));
                      }

                      tmp1 = AllocateNode();
                      tmp1->type = OpConst|OpLftProp|OpRgtVal;
                      tmp1->lft.val = pred;
                      FetchToken();
                      return( tmp1 );

                  } else if( IsPropTok(CurToken) )
                  {   /* Avoid compiler warning! */
                      pred = 0;

                      tmp1 = AllocateNode();
                      tmp1->type = OpLftProp|OpRgtVal;
                      switch( CurToken )
                      {   case(TemperatureTok): pred = PropTemp;    break;
                          case(RadiusTok):      pred = PropRad;     break;
                          case(AtomNoTok):      pred = PropIdent;   break;
                          case(ElemNoTok):      pred = PropElemNo;  break;
                          case(ResNoTok):       pred = PropResId;   break;
                          case(ModelTok):       pred = PropModel;   break;
                      }
                      tmp1->lft.val = pred;

                      FetchToken();
                      if( CurToken=='=' )
                      {   tmp1->type |= OpEqual;
                          if( FetchToken()=='=' )
                              FetchToken();
                      } else if( CurToken=='<' )
                      {   FetchToken();
                          if( CurToken=='>' )
                          {   tmp1->type |= OpNotEq;
                              FetchToken();
                          } else if( CurToken=='=' )
                          {   tmp1->type |= OpLessEq;
                              FetchToken();
                          } else tmp1->type |= OpLess;
                      } else if( CurToken=='>' )
                      {   if( FetchToken()=='=' )
                          {   tmp1->type |= OpMoreEq;
                              FetchToken();
                          } else tmp1->type |= OpMore;
                      } else if( (CurToken=='!') || (CurToken=='/') )
                      {   if( NextIf('=',ErrBadExpr) )
                          {   DeAllocateExpr( tmp1 );
                              return( (Expr*)0 );
                          } else tmp1->type |= OpNotEq;
                          FetchToken();
                      } else
                      {   CommandError(ErrorMsg[ErrBadExpr]);
                          DeAllocateExpr( tmp1 );
                          return( (Expr*)0 );
                      }


                      if( CurToken == '-' )
                      {   FetchToken();
                          neg = True;
                      } else neg = False;

                      if( CurToken!=NumberTok )
                      {   CommandError(ErrorMsg[ErrNotNum]);
                          DeAllocateExpr( tmp1 );
                          return( (Expr*)0 );
                      } 

                      if( neg )
                      {     tmp1->rgt.val = -(int)TokenValue; 
                      } else tmp1->rgt.val = (int)TokenValue;
                      FetchToken();
                      return( tmp1 );
                      
                  } else switch( CurToken )
                  {   case('('):    FetchToken();
                                    if( !(tmp1=ParseExpression(0)) )
                                        return( (Expr*)0 );

                                    if( CurToken!=')' )
                                    {   CommandError(ErrorMsg[ErrParen]);
                                        DeAllocateExpr( tmp1 );
                                        return( (Expr*)0 );
                                    }
                                    FetchToken();
                                    return(tmp1);

                      case('!'): case('~'):
                      case(NotTok): FetchToken();
                                    if( !(tmp1=ParseExpression(2)) )
                                        return( (Expr*)0 );

                                    tmp2 = AllocateNode();
                                    tmp2->type = OpNot | OpRgtVal;
                                    tmp2->lft.ptr = tmp1;
                                    return( tmp2 );

                      case('-'):    if( NextIf(NumberTok,ErrNotNum) )
                                        return( (Expr*)0 );
                                    return( ParseRange(True) );

                      case(NumberTok):
                                    return( ParseRange(False) );

                      case(WithinTok):
                                    if( NextIf('(',ErrFunc) )
                                        return( (Expr*)0 );

                                    FetchToken();
                                    if( CurToken==NumberTok )
                                    {   if( *TokenPtr=='.' )
                                        {   TokenPtr++;
                                            FetchFloat(TokenValue,250);
                                        }
                                    } else if( CurToken!='.' )
                                    {   CommandError(ErrorMsg[ErrNotNum]);
                                        return( (Expr*)0 );
                                    } else FetchFloat(0,250);

                                    if( TokenValue>10000 )
                                    {   CommandError(ErrorMsg[ErrBigNum]);
                                        return( (Expr*)0 );
                                    }
                                    pred = (int)TokenValue;
                                    if( NextIf(',',ErrNotSep) )
                                        return( (Expr*)0 );

                                    FetchToken();
                                    if( !(tmp1=ParseExpression(0)) )
                                        return( (Expr*)0 );

                                    if( CurToken!=')' )
                                    {   CommandError(ErrorMsg[ErrParen]);
                                        DeAllocateExpr( tmp1 );
                                        return( (Expr*)0 );
                                    }

                                    FetchToken();
                                    if( !pred )
                                        return( tmp1 );

                                    tmp2 = AllocateNode();
                                    tmp2->type = OpWithin;
                                    tmp2->lft.limit = (Long)pred*pred;
                                    tmp2->rgt.set = BuildAtomSet(tmp1);
                                    DeAllocateExpr(tmp1);
                                    return( tmp2 );

                      default:      if( CurToken==IdentTok )
                                    {   tmp1 = LookUpSetExpr(TokenIdent);
                                        if( !tmp1 ) 
                                            tmp1 = LookUpElement(TokenIdent);

                                        if( tmp1 )
                                        {   FetchToken();
                                            return(tmp1);
                                        }
                                    }

                                    TokenPtr = TokenStart;
                                    done = ParsePrimitiveExpr(&TokenPtr);
                                    FetchToken();

                                    if( !done )
                                    {   CommandError(ErrorMsg[ErrBadExpr]);
                                        DeAllocateExpr( QueryExpr );
                                        return( (Expr*)0 );
                                    } else return( QueryExpr );
                  }
    }
    return (Expr*)0;
}



/*======================================*/
/*  RasMol Command Parsing & Execution  */
/*  Commands listed alphabetically      */
/*======================================*/

static void ExecuteAxesCommand( void )
{
    FetchToken();
    if( !CurToken || (CurToken==FalseTok) )
    {   ReDrawFlag |= RFRefresh;
        DrawAxes = False;
    } else if( CurToken == TrueTok )
    {   ReDrawFlag |= RFRefresh;
        DrawAxes = True;
    } else CommandError(ErrorMsg[ErrBadOpt]);
}


static void ExecuteBoundBoxCommand( void )
{
    FetchToken();
    if( !CurToken || (CurToken==FalseTok) )
    {   ReDrawFlag |= RFRefresh;
        DrawBoundBox = False;
    } else if( CurToken == TrueTok )
    {   ReDrawFlag |= RFRefresh;
        DrawBoundBox = True;
    } else CommandError(ErrorMsg[ErrBadOpt]);
}


static void ExecuteCentreCommand( void )
{
    register Real x, y, z;
    register Long count;

    FetchToken();
    if( !CurToken || (CurToken==AllTok) )
    {   CenX = CenY = CenZ = 0;
        return;
    }

    QueryExpr = ParseExpression(0);
    if( !QueryExpr ) return;

    if( CurToken )
    {   CommandError(ErrorMsg[ErrSyntax]);
        DeAllocateExpr(QueryExpr);
        return;
    }

    /* CentreZoneExpr(QueryExpr); */
    if( !Database ) return;

    count = 0;
    x = y = z = 0.0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
        for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
            for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
                if( EvaluateExpr(QueryExpr) )
                {   x += (Real)QAtom->xorg;
                    y += (Real)QAtom->yorg;
                    z += (Real)QAtom->zorg;
                    count++;
                }

    if( count )
    {   CenX = (Long)(x/count);
        CenY = (Long)(y/count);
        CenZ = (Long)(z/count);
        /* TranslateToCentre(CenX,CenY,CenZ); */
    } else
    {   InvalidateCmndLine();
        WriteString("No Atoms to Centre!\n");
    }
    DeAllocateExpr(QueryExpr);
}


static void ExecuteClipboardCommand( void )
{
    if( !ClipboardImage() )
    {   InvalidateCmndLine();
        WriteString("Unable to copy to clipboard!\n");
    }
}


static void ExecuteLoadCommand( void )
{
    register int format;
    register int info;
    register FILE *fp;

    FetchToken();
    format = FormatPDB;
    if( !*TokenPtr || *TokenPtr==' ' )
    {   if( IsMoleculeToken(CurToken) )
        {   format = Tok2Format(CurToken);
            FetchToken();
        } else if( CurToken == DotsTok )
        {   format = FormatDots;
            FetchToken();
        } else if( CurToken == GraspTok )
        {   format = FormatGrasp;
            FetchToken();
        }
    }

    if( !CurToken )
    {   CommandError(ErrorMsg[ErrFilNam]);
        return;
    }

#ifdef STRICT
    if( (CurToken!=StringTok) && (CurToken!=IdentTok) )
    {   CommandError(ErrorMsg[ErrFilNam]);
        return;
    }
#endif

    info = (FileDepth == -1);
    if( IsMoleculeFormat(format) )
    {   if( Database )
        {   CommandError(ErrorMsg[ErrBadLoad]);
            return;
        }

        if( CurToken==InLineTok )
        {   if( (FileDepth!=-1) && LineStack[FileDepth] )
            {   FetchFile(format,info,NameStack[FileDepth]);
            } else CommandError(ErrorMsg[ErrOutScrpt]);
        } else if( CurToken==StringTok )
        {      FetchFile(format,info,TokenIdent);
        } else FetchFile(format,info,TokenStart);
        DefaultRepresentation();
    } else if( format == FormatDots )
    {   if( !Database )
        {   CommandError(ErrorMsg[ErrBadMolDB]);
            return;
        }

        if( CurToken==StringTok )
        {      ProcessFileName(TokenIdent);
        } else ProcessFileName(TokenStart);

        if( !(fp=fopen(DataFileName,"rb")) )
        {   CommandError( (char*)0 );
            WriteString("Cannot open dots file '");
            WriteString(DataFileName);
            WriteString("'\n");
            return;
        } else
        {   LoadDotsFile(fp,info);
            fclose(fp);
        }
    } else if( format == FormatGrasp )
    {   if( !Database )
        {   CommandError(ErrorMsg[ErrBadMolDB]);
            return;
        }

        if( CurToken==StringTok )
        {      ProcessFileName(TokenIdent);
        } else ProcessFileName(TokenStart);

        if( !(fp=fopen(DataFileName,"rb")) )
        {   CommandError( (char*)0 );
            WriteString("Cannot open grasp file '");
            WriteString(DataFileName);
            WriteString("'\n");
            return;
        } else
        {   LoadGraspSurface(fp,info);
            fclose(fp);
        }
    }
    CurToken = 0;
}


static void ExecutePauseCommand( void )
{
    if( FileDepth == -1 )
    {   CommandError(ErrorMsg[ErrOutScrpt]);
        return;
    }

    /* Ignore Pause Commands via IPC! */
    if( LineStack[FileDepth] )
    {   CommandActive = True;
        IsPaused = True;

#ifdef MSWIN
        /* Disable Drag & Drop! */
        DragAcceptFiles(CanvWin,FALSE);
#endif
    }
}


static void ExecutePickingCommand( void )
{
    switch( FetchToken() )
    {   case(TrueTok):     case(0):
        case(IdentifyTok): SetPickMode(PickIdent); break;
        case(FalseTok):
        case(NoneTok):     SetPickMode(PickNone);  break;
        case(LabelTok):    SetPickMode(PickLabel); break;
        case(DistanceTok): SetPickMode(PickDist);  break;
        case(AngleTok):    SetPickMode(PickAngle); break;
        case(TorsionTok):  SetPickMode(PickTorsn); break;
        case(MonitorTok):  SetPickMode(PickMonit); break;
        case(CentreTok):   SetPickMode(PickCentr); break;
        case(OriginTok):   SetPickMode(PickOrign); break;
        default:           CommandError(ErrorMsg[ErrBadOpt]);
    }
}


static void ExecutePrintCommand( void )
{
    if( !PrintImage() )
    {   InvalidateCmndLine();
        WriteString("Warning: No suitable printer!\n");
    }
}


static void ExecuteTitleCommand( void )
{
    FetchToken();
    if( !CurToken )
    {   SetCanvasTitle("RasMol Version 2.6");
    } else if( CurToken == StringTok )
    {      SetCanvasTitle(TokenIdent);
    } else SetCanvasTitle(TokenStart);
    CurToken = 0;
}


static void ExecuteUnitCellCommand( void )
{
    FetchToken();
    if( !CurToken || (CurToken==FalseTok) )
    {   ReDrawFlag |= RFRefresh;
        DrawUnitCell = False;
    } else if( CurToken == TrueTok )
    {   ReDrawFlag |= RFRefresh;
        DrawUnitCell = True;
    } else CommandError(ErrorMsg[ErrBadOpt]);
}



/*=======================================*/
/*  Generic Command Parsing & Execution  */
/*=======================================*/

static void ExecuteSetCommand( void )
{
    register int option;

    switch( FetchToken() )
    {   case(SlabTok):
        case(SlabModeTok):
            option = -1;
            FetchToken();
            if( CurToken==RejectTok )
            {   option = SlabReject;
            } else if( CurToken==HalfTok )
            {   option = SlabHalf;
            } else if( CurToken==HollowTok )
            {   option = SlabHollow;
            } else if( CurToken==SolidTok )
            {   option = SlabClose;
            } else if( CurToken==SectionTok )
                option = SlabSection;

            if( option != -1 )
            {   if( UseSlabPlane && (SlabMode!=option) )
                    ReDrawFlag |= RFRefresh;
                SlabMode = option;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(ShadowTok):
            FetchToken();
            if( CurToken==TrueTok )
            {   UseShadow = True;
                ReviseInvMatrix();
                VoxelsClean = False;
                UseSlabPlane = False;
                ReDrawFlag |= RFRefresh;
                ReAllocBuffers();
            } else if( CurToken==FalseTok )
            {   ReDrawFlag |= RFRefresh;
                UseShadow = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;
                                  
        case(SpecularTok):
            FetchToken();
            if( CurToken==TrueTok )
            {   FakeSpecular = True;
                ReDrawFlag |= RFColour;
            } else if( CurToken==FalseTok )
            {   FakeSpecular = False;
                ReDrawFlag |= RFColour;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(SpecPowerTok):
            FetchToken();
            if( !CurToken )
            {   SpecPower = 8;
                ReDrawFlag |= RFColour;
            } else if( CurToken==NumberTok )
            {   if( TokenValue<=100 )
                {   ReDrawFlag |= RFColour;
                    SpecPower = (int)TokenValue;
                } else 
                    CommandError(ErrorMsg[ErrBigNum]);
            } else CommandError(ErrorMsg[ErrNotNum]);
            break;

        case(AmbientTok):
            FetchToken();
            if( !CurToken )
            {   ReDrawFlag |= RFColour;
                Ambient = DefaultAmbient;
            } else if( CurToken==NumberTok )
            {   if( TokenValue<=100 )
                {   Ambient = TokenValue/100.0;
                    ReDrawFlag |= RFColour;
                } else
                    CommandError(ErrorMsg[ErrBigNum]); 
            } else CommandError(ErrorMsg[ErrNotNum]);
            break;

        case(HeteroTok):
            FetchToken();
            if( CurToken==TrueTok )
            {   HetaGroups = True;
            } else if( CurToken==FalseTok )
            {   HetaGroups = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;
                                  
        case(HydrogenTok):
            FetchToken();
            if( CurToken==TrueTok )
            {   Hydrogens = True;
            } else if( CurToken==FalseTok )
            {   Hydrogens = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;
                                  

        case(BackgroundTok):
            FetchToken();
            if( CurToken == TransparentTok )
            {   UseTransparent = True;
            } else if( CurToken == NormalTok )
            {   UseTransparent = False;
            } else if( ParseColour() )
            {   ReDrawFlag |= RFColour;
                BackR = RVal;
                BackG = GVal;
                BackB = BVal;
#ifndef IBMPC
                FBClear = False;
#endif
            } else if( CurToken )
            {   CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(BondModeTok):
            FetchToken();
            if( !CurToken || (CurToken==AndTok) )
            {   ZoneBoth = True;
            } else if( CurToken==OrTok )
            {   ZoneBoth = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;
            
        case(HBondTok):
            FetchToken();
            if( (CurToken==BackboneTok) || (CurToken==MainChainTok) )
            {   ReDrawFlag |= RFRefresh;
                HBondMode = True;
            } else if( !CurToken || (CurToken==SidechainTok) )
            {   ReDrawFlag |= RFRefresh;
                HBondMode = False;
            } else if( CurToken == ChainTok )
            {   FetchToken();
                if( !CurToken || (CurToken==TrueTok) )
                {   if( !HBondChainsFlag && (Info.hbondcount>=0) )
                    {   ReDrawFlag |= RFRefresh;
                        HBondChainsFlag = True;
                        CalcHydrogenBonds();
                    }
                } else if( CurToken == FalseTok )
                {   if( HBondChainsFlag && (Info.hbondcount>=0) )
                    {   ReDrawFlag |= RFRefresh;
                        HBondChainsFlag = False;
                        CalcHydrogenBonds();
                    }
                } else CommandError(ErrorMsg[ErrBadOpt]);
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(SSBondTok):
            FetchToken();
            if( (CurToken==BackboneTok) || (CurToken==MainChainTok) )
            {   ReDrawFlag |= RFRefresh;
                SSBondMode = True;
            } else if( !CurToken || (CurToken==SidechainTok) )
            {   ReDrawFlag |= RFRefresh;
                SSBondMode = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(HourGlassTok):
            FetchToken();
            if( CurToken==TrueTok )
            {   UseHourGlass = True;
            } else if( CurToken==FalseTok )
            {   UseHourGlass = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(StrandsTok):
            FetchToken();
            if( !CurToken )
            {   ReDrawFlag |= RFRefresh;
                SplineCount = 5;
            } else if( CurToken==NumberTok )
            {   if( (TokenValue>0) && (TokenValue<=5) )
                {   SplineCount = (int)TokenValue;
                    ReDrawFlag |= RFRefresh;
                } else if( TokenValue==9 )
                {   ReDrawFlag |= RFRefresh;
                    SplineCount = 9;
                } else CommandError(ErrorMsg[ErrBadOpt]);
            } else CommandError(ErrorMsg[ErrNotNum]);
            break;

        case(MouseTok):
            FetchToken();
            if( !CurToken || (CurToken==RasMolTok) )
            {   if( Interactive )
                    SetMouseMode( MMRasMol );
            } else if( CurToken==InsightTok )
            {   if( Interactive )
                    SetMouseMode( MMInsight );
            } else if( CurToken==QuantaTok )
            {   if( Interactive )
                    SetMouseMode( MMQuanta );
            } else if( CurToken==SybylTok )
            {   if( Interactive )
                    SetMouseMode( MMSybyl );
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(DisplayTok):
            FetchToken();
            /* Affect StereoMode Parameters?? */
            if( !CurToken || (CurToken==NormalTok) )
            {   ReDrawFlag |= RFRefresh | RFColour;
                DisplayMode = 0;
            } else if( CurToken==SelectedTok )
            {   ReDrawFlag |= RFRefresh | RFColour;
                DisplayMode = 1;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(VectPSTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   UseOutLine = False;
            } else if( CurToken == TrueTok )
            {   UseOutLine = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(KinemageTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   KinemageFlag = False;
            } else if( CurToken == TrueTok )
            {   KinemageFlag = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(MenusTok):
            FetchToken();
            if( !CurToken || (CurToken==TrueTok) )
            {   EnableMenus(True);
            } else if( CurToken == FalseTok )
            {   EnableMenus(False);
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(RadiusTok):
            FetchToken();
            if( !CurToken )
            {   ProbeRadius = SolventDots? 300 : 0;
            } else if( CurToken==NumberTok )
            {   if( *TokenPtr=='.' )
                {   TokenPtr++;
                    FetchFloat(TokenValue,250);
                }

                if( TokenValue>750 )
                {   CommandError(ErrorMsg[ErrBigNum]);
                } else ProbeRadius = (int)TokenValue;
            } else if( CurToken=='.' )
            {   FetchFloat(0,250);
                if( TokenValue>750 )
                {   CommandError(ErrorMsg[ErrBigNum]);
                } else ProbeRadius = (int)TokenValue;

            } else CommandError(ErrorMsg[ErrNotNum]);
            break;

        case(SolventTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   SolventDots = False;
                ProbeRadius = 0;
            } else if( CurToken == TrueTok )
            {   SolventDots = True;
                ProbeRadius = 300;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(FontSizeTok):
            FetchToken();
            if( CurToken==NumberTok )
            {   if( TokenValue<=32 )
                {   if( LabelList || (MonitList && DrawMonitDistance) )
                        ReDrawFlag |= RFRefresh;
                    SetFontSize((int)TokenValue);
                } else CommandError(ErrorMsg[ErrBigNum]);
            } else if( !CurToken )
            {   if( LabelList )
                    ReDrawFlag |= RFRefresh;
                SetFontSize(8);
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(WriteTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   AllowWrite = False;
            } else if( CurToken == TrueTok )
            {   if( (FileDepth!=-1) && LineStack[FileDepth] )
                {   CommandError(ErrorMsg[ErrInScrpt]);
                } else AllowWrite = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(StereoTok):  
            FetchToken();
            if( !CurToken )
            {   SetStereoMode(False);
                StereoAngle = 6.0;
            } else if( CurToken==TrueTok )
            {   SetStereoMode(True);
            } else if( CurToken==FalseTok )
            {   SetStereoMode(False);
            } else if( CurToken == '-' )
            {   if( !NextIf(NumberTok,ErrNotNum) )
                {   StereoAngle = -TokenValue;
                    SetStereoMode(True);
                }
            } else if( CurToken == '+' )
            {   if( !NextIf(NumberTok,ErrNotNum) )
                {   StereoAngle = TokenValue;
                    SetStereoMode(True);
                }
            } else if( CurToken==NumberTok )
            {   StereoAngle = TokenValue;
                SetStereoMode(True);
            } else CommandError(ErrorMsg[ErrSyntax]);
            break;

        case(BondTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   ReDrawFlag |= RFRefresh;
                DrawDoubleBonds = False;
            } else if( CurToken == TrueTok )
            {   ReDrawFlag |= RFRefresh;
                DrawDoubleBonds = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(MonitorTok):
            FetchToken();
            if( !CurToken )
            {   ReDrawFlag |= RFRefresh;
                DrawMonitDistance = True;
                MonitRadius = 0;
            } else if( CurToken == TrueTok )
            {   ReDrawFlag |= RFRefresh;
                DrawMonitDistance = True;
            } else if( CurToken == FalseTok )
            {   ReDrawFlag |= RFRefresh;
                DrawMonitDistance = False;
            } else if( CurToken==NumberTok )
            {   if( *TokenPtr=='.' )
                {   TokenPtr++;
                    FetchFloat(TokenValue,250);
                }

                if( TokenValue <= 750 )
                {   MonitRadius = (int)TokenValue;
                    ReDrawFlag |= RFRefresh;
                } else CommandError(ErrorMsg[ErrBigNum]);
            } else if( CurToken=='.' )
            {   FetchFloat(0,250);
                if( TokenValue <= 750 )
                {   MonitRadius = (int)TokenValue;
                    ReDrawFlag |= RFRefresh;
                } else CommandError(ErrorMsg[ErrBigNum]);
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(CartoonTok):
            FetchToken();
            if( !CurToken )
            {   ReDrawFlag |= RFRefresh;
                DrawBetaArrows = True;
                CartoonHeight = 120;
            } else if( CurToken==TrueTok )
            {   ReDrawFlag |= RFRefresh;
                DrawBetaArrows = True;
            } else if( CurToken==FalseTok )
            {   ReDrawFlag |= RFRefresh;
                DrawBetaArrows = False;
            } else if( CurToken==NumberTok )
            {   if( *TokenPtr=='.' )
                {   TokenPtr++;
                    FetchFloat(TokenValue,250);
                }

                if( TokenValue <= 500 )
                {   CartoonHeight = (int)TokenValue;
                    ReDrawFlag |= RFRefresh;
                } else CommandError(ErrorMsg[ErrBigNum]);
            } else if( CurToken=='.' )
            {   FetchFloat(0,250);
                if( TokenValue <= 500 )
                {   CartoonHeight = (int)TokenValue;
                    ReDrawFlag |= RFRefresh;
                } else CommandError(ErrorMsg[ErrBigNum]);

            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(BackFadeTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   ReDrawFlag |= RFColour;
                UseBackFade = False;
            } else if( CurToken == TrueTok )
            {   ReDrawFlag |= RFColour;
                UseBackFade = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(TransparentTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   UseTransparent = False;
            } else if( CurToken == TrueTok )
            {   UseTransparent = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(DepthCueTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   ReDrawFlag |= RFColour;
                UseDepthCue = False;
            } else if( CurToken == TrueTok )
            {   ReDrawFlag |= RFColour;
                UseDepthCue = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(SequenceTok):
            FetchToken();
            if( !CurToken || (CurToken==FalseTok) )
            {   SeqFormat = False;
            } else if( CurToken == TrueTok )
            {   SeqFormat = True;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(ConnectTok):
            FetchToken();
            if( !CurToken || (CurToken==TrueTok) )
            {   CalcBondsFlag = True;
            } else if( CurToken == FalseTok )
            {   CalcBondsFlag = False;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(AxesTok):     ExecuteAxesCommand();     break;
        case(BoundBoxTok): ExecuteBoundBoxCommand(); break;
        case(PickingTok):  ExecutePickingCommand();  break;
        case(TitleTok):    ExecuteTitleCommand();    break;
        case(UnitCellTok): ExecuteUnitCellCommand(); break;

        case(DotsTok):
            FetchToken();
            if( !CurToken )
            {   SurfaceChainsFlag = False;
                SolventDots = False;
                ProbeRadius = 0;
            } else if( CurToken == SolventTok )
            {   FetchToken();
                if( !CurToken || (CurToken==FalseTok) )
                {   SolventDots = False;
                    ProbeRadius = 0;
                } else if( CurToken == TrueTok )
                {   SolventDots = True;
                    ProbeRadius = 300;
                } else CommandError(ErrorMsg[ErrBadOpt]);
            } else if( CurToken == ChainTok )
            {   FetchToken();
                if( !CurToken || (CurToken==TrueTok) )
                {   SurfaceChainsFlag = True;
                } else if( CurToken == FalseTok )
                {   SurfaceChainsFlag = False;
                } else CommandError(ErrorMsg[ErrBadOpt]);
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(STLTok):
        case(DXFTok):
        case(PLYTok):
            FetchToken();
            if( !CurToken )
            {   ResolutionFlag = ResDefault;
                BinaryFlag = True;
            } else if( CurToken == TrueTok )
            {   BinaryFlag = True;
            } else if( CurToken == FalseTok )
            {   BinaryFlag = False;
            } else if( CurToken == HighTok )
            {   ResolutionFlag = ResHigh;
            } else if( CurToken == MediumTok )
            {   ResolutionFlag = ResMedium;
            } else if( CurToken == LowTok )
            {   ResolutionFlag = ResLow;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        case(SurfaceTok):
            FetchToken();
            if( CurToken == NoneTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfNone;
            } else if( CurToken == DotsTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfDots;
            } else if( CurToken == LineTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfLines;
            } else if( CurToken == VectorTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfVectors;
            } else if( CurToken == NormalTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfNormal;
            } else if( CurToken == SolidTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfSolid;
            } else if( CurToken == TransparentTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfTransparent;
            } else if( CurToken == TranslucentTok ) {
                ReDrawFlag |= RFRefresh;
                SurfaceMode = SurfTranslucent;
            } else CommandError(ErrorMsg[ErrBadOpt]);
            break;

        default:
            CommandError(ErrorMsg[ErrParam]);
    }
}


static void ExecuteColourCommand( void )
{
    register int flag;

    flag = 0;
    switch( FetchToken() )
    {   case(AtomTok):
            FetchToken();
        default:
            switch( CurToken )
            {   case(CPKTok):         CPKColourAttrib(); 
                                      ReDrawFlag |= RFColour; break;

                case(AminoTok):       AminoColourAttrib();
                                      ReDrawFlag |= RFColour; break;

                case(ShapelyTok):     ShapelyColourAttrib();
                                      ReDrawFlag |= RFColour; break;
                
                case(UserTok):        UserMaskAttrib(MaskColourFlag);
                                      ReDrawFlag |= RFColour; break;

                case(GroupTok):       ScaleColourAttrib(GroupAttr);
                                      ReDrawFlag |= RFColour; break;

                case(ChainTok):       ScaleColourAttrib(ChainAttr);
                                      ReDrawFlag |= RFColour; break;

                case(ChargeTok):      ScaleColourAttrib(ChargeAttr);
                                      ReDrawFlag |= RFColour; break;

                case(TemperatureTok): ScaleColourAttrib(TempAttr);
                                      ReDrawFlag |= RFColour; break;

                case(StructureTok):   StructColourAttrib();
                                      ReDrawFlag |= RFColour; break;

                default:  if( ParseColour() )
                          {   MonoColourAttrib(RVal,GVal,BVal);
                              ReDrawFlag |= RFColour;
                          } else if( CurToken )
                          {      CommandError(ErrorMsg[ErrColour]);
                          } else CommandError(ErrorMsg[ErrNoCol]);
            }
            break;

        case(BondTok):    
        case(DashTok):
            FetchToken();
            if( CurToken==NoneTok )
            {   ColourBondNone();
                ReDrawFlag |= RFColour;
            } else if( ParseColour() )
            {   ColourBondAttrib(RVal,GVal,BVal);
                ReDrawFlag |= RFColour;
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(BackboneTok):
            FetchToken();
            if( CurToken==NoneTok )
            {   ColourBackNone();
                ReDrawFlag |= RFColour;
            } else if( ParseColour() )
            {   ColourBackAttrib(RVal,GVal,BVal);
                ReDrawFlag |= RFColour;
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(SSBondTok):
            FetchToken();
            if( CurToken==NoneTok )
            {   ReDrawFlag |= RFColour;
                ColourHBondNone( False );
            } else if( ParseColour() )
            {   ReDrawFlag |= RFColour;
                ColourHBondAttrib(False,RVal,GVal,BVal);
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(HBondTok):
            FetchToken();
            if( CurToken==NoneTok )
            {   ReDrawFlag |= RFColour;
                ColourHBondNone( True );
            } else if( CurToken==TypeTok )
            {   ReDrawFlag |= RFColour;
                ColourHBondType();
            } else if( ParseColour() )
            {   ReDrawFlag |= RFColour;
                ColourHBondAttrib(True,RVal,GVal,BVal);
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(DotsTok):
            FetchToken();
            if( CurToken==PotentialTok )
            {   ReDrawFlag |= RFColour;
                ColourDotsPotential();
            } else if( ParseColour() )
            {   ReDrawFlag |= RFColour;
                ColourDotsAttrib(RVal,GVal,BVal);
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(SurfaceTok):
            FetchToken();
            if( CurToken == PotentialTok )
            {   ReDrawFlag |= RFColour;
                ColourSurfacePotential();
            } else if( ParseColour() )
            {   ReDrawFlag |= RFColour;
                ColourSurfaceAttrib(RVal,GVal,BVal);
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(MonitorTok):
            FetchToken();
            if( CurToken == NoneTok )
            {   ColourMonitNone();
            } else if( ParseColour() )
            {   ReDrawFlag |= RFColour;
                ColourMonitAttrib(RVal,GVal,BVal);
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(AxesTok):
        case(BoundBoxTok):
        case(UnitCellTok):
            FetchToken();
            if( ParseColour() )
            {   BoxR = RVal;  BoxG = GVal;  BoxB = BVal;
                ReDrawFlag |= RFColour;
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(LabelTok):
            FetchToken();
            if( CurToken==NoneTok )
            {   ReDrawFlag |= RFColour;
                UseLabelCol = False;
            } else if( ParseColour() )
            {   LabR = RVal;  LabG = GVal;  LabB = BVal;
                ReDrawFlag |= RFColour;
                UseLabelCol = True;
            } else if( CurToken )
            {      CommandError(ErrorMsg[ErrColour]);
            } else CommandError(ErrorMsg[ErrNoCol]);
            break;

        case(TraceTok): 
        case(RibbonTok):
        case(CartoonTok):  flag = RibColBoth;     break;
        case(Ribbon1Tok):  flag = RibColInside;   break;
        case(Ribbon2Tok):  flag = RibColOutside;  break;
    }

    if( flag )
    {   FetchToken();
        if( CurToken==NoneTok )
        {   ReDrawFlag |= RFColour;
            ColourRibbonNone(flag);
        } else if( ParseColour() )
        {   ReDrawFlag |= RFColour;
            ColourRibbonAttrib(flag,RVal,GVal,BVal);
        } else if( CurToken )
        {      CommandError(ErrorMsg[ErrColour]);
        } else CommandError(ErrorMsg[ErrNoCol]);
    }
}


static void DescribeSequence( void )
{
    register Chain __far *chn;
    register Group __far *grp;
    register int chain,count;
    register char *str;
    char buffer[40];

    InvalidateCmndLine();
    if( !Database )
        return;

    for( chn=Database->clist; chn; chn=chn->cnext )
    {   chain = (Info.chaincount<2);  count = 0;
        for( grp=chn->glist; grp; grp=grp->gnext )
            if( grp->alist && !(grp->alist->flag&HeteroFlag) )
            {   if( !chain )
                {   WriteString("Chain ");
                    WriteChar(chn->ident);
                    WriteString(":\n");
                    chain = True;
                }

                if( !SeqFormat )
                {   if( count == 10 )
                    {   WriteChar('\n');
                        count = 1;
                    } else count++;

                    str = Residue[grp->refno];
                    WriteChar(str[0]);
                    WriteChar(str[1]);
                    WriteChar(str[2]);

                    sprintf(buffer,"%-3d ",grp->serno);
                    WriteString(buffer);
                } else
                {   if( count == 60 )
                    {   WriteChar('\n');
                        count = 1;
                    } else count++;

                    if( grp->refno < 29 )
                    {   WriteChar(ResidueChar[grp->refno]);
                    } else WriteChar('*');
                }
            }
        WriteChar('\n');
    }
    WriteChar('\n');
}


static void ExecuteShowCommand( void )
{
    register Real temp;
    char buffer[40];

    switch( FetchToken() )
    {   case(InfoTok):
                DescribeMolecule();
                break;

        case(SequenceTok):
                DescribeSequence();
                break;

        case(SymmetryTok):
                InvalidateCmndLine();

                if( *Info.spacegroup )
                {   sprintf(buffer,"Space Group ...... %s\n",Info.spacegroup);
                    WriteString(buffer);

                    sprintf(buffer,"Unit cell A ...... %g\n",Info.cella);
                    WriteString(buffer);
                    sprintf(buffer,"Unit cell B ...... %g\n",Info.cellb);
                    WriteString(buffer);
                    sprintf(buffer,"Unit cell C ...... %g\n",Info.cellc);
                    WriteString(buffer);

                    temp = Rad2Deg*Info.cellalpha;
                    sprintf(buffer,"Unit cell alpha .. %g\n",temp);
                    WriteString(buffer);
                    temp = Rad2Deg*Info.cellbeta;
                    sprintf(buffer,"Unit cell beta ... %g\n",temp);
                    WriteString(buffer);
                    temp = Rad2Deg*Info.cellgamma;
                    sprintf(buffer,"Unit cell gamma .. %g\n",temp);
                    WriteString(buffer);

                } else WriteString("No crystal symmetry data!\n");
                WriteChar('\n');
                break;

        default:
            CommandError(ErrorMsg[ErrBadArg]);
    }
}


void ZapDatabase( void )
{
    register int i;

    for( i=0; i<8; i++ )
        DialValue[i] = 0.0;
    SelectCount = 0;

    DestroyDatabase();
    ResetSymbolTable();
    ResetTransform();
    ResetRenderer();
    ResetRepres();

    ZoneBoth = True;
    HetaGroups = True;    
    Hydrogens = True;

    BackR = BackG = BackB = 0;
#ifndef IBMPC
    FBClear = False;
#endif

    ResetColourMap();
    DefineColourMap();
    ClearBuffers();
    ReDrawFlag = 0;

    if( Interactive )
    {   UpdateScrollBars();
        ClearImage();
    }
    AdviseUpdate(AdvName);
    AdviseUpdate(AdvClass);
    AdviseUpdate(AdvIdent);
}


static void WriteImageFile( char *name, int type )
{
    if( !type )
#ifdef EIGHTBIT
        type = GIFTok;
#else
        type = PPMTok;
#endif

    switch( type )
    {   case(GIFTok):     WriteGIFFile(name);             break;
        case(BMPTok):     WriteBMPFile(name);             break;
        case(PPMTok):     WritePPMFile(name,True);        break;
        case(SUNTok):     WriteRastFile(name,False);      break;
        case(SUNRLETok):  WriteRastFile(name,True);       break;
        case(PICTTok):    WritePICTFile(name);            break;
        case(IRISTok):    WriteIRISFile(name);            break;
        case(EPSFTok):    WriteEPSFFile(name,True,True);  break;
        case(MonoPSTok):  WriteEPSFFile(name,False,True); break;
        case(VectPSTok):  WriteVectPSFile(name);          break;

        case(RasMolTok):
        case(ScriptTok):     WriteScriptFile(name);     break;
        case(KinemageTok):   WriteKinemageFile(name);   break;
        case(MolScriptTok):  WriteMolScriptFile(name);  break;
        case(POVRayTok):     WritePOVRayFile(name);     break;
        case(GraspTok):      WriteGraspFile(name);      break;
        case(VRMLTok):       WriteVRMLFile(name);       break;
        case(STLTok):        WriteSTLFile(name);        break;
        case(DXFTok):        WriteDXFFile(name);        break;
        case(PLYTok):        WritePLYFile(name);        break;
    }
}


void ResumePauseCommand( void )
{
    register int ch,len;
    register FILE *fp;
    register int stat;

    CommandActive = False;
    IsPaused = False;

#ifdef MSWIN
    /* Re-enable Drag & Drop! */
    DragAcceptFiles(CanvWin,TRUE);
#endif

    while( FileDepth >= 0 )
    {   fp = FileStack[FileDepth];
        do {
            len = 0;
            ch = getc(fp);
            while( (ch!='\n') && (ch!='\r') && (ch!=EOF) )
            {   if( len<MAXBUFFLEN )
                    CurLine[len++] = ch;
                ch = getc(fp);
            }

            if( ch == '\r' )
            {   ch = getc(fp);
                if( ch != '\n' )
                    ungetc(ch,fp);
            }

            LineStack[FileDepth]++;
            if( len<MAXBUFFLEN )
            {   CurLine[len] = '\0';
                stat = ExecuteCommand();
                if( stat )
                {   if( stat == QuitTok )
                    {   while( FileDepth >= 0 )
                        {   fclose(FileStack[FileDepth]);
                            free(NameStack[FileDepth]);
                            FileDepth--;
                        }
                        RasMolExit();
                    } else /* ExitTok */
                        break;
                } else if( IsPaused )
                    return;
            } else CommandError("Script command line too long");
        } while( ch != EOF );
        free(NameStack[FileDepth]);
        fclose( fp );
        FileDepth--;
    }
}


void InterruptPauseCommand( void )
{
    WriteString("*** RasMol script interrupted! ***\n\n");
    CommandActive = False;
    IsPaused = False;

#ifdef MSWIN
    /* Re-enable Drag & Drop! */
    DragAcceptFiles(CanvWin,TRUE);
#endif

    while( FileDepth >= 0 )
    {   fclose(FileStack[FileDepth]);
        free(NameStack[FileDepth]);
        FileDepth--;
    }
}


static void ExecuteConnectCommand( void )
{
    register Bond __far *bptr;
    register int info,flag;

    FetchToken();
    if( !CurToken )
    {   flag = (MainAtomCount+HetaAtomCount<256);
    } else if( CurToken == TrueTok )
    {   flag = True;
    } else if( CurToken == FalseTok )
    {   flag = False;
    } else 
    {   CommandError(ErrorMsg[ErrSyntax]);
        return;
    }

    if( Database )
    {   ForEachBond
            if( bptr->col )
                Shade[Colour2Shade(bptr->col)].refcount--;
        info = (FileDepth == -1);
        CreateMoleculeBonds(info,flag);
        ReDrawFlag |= RFRefresh|RFColour;
        EnableWireframe(WireFlag,0);
    }
}


int ExecuteCommand( void )
{
    register char *param;
    register int option;
    register int i,done;
    register Long temp;
    FILE *script;

    TokenPtr = CurLine;
    if( !FetchToken() )
    {   TokenPtr = (char*)0;
        return( False );
    }

    switch( CurToken )
    {   case(AxesTok):       ExecuteAxesCommand();       break;
        case(BoundBoxTok):   ExecuteBoundBoxCommand();   break;
        case(CentreTok):     ExecuteCentreCommand();     break;
        case(ClipboardTok):  ExecuteClipboardCommand();  break;
        case(ColourTok):     ExecuteColourCommand();     break;
        case(ConnectTok):    ExecuteConnectCommand();    break;
        case(LoadTok):       ExecuteLoadCommand();       break;
        case(WaitTok):       ExecutePauseCommand();      break;
        case(PickingTok):    ExecutePickingCommand();    break;
        case(PrintTok):      ExecutePrintCommand();      break;
        case(SetTok):        ExecuteSetCommand();        break;
        case(ShowTok):       ExecuteShowCommand();       break;
        case(TitleTok):      ExecuteTitleCommand();      break;
        case(UnitCellTok):   ExecuteUnitCellCommand();   break;

        case(RefreshTok):    RefreshScreen();            break;
        case(ZapTok):        ZapDatabase();              break;


        case(SelectTok):  FetchToken();
                          if( !CurToken )
                          {   option = NormAtomFlag;
                              if( HetaGroups ) option |= HeteroFlag;
                              if( Hydrogens )  option |= HydrogenFlag;
                              SelectZone(option);
                          } else if( CurToken==AllTok )
                          {   SelectZone(AllAtomFlag);
                          } else if( CurToken==NoneTok )
                          {   SelectZone(0x00);
                          } else
                          {   QueryExpr = ParseExpression(0);
                              if( QueryExpr )
                              {   if( !CurToken )
                                  {   SelectZoneExpr(QueryExpr);
                                  } else CommandError(ErrorMsg[ErrSyntax]);
                                  DeAllocateExpr(QueryExpr);
                              }
                          }
                          break;

        case(RestrictTok):
                          FetchToken();
                          if( !CurToken )
                          {   option = NormAtomFlag;
                              if( HetaGroups ) option |= HeteroFlag;
                              if( Hydrogens )  option |= HydrogenFlag;
                              RestrictZone(option);
                              ReDrawFlag |= RFRefresh;
                          } else if( CurToken==AllTok )
                          {   RestrictZone(AllAtomFlag);
                              ReDrawFlag |= RFRefresh;
                          } else if( CurToken==NoneTok )
                          {   RestrictZone(0x00);
                              ReDrawFlag |= RFRefresh;
                          } else
                          {   QueryExpr = ParseExpression(0);
                              if( QueryExpr )
                              {   if( !CurToken )
                                  {   RestrictZoneExpr(QueryExpr);
                                      ReDrawFlag |= RFRefresh;
                                  } else CommandError(ErrorMsg[ErrSyntax]);
                                  DeAllocateExpr(QueryExpr);
                              }
                          }
                          break;



        case(WireframeTok):
                          FetchToken();
                          if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DisableWireframe();
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              EnableWireframe(WireFlag,0);
                          } else if( CurToken==DashTok )
                          {   ReDrawFlag |= RFRefresh;
                              EnableWireframe(DashFlag,0);
                          } else if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=500 )
                              {   EnableWireframe(CylinderFlag,
                                                  (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=500 )
                              {   EnableWireframe(CylinderFlag,
                                                  (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(BackboneTok):
                          FetchToken();
                          if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DisableBackbone();
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              EnableBackbone(WireFlag,0);
                          } else if( CurToken==DashTok )
                          {   ReDrawFlag |= RFRefresh;
                              EnableBackbone(DashFlag,0);
                          } else if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }


                              if( TokenValue<=500 )
                              {   EnableBackbone(CylinderFlag,
                                                 (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=500 )
                              {   EnableBackbone(CylinderFlag,
                                                 (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(CPKTok):
        case(SpacefillTok):
                          FetchToken();
                          if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DisableSpacefill();
                          } else if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=750 )
                              {   SetRadiusValue(MaxFun((int)TokenValue,1));
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=750 )
                              {   SetRadiusValue(MaxFun((int)TokenValue,1));
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken==UserTok )
                          {   UserMaskAttrib(MaskRadiusFlag);
                              ReDrawFlag |= RFRefresh;
                          } else if( CurToken==TemperatureTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetRadiusTemperature();
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetVanWaalRadius();
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(DashTok):    FetchToken();
                          if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DisableWireframe();
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              EnableWireframe(DashFlag,0);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(SSBondTok):  FetchToken();
                          if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=500 )
                              {   SetHBondStatus(False,True,(int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=500 )
                              {   SetHBondStatus(False,True,(int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetHBondStatus(False,False,0);
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetHBondStatus(False,True,0);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(HBondTok):   FetchToken();
                          if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=500 )
                              {   SetHBondStatus(True,True,(int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=500 )
                              {   SetHBondStatus(True,True,(int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetHBondStatus(True,False,0);
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetHBondStatus(True,True,0);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(RibbonTok):  FetchToken();
                          if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=1000 )
                              {   SetRibbonStatus(True,RibbonFlag,
                                                  (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=1000 )
                              {   SetRibbonStatus(True,RibbonFlag,
                                                  (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(False,RibbonFlag,0);
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,RibbonFlag,0);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(StrandsTok): FetchToken();
                          if( CurToken == DashTok )
                          {   option = DashStrandFlag;
                              FetchToken();
                          } else option = StrandFlag;

                          if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=1000 )
                              {   SetRibbonStatus(True,option,(int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=1000 )
                              {   SetRibbonStatus(True,option,(int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(False,option,0);
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,option,0);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(TraceTok):   FetchToken();
                          if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(False,TraceFlag,80);
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,TraceFlag,80);
                          } else if( CurToken==TemperatureTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetTraceTemperature();
                          } else if( CurToken==DotsTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,DotsFlag,80);
                          } else if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=500 )
                              {   SetRibbonStatus(True,TraceFlag,
                                                 (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=500 )
                              {   SetRibbonStatus(True,TraceFlag,
                                                 (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(CartoonTok): FetchToken();
                          if( CurToken==NumberTok )
                          {   if( *TokenPtr=='.' )
                              {   TokenPtr++;
                                  FetchFloat(TokenValue,250);
                              }

                              if( TokenValue<=1000 )
                              {   SetRibbonStatus(True,CartoonFlag,
                                                  (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken=='.' )
                          {   FetchFloat(0,250);
                              if( TokenValue<=1000 )
                              {   SetRibbonStatus(True,CartoonFlag,
                                                  (int)TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken==FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(False,CartoonFlag,0);
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,CartoonFlag,0);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(DotsTok):    FetchToken();
                          if( CurToken == NumberTok )
                          {   if( TokenValue<=1000 )
                              {   if( TokenValue )
                                  {   CalculateDots((int)TokenValue);
                                  } else CalculateDots(1);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else if( CurToken == FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DeleteDots();
                          } else if( (CurToken==TrueTok) || !CurToken )
                          {   ReDrawFlag |= RFRefresh;
                              CalculateDots(100);
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(SurfaceTok): FetchToken();
                          if( CurToken == FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DeleteSurface();
                          } else if( CurToken == RibbonTok )
                          {   ReDrawFlag |= RFRefresh;
                              CreateRibbons();
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(MonitorTok): FetchToken();
                          if( CurToken == NumberTok )
                          {   temp = TokenValue;
                              FetchToken();
                              if( CurToken == ',' )
                                  FetchToken();

                              if( CurToken == NumberTok )
                              {   CreateMonitor(temp,TokenValue);
                                  ReDrawFlag |= RFRefresh;
                              } else CommandError(ErrorMsg[ErrNotNum]);
                          } else if( CurToken == FalseTok )
                          {   ReDrawFlag |= RFRefresh;
                              DeleteMonitors();
                          } else CommandError(ErrorMsg[ErrBadArg]);
                          break;

        case(SlabTok):    FetchToken();
                          if( (CurToken==NumberTok) || (CurToken=='.') )
                          {   if( CurToken==NumberTok )
                              {   if( *TokenPtr=='.' )
                                  {   TokenPtr++;
                                      FetchFloat(TokenValue,100);
                                  } else TokenValue *= 100;
                              } else FetchFloat(0,100);

                              if( TokenValue<=10000 )
                              {   DialValue[7] = (TokenValue-5000)/5000.0;
                                  /* UpdateScrollBars(); */
                                  ReDrawFlag |= RFSlab;
                                  UseSlabPlane = True;
                                  UseShadow = False;
                              } else CommandError(ErrorMsg[ErrBigNum]);

                          } else if( CurToken==FalseTok )
                          {   if( UseSlabPlane )
                              {   ReDrawFlag |= RFRefresh;
                                  UseSlabPlane = False;
                              }
                          } else if( !CurToken || (CurToken==TrueTok) )
                          {   if( !UseSlabPlane )
                              {   ReDrawFlag |= RFRefresh;
                                  UseSlabPlane = True;
                                  UseShadow = False;
                              }
                          } else CommandError(ErrorMsg[ErrSyntax]);
                          break;

        case(ZoomTok):    FetchToken();
                          if( (CurToken==NumberTok) || (CurToken=='.') )
                          {   if( CurToken==NumberTok )
                              {   if( *TokenPtr=='.' )
                                  {   TokenPtr++;
                                      FetchFloat(TokenValue,100);
                                  } else TokenValue *= 100;
                              } else FetchFloat(0,100);

                              if( TokenValue<=10000 )
                              {   DialValue[3] = (TokenValue-10000)/10000.0;
                                  ReDrawFlag |= RFZoom;
                              } else if( Database )
                              {   /* Magnification */
                                  TokenValue -= 10000;
                                  temp = (Long)(MaxZoom*10000);
                                  if( TokenValue<=temp )
                                  {   DialValue[3] = (Real)TokenValue/temp;
                                      ReDrawFlag |= RFZoom;
                                  } else CommandError(ErrorMsg[ErrBigNum]);
                              }
                          } else if( CurToken==TrueTok )
                          {   ReDrawFlag |= RFZoom;
                              DialValue[3] = 0.5;
                          } else if( !CurToken || (CurToken==FalseTok) )
                          {   ReDrawFlag |= RFZoom;
                              DialValue[3] = 0.0;
                          } else CommandError(ErrorMsg[ErrSyntax]);
                          /* UpdateScrollBars(); */
                          break;

        case(RotateTok):  FetchToken();
                          if( CurToken==XTok )
                          {   option = 0;
                          } else if( CurToken==YTok )
                          {   option = 1;
                          } else if( CurToken==ZTok )
                          {   option = 2;
                          } else
                          {   CommandError(ErrorMsg[ErrSyntax]);
                              break;
                          }

                          FetchToken();
                          if( CurToken == '-' )
                          {   FetchToken();
                              done = True;
                          } else done = False;
#ifdef INVERT
                          if( option != 1 )
                              done = !done;
#endif
                          if( (CurToken==NumberTok) || (CurToken=='.') )
                          {   if( CurToken==NumberTok )
                              {   if( *TokenPtr=='.' )
                                  {   TokenPtr++;
                                      FetchFloat(TokenValue,100);
                                  } else TokenValue *= 100;
                              } else FetchFloat(0,100);

                              if( TokenValue )
                              {   if( ReDrawFlag & RFRotate )
                                      PrepareTransform();

                                  ReDrawFlag |= (1<<option);
                                  if( done ) TokenValue = -TokenValue;
                                  DialValue[option] += TokenValue/18000.0;

                                  while( DialValue[option]<-1.0 )
                                      DialValue[option] += 2.0;
                                  while( DialValue[option]>1.0 )
                                      DialValue[option] -= 2.0;
                                  if( Interactive )
                                      UpdateScrollBars();
                              }
                          } else CommandError(ErrorMsg[ErrNotNum]);
                          break;

        case(TranslateTok):
                          FetchToken();
                          if( CurToken==XTok )
                          {   option = 4;
                          } else if( CurToken==YTok )
                          {   option = 5;
                          } else if( CurToken==ZTok )
                          {   option = 6;
                          } else
                          {   CommandError(ErrorMsg[ErrSyntax]);
                              break;
                          }

                          FetchToken();
                          if( CurToken == '-' )
                          {   FetchToken();
                              done = True;
                          } else done = False;
#ifdef INVERT
                          if( option == 5 )
                              done = !done;
#endif

                          if( (CurToken==NumberTok) || (CurToken=='.') )
                          {   if( CurToken==NumberTok )
                              {   if( *TokenPtr=='.' )
                                  {   TokenPtr++;
                                      FetchFloat(TokenValue,100);
                                  } else TokenValue *= 100;
                              } else FetchFloat(0,100);

                              if( TokenValue<=10000 )
                              {   ReDrawFlag |= (1<<option);
                                  if( done ) TokenValue = -TokenValue;
                                  DialValue[option] = TokenValue/10000.0;
                                  /* UpdateScrollBars(); */
                              } else CommandError(ErrorMsg[ErrBigNum]);
                          } else CommandError(ErrorMsg[ErrNotNum]);
                          break;

        case(StereoTok):  FetchToken();
                          if( !CurToken || (CurToken==TrueTok) )
                          {   SetStereoMode(True);
                          } else if( CurToken==FalseTok )
                          {   SetStereoMode(False);
                          } else if( CurToken == '-' )
                          {   if( !NextIf(NumberTok,ErrNotNum) )
                              {   StereoAngle = -TokenValue;
                                  SetStereoMode(True);
                              }
                          } else if( CurToken==NumberTok )
                          {   StereoAngle = TokenValue;
                              SetStereoMode(True);
                          } else CommandError(ErrorMsg[ErrSyntax]);
                          break;

        case(ResizeTok):  FetchToken();
                          break;

        case(ResetTok):   for( i=0; i<8; i++ )
                              DialValue[i] = 0.0;
                          ReDrawFlag |= RFDials;
                          ResetTransform();

                          /* ReDrawFlag |= RFRefresh|RFColour; */
                          /* DisplayMode = 0;                  */

                          if( Interactive )
                              UpdateScrollBars();
                          break;

        case('?'):
        case(HelpTok):    if( !HelpFileName )
                              InitHelpFile();
                          if( HelpInfo )
                              FindHelpInfo();
                          CurToken=0;
                          break;

        case(LabelTok):   FetchToken();
                          if( !CurToken || (CurToken==TrueTok) )
                          {   if( Info.chaincount>1 )
                              {   DefineLabels("%n%r:%c.%a");
                              } else if( MainGroupCount>1 )
                              {   DefineLabels("%n%r.%a");
                              } else DefineLabels("%e%i");
                          } else if( CurToken==FalseTok )
                          {   DeleteLabels();
                          } else if( CurToken!=StringTok )
                          {   DefineLabels(TokenStart);
                              CurToken = 0;
                          } else DefineLabels(TokenIdent);
                          ReDrawFlag |= RFRefresh;
                          break;

        case(EchoTok):    FetchToken();
                          InvalidateCmndLine();
                          if( CurToken==StringTok )
                          {   WriteString(TokenIdent);
                          } else if( CurToken )
                              WriteString(TokenStart);
                          WriteChar('\n');
                          CurToken = 0;
                          break;

        case(DefineTok):  FetchToken();
                          if( CurToken != IdentTok ) 
                          {   CommandError(ErrorMsg[ErrSetName]);
                              break;
                          }

                          param = (char*)malloc(TokenLength+1);
                          if( param )
                          {   memcpy(param,TokenIdent,TokenLength+1);
                              if( FetchToken() )
                              {   QueryExpr = ParseExpression(0);
                                  if( QueryExpr )
                                  {   done = DefineSetExpr(param,QueryExpr);
                                  } else done = True;
                              } else done = DefineSetExpr(param,(Expr*)0);
                          } else done = False;

                          if( !done )
                              CommandError(ErrorMsg[ErrBadSet]);
                          break;

        case(BackgroundTok):
                          FetchToken();
                          if( CurToken == TransparentTok )
                          {   UseTransparent = True;
                          } else if( CurToken == NormalTok )
                          {   UseTransparent = False;
                          } else if( ParseColour() )
                          {   ReDrawFlag |= RFColour;
                              BackR = RVal;
                              BackG = GVal;
                              BackB = BVal;
#ifndef IBMPC
                              FBClear = False;
#endif
                          } else if( CurToken )
                          {      CommandError(ErrorMsg[ErrColour]);
                          } else CommandError(ErrorMsg[ErrNoCol]);
                          break;

        case(WriteTok):
        case(SaveTok):    i = CurToken; /* Save keyword! */
                          if( !AllowWrite )
                              if( (FileDepth!=-1) && LineStack[FileDepth] )
                              {   CommandError(ErrorMsg[ErrInScrpt]);
                                  break;
                              }

                          option = FetchToken();
                          if( (option==RasMolTok)
                              || (option==ScriptTok)
                              || (option==GraspTok)
                              || IsMoleculeToken(option)
                              || IsImageToken(option) )
                          {   if( !*TokenPtr || *TokenPtr==' ' )
                                  FetchToken();
                          } else if( i==SaveTok )
                          {   option = PDBTok;
                          } else option = 0;

                          if( !CurToken )
                          {   CommandError(ErrorMsg[ErrFilNam]);
                              break;
                          } else if( CurToken==StringTok )
                          {      ProcessFileName(TokenIdent);
                          } else ProcessFileName(TokenStart);
                          param = DataFileName;
                          CurToken = 0;

                          if( !IsMoleculeToken(option) )
                          {   if( ReDrawFlag ) RefreshScreen();
                              WriteImageFile( param, option );

                          } else switch(option)
                          {   case(NMRPDBTok):
                              case(PDBTok):  SavePDBMolecule(param);  break;
                              case(MDLTok):  SaveMDLMolecule(param);  break;
                              case(XYZTok):  SaveXYZMolecule(param);  break;
                              case(CIFTok):  SaveCIFMolecule(param);  break;
                              case(Mol2Tok): SaveMol2Molecule(param); break;
                              case(AlchemyTok): SaveAlchemyMolecule(param);
                                                break;
#ifdef CEXIOLIB
                              case(CEXTok):  SaveCEXMolecule(param); break;
#endif
                          } break;

        case(SourceTok):
        case(ScriptTok):  FetchToken();
                          if( FileDepth<STACKSIZE )
                          {   if( !CurToken )
                              {   CommandError(ErrorMsg[ErrFilNam]);
                                  break;
                              } else if( CurToken==StringTok )
                              {      ProcessFileName(TokenIdent);
                              } else ProcessFileName(TokenStart);
                              CurToken = 0;

                              script = fopen(DataFileName,"rb");
                              LoadScriptFile(script,DataFileName);
                          } else CommandError(ErrorMsg[ErrScript]);
                          break;

        case(RenumTok):   FetchToken();
                          if( CurToken == ChainTok ) {
                              RenumberChains();
                          } else if( CurToken )
                          {   if( CurToken == '-' )
                              {   FetchToken();
                                  done = True;
                              } else done = False;

                              if( CurToken == NumberTok )
                              {   if( done )
                                  {     RenumberAtoms(-(int)TokenValue);
                                  } else RenumberAtoms((int)TokenValue); 
                              } else CommandError(ErrorMsg[ErrNotNum]);
                          } else RenumberAtoms(1);
                          break;

        case(StructureTok):
                          FetchToken();
                          if( !CurToken || (CurToken==FalseTok) )
                          {   DetermineStructure(False);
                          } else if( CurToken==TrueTok )
                          {   DetermineStructure(True);
                          } else CommandError(ErrorMsg[ErrSyntax]);
                          break;

        case(ExitTok):    return( ExitTok );
        case(QuitTok):    return( QuitTok );
        default:          CommandError("Unrecognised command");
                          break;
    }

    if( CurToken )
        if( FetchToken() )
            CommandError("Warning: Ignoring rest of command");
    TokenPtr = (char*)0;
    return False;
}


int ExecuteIPCCommand( char __huge *ptr )
{
    auto char buffer[256];
    register int stat,result;
    register int len,depth;
    register char *dst;

    /* Ignore IPC commands while paused! */
    if( IsPaused ) return 0;

    FileDepth = 0;
    *LineStack = 0;
    result = IPC_Ok;
    *NameStack = "IPC Error";

    /* Save command line */
    strcpy(buffer,CurLine);

    while( *ptr && (*ptr==' ') )
        ptr++;

    if( *ptr == '[' )
    {   while( *ptr++ == '[' )
        {   dst = CurLine;
            depth = 0;  
            len = 0;

            while( *ptr )
            {   if( *ptr == ']' )
                {   if( !depth )
                    {   ptr++; 
                        break;
                    } else depth--;
                } else if( *ptr=='[' )
                    depth++;

                if( len < MAXBUFFLEN-1 )
                {   *dst++ = *ptr++;
                    len++;
                } else ptr++;
            }
            *dst = '\0';

            if( len < MAXBUFFLEN-1 )
            {   stat = ExecuteCommand();
                if( stat == QuitTok )
                {  result = IPC_Quit;
                } else if( stat )
                   result = IPC_Exit;
            } else
            {   InvalidateCmndLine();
                WriteString("Warning: Remote command too long!\n");
            }

            while( *ptr && ((*ptr==' ')||(*ptr==';')) )
                ptr++;
        }
    } else if( *ptr )
    {   dst = CurLine;
        len = 0;
        while( *ptr )
        {   if( len < MAXBUFFLEN-1 )
            {   *dst++ = *ptr++;
                len++;
            } else break;
        }
        *dst = '\0';

        if( len < MAXBUFFLEN-1 )
        {   stat = ExecuteCommand();
            if( stat == QuitTok )
            {  result = IPC_Quit;
            } else if( stat )
               result = IPC_Exit;
        } else
        {   InvalidateCmndLine();
            WriteString("Warning: Remote command too long!\n");
        }
    }

    FileDepth = -1;
    if( CommandActive )
    {   strcpy(CurLine,buffer);
        if( !result ) result = IPC_Error;
    }
    return result;
}


void InitialiseCommand( void )
{
    HelpFileName = (char*)0;
    FreeInfo = (HlpEntry __far*)0;
    HelpInfo = (HlpEntry __far*)0;

    TokenPtr = (char*)0;
    SelectCount = 0;
    FileDepth = -1;

    SeqFormat = False;
    IsPaused = False;
}

