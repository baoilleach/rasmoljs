/* command.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

/* Format values are related to Tokens */
#define Tok2Format(x) ((x)-359)
#define Format2Tok(x) ((x)+359)

#define IsMoleculeFormat(x)  ((x)<=16)

#define FormatPDB        1
#define FormatMacroMod   2
#define FormatGaussian   3
#define FormatAlchemy    4
#define FormatNMRPDB     5
#define FormatCharmm     6
#define FormatBiosym     7
#define FormatMOPAC      8
#define FormatSHELX      9
#define FormatMol2      10
#define FormatFDAT      11
#define FormatMMDB      12
#define FormatMDL       13
#define FormatXYZ       14
#define FormatCIF       15
#define FormatCEX       16

#define FormatDots      20


#define IPC_Ok      0
#define IPC_Error   1
#define IPC_Exit    2
#define IPC_Quit    3

#ifdef COMMAND
int DataFileFormat;
char DataFileName[1024];
Long SelectCount;
int Interactive;
int FileDepth;
int IsPaused;

int CalcBondsFlag;
int AllowWrite;

#else
extern int DataFileFormat;
extern char DataFileName[1024];
extern Long SelectCount;
extern int Interactive;
extern int FileDepth;
extern int IsPaused;

extern int CalcBondsFlag;
extern int AllowWrite;
#endif


int FetchFile( int, int, char* );
int ProcessFile( int, int, FILE* );
void LoadScriptFile( FILE*, char* );

void InterruptPauseCommand( void );
void ResumePauseCommand( void );
int ExecuteIPCCommand( char __huge* );
int ExecuteCommand( void );

void InitialiseCommand( void );
void ZapDatabase( void );

