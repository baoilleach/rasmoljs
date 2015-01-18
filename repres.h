/* repres.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#define DotMax    100
typedef struct _DotStruct {
        struct _DotStruct __far *next;
        short col[DotMax];
        Long xpos[DotMax];
        Long ypos[DotMax];
        Long zpos[DotMax];
        int count;
    } DotStruct;


typedef struct _Monitor {
        struct _Monitor *next;
        Atom __far *src;
        Atom __far *dst;
        unsigned short dist;
        short col;
    } Monitor;


typedef struct _Label {
        struct _Label *next;
        Long  refcount;
        char *label;
    } Label;



#ifdef REPRES
DotStruct __far *DotPtr;
Monitor *MonitList;
Label *LabelList;

int CartoonHeight;
int SolventDots;
int ProbeRadius;

int SurfaceChainsFlag;
int DrawMonitDistance;
int DrawBetaArrows;

#else
extern DotStruct __far *DotPtr;
extern Monitor *MonitList;
extern Label *LabelList;

extern int CartoonHeight;
extern int ProbeRadius;
extern int SolventDots;

extern int SurfaceChainsFlag;
extern int DrawMonitDistance;
extern int DrawBetaArrows;
#endif


int DeleteLabels( void );
void DeleteLabel( Label* );
Label *CreateLabel( char*, int );
void LabelTerminii( int );
void DefaultLabels( int );
void DefineLabels( char* );
void DisplayLabels( void );

void DeleteMonitors( void );
void AddMonitors( Atom __far*, Atom __far* );
void CreateMonitor( Long, Long );
void DisplayMonitors( void );

void DeleteSurface( void );
void CalculateSurface( int );
void DisplaySurface( void );
void LoadDotsFile( FILE*, int );

void DisplayRibbon( Chain __far* );

void ResetRepres( void );
void InitialiseRepres( void );

