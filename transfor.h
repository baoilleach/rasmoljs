/* transfor.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#define GroupAttr       0x00
#define ChainAttr       0x01
#define TempAttr	0x02
#define ChargeAttr      0x03

#define MaskColourFlag  0x01
#define MaskRadiusFlag  0x02
#define MaskBothFlag    0x03

#define RibColInside    0x01
#define RibColOutside   0x02
#define RibColBoth      0x03

#ifdef EIGHTBIT
#define ColourDepth       16
#define ColourMask        15
#ifdef APPLEMAC
#define LastShade         14
#else
#define LastShade         15
#endif
#else
#define ColourDepth       32
#define ColourMask        31
#define LastShade         31
#endif


#ifdef __esv
/* Evans & Sutherland Gamma! */
#define DefaultAmbient    0.05
#else
#define DefaultAmbient    0.4
#endif


typedef struct { 
        Long refcount;
        unsigned char r;
        unsigned char g;
        unsigned char b;
    } ShadeDesc;


#ifdef IBMPC
#define BackCol    0
#define BoxCol     1
#define LabelCol   2
#define FirstCol   3
#endif

#ifdef APPLEMAC
#define BackCol    1
#define BoxCol     2
#define LabelCol   3
#define FirstCol   4
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
#define BackCol    5
#define BoxCol     6
#define LabelCol   7
#define FirstCol   8
#endif

#define Colour2Shade(x)  ((int)((x)-FirstCol)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+FirstCol)


#ifdef TRANSFORM
ShadeDesc Shade[LastShade];
Real RotX[3],RotY[3],RotZ[3];
Real MatX[3],MatY[3],MatZ[3];
Real InvX[3],InvY[3],InvZ[3];
Long OrigCX,OrigCY,OrigCZ;
Long CenX, CenY, CenZ;

int FakeSpecular,SpecPower;
int BackR,BackG,BackB;
int LabR,LabG,LabB;
int BoxR,BoxG,BoxB;
int UseLabelCol;
int UseBackFade;
Real Ambient;

Real Scale,MaxZoom;
Long SideLen,Offset;
Card WorldRadius,WorldSize;
int XOffset,YOffset,ZOffset;
int UseScreenClip;
int ZoomRange;

int Hydrogens,HetaGroups;
int DrawAtoms,MaxAtomRadius;
int DrawBonds,MaxBondRadius;
int DrawRibbon;
int ZoneBoth;

#else
extern ShadeDesc Shade[LastShade];
extern Real RotX[3],RotY[3],RotZ[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Long OrigCX, OrigCY, OrigCZ;
extern Long CenX, CenY, CenZ;

extern int FakeSpecular,SpecPower;
extern int BackR,BackG,BackB;
extern int LabR,LabG,LabB;
extern int BoxR,BoxG,BoxB;
extern int UseLabelCol;
extern int UseBackFade;
extern Real Ambient;

extern Real Scale,MaxZoom;
extern Long SideLen,Offset;
extern Card WorldRadius,WorldSize;
extern int XOffset,YOffset,ZOffset;
extern int UseScreenClip;
extern int ZoomRange;

extern int Hydrogens,HetaGroups;
extern int DrawAtoms,MaxAtomRadius;
extern int DrawBonds,MaxBondRadius;
extern int DrawRibbon;
extern int ZoneBoth;
#endif


void SetRadiusValue( int );
void SetRadiusTemperature( void );
void SetVanWaalRadius( void );
void DisableSpacefill( void );
void SetHBondStatus( int, int, int );
void SetRibbonStatus( int, int, int );
void SetRibbonCartoons( void );
void SetTraceTemperature( void );
void EnableWireframe( int, int );
void EnableBackbone( int, int );
void DisableWireframe( void );
void DisableBackbone( void );

void SelectZoneExpr( Expr* );
void RestrictZoneExpr( Expr* );
void RestrictZone( int );
void SelectZone( int );

int IsCPKColour( Atom __far * );
int IsVDWRadius( Atom __far * );

void DefineColourMap( void );
void ResetColourMap( void );
int DefineShade( int, int, int );

void ColourBackNone( void );
void ColourBondNone( void );
void ColourHBondType( void );
void ColourHBondNone( int );
void ColourRibbonNone( int );
void ColourMonitNone( void );
void ColourBackAttrib( int, int, int );
void ColourBondAttrib( int, int, int );
void ColourHBondAttrib( int, int, int, int );
void ColourRibbonAttrib( int, int, int, int );
void ColourMonitAttrib( int, int, int );
void ColourDotsAttrib( int, int, int );
void ColourDotsPotential( void );
void MonoColourAttrib( int, int, int );
void ScaleColourAttrib( int );
void CPKColourAttrib( void );
void AminoColourAttrib( void );
void ShapelyColourAttrib( void );
void StructColourAttrib( void );
void UserMaskAttrib( int );

void DefaultRepresentation( void );

void DetermineClipping( void );
void InitialiseTransform( void );
void InitialTransform( void );
void PrepareTransform( void );
void ReviseInvMatrix( void );
void ApplyTransform( void );
void ResetTransform( void );
void TranslateToCentre( Long, Long, Long );

