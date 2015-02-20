/* molecule.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#define MAXMASK 40
#define MAXELEM 1024
#define MINELEM 29
#define MAXRES  100
#define MINRES  54


#define IsAmino(x)       ((x)<=23)
#define IsAminoNucleo(x) ((x)<=42)
#define IsNucleo(x)      (((x)>=24) && ((x)<=42))
#define IsProtein(x)     (((x)<=23) || (((x)>=43) && ((x)<=45)))
#define IsDNA(x)         (((x)>=24) && ((x)<=27))
#define IsSolvent(x)     (((x)>=46) && ((x)<=49))
#define IsWater(x)       (((x)==46) || ((x)==47))
#define IsIon(x)         (((x)==48) || ((x)==49))

#define IsPyrimidine(x)  (IsCytosine(x) || IsThymine(x))
#define IsPurine(x)      (IsAdenine(x) || IsGuanine(x))
#define IsRNA(x)         (IsNucleo(x) && !IsThymine(x))
#define NucleicCompl(x)  ((x)^3)


#define IsProline(x)     ((x)==11)
#define IsCysteine(x)    ((x)==17)
#define IsAdenine(x)     ((x)==24)
#define IsCytosine(x)    ((x)==25)
#define IsGuanine(x)     ((x)==26)
#define IsThymine(x)     ((x)==27)



#define IsAlphaCarbon(x)     ((x)==1)
#define IsSugarPhosphate(x)  ((x)==7)
#define IsAminoBackbone(x)   ((x)<=3)
#define IsShapelyBackbone(x) ((x)<=7)
#define IsNucleicBackbone(x) (((x)>=7) && ((x)<=18))
#define IsShapelySpecial(x)  ((x)==19)
#define IsCysteineSulphur(x) ((x)==20)
#define IsCoenzyme(x)        (((x)>=50) && ((x)<=53))


/*=================*/
/*  Database Flags */
/*=================*/

#define SelectFlag      0x01
#define DrawBondFlag    0x0e
#define AllAtomFlag     0x1c
#define HelixFlag       0x03
#define DrawKnotFlag    0x7e
#define WideKnotFlag    0x0e

/* Atom Flags */
#define SphereFlag      0x02     /* Sphere representation */
#define HeteroFlag      0x04     /* HETATM record         */
#define HydrogenFlag    0x08     /* Hydrogen atom         */
#define NormAtomFlag    0x10
#define NonBondFlag     0x20
#define BreakFlag       0x40     /* Break in backbone     */

/* Bond Flags */
#define WireFlag        0x02     /* Depth-cued wireframe         */
#define DashFlag        0x04     /* Dashed Depth-cued wireframe  */
#define CylinderFlag    0x08     /* Line/Cylinder representation */

#define HydrBondFlag    0x00     /* Hydrogen bond [place keeper] */
#define NormBondFlag    0x10
#define DoubBondFlag    0x20
#define TripBondFlag    0x40
#define AromBondFlag    0x80


/* Group Flags */
#define CystineFlag     0x01     /* Disulphide bonded cysteine  */
#define StrandFlag      0x02     /* Strands representation      */
#define DashStrandFlag  0x04     /* Dash Strands representation */
#define RibbonFlag      0x08     /* Solid Ribbon representation */
#define TraceFlag       0x10     /* Smooth trace representation */
#define CartoonFlag     0x20     /* Richardson protein cartoon  */
#define DotsFlag        0x40     /* Dotted trace representation */

/* Structure Flags */
#define Helix3Flag      0x01     /* 3,10-Helix structure       */
#define Helix4Flag      0x02     /* Alpha Helix structure      */
#define Helix5Flag      0x03     /* 5-Helix structure          */
#define SheetFlag       0x04     /* Beta Sheet structure       */
#define TurnFlag        0x08     /* Turn Secondary structure   */


/*=====================*/
/*  Molecule Database  */
/*=====================*/

typedef struct _Atom {
        struct _Atom __far *anext;        /* Linked list of atoms  */
        struct _Atom __far *bucket;       /* Sphere Y-Bucket       */
        struct _Atom __far *next;         /* Active Object List    */
        Long   xorg, yorg, zorg;          /* World Co-ordinates    */
        short  x, y, z;                   /* Image Co-ordinates    */
        short  radius;                    /* World Radius          */
        short  temp;                      /* Temperature Factor    */
        short  col;                       /* Atom Colour           */
        Long   serno;                     /* Atom Serial Number    */
        void   *label;                    /* Atom Label Structure  */
        unsigned char elemno;             /* Atomic Number         */
        unsigned char refno;              /* ElemDesc index number */
        unsigned char flag;               /* Database flags        */
        char   altl;                      /* Alternate Location    */
        short  irad;                      /* Image Radius          */
        short  mbox;                      /* Shadow Casting NOnce  */
    } Atom;


typedef struct _Bond {
        struct _Bond __far *bnext;       /* Linked list of bonds  */
        Atom __far *srcatom;             /* Source Atom Ptr       */
        Atom __far *dstatom;             /* Destination Atom Ptr  */
        short radius;                    /* World Radius          */
        short irad;                      /* Image Radius          */
        short col;                       /* Bond Colour           */
        unsigned char flag;              /* Database flags        */
    } Bond;

typedef struct _Group {
        struct _Group __far *gnext;       /* Linked list of groups */
        Atom __far *alist;                /* Linked list of atoms  */
        short serno;                      /* Group serial number   */
        short width;                      /* Ribbon Width          */
        short col1;                       /* Ribbon Colour #1      */
        short col2;                       /* Ribbon Colour #2      */
        char  insert;                     /* PDB insertion code    */
        unsigned char refno;              /* Residue index number  */
        unsigned char struc;              /* Secondary Structure   */
        unsigned char flag;               /* Database flags        */
    } Group;
 
#ifdef APPLEMAC
/* Avoid Name Clash! */
#define Chain ChainSeg
#endif

typedef struct _ChainSeg {
        struct _ChainSeg __far *cnext;    /* Linked list of chains     */
        Group __far *glist;               /* Linked list of groups     */
        Bond __far *blist;                /* Linked list of back bonds */
        short model;                      /* NMR Model / Symmetry      */
        char ident;                       /* Chain identifier          */
    } Chain;

typedef struct _AtomRef {
        Chain __far *chn;
        Group __far *grp;
        Atom  __far *atm;
    } AtomRef;

typedef struct _HBond {
        struct _HBond __far *hnext;       /* Ordered list of hbonds   */
        Atom __far *srcCA;                /* Source Alpha Carbon      */
        Atom __far *dstCA;                /* Destination Alpha Carbon */
        Atom __far *dst;                  /* Acceptor [=CO] Atom Ptr  */
        Atom __far *src;                  /* Donor [=NH] Atom Ptr     */
        short energy;                     /* Hydrogen bond energy     */
        short radius;                     /* World Radius             */
        short irad;                       /* Image Radius             */
        Char offset;                      /* Signed Offset            */
        unsigned char flag;               /* Database flags           */
        unsigned char col;                /* Hydrogen bond colour     */
    } HBond;

typedef struct _Molecule {
        HBond __far *slist;               /* Linked list of SS bonds  */
        HBond __far *hlist;               /* Linked list of hbonds    */
        Chain __far *clist;               /* Linked list of chains    */
        Bond __far *blist;                /* Linked list of bonds     */
    } Molecule;



/*========================*/
/* Other Consts & Structs */
/*========================*/

#define SourceNone   0
#define SourcePDB    1
#define SourceCalc   2
 
#define SerNoFlag 0x01
#define ResNoFlag 0x02

typedef struct {
        short radius;
        char  mask[19];
        unsigned char flags;
        unsigned char r;
        unsigned char g;
        unsigned char b;
    } MaskDesc;

typedef struct _IntCoord {
        struct _IntCoord __far *inext;
        short na,nb,nc;
        short refno;
        Real dihed;
        Real angle;
        Real dist;
    } IntCoord;

typedef struct _InfoStruct {
        char filename[1024];
        char moleculename[80];
        char classification[42];
        char identcode[6];

        char spacegroup[11];
        Real cellalpha, cellbeta, cellgamma;
        Real cella, cellb, cellc;

        Long bondcount;
        int chaincount;
        int ssbondcount;
        int hbondcount;

        int structsource;
        int laddercount;
        int helixcount;
        int turncount;
    } InfoStruct;



#ifdef MOLECULE
/* Avoid SGI Compiler Warnings! */
char Residue[MAXRES][4] = {
    /*===============*/
    /*  Amino Acids  */
    /*===============*/

/* Ordered by Cumulative Frequency in Brookhaven *
 * Protein Databank, December 1991               */

          "ALA", /* 8.4% */     "GLY", /* 8.3% */
          "LEU", /* 8.0% */     "SER", /* 7.5% */
          "VAL", /* 7.1% */     "THR", /* 6.4% */
          "LYS", /* 5.8% */     "ASP", /* 5.5% */
          "ILE", /* 5.2% */     "ASN", /* 4.9% */
          "GLU", /* 4.9% */     "PRO", /* 4.4% */
          "ARG", /* 3.8% */     "PHE", /* 3.7% */
          "GLN", /* 3.5% */     "TYR", /* 3.5% */
          "HIS", /* 2.3% */     "CYS", /* 2.0% */
          "MET", /* 1.8% */     "TRP", /* 1.4% */

          "ASX", "GLX", "PCA", "HYP",

    /*===================*/
    /*  DNA Nucleotides  */
    /*===================*/
          "  A", "  C", "  G", "  T",

    /*===================*/
    /*  RNA Nucleotides  */
    /*===================*/
          "  U", " +U", "  I", "1MA", 
          "5MC", "OMC", "1MG", "2MG", 
          "M2G", "7MG", "OMG", " YG", 
          "H2U", "5MU", "PSU",

    /*=================*/
    /*  Miscellaneous  */ 
    /*=================*/
          "UNK", "ACE", "FOR", "HOH",
          "DOD", "SO4", "PO4", "NAD",
          "COA", "NAP", "NDP"  };


/* Avoid SGI Compiler Warnings! */
char ElemDesc[MAXELEM][4] = {
    { ' ', 'N', ' ', ' ' },  /* 0*/
    { ' ', 'C', 'A', ' ' },  /* 1*/
    { ' ', 'C', ' ', ' ' },  /* 2*/
    { ' ', 'O', ' ', ' ' },  /* 3*/   /* 0-3   Amino Acid Backbone    */
    { ' ', 'C', '\'', ' ' }, /* 4*/
    { ' ', 'O', 'T', ' ' },  /* 5*/
    { ' ', 'S', ' ', ' ' },  /* 6*/
    { ' ', 'P', ' ', ' ' },  /* 7*/   /* 4-7   Shapely Amino Backbone */
    { ' ', 'O', '1', 'P' },  /* 8*/
    { ' ', 'O', '2', 'P' },  /* 9*/
    { ' ', 'O', '5', '*' },  /*10*/
    { ' ', 'C', '5', '*' },  /*11*/
    { ' ', 'C', '4', '*' },  /*12*/
    { ' ', 'O', '4', '*' },  /*13*/
    { ' ', 'C', '3', '*' },  /*14*/
    { ' ', 'O', '3', '*' },  /*15*/
    { ' ', 'C', '2', '*' },  /*16*/
    { ' ', 'O', '2', '*' },  /*17*/
    { ' ', 'C', '1', '*' },  /*18*/   /* 7-18  Nucleic Acid Backbone  */
    { ' ', 'C', 'A', '2' },  /*19*/   /* 19    Shapely Special        */
    { ' ', 'S', 'G', ' ' },  /*20*/   /* 20    Cysteine Sulphur       */
    { ' ', 'N', '1', ' ' },  /*21*/
    { ' ', 'N', '2', ' ' },  /*22*/
    { ' ', 'N', '3', ' ' },  /*23*/
    { ' ', 'N', '4', ' ' },  /*24*/
    { ' ', 'N', '6', ' ' },  /*25*/
    { ' ', 'O', '2', ' ' },  /*26*/
    { ' ', 'O', '4', ' ' },  /*27*/
    { ' ', 'O', '6', ' ' }   /*28*/   /* 21-28 Nucleic Acid H-Bonding */
    };


InfoStruct Info;
int MainGroupCount;
int HetaGroupCount;
Long MainAtomCount; 
Long HetaAtomCount;

Long MinX, MinY, MinZ;
Long MaxX, MaxY, MaxZ;

int HMinMaxFlag, MMinMaxFlag;
int MinMainTemp, MaxMainTemp;
int MinHetaTemp, MaxHetaTemp;
int MinMainRes,  MaxMainRes;
int MinHetaRes,  MaxHetaRes;

Molecule __far *CurMolecule;
Chain __far *CurChain;
Group __far *CurGroup;
Atom __far *CurAtom;

IntCoord __far *IntList;
Molecule __far *Database;
MaskDesc UserMask[MAXMASK];
Long MinHBondDist, MaxHBondDist;
Long MinBondDist,  MaxBondDist;
int ElemNo,ResNo;
int HasHydrogen;
int MaskCount;
int NMRModel;

int HBondChainsFlag;

#else
extern char Residue[MAXRES][4];
extern char ElemDesc[MAXELEM][4];
extern InfoStruct Info;

extern int MainGroupCount;
extern int HetaGroupCount;
extern Long MainAtomCount;
extern Long HetaAtomCount;

extern Long MinX, MinY, MinZ;
extern Long MaxX, MaxY, MaxZ;

extern int HMinMaxFlag, MMinMaxFlag;
extern int MinMainTemp, MaxMainTemp;
extern int MinHetaTemp, MaxHetaTemp;
extern int MinMainRes,  MaxMainRes;
extern int MinHetaRes,  MaxHetaRes;

extern Molecule __far *CurMolecule;
extern Chain __far *CurChain;
extern Group __far *CurGroup;
extern Atom __far *CurAtom;

extern IntCoord __far *IntList;
extern Molecule __far *Database;
extern MaskDesc UserMask[MAXMASK];
extern Long MinHBondDist, MaxHBondDist;
extern Long MinBondDist,  MaxBondDist;
extern int ElemNo,ResNo;
extern int HasHydrogen;
extern int MaskCount;
extern int NMRModel;

extern int HBondChainsFlag;
#endif


void CreateChain( int );
void CreateGroup( int );
void ProcessGroup( int );
void CreateMolGroup( void );
void CreateNextMolGroup( void );
int FindResNo( char* );

Atom __far *CreateAtom( void );
Atom __far *FindGroupAtom( Group __far*, int );
void ProcessAtom( Atom __far* );

int NewAtomType( char* );
int SimpleAtomType( char* );
int ComplexAtomType( char* );

Bond __far *ProcessBond( Atom __far*, Atom __far*, int );
void CreateBond( Long, Long, int );
void CreateBondOrder( Long, Long );
void CreateMoleculeBonds( int, int );
void FindDisulphideBridges( void );
void CalcHydrogenBonds( void );

void InitInternalCoords( void );
IntCoord __far* AllocInternalCoord( void );
int ConvertInternal2Cartesian( void );
void FreeInternalCoords( void );

void DetermineStructure( int );
void RenumberChains( void );
void RenumberAtoms( int );

void InitialiseDatabase( void );
void DescribeMolecule( void );
void DestroyDatabase( void );
void PurgeDatabase( void );

