/* cexio.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include <stdlib.h>
#include <stdio.h>

#include "rasmol.h"
#include "molecule.h"
#include "command.h"
#include "graphics.h"

#include "cx.h"
#include "cx_molecule.h"


#define CEXMODELS
#define CEXTITLES


#ifdef CEXTITLES
#include "cx_message.h"
#endif


#define PROP_ATSYM   "atomic symbol"
#define PROP_BONDO   "bond order"

#define PROP_MNAME   "molname"
#define PROP_ALABS   "atomname"
#define PROP_RESID   "resname"
#define PROP_CHAIN   "chain"
#define PROP_RESNO   "resno"
#define PROP_COORD   "coordinates"
#define PROP_BVALU   "bvalue"
#define PROP_CHARG   "charge"
#define PROP_RADIU   "radius"


static cx_String prop_alabs;
static cx_String prop_resid;
static cx_String prop_chain;
static cx_String prop_resno;
static cx_String prop_coord;
static cx_String prop_tempr;

static Atom __far *CEXConnectAtom;
static int CEXSerialNumber;


static cx_String get_atom_tuple( cx_Object mol, cx_String prop )
{
    cx_Object tupleseq, tuple;

    if( (tupleseq = cx_prefix2atuples(mol,prop)) )
    {   tuple = cx_next(tupleseq);
    } else tuple = (cx_Object)0;
 
    if( !tuple )
    {   /* No such property! */
        cx_destroy( tupleseq );
        return (cx_String)0;
    }
 
    /* Ambiguous atom property prefix! */
    /* cx_next(tupleseq) or !cx_atend(tupleseq) */

    cx_destroy(tupleseq);
    return cx_atomtuple_name(tuple);
}
 

static void padfield( char *dst, char *src, int len )
{
    register int i;

    for( i=0; i<len; i++ )
    {   if( *src )
        {   *dst++ = *src++;
        } else *dst++ = ' ';
    }
    *dst = '\0';
}


static int process_cx_atom( cx_Object atom, int serno )
{
    auto char buffer[5];
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *ptr;
    register int refno,resno;
    register Long dx,dy,dz;
    register int chain;
    register char *prop;
    auto double x,y,z,t;

    /* Test for valid 3D co-ordinates */
    prop = cx_sprop(atom,prop_coord);
    if( !prop ) return False;
    if( sscanf(prop,"%lg,%lg,%lg",&x,&y,&z) != 3 )
        return False;

    /* Determine Chain */
    prop = cx_sprop(atom,prop_chain);
    if( prop && *prop && (*prop!=',') ) 
    {   chain = *prop;
    } else chain = ' ';

    if( CurChain )  
    {   if( (CurChain->ident!=chain) || (CurChain->model!=NMRModel) )
        {   CEXConnectAtom = (Atom __far*)0;
            CurChain = CurMolecule->clist;
            while( CurChain )
                if( (CurChain->ident!=chain) || (CurChain->model!=NMRModel) )
                {   CurChain = CurChain->cnext;
                } else break;
        }
    }

    if( !CurChain ) {
        CEXConnectAtom = (Atom __far*)0;
        CreateChain( chain );
    }

    /* Determine Residue */
    resno = cx_iprop(atom,prop_resno);
    prop = cx_sprop(atom,prop_resid);
    if( !prop || !*prop ) 
    {   refno = FindResNo("MOL");
    } else refno = FindResNo(prop);

    if( !CurGroup || CurGroup->refno!=refno || CurGroup->serno!=resno )
    {   group = CurChain->glist;
        CurGroup = (Group __far*)0;
        CurAtom = (Atom __far*)0;

        while( group && group->serno <= resno )
        {   CurGroup = group;
            if( (group->serno==resno) && (group->refno==refno) ) break;
            group = group->gnext;
        }

        if( !group || (group->serno>resno) ) 
        {   CreateGroup( 8 );
            CurGroup->serno = resno;
            CurGroup->refno = refno;
            ProcessGroup( False );
        }
    }

    /* Process Atom */
    ptr = CreateAtom();
    if( !IsProtein(refno) && !IsNucleo(refno) )
        ptr->flag |= HeteroFlag;
    ptr->serno = serno;

    prop = cx_sprop(atom,prop_alabs);
    if( prop )
    {   padfield(buffer,prop,4);
        ptr->refno = ComplexAtomType(buffer);
    } else
    {   prop = cx_sprop(atom,PROP_ATSYM);
        if( prop )
        {   padfield(buffer,prop,4);
            ptr->refno = SimpleAtomType(buffer);
        } else ptr->refno = SimpleAtomType("Du  ");
    }

    t = cx_rprop(atom,prop_tempr);
    ptr->temp = (short)(t*100.0);

    ptr->xorg =  (Long)(x*250.0);
    ptr->yorg =  (Long)(y*250.0);
    ptr->zorg = -(Long)(z*250.0);
    ProcessAtom( ptr );

    /* Create biopolymer Backbone */
    if( IsAlphaCarbon(ptr->refno) && IsProtein(CurGroup->refno) )
    {   if( CEXConnectAtom )
        {   dx = CEXConnectAtom->xorg - ptr->xorg;
            dy = CEXConnectAtom->yorg - ptr->yorg;
            dz = CEXConnectAtom->zorg - ptr->zorg;

            /* Break backbone if CA-CA > 4.20A */
            if( dx*dx+dy*dy+dz*dz < (Long)1050*1050 )
            {   bptr = ProcessBond(ptr,CEXConnectAtom,NormBondFlag);
                bptr->bnext = CurChain->blist;
                CurChain->blist = bptr;
            } else ptr->flag |= BreakFlag;
        }
        CEXConnectAtom = ptr;
    } else if( IsSugarPhosphate(ptr->refno) && IsNucleo(CurGroup->refno) )
    {   if( CEXConnectAtom )
        {   bptr = ProcessBond(CEXConnectAtom,ptr,NormBondFlag);
            bptr->bnext = CurChain->blist;
            CurChain->blist = bptr;
        }
        CEXConnectAtom = ptr;
    }
    return True;
}


static void process_cx_bond( cx_Object bond )
{
    cx_Object atoms,atom;
    register int src,dst;
    register int i,flag;

    i = cx_iprop(bond,PROP_BONDO);
    if( i==2 )
    {   flag = DoubBondFlag;
    } else if( i == 3 )
    {   flag = TripBondFlag;
    } else flag = NormBondFlag;

    atoms = cx_stream(bond,CX_OB_ATOM);
    if( (atom=cx_next(atoms)) != (cx_Object)0 )
        if( (src=cx_iprop(atom,"_serno")) != 0 )
            if( (atom=cx_next(atoms)) != (cx_Object)0 )
                if( (dst=cx_iprop(atom,"_serno")) != 0 )
                    CreateBond(src,dst,flag);
    cx_destroy(atoms);
}


static void load_cx_molecule( cx_Object mol )
{
    cx_String molname;
    cx_Object atoms,atom;
    cx_Object bonds,bond;
    register int i;

    if( !Database )
    {   if( (molname = cx_sprop(mol,PROP_MNAME)) )
        {   for( i=0; molname[i] && (i<62); i++ )
                Info.moleculename[i] = molname[i];
            while( i && molname[i-1]==' ' ) i--;
            Info.moleculename[i] = '\0';
        }
    }

    /* Molecule must have 3d co-ordinates */
    prop_coord = get_atom_tuple(mol,PROP_COORD);
    if( !prop_coord ) return;

    prop_alabs = get_atom_tuple(mol,PROP_ALABS);
    prop_resid = get_atom_tuple(mol,PROP_RESID);

    prop_tempr=get_atom_tuple(mol,PROP_BVALU);
    if( !prop_tempr ) prop_tempr = get_atom_tuple(mol,PROP_CHARG);
    if( !prop_tempr ) prop_tempr = get_atom_tuple(mol,PROP_RADIU);

    prop_chain = get_atom_tuple(mol,PROP_CHAIN);
    prop_resno = get_atom_tuple(mol,PROP_RESNO);

    CEXConnectAtom = (Atom __far*)0;
    atoms = cx_stream(mol,CX_OB_ATOM);
    while( (atom=cx_next(atoms)) != (cx_Object)0 )
    {   if( process_cx_atom(atom,CEXSerialNumber) )
        {   cx_set_iprop(atom,"_serno",CEXSerialNumber++);
        } else cx_set_iprop(atom,"_serno",0);
    }
    cx_destroy(atoms);

    bonds = cx_stream(mol,CX_OB_BOND);
    while( (bond=cx_next(bonds)) != (cx_Object)0 )
        process_cx_bond(bond);
    cx_destroy(bonds);
}


static void InitCEXMolecule( void )
{
    CEXSerialNumber = 1;
    NMRModel = 0;

    cx_molecule_pkg();
#ifdef CEXTITLES
    cx_message_pkg();
#endif
}


static void ProcessCEXObject( cx_Object obj )
{
#ifdef CEXTITLES
    register cx_Object sos;
    register cx_Object str;
    register char *ptr;
#endif

    if( cx_type(obj) == CX_OB_MOLECULE )
    {   
#ifdef CEXMODELS
        if( Database )
            NMRModel++;
        load_cx_molecule(obj);
#else
        if( !Database )
            load_cx_molecule(obj);
#endif
#ifdef CEXTITLES
    } else if( cx_type(obj) == CX_OB_MESSAGE )
    {   if( Interactive )
        {   ptr = cx_sprop(obj,"name");
            if( ptr && !strcmp(ptr,"title") )
            {   sos = cx_stream(obj,CX_OB_STRING);
                str = cx_next(sos);
                ptr = cx_stringvalue(str);
                if( ptr && *ptr )
                    SetCanvasTitle(ptr);
                cx_destroy(sos);
            }
        }
#endif
    }
}


int LoadCEXMolecule( FILE *fp )
{
#ifdef CX_OB_IOSTREAM
    register cx_Object ins;
    register cx_Object obj;

    UnusedArgument(fp);

    InitCEXMolecule();
    ins = cx_create_iostream(DataFileName,CX_IO_READ);
    while( (obj = cx_next(ins)) ) 
    {   ProcessCEXObject(obj);
        cx_destroy(obj);
    }
    cx_destroy(ins);
    /* cx_cleanup(); */
    return True;
    
#else
    register cx_Object obj;

    if( !fp ) return False;

    InitCEXMolecule();
    while( (obj = cx_receive((cx_Object)0,fp,(cx_Object)0)) ) 
    {   ProcessCEXObject(obj);
        cx_destroy(obj);
    }
    /* cx_cleanup(); */
    return True;
#endif
}
 
 
int SaveCEXMolecule( char *filename )
{
    UnusedArgument(filename);

    if( !Database )
        return False;
    return True;
}
 
