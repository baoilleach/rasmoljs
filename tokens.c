/* tokens.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include <string.h>
#include <ctype.h>

#include "tokens.h"

int LookUpKeyword( char *ptr )
{
    switch(*ptr++) {
        case('A'):
            switch(*ptr++) {
                case('C'):
                    if( !strcmp(ptr,"IDIC") ) {
                        return( AcidicTok );
                    } else if( !strcmp(ptr,"YCLIC") ) {
                        return( AcyclicTok );
                    }
                    break;

                case('L'):
                    if( !strcmp(ptr,"CHEMY") ) {
                        return( AlchemyTok );
                    } else if( !strcmp(ptr,"IPHATIC") ) {
                        return( AliphaticTok );
                    } else if( (*ptr=='L') && !ptr[1] ) {
                        return( AllTok );
                    } else if( !strcmp(ptr,"PHA") ) {
                        return( AlphaTok );
                    }
                    break;

                case('M'):
                    if( !strcmp(ptr,"BIENT") ) {
                        return( AmbientTok );
                    } else if( !strcmp(ptr,"INO") ) {
                        return( AminoTok );
                    }
                    break;

                case('N'):
                    if( (*ptr=='D') && !ptr[1] ) {
                        return( AndTok );
                    } else if( !strcmp(ptr,"GLE") ) {
                        return( AngleTok );
                    } else if( !strcmp(ptr,"GLES") ) {
                        return( AngleTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"OMATIC") ) {
                        return( AromaticTok );
                    }
                    break;

                case('T'):
                    if( !*ptr ) {
                        return( ATTok );
                    } else if( !strcmp(ptr,"OM") ) {
                        return( AtomTok );
                    } else if( !strcmp(ptr,"OMNO") ) {
                        return( AtomNoTok );
                    } else if( !strcmp(ptr,"OMS") ) {
                        return( AtomTok );
                    }
                    break;

                case('X'):
                    if( !strcmp(ptr,"ES") ) {
                        return( AxesTok );
                    } else if( !strcmp(ptr,"IS") ) {
                        return( AxesTok );
                    }
                    break;

            }
            break;

        case('B'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"CKBONE") ) {
                        return( BackboneTok );
                    } else if( !strcmp(ptr,"CKFADE") ) {
                        return( BackFadeTok );
                    } else if( !strcmp(ptr,"CKGROUND") ) {
                        return( BackgroundTok );
                    } else if( !strcmp(ptr,"SIC") ) {
                        return( BasicTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"OSYM") ) {
                        return( BiosymTok );
                    }
                    break;

                case('L'):
                    if( !strcmp(ptr,"ACK") ) {
                        return( BlackTok );
                    } else if( !strcmp(ptr,"UE") ) {
                        return( BlueTok );
                    }
                    break;

                case('M'):
                    if( (*ptr=='P') && !ptr[1] ) {
                        return( BMPTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"ND") ) {
                        return( BondTok );
                    } else if( !strcmp(ptr,"NDED") ) {
                        return( BondedTok );
                    } else if( !strcmp(ptr,"NDMODE") ) {
                        return( BondModeTok );
                    } else if( !strcmp(ptr,"NDS") ) {
                        return( BondTok );
                    } else if( !strcmp(ptr,"UNDBOX") ) {
                        return( BoundBoxTok );
                    } else if( !strcmp(ptr,"UNDINGBOX") ) {
                        return( BoundBoxTok );
                    }
                    break;

                case('U'):
                    if( !strcmp(ptr,"RIED") ) {
                        return( BuriedTok );
                    }
                    break;

            }
            break;

        case('C'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"RTOON") ) {
                        return( CartoonTok );
                    } else if( !strcmp(ptr,"RTOONS") ) {
                        return( CartoonTok );
                    }
                    break;

                case('E'):
                    if( !strcmp(ptr,"NTER") ) {
                        return( CentreTok );
                    } else if( !strcmp(ptr,"NTRE") ) {
                        return( CentreTok );
                    } else if( (*ptr=='X') && !ptr[1] ) {
                        return( CEXTok );
                    }
                    break;

                case('G'):
                    if( !*ptr ) {
                        return( CGTok );
                    }
                    break;

                case('H'):
                    if( !strcmp(ptr,"AIN") ) {
                        return( ChainTok );
                    } else if( !strcmp(ptr,"AINS") ) {
                        return( ChainTok );
                    } else if( !strcmp(ptr,"ARGE") ) {
                        return( ChargeTok );
                    } else if( !strcmp(ptr,"ARGED") ) {
                        return( ChargedTok );
                    } else if( !strcmp(ptr,"ARGES") ) {
                        return( ChargeTok );
                    } else if( !strcmp(ptr,"ARMM") ) {
                        return( CharmmTok );
                    }
                    break;

                case('I'):
                    if( (*ptr=='F') && !ptr[1] ) {
                        return( CIFTok );
                    }
                    break;

                case('L'):
                    if( !strcmp(ptr,"IPBOARD") ) {
                        return( ClipboardTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"LOR") ) {
                        return( ColourTok );
                    } else if( !strcmp(ptr,"LORS") ) {
                        return( ColourTok );
                    } else if( !strcmp(ptr,"LOUR") ) {
                        return( ColourTok );
                    } else if( !strcmp(ptr,"LOURS") ) {
                        return( ColourTok );
                    } else if( !strcmp(ptr,"NNECT") ) {
                        return( ConnectTok );
                    }
                    break;

                case('P'):
                    if( (*ptr=='K') && !ptr[1] ) {
                        return( CPKTok );
                    }
                    break;

                case('Y'):
                    if( !strcmp(ptr,"AN") ) {
                        return( CyanTok );
                    } else if( !strcmp(ptr,"CLIC") ) {
                        return( CyclicTok );
                    } else if( !strcmp(ptr,"STINE") ) {
                        return( CystineTok );
                    }
                    break;

            }
            break;

        case('D'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"SH") ) {
                        return( DashTok );
                    } else if( !strcmp(ptr,"SHES") ) {
                        return( DashTok );
                    }
                    break;

                case('E'):
                    if( !strcmp(ptr,"FINE") ) {
                        return( DefineTok );
                    } else if( !strcmp(ptr,"PTHCUE") ) {
                        return( DepthCueTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"SPLAY") ) {
                        return( DisplayTok );
                    } else if( !strcmp(ptr,"STANCE") ) {
                        return( DistanceTok );
                    } else if( !strcmp(ptr,"STANCES") ) {
                        return( DistanceTok );
                    }
                    break;

                case('N'):
                    if( (*ptr=='A') && !ptr[1] ) {
                        return( DNATok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"TS") ) {
                        return( DotsTok );
                    }
                    break;

                case('X'):
                    if( (*ptr=='F') && !ptr[1] ) {
                        return( DXFTok );
                    }
                    break;
            }
            break;

        case('E'):
            switch(*ptr++) {
                case('C'):
                    if( !strcmp(ptr,"HO") ) {
                        return( EchoTok );
                    }
                    break;

                case('L'):
                    if( !strcmp(ptr,"EMNO") ) {
                        return( ElemNoTok );
                    }
                    break;

                case('P'):
                    if( !strcmp(ptr,"SF") ) {
                        return( EPSFTok );
                    }
                    break;

                case('X'):
                    if( !strcmp(ptr,"IT") ) {
                        return( ExitTok );
                    }
                    break;

            }
            break;

        case('F'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"LSE") ) {
                        return( FalseTok );
                    }
                    break;

                case('D'):
                    if( !strcmp(ptr,"AT") ) {
                        return( FDATTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"NTSIZE") ) {
                        return( FontSizeTok );
                    }
                    break;

            }
            break;

        case('G'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"USSIAN") ) {
                        return( GaussianTok );
                    }
                    break;

                case('I'):
                    if( (*ptr=='F') && !ptr[1] ) {
                        return( GIFTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"ASP") ) {
                        return( GraspTok );
                    } else if( !strcmp(ptr,"EEN") ) {
                        return( GreenTok );
                    } else if( !strcmp(ptr,"EENBLUE") ) {
                        return( GreenBlueTok );
                    } else if( !strcmp(ptr,"OUP") ) {
                        return( GroupTok );
                    }
                    break;

            }
            break;

        case('H'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"LF") ) {
                        return( HalfTok );
                    } else if( !strcmp(ptr,"RDWARE") ) {
                        return( HardwareTok );
                    }
                    break;

                case('B'):
                    if( !strcmp(ptr,"OND") ) {
                        return( HBondTok );
                    } else if( !strcmp(ptr,"ONDS") ) {
                        return( HBondTok );
                    }
                    break;

                case('E'):
                    if( !strcmp(ptr,"LICES") ) {
                        return( HelixTok );
                    } else if( !strcmp(ptr,"LIX") ) {
                        return( HelixTok );
                    } else if( !strcmp(ptr,"LP") ) {
                        return( HelpTok );
                    } else if( !strcmp(ptr,"TERO") ) {
                        return( HeteroTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"GH") ) {
                        return( HighTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"LLOW") ) {
                        return( HollowTok );
                    } else if( !strcmp(ptr,"URGLASS") ) {
                        return( HourGlassTok );
                    }
                    break;

                case('Y'):
                    if( !strcmp(ptr,"DROGEN") ) {
                        return( HydrogenTok );
                    } else if( !strcmp(ptr,"DROPHOBIC") ) {
                        return( HydrophobicTok );
                    }
                    break;

            }
            break;

        case('I'):
            switch(*ptr++) {
                case('D'):
                    if( !strcmp(ptr,"ENT") ) {
                        return( IdentifyTok );
                    } else if( !strcmp(ptr,"ENTIFY") ) {
                        return( IdentifyTok );
                    }
                    break;

                case('N'):
                    if( !strcmp(ptr,"FO") ) {
                        return( InfoTok );
                    } else if( !strcmp(ptr,"FORMATION") ) {
                        return( InfoTok );
                    } else if( !strcmp(ptr,"LINE") ) {
                        return( InLineTok );
                    } else if( !strcmp(ptr,"SIGHT") ) {
                        return( InsightTok );
                    }
                    break;

                case('O'):
                    if( (*ptr=='N') && !ptr[1] ) {
                        return( IonTok );
                    } else if( !strcmp(ptr,"NS") ) {
                        return( IonTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"IS") ) {
                        return( IRISTok );
                    }
                    break;

            }
            break;

        case('J'):
            if( !strcmp(ptr,"PEG") ) {
                return( JPEGTok );
            }
            break;

        case('K'):
            if( !strcmp(ptr,"INEMAGE") ) {
                return( KinemageTok );
            }
            break;

        case('L'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"BEL") ) {
                        return( LabelTok );
                    } else if( !strcmp(ptr,"BELS") ) {
                        return( LabelTok );
                    } else if( !strcmp(ptr,"RGE") ) {
                        return( LargeTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"GAND") ) {
                        return( LigandTok );
                    } else if( !strcmp(ptr,"GANDS") ) {
                        return( LigandTok );
                    } else if( !strcmp(ptr,"NE") ) {
                        return( LineTok );
                    } else if( !strcmp(ptr,"NES") ) {
                        return( LineTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"AD") ) {
                        return( LoadTok );
                    } else if( (*ptr=='W') && !ptr[1] ) {
                        return( LowTok );
                    }
                    break;

            }
            break;

        case('M'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"CROMODEL") ) {
                        return( MacroModelTok );
                    } else if( !strcmp(ptr,"GENTA") ) {
                        return( MagentaTok );
                    } else if( !strcmp(ptr,"INCHAIN") ) {
                        return( MainChainTok );
                    }
                    break;

                case('D'):
                    if( (*ptr=='L') && !ptr[1] ) {
                        return( MDLTok );
                    }
                    break;

                case('E'):
                    if( !strcmp(ptr,"DIUM") ) {
                        return( MediumTok );
                    } else if( !strcmp(ptr,"NUS") ) {
                        return( MenusTok );
                    }
                    break;

                case('M'):
                    if( !strcmp(ptr,"DB") ) {
                        return( MMDBTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"DEL") ) {
                        return( ModelTok );
                    } else if( !strcmp(ptr,"L2") ) {
                        return( Mol2Tok );
                    } else if( !strcmp(ptr,"LSCRIPT") ) {
                        return( MolScriptTok );
                    } else if( !strcmp(ptr,"NITOR") ) {
                        return( MonitorTok );
                    } else if( !strcmp(ptr,"NITORS") ) {
                        return( MonitorTok );
                    } else if( !strcmp(ptr,"NO") ) {
                        return( MonoTok );
                    } else if( !strcmp(ptr,"NOCHROME") ) {
                        return( MonoTok );
                    } else if( !strcmp(ptr,"NOPS") ) {
                        return( MonoPSTok );
                    } else if( !strcmp(ptr,"PAC") ) {
                        return( MOPACTok );
                    } else if( !strcmp(ptr,"USE") ) {
                        return( MouseTok );
                    } else if( !strcmp(ptr,"USEMODE") ) {
                        return( MouseTok );
                    }
                    break;

            }
            break;

        case('N'):
            switch(*ptr++) {
                case('E'):
                    if( !strcmp(ptr,"GATIVE") ) {
                        return( AcidicTok );
                    } else if( !strcmp(ptr,"UTRAL") ) {
                        return( NeutralTok );
                    }
                    break;

                case('M'):
                    if( !strcmp(ptr,"RPDB") ) {
                        return( NMRPDBTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"NE") ) {
                        return( NoneTok );
                    } else if( !strcmp(ptr,"RMAL") ) {
                        return( NormalTok );
                    } else if( (*ptr=='T') && !ptr[1] ) {
                        return( NotTok );
                    }
                    break;

                case('U'):
                    if( !strcmp(ptr,"CLEIC") ) {
                        return( NucleicTok );
                    }
                    break;

            }
            break;

        case('O'):
            switch(*ptr++) {
                case('F'):
                    if( (*ptr=='F') && !ptr[1] ) {
                        return( FalseTok );
                    }
                    break;

                case('N'):
                    if( !*ptr ) {
                        return( TrueTok );
                    }
                    break;

                case('R'):
                    if( !*ptr ) {
                        return( OrTok );
                    } else if( !strcmp(ptr,"ANGE") ) {
                        return( OrangeTok );
                    } else if( !strcmp(ptr,"IGIN") ) {
                        return( OriginTok );
                    }
                    break;

            }
            break;

        case('P'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"USE") ) {
                        return( WaitTok );
                    }
                    break;

                case('D'):
                    if( (*ptr=='B') && !ptr[1] ) {
                        return( PDBTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"CKING") ) {
                        return( PickingTok );
                    } else if( !strcmp(ptr,"CT") ) {
                        return( PICTTok );
                    }
                    break;

                case('L'):
                    if( (*ptr=='Y') && !ptr[1] ) {
                        return( PLYTok );
                    }
                    break;

                case('N'):
                    if( (*ptr=='G') && !ptr[1] ) {
                        return( PNGTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"LAR") ) {
                        return( PolarTok );
                    } else if( !strcmp(ptr,"SITIVE") ) {
                        return( BasicTok );
                    } else if( !strcmp(ptr,"TENTIAL") ) {
                        return( PotentialTok );
                    } else if( !strcmp(ptr,"VRAY") ) {
                        return( POVRayTok );
                    }
                    break;

                case('P'):
                    if( (*ptr=='M') && !ptr[1] ) {
                        return( PPMTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"INT") ) {
                        return( PrintTok );
                    } else if( !strcmp(ptr,"OTEIN") ) {
                        return( ProteinTok );
                    }
                    break;

                case('S'):
                    if( !*ptr ) {
                        return( EPSFTok );
                    }
                    break;

                case('U'):
                    if( !strcmp(ptr,"RINE") ) {
                        return( PurineTok );
                    } else if( !strcmp(ptr,"RINES") ) {
                        return( PurineTok );
                    } else if( !strcmp(ptr,"RPLE") ) {
                        return( PurpleTok );
                    }
                    break;

                case('Y'):
                    if( !strcmp(ptr,"RIMIDINE") ) {
                        return( PyrimidineTok );
                    } else if( !strcmp(ptr,"RIMIDINES") ) {
                        return( PyrimidineTok );
                    }
                    break;

            }
            break;

        case('Q'):
            if( *ptr++ == 'U' ) {
                if( !strcmp(ptr,"ANTA") ) {
                    return( QuantaTok );
                } else if( !strcmp(ptr,"IT") ) {
                    return( QuitTok );
                }
            }
            break;

        case('R'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"DIUS") ) {
                        return( RadiusTok );
                    } else if( !strcmp(ptr,"SMOL") ) {
                        return( RasMolTok );
                    } else if( !strcmp(ptr,"SWIN") ) {
                        return( RasMolTok );
                    }
                    break;

                case('E'):
                    if( (*ptr=='D') && !ptr[1] ) {
                        return( RedTok );
                    } else if( !strcmp(ptr,"DORANGE") ) {
                        return( RedOrangeTok );
                    } else if( !strcmp(ptr,"FRESH") ) {
                        return( RefreshTok );
                    } else if( !strcmp(ptr,"JECT") ) {
                        return( RejectTok );
                    } else if( !strcmp(ptr,"NUM") ) {
                        return( RenumTok );
                    } else if( !strcmp(ptr,"NUMBER") ) {
                        return( RenumTok );
                    } else if( !strcmp(ptr,"SET") ) {
                        return( ResetTok );
                    } else if( !strcmp(ptr,"SIDUE") ) {
                        return( ResidueTok );
                    } else if( !strcmp(ptr,"SIZE") ) {
                        return( ResizeTok );
                    } else if( !strcmp(ptr,"SNO") ) {
                        return( ResNoTok );
                    } else if( !strcmp(ptr,"STRICT") ) {
                        return( RestrictTok );
                    }
                    break;

                case('G'):
                    if( (*ptr=='B') && !ptr[1] ) {
                        return( IRISTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"BBON") ) {
                        return( RibbonTok );
                    } else if( !strcmp(ptr,"BBON1") ) {
                        return( Ribbon1Tok );
                    } else if( !strcmp(ptr,"BBON2") ) {
                        return( Ribbon2Tok );
                    } else if( !strcmp(ptr,"BBONS") ) {
                        return( RibbonTok );
                    } else if( !strcmp(ptr,"BBONS1") ) {
                        return( Ribbon1Tok );
                    } else if( !strcmp(ptr,"BBONS2") ) {
                        return( Ribbon2Tok );
                    }
                    break;

                case('N'):
                    if( (*ptr=='A') && !ptr[1] ) {
                        return( RNATok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"TATE") ) {
                        return( RotateTok );
                    }
                    break;

            }
            break;

        case('S'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"VE") ) {
                        return( SaveTok );
                    }
                    break;

                case('C'):
                    if( !strcmp(ptr,"RIPT") ) {
                        return( ScriptTok );
                    }
                    break;

                case('E'):
                    if( !strcmp(ptr,"CTION") ) {
                        return( SectionTok );
                    } else if( !strcmp(ptr,"LECT") ) {
                        return( SelectTok );
                    } else if( !strcmp(ptr,"LECTED") ) {
                        return( SelectedTok );
                    } else if( !strcmp(ptr,"QUENCE") ) {
                        return( SequenceTok );
                    } else if( (*ptr=='T') && !ptr[1] ) {
                        return( SetTok );
                    }
                    break;

                case('H'):
                    if( !strcmp(ptr,"ADOW") ) {
                        return( ShadowTok );
                    } else if( !strcmp(ptr,"ADOWS") ) {
                        return( ShadowTok );
                    } else if( !strcmp(ptr,"APELY") ) {
                        return( ShapelyTok );
                    } else if( !strcmp(ptr,"EET") ) {
                        return( SheetTok );
                    } else if( !strcmp(ptr,"EETS") ) {
                        return( SheetTok );
                    } else if( !strcmp(ptr,"ELX") ) {
                        return( SHELXTok );
                    } else if( !strcmp(ptr,"OW") ) {
                        return( ShowTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"DECHAIN") ) {
                        return( SidechainTok );
                    }
                    break;

                case('L'):
                    if( !strcmp(ptr,"AB") ) {
                        return( SlabTok );
                    } else if( !strcmp(ptr,"ABMODE") ) {
                        return( SlabModeTok );
                    }
                    break;

                case('M'):
                    if( !strcmp(ptr,"ALL") ) {
                        return( SmallTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"LID") ) {
                        return( SolidTok );
                    } else if( !strcmp(ptr,"LVENT") ) {
                        return( SolventTok );
                    } else if( !strcmp(ptr,"LVENTS") ) {
                        return( SolventTok );
                    } else if( !strcmp(ptr,"URCE") ) {
                        return( SourceTok );
                    }
                    break;

                case('P'):
                    if( !strcmp(ptr,"ACEFILL") ) {
                        return( SpacefillTok );
                    } else if( !strcmp(ptr,"ECPOWER") ) {
                        return( SpecPowerTok );
                    } else if( !strcmp(ptr,"ECULAR") ) {
                        return( SpecularTok );
                    }
                    break;

                case('S'):
                    if( !strcmp(ptr,"BOND") ) {
                        return( SSBondTok );
                    } else if( !strcmp(ptr,"BONDS") ) {
                        return( SSBondTok );
                    }
                    break;

                case('T'):
                    if( !strcmp(ptr,"EREO") ) {
                        return( StereoTok );
                    } else if( (*ptr=='L') && !ptr[1] ) {
                        return( STLTok );
                    } else if( !strcmp(ptr,"RANDS") ) {
                        return( StrandsTok );
                    } else if( !strcmp(ptr,"RUCTURE") ) {
                        return( StructureTok );
                    }
                    break;

                case('U'):
                    if( (*ptr=='N') && !ptr[1] ) {
                        return( SUNTok );
                    } else if( !strcmp(ptr,"NRLE") ) {
                        return( SUNRLETok );
                    } else if( !strcmp(ptr,"RFACE") ) {
                        return( SurfaceTok );
                    }
                    break;

                case('Y'):
                    if( !strcmp(ptr,"BYL") ) {
                        return( SybylTok );
                    } else if( !strcmp(ptr,"MMETRY") ) {
                        return( SymmetryTok );
                    }
                    break;

            }
            break;

        case('T'):
            switch(*ptr++) {
                case('E'):
                    if( !strcmp(ptr,"MPERATURE") ) {
                        return( TemperatureTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"TLE") ) {
                        return( TitleTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"RSION") ) {
                        return( TorsionTok );
                    } else if( !strcmp(ptr,"RSIONS") ) {
                        return( TorsionTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"ACE") ) {
                        return( TraceTok );
                    } else if( !strcmp(ptr,"ANSLATE") ) {
                        return( TranslateTok );
                    } else if( !strcmp(ptr,"ANSLUCENT") ) {
                        return( TranslucentTok );
                    } else if( !strcmp(ptr,"ANSPARENT") ) {
                        return( TransparentTok );
                    } else if( !strcmp(ptr,"UE") ) {
                        return( TrueTok );
                    }
                    break;

                case('U'):
                    if( !strcmp(ptr,"RN") ) {
                        return( TurnTok );
                    } else if( !strcmp(ptr,"RNS") ) {
                        return( TurnTok );
                    }
                    break;

                case('Y'):
                    if( !strcmp(ptr,"PE") ) {
                        return( TypeTok );
                    }
                    break;

            }
            break;

        case('U'):
            switch(*ptr++) {
                case('N'):
                    if( !strcmp(ptr,"ITCELL") ) {
                        return( UnitCellTok );
                    }
                    break;

                case('S'):
                    if( !strcmp(ptr,"ER") ) {
                        return( UserTok );
                    }
                    break;

            }
            break;

        case('V'):
            switch(*ptr++) {
                case('D'):
                    if( (*ptr=='W') && !ptr[1] ) {
                        return( VDWTok );
                    }
                    break;

                case('E'):
                    if( !strcmp(ptr,"CTPS") ) {
                        return( VectPSTok );
                    } else if( !strcmp(ptr,"CTOR")) {
                        return( VectorTok );
                    } else if( !strcmp(ptr,"CTORS")) {
                        return( VectorTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"OLET") ) {
                        return( VioletTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"ML") ) {
                        return( VRMLTok );
                    }
                    break;

            }
            break;

        case('W'):
            switch(*ptr++) {
                case('A'):
                    if( !strcmp(ptr,"IT") ) {
                        return( WaitTok );
                    } else if( !strcmp(ptr,"TER") ) {
                        return( WaterTok );
                    } else if( !strcmp(ptr,"TERS") ) {
                        return( WaterTok );
                    }
                    break;

                case('H'):
                    if( !strcmp(ptr,"ITE") ) {
                        return( WhiteTok );
                    }
                    break;

                case('I'):
                    if( !strcmp(ptr,"REFRAME") ) {
                        return( WireframeTok );
                    } else if( !strcmp(ptr,"THIN") ) {
                        return( WithinTok );
                    }
                    break;

                case('R'):
                    if( !strcmp(ptr,"ITE") ) {
                        return( WriteTok );
                    }
                    break;

            }
            break;

        case('X'):
            switch(*ptr++) {
                case(0):
                    return( XTok );

                case('Y'):
                    if( (*ptr=='Z') && !ptr[1] ) {
                        return( XYZTok );
                    }
                    break;

            }
            break;

        case('Y'):
            switch(*ptr++) {
                case(0):
                    return( YTok );

                case('E'):
                    if( !strcmp(ptr,"LLOW") ) {
                        return( YellowTok );
                    }
                    break;

            }
            break;

        case('Z'):
            switch(*ptr++) {
                case(0):
                    return( ZTok );

                case('A'):
                    if( (*ptr=='P') && !ptr[1] ) {
                        return( ZapTok );
                    }
                    break;

                case('O'):
                    if( !strcmp(ptr,"OM") ) {
                        return( ZoomTok );
                    }
                    break;

            }
            break;

    }
    return( IdentTok );
}

