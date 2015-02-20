/* scripts.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#define ResDefault  0x00
#define ResLow      0x01
#define ResMedium   0x02
#define ResHigh     0x03

#ifdef SCRIPTS
int ResolutionFlag;
int KinemageFlag;
int BinaryFlag;

#else
extern int ResolutionFlag;
extern int KinemageFlag;
extern int BinaryFlag;
#endif


int WriteMolScriptFile( char* );
int WriteKinemageFile( char* );
int WriteScriptFile( char* );
int WritePOVRayFile( char* );
int WriteGraspFile( char* );
int WriteVRMLFile( char* );
int WriteSTLFile( char* );
int WriteDXFFile( char* );
int WritePLYFile( char* );

