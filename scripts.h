/* scripts.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#ifdef SCRIPTS
int KinemageFlag;

#else
extern int KinemageFlag;
#endif


int WriteMolScriptFile( char* );
int WriteKinemageFile( char* );
int WriteScriptFile( char* );
int WritePOVRayFile( char* );
int WriteVRMLFile( char* );

