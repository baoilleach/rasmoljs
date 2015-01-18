/* outfile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#ifdef OUTFILE
int UseTransparent;
int UseOutLine;

#else
extern int UseTransparent;
extern int UseOutLine;
#endif


int WriteVectPSFile( char* );
int WriteEPSFFile( char*, int, int );
int WriteRastFile( char*, int );
int WritePICTFile( char* );
int WriteIRISFile( char* );
int WritePPMFile( char*, int );
int WriteGIFFile( char* );
int WriteBMPFile( char* );
void InitialiseOutFile( void );

