/* infile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

int FetchRecord( FILE*, char* );

int LoadPDBMolecule( FILE*, int );
int LoadMacroModelMolecule( FILE* );
int LoadAlchemyMolecule( FILE* );
int LoadCharmmMolecule( FILE* );
int LoadBiosymMolecule( FILE* );
int LoadMOPACMolecule( FILE* );
int LoadSHELXMolecule( FILE* );
int LoadMol2Molecule( FILE* );
int LoadFDATMolecule( FILE* );
int LoadMDLMolecule( FILE* );
int LoadXYZMolecule( FILE* );
int LoadCEXMolecule( FILE* );

int SaveAlchemyMolecule( char* );
int SavePDBMolecule( char* );
int SaveMDLMolecule( char* );
int SaveXYZMolecule( char* );
int SaveCIFMolecule( char* );
int SaveCEXMolecule( char* );

