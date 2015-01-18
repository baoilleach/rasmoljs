/* prepdoc.c
 * Document Preparation System
 * Roger Sayle, January 1994
 * Version 1.1
 */

#include <string.h>
#include <stdio.h>
#include <ctype.h>

#ifndef True
#define True  1
#define False 0
#endif

#define RTFForm     0x00
#define LaTeXForm   0x01
#define HTMLForm    0x02
#define TextForm    0x03
#define HelpForm    0x04
#define ManForm     0x05
#define VaxForm     0x06



static char Buffer[514];
static FILE *OutFile;
static FILE *InFile;
static int SplitFlag;
static int Format;


static int ReadLine()
{
    register char *ptr;
    register int len;
    register int ch;

    if( feof(InFile) )
    {   *Buffer = 0;
        return( False );
    }

    ptr = Buffer;
    do {
        ch = getc(InFile);
        if( (ch=='\n') || (ch=='\r') )
        {   while( ptr != Buffer )
                if( *(ptr-1)!=' ' )
                {   *ptr = 0;
                    return( True );
                } else ptr--;
        } else if( ch==EOF )
        {   *ptr = 0;
            return( True );
        } else *ptr++ = ch;
    } while( ptr < Buffer+512 );
    *ptr = 0;

    /* skip to the end of the line! */
    do { ch = getc(InFile);
    } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );
    return( True );
}


static int TextFlag;
static int TextCol;

static void DisplayLowerCase( ptr )
    char *ptr;
{
    register char ch;

    while( (ch = *ptr++) )
        if( isupper(ch) )
        {   putc(tolower(ch),OutFile);
        } else putc(ch,OutFile);
}

static void OpenLowerCase( ptr )
    char *ptr;
{
    register char *dst;
    register char ch;
    char buffer[256];

    fputs("\n<p>\n",OutFile);
    if( OutFile != stdout )
        fclose(OutFile);

    dst = buffer;
    while( (ch = *ptr++) )
        if( isupper(ch) )
        {   *dst++ = tolower(ch);
        } else if( ch != ' ' )
            *dst++ = ch;
    strcpy(dst,".html");

    if( !(OutFile = fopen(buffer,"w")) )
    {   fprintf(stderr,"Warning: Unable to create file `%s'!\n",buffer);
        OutFile = stdout;
    }
}


static void DisplayOnlyLowerCase( ptr )
    char *ptr;
{
    register char ch;
    
    while( (ch = *ptr++) )
        if( isupper(ch) )
        {   putc(tolower(ch),OutFile);
        } else if( ch != ' ' )
            putc(ch,OutFile);
}
    

static void DisplayText( src )
    char *src;
{
    register char *ptr;
    register char *dst;
    register int max;

    max = (Format==VaxForm)? 68 : 76;

    while( *src )
    {   dst = src;
        while( *dst && *dst!=' ' )
            dst++;

        if( TextCol+(dst-src) > max )
        {   putc('\n',OutFile);
            TextCol = 0;
        }

        if( !TextCol && (Format==VaxForm) )
            putc(' ',OutFile);

        while( src != dst )
        {   putc( *src++, OutFile );
            TextCol++;
        }

        while( *src==' ' )
            src++;
        putc(' ',OutFile);
        TextCol++;
    }
}

static char *VAXSection( ptr )
    char *ptr;
{
    register char *tmp;

    if( !strncmp(ptr,"Set ",4) )
        ptr += 4;

    for( tmp=ptr; *tmp; tmp++ )
        if( *tmp==' ' )
            *tmp = '_';
    return( ptr );
}


static void SplitRawHTML( ptr )
    char *ptr;
{
    register char *src;
    register char *dst;
    register int state;
    register int flag;
    char buffer[80];

    flag = True;
    while( *ptr )
        if( *ptr == '<' )
        {   if( ptr[1]=='a' )
            {   /* Anchor Format Problems? */
                if( !strncmp(ptr,"<a name=\"",9) )
                {   dst = buffer;
                    for( src = ptr+9; *src && (*src!='"'); src++ )
                        *dst++ = *src;
                    *dst = '\0';

                    OpenLowerCase(buffer);
                    fputs("<title>RasMol Manual</title>",OutFile);

                    if( *src )
                    {   ptr = src+1;
                        while( *ptr && (*ptr!='>') )
                            ptr++; 
                        if( *ptr ) ptr++;
                    } else ptr = src;
                    flag = False;
                } else if( !strncmp(ptr,"<a href=\"#",10) )
                {   state = False;
                    fputs("<a href=\"",OutFile);
                    for( ptr+=10; *ptr && (*ptr!='"'); ptr++ )
                        putc(*ptr,OutFile);
                    fputs(".html\"",OutFile);
                    ptr++;
                    flag = True;

                } else /* Unrecognised anchor */
                {   putc(*ptr,OutFile);
                   ptr++;
                }
            } else if( (ptr[1]=='/') && (ptr[2]=='a') && (ptr[3]=='>') )
            {   if( flag ) fputs("</a>",OutFile);
                ptr += 4;
            } else /* Unrecognised code */
            {   putc(*ptr,OutFile);
                ptr++;
            }
        } else  /* Normal Character! */
        {   putc(*ptr,OutFile);
            ptr++;
        }

    /* EndOfLine Char */
    putc('\n',OutFile);
}

static void FilterHTMLChar( ch )
    char ch;
{
    if( ch == '&' )
    {   fprintf(OutFile,"&amp;");
    } else if( ch == '<' )
    {   fprintf(OutFile,"&lt;");
    } else if( ch == '>' )
    {   fprintf(OutFile,"&gt;");
    } else putc(ch,OutFile);
}


static void ProcessCommand()
{
    static char buffer[80];
    register char *ptr;
    register int i,len;

    ptr = Buffer+2;
    switch( Buffer[1] )
    {   case('R'):  if( TextCol )
                        putc('\n',OutFile);
                    if( SplitFlag )
                    {   SplitRawHTML(ptr);
                    } else fprintf(OutFile,"%s\n",ptr);
                    TextFlag = True;
                    TextCol = 0;
                    break;

        case('T'):  if( Format == RTFForm )
                    {   fprintf(OutFile,"%s {}\n",ptr);
                        TextCol = 0;
                    } else if( Format == HTMLForm )
                    {   while( *ptr )
                            FilterHTMLChar(*ptr++);
                        putc('\n',OutFile);
                        TextCol = 0;
                    } else if( (Format!=TextForm) &&
                               (Format!=HelpForm) &&
                               (Format!=VaxForm) )
                    {   fprintf(OutFile,"%s\n",ptr);
                        TextCol = 0;
                    } else DisplayText(ptr);
                    TextFlag = True;
                    break;

        case('P'):  if( TextFlag )
                    {   if( Format == HTMLForm )
                        {   fputs("<p>\n",OutFile);
                        } else if( Format==RTFForm )
                        {   fputs("\\par\\par\n",OutFile);
                        } else if( (Format==TextForm) || 
                                   (Format==HelpForm) ||
                                   (Format==VaxForm) )
                        {   fputs("\n\n",OutFile);
                        } else putc('\n',OutFile);
                        TextFlag = False;
                        TextCol = 0;
                    }
                    break;

        case('B'):  if( TextFlag )
                    {   if( Format==HTMLForm )
                        {   if( !SplitFlag )
                                fputs("<p><hr>\n",OutFile);
                        } else if( Format==RTFForm )
                        {   fputs("\\par\\page\\par\n",OutFile);
                            fputs("+{\\footnote doc}\n",OutFile);
                        } else if( (Format==TextForm) || 
                                   (Format==HelpForm) ||
                                   (Format==VaxForm) )
                        {   fputs("\n\n",OutFile);
                        } else putc('\n',OutFile);
                        TextFlag = False;
                        TextCol = 0;
                    }
                    break;

        case('S'):  if( Format==TextForm )
                    {   len = strlen(ptr);
                        fprintf(OutFile,"%s\n",ptr);
                        for( i=0; i<len; i++ )
                            putc('-',OutFile);
                        fputs("\n\n",OutFile);
                    } else if( Format==HTMLForm )
                    {   if( SplitFlag )
                        {   OpenLowerCase(ptr);
                            fputs("<title>RasMol Manual: ",OutFile);
                            fprintf(OutFile,"%s</title>\n",ptr);
                            fprintf(OutFile,"<h1>%s</h1><p>\n",ptr);
                        } else
                        {   fputs("<a name=\"",OutFile);
                            DisplayOnlyLowerCase(ptr);
                            fprintf(OutFile,"\"><h3>%s</h3></a><p>\n",ptr);
                        }
                    } else if( Format==RTFForm )
                    {   fputs("#{\\footnote ",OutFile);
                        DisplayOnlyLowerCase(ptr);
                        fprintf(OutFile,"}\n${\\footnote %s}\n",ptr);
                        fputs("K{\\footnote ",OutFile);
                        DisplayLowerCase(ptr);
                        fprintf(OutFile,"}\n{\\b %s}\\par\\par\n",ptr);
                    } else if( Format==HelpForm )
                    {   putc('?',OutFile);
                        DisplayLowerCase(ptr);
                        fprintf(OutFile,"\n%s\n",ptr);
                    } else if( Format==ManForm )
                    {   fprintf(OutFile,".TP\n.B %s\n",ptr);
                    } else if( Format==VaxForm )
                        fprintf(OutFile,"3 %s\n",VAXSection(ptr));

                    TextCol = 0;
                    break;

        case('X'):  while( *ptr!=' ' )
                        ptr++;
                    *ptr++ = 0;
                    
                    if( Format == RTFForm )
                    {   fprintf(OutFile,"{\\uldb %s}",ptr);
                        fprintf(OutFile,"{\\v %s} {}\n",Buffer+2);
                        break;
                    } else if( Format == HTMLForm )
                    {   if( SplitFlag )
                        {   fprintf(OutFile,"<a href=\"%s.html\">",Buffer+2);
                        } else fprintf(OutFile,"<a href=\"#%s\">",Buffer+2);
                        fprintf(OutFile,"<tt><b>%s</b></tt></a>\n",ptr);
                        break;
                    }
                    
        case('C'):  if( Format == RTFForm )
                    {   if( *ptr=='"' )
                        {   fprintf(OutFile,"\"{\\f2\\b %s}\" {}\n",ptr+1);
                        } else fprintf(OutFile,"{\\f2\\b %s} {}\n",ptr);
                        TextCol = 0;
                    } else if( Format == HTMLForm )
                    {   if( *ptr=='"' )
                        {   fprintf(OutFile,"\"<tt><b>%s</b></tt>\"\n",ptr+1);
                        } else fprintf(OutFile,"<tt><b>%s</b></tt>\n",ptr);
                    } else if( Format == ManForm )
                    {   if( *ptr=='"') ptr++;
                        fprintf(OutFile,".B %s\n",ptr);
                    } else if( (Format!=TextForm) &&
                               (Format!=HelpForm) &&
                               (Format!=VaxForm) )
                    {   if( *ptr=='*' )
                        {   fprintf(OutFile,"\"%s\"\n",ptr);
                        } else fprintf(OutFile,"`%s'\n",ptr);
                        TextCol = 0;
                    } else /* DisplayText! */
                    {   if( *ptr=='"' )
                        {   sprintf(buffer,"\"%s\"",ptr+1);
                        } else sprintf(buffer,"`%s'",ptr);
                        DisplayText(buffer);
                    }
                    TextFlag = True;
                    break;
    }
}


int main( argc, argv )
    int argc;  char *argv[];
{
    register char *fname;
    register char *ptr;
    register int flag;

    fputs("Document Preparation System\n",stderr);
    fputs("Roger Sayle, January 1994\n",stderr);
    fputs("Version 1.1\n\n",stderr);

    Format = TextForm;
    SplitFlag = False;
    OutFile = stdout;

    if( argc==2 ) 
    {   fname = argv[1];
    } else if( argc==3 )
    {   fname = argv[2];
        ptr = argv[1];

        if( *ptr=='-' )
            ptr++;

        if( !strcmp(ptr,"latex") )
        {   Format = LaTeXForm;
        } else if( !strcmp(ptr,"help") )
        {   Format = HelpForm;
        } else if( !strcmp(ptr,"html") )
        {   Format = HTMLForm;
        } else if( !strcmp(ptr,"splithtml") )
        {   Format = HTMLForm; SplitFlag = True;
        } else if( !strcmp(ptr,"rtf") || !strcmp(ptr,"mshelp") )
        {   Format = RTFForm;
        } else if( !strcmp(ptr,"text") || !strcmp(ptr,"ascii") )
        {   Format = TextForm;
        } else if( !strcmp(ptr,"man") || !strcmp(ptr,"troff") )
        {   Format = ManForm;
        } else if( !strcmp(ptr,"vax") || !strcmp(ptr,"vms") )
        {   Format = VaxForm;
        } else
        {   fputs("Formats:  -latex  LaTeX .tex file\n",stderr);
            fputs("          -troff  UNIX man(1) pages\n",stderr);
            fputs("          -html   HyperText metalanguage\n",stderr);
            fputs("          -help   RasMol on-line help file\n",stderr);
            fputs("          -rtf    Microsoft Help (Rich Text)\n",stderr);
            fputs("          -text   Standard ASCII text\n\n",stderr);
            fputs("          -vax    VAX VMS Help file\n\n",stderr);
            exit(1);
        }

    } else /* DisplayUsage */
    {   fputs("Usage: prepdoc [format] <filename>\n",stderr);
        exit(1);
    }

    if( !(InFile=fopen(fname,"r")) )
    {   fputs("Error: Unable to open input file!\n",stderr);
        exit(1);
    }

    TextFlag = False;
    TextCol = 0;

    while( !feof(InFile) )
    {   ReadLine();
        switch( *Buffer )
        {   case('V'):  flag = (Format==LaTeXForm) ||
                               (Format==HTMLForm)  ||
                               (Format==RTFForm)   ||
                               (Format==TextForm);   break;

            case('D'):  flag = (Format==TextForm) ||
                               (Format==HelpForm) ||
                               (Format==ManForm);    break;

            case('S'):  flag = (Format==HelpForm) ||
                               (Format==ManForm);    break;

            case('N'):  flag = (Format==TextForm) ||
                               (Format==HelpForm);   break;

            case('L'):  flag = (Format==LaTeXForm);  break;
            case('H'):  flag = (Format==HTMLForm);   break;
            case('P'):  flag = (Format==HelpForm);   break;
            case('T'):  flag = (Format==TextForm);   break;
            case('M'):  flag = (Format==ManForm);    break;
            case('R'):  flag = (Format==RTFForm);    break;
            case('X'):  flag = (Format==VaxForm);    break;
            case('A'):  flag = True;                 break;
            default:    flag = False;
        }

        if( flag )
            ProcessCommand();
    }
    if( OutFile != stdout )
    {   fputs("\n<p>\n",OutFile);
        fclose(OutFile);
    }
    fclose(InFile);
    exit(0);
}

