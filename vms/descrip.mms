CC = CC
LINK = LINK

!LIBS = sys$share:decw$xlibshr/share
LIBS = 

!CFLAGS = /standard=vaxc/debug/noopt
CFLAGS= /optimize/standard=vaxc
LFLAGS = 
       
 
rasmol : rasmol.obj, molecule.obj, infile.obj, transfor.obj, command.obj, -
         abstree.obj, render.obj, repres.obj, x11win.obj, pixutils.obj, -
         outfile.obj, script.obj
         $(LINK) /exec=rasmol $(LFLAGS) rasmol.obj, molecule.obj, -
                infile.obj, transfor.obj, command.obj, abstree.obj, -
                render.obj, repres.obj, x11win.obj, pixutils.obj, -
                outfile.obj, script.obj, rasmol/opt

rasmol.obj :   rasmol.c, rasmol.h, molecule.h, transfor.h, command.h, -
               abstree.h, render.h, graphics.h, pixutils.h, outfile.h
               $(CC) $(CFLAGS) rasmol.c

molecule.obj : molecule.c, molecule.h, rasmol.h, abstree.h, transfor.h, -
               render.h
               $(CC) $(CFLAGS) molecule.c

infile.obj:    infile.c, infile.h, rasmol.h, molecule.h
               $(CC) $(CFLAGS) infile.c

transfor.obj : transfor.c, transfor.h, rasmol.h, molecule.h, command.h, -
               render.h, graphics.h
               $(CC) $(CFLAGS) transfor.c

command.obj :  command.c, command.h, rasmol.h, tokens.h, abstree.h, -
               molecule.h, transfor.h, render.h, graphics.h, pixutils.h, -
               outfile.h, script.h
               $(CC) $(CFLAGS) command.c

abstree.obj :  abstree.c, abstree.h, rasmol.h, molecule.h
               $(CC) $(CFLAGS) abstree.c

render.obj :   render.c, render.h, rasmol.h, molecule.h, transfor.h, -
               command.h, graphics.h, pixutils.h
               $(CC) $(CFLAGS) render.c

repres.obj:    repres.c, repres.h, rasmol.h
               $(CC) $(CFLAGS) repres.c

x11win.obj :   x11win.c, graphics.h, rasmol.h, bitmaps.h
               $(CC) $(CFLAGS) x11win.c

pixutils.obj : pixutils.c, pixutils.h, rasmol.h, render.h, graphics.h
               $(CC) $(CFLAGS) pixutils.c

outfile.obj :  outfile.c, outfile.h, rasmol.h, molecule.h, command.h, -
               transfor.h, render.h, graphics.h, pixutils.h, script.h
               $(CC) $(CFLAGS) outfile.c

script.obj:    script.c, script.h, rasmol.h, molecule.h, command.h, -
               abstree.h, transfor.h, render.h, graphics.h, pixutils.h
               $(CC) $(CFLAGS) script.c
