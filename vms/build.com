$! BUILD.COM
$! VMS Compilation Script
$! RasMol v2.6
$!
$ If F$Trnlnm ("X11") .eqs. "" Then Define X11 DECW$Include
$ If F$Trnlnm ("SYS") .eqs. "" Then Define Sys Sys$Share
$!
$ cc/optimize rasmol.c
$ cc/optimize molecule.c
$ cc/optimize infile.c
$ cc/optimize transfor.c
$ cc/optimize command.c
$ cc/optimize abstree.c
$ cc/optimize render.c
$ cc/optimize repres.c
$ cc/optimize x11win.c
$ cc/optimize pixutils.c
$ cc/optimize outfile.c
$ cc/optimize script.c
$!
$ link /exec=rasmol rasmol.obj, molecule.obj, infile.obj, transfor.obj, -
    command.obj, abstree.obj, render.obj, repres.obj, x11win.obj, -
    pixutils.obj, outfile.obj, script.obj, rasmol/opt
$!
$ set prot=w:re rasmol.exe
$ set prot=w:r rasmol.hlp
