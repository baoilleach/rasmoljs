cmake_minimum_required(VERSION 2.8)


set(headers abstree.h cmndline.h command.h font.h graphics.h infile.h
       molecule.h outfile.h pixutils.h rasmol.h render.h repres.h
       scripts.h tmesh.h tokens.h transfor.h)
set(srcs abstree.c cad.c cmndline.c command.c infile.c molecule.c outfile.c
       pixutils.c render.c repres.c scripts.c tmesh.c transfor.c tokens.c)
add_library(raslib ${srcs} ${headers})

set(headersx11 bitmaps.h)
set(srcsx11 rasmol.c x11win.c)

add_definitions(-DTHIRTYTWOBIT)

if (CMAKE_COMPILER_IS_GNUCC)
  add_executable(rasmol ${srcsx11} ${headersx11})
  target_link_libraries(rasmol raslib m Xext X11 Xi)

  add_executable(rastxt rastxt.c)
  target_link_libraries(rastxt raslib m)

  find_package(SDL)
  add_executable(rassdl rassdl.c)
  target_link_libraries(rassdl raslib m ${SDL_LIBRARY})

else()
  add_executable(rassdl rassdl.c)
  target_link_libraries(rassdl raslib)

endif()
