##############################################################################
#
# Makefile
#
# This is a makefile for use with a 'c' version of the healpix toolkit. 
#
# Usage:
#
# To build the library, type
#
# % make [CFITSIO_INCDIR=<cfitsio_incdir>] [CFITSIO_LIBDIR=<cfitsio_libdir>]
#
# where
#  CFITSIO_INCDIR is where the header files for the 'cfitsio' library is kept
#
# For example:
#  % make CFITSIO_INCDIR=/scr/kmg/include
#
# To install the files, type
#
# % make install [INCDIR=<incdir>] [LIBDIR=<libdir>] [RANLIB=<ranlib>] 
#
# For example:
#  % make install LIBDIR=/home/kmg/lib INCDIR=/home/kmg/include RANLIB="ar -ts"
#
# where
#  INCDIR is the directory in which you want to install the header files
#  LIBDIR is the directory in which you want to install the library files
#  RANLIB is the 'ranlib' command. It defaults to 'ranlib', but if your
#         system does not have or need this, you should use something like
#         RANLIB="ar -ts" 
#
# You could also build the package by directly editing the Makefile
#  (not recommended unless you're going to do it a lot and don't think
#   you'll be keeping up with changes in the distribution).
# To build the library:
#  1) Edit the 'CC', 'LIBDIR' and 'INCDIR' lines in 'Makefile' to reflect your
#     system.
#  2) Type 'make'
#
# To install the library and header file:
#  1) Type 'make install'
#
##############################################################################
#
# Compiler Options
CC = gcc -D__USE_FIXED_PROTOTYPES__ -Wall
#
# Where you want to install the library and header file
LIBDIR = /home/efh/Healpix/lib2
INCDIR = /home/efh/Healpix/include2
RANLIB = ranlib
#
# Where you have the cfitsio installation
CFITSIO_INCDIR = /usr/local/include
CFITSIO_LIBDIR = /usr/local/lib
#
# You should not have to edit below this line ################################
#
# The sources in the package
MODS = ang2pix_nest ang2pix_ring \
	pix2ang_nest pix2ang_ring nest2ring ring2nest \
	mk_pix2xy mk_xy2pix nside2npix  npix2nside \
	ang2vec vec2ang pix2vec_nest pix2vec_ring vec2pix_nest vec2pix_ring \
	printerror get_fits_size write_healpix_map read_healpix_map

OBJS = ${MODS:%=%.o}

#
# Make the library itself
libchealpix.a : ${OBJS}
	ar -r $@ ${OBJS}
#
# Install the library
install : libchealpix.a chealpix.h
	cp libchealpix.a ${LIBDIR}/.
	${RANLIB}  ${LIBDIR}/libchealpix.a
	cp chealpix.h ${INCDIR}/.
#
# Make the programs to test the package.
test_chealpix : test_chealpix.o libchealpix.a chealpix.h
	${CC} -o $@ test_chealpix.o  -L. -lchealpix -lm

test_chealpix2 : test_chealpix2.o libchealpix.a chealpix.h
	${CC} -o $@ test_chealpix2.o -L. -lchealpix -lm

test_cio : test_cio.o libchealpix.a chealpix.h
	${CC} -o $@ test_cio.o       -L. -lchealpix -lm -L${CFITSIO_LIBDIR} -lcfitsio

tests : test_chealpix test_chealpix2 test_cio

#
# General compilation rules
%.o : %.c
	${CC} -c -o $@ $< -I${CFITSIO_INCDIR}

# Tidy
clean :
	rm -f *~ *.o core libchealpix.a test_chealpix test_chealpix2 test_cio \
		test_chealpix.dat test_chealpix.out test_output.fits
#
# Modification History:
#
#    Date     Programmer Remarks
# ----------- ---------- -------
# 06-Jan-2000 K. Ganga   First version. 
# 16-Nov-2001 K. Ganga   Added documentation
# 27-Nov-2001 S. Du      Introduced RANLIB parameter so this works on 
#			systems without 'ranlib'
#	      K. Ganga   Added some docs. 
# 30-Jan-2003 E. Hivon   Edited for distribution

