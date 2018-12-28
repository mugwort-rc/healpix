# -*- Makefile -*-

# VPATH	= $(LIBDIR)
srcf90	= \
	healpix_types.f90 \
	extension.F90 \
	misc_utils.f90 \
	utilities.f90 \
	healpix_fft.F90 \
	num_rec.f90 \
	indmed.f90 \
	statistics.90 \
	head_fits.F90 \
	bit_manipulation.f90\
	pix_tools.f90 \
	fitstools.f90 \
	rngmod.f90 \
	alm_tools.f90 \
	udgrade_nr.f90 \
	paramfile_io.f90 \
	obsolete.f90 \
	coord_v_convert.f90 \
	ran_tools_dist.f90 \
	healpix_modules.f90

incalm  = alm_map_ss_inc.f90 alm_map_dd_inc.f90

incpix = convert_inplace_1d_inc.f90 \
	convert_inplace_nd_inc.f90 \
	convert_nest2ring_nd_inc.f90 \
	convert_ring2nest_nd_inc.f90 \
	remove_dipole_inc.f90

incfits = fits_s_inc.f90 fits_d_inc.f90

incudgrade = udgrade_s_inc.f90 udgrade_d_inc.f90

srcf	=
srcc	= cgetEnvironment.c

obj	= healpix_types.o extension.o misc_utils.o utilities.o healpix_fft.o \
	num_rec.o indmed.o statistics.o \
	head_fits.o \
	bit_manipulation.o pix_tools.o cgetEnvironment.o \
	fitstools.o rngmod.o alm_tools.o \
	udgrade_nr.o paramfile_io.o obsolete.o coord_v_convert.o ran_tools_dist.o healpix_modules.o

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.F90
	case $(OS) in \
	Darwin ) \
	ln -s -f $< $*.junk.F90; \
	$(FC) $(FFLAGS) $(ADDUS) $(PPFLAGS) -c $*.junk.F90 -o $@; \
	rm $*.junk.F90;; \
	OSF ) \
	ln -s -f $< $*.f90; \
	$(FC) $(FFLAGS) $(ADDUS) -cpp -c $(basename $<).f90 -o $@;; \
	* ) \
	$(FC) $(FFLAGS) $(ADDUS) $(PPFLAGS) -c $< -o $@;; \
	esac


default: libhealpix.a healpix_fft.o

libhealpix.a: $(obj) $(incpix) $(incalm) $(incfits) $(incudgrade)
	$(AR) libhealpix.a $(obj)
	@cp -p libhealpix.a $(LIBDIR)/
	@for p in *.$(MOD); do \
	cp $$p $(INCDIR)/; \
	done

alm_map_ss_inc.f90: alm_map_template.f90 gen_alm_code
	./gen_alm_code

alm_map_dd_inc.f90: alm_map_template.f90 gen_alm_code
	./gen_alm_code

fits_s_inc.f90: fits_template.f90 gen_fits_code
	./gen_fits_code

fits_d_inc.f90: fits_template.f90 gen_fits_code
	./gen_fits_code

udgrade_s_inc.f90: udgrade_template.f90 gen_udgrade_code
	./gen_udgrade_code

udgrade_d_inc.f90: udgrade_template.f90 gen_udgrade_code
	./gen_udgrade_code

pix_tools.o:pix_tools.f90  $(incpix)
	$(FC) $(FFLAGS) -c -o $@ $<

alm_tools.o:alm_tools.f90  $(incalm)
	$(FC) $(FFLAGS) -c -o $@ $<

fitstools.o: fitstools.f90 $(incfits)
	$(FC) $(FFLAGS) -c -o $@ $<

udgrade_nr.o: udgrade_nr.f90 $(incudgrade)
	$(FC) $(FFLAGS) -c -o $@ $<
