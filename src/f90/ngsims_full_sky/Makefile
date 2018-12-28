# -*- Makefile -*-

VPATH	= $(BINDIR)
exec1	= $(BINDIR)/sky_ng_sim
src1	= sky_ng_sim
obj1 	= sky_sub.o  sub_ngpdf_powergauss.o  sub_ngpdf_sho.o  sky_ng_sim.o

# exec2	= 
# src2	= 
# obj2 	= 

%.o: %.f90
	$(FC) $(FFLAGS)  -c $< -o $@

%.o: %.F90
	case $(OS) in \
	Darwin ) \
	ln -s -f $< $*.junk.F90; \
	$(FC) $(FFLAGS)  $(ADDUS) $(PPFLAGS) $(PGFLAG) -c $*.junk.F90 -o $@; \
	rm $*.junk.F90;; \
	OSF ) \
	ln -s -f $< $*.f90; \
	$(FC) $(FFLAGS)  $(ADDUS) -cpp       $(PGFLAG) -c $(basename $<).f90 -o $@;; \
	* ) \
	$(FC) $(FFLAGS)  $(ADDUS) $(PPFLAGS) $(PGFLAG) -c $< -o $@;; \
	esac

default: $(exec1)

$(exec1): $(obj1) $(LIBDIR)/libhealpix.a
	$(FC)  $(FFLAGS) -o $(exec1) $(obj1) $(PGLIBS) $(LDFLAGS)
