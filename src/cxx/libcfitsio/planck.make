include $(PARAMFILE)

all: fitslib

fitslib: $(LIBDIR)/libcfitsio.a

$(LIBDIR)/libcfitsio.a: $(SRCROOT)/libcfitsio/cfitsio2510.tar.gz
	rm -rf cfitsio
	gunzip -c $(SRCROOT)/libcfitsio/cfitsio2510.tar.gz | tar xvf -
	cd cfitsio; \
	mv compress_alternate.c compress.c; \
	MAKE="$(MAKE)" FC="$(FC)" CC="$(CC)" CFLAGS="$(CCFLAGS_NO_C)" \
	  ./configure --includedir="$(INCDIR)" --libdir="$(LIBDIR)"; \
	$(MAKE)
	cp cfitsio/libcfitsio.a $(LIBDIR)
	cp cfitsio/*.h $(INCDIR)

clean:
	if test \( -d cfitsio \); then rm -rf cfitsio; fi
