NETCDF = /home/decker/netcdf/netcdf-4.2

FC = gfortran
#FC = pgf90
FFLAGS = -std=f95 -pedantic -Wall -Wextra -march=native -ffast-math -funroll-loops -O3 -I/${NETCDF}/include
#FFLAGS = -fast -Mscalarsse -Mvect=sse -Mflushz -Mstandard -g77libs
#FFLAGS = -g -Mstandard -g77libs
LIBS = ${NETCDF}/lib/libnetcdff.a ${NETCDF}/lib/libnetcdf.a ${GEMLIB}/gemlib.a ${GEMLIB}/textlib.a ${GEMLIB}/libxml2.a ${GEMLIB}/libz.a ${GEMLIB}/libxslt.a ${GEMLIB}/cgemlib.a

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $<

FOBJ = dateutil.o diagnostics.o gempak.o registry.o wrf2gem.o wrf2gemsubs.o 

wrf2gem: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(FOBJ) $(LIBS)

dateutil.o: dateutil.f90
diagnostics.o: diagnostics.f90 dateutil.o
gempak.o: gempak.f90 wrf2gemsubs.o 
registry.o: registry.f90 diagnostics.o gempak.o wrf2gemsubs.o 
wrf2gem.o: wrf2gem.f90 registry.o gempak.o wrf2gemsubs.o 
wrf2gemsubs.o: wrf2gemsubs.f90 

clean:
	rm -f *.o *.mod
