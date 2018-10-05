# ============================================================================
# Name        : Makefile
# Author      : Roy Walters
# Version     :
# Copyright   : Copyright (C) 2012 under GNU General Public License
# Description : Makefile for lanslide geometry versions
# ============================================================================

.PHONY: all clean

# Change these lines if you are using different
# Fortran or compiler options.
# Portland Group fortran
FC = pgfortran
F_FLAGS = -O2 -fast
DBG_FLAGS = -C -g -Ktrap=fp -traceback
OMP_FLAGS = -mp
# intel fortran
#FC = ifort
#F_FLAGS = -fast
#DBG_FLAGS = -fpe0 -check bounds -check uninit -traceback
#OMP_FLAGS = -qopenmp
# gfortran
#FC = gfortran
#F_FLAGS = -O2
#DBG_FLAGS = -g -O2 -fcheck=bounds -ffpe-trap='invalid','zero','overflow' -fbacktrace
#OMP_FLAGS = -mp
# gfortran - needed for netCDF for library compatibility

lg82:
	$(FC) $(F_FLAGS) -c \
		Landslide_geometry82.f90
	$(FC) *.o -o lg82
	rm *.o
	rm *.mod

lg82_dbg:
	$(FC) $(DBG_FLAGS) -c \
		Landslide_geometry82.f90
	$(FC) *.o -o lg82_dbg
	rm *.o
	rm *.mod

lg83:
	$(FC) $(F_FLAGS) -c \
		Landslide_geometry83.f90
	$(FC) *.o -o lg83
	rm *.o
	rm *.mod

lg83_dbg:
	$(FC) $(DBG_FLAGS) -c \
		Landslide_geometry83.f90
	$(FC) *.o -o lg83_dbg
	rm *.o
	rm *.mod

clean:
	rm -f *.o *.mod
