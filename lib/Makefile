# Makefile for Modules

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../UtilXlib $(MOD_FLAG).

MODULES = \
asy.o \
debye_module.o \
bz_2d_form.o \
bz_asy.o \
elastic_const.o \
eldos_module.o \
isoentropic.o \
latgen_2d.o \
gnuplot.o \
gnuplot_color.o \
lattices.o \
linear_solvers.o \
magnetic_point_group.o \
mp_asyn.o \
nye.o \
phdos_module.o \
ph_freq_module.o \
piezoelectric_tensor.o \
point_group.o \
quadratic_surfaces.o \
quartic_surfaces.o \
rotate.o \
space_groups.o \
strain.o \
voigt.o \
xrdp.o 

all : thermo_lib.a

thermo_lib.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L 

include make.depend
