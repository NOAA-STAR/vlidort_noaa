
#  5/5/22. R. Spurr.
#    --Makefile for testing the new experimental BRDF system

# Define some variables

TEST_PATH     = .
MOD_PATH = mod
OBJ_PATH = obj
MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

SHELL = /bin/bash
FC = gfortran
FFLAGS = -c -I$(MOD_PATH) -J$(MOD_PATH)
FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow

ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif

#  For testing individual kernel output (selected geometry)

SOURCES_BRDF_EXPTESTER_KERNELS_ONLY += \
   $(TEST_PATH)/vlidort_experimental_local_pars.f90 \
   $(TEST_PATH)/vbrdf_sup_inputs_def.f90 \
   $(TEST_PATH)/vbrdf_sup_outputs_def.f90\
   $(TEST_PATH)/vbrdf_sup_aux.f90	  \
   $(TEST_PATH)/vbrdf_sup_kernels_old_reduced.f90	  \
   $(TEST_PATH)/vbrdf_sup_kernels_experimental.f90 \
   $(TEST_PATH)/V2p8p5_VBRDF_exptester_Konly.f90

#  For testing full BRDF output (DBbounce and Fouriers, 36 geometries)

SOURCES_BRDF_EXPTESTER_FULL += \
   $(TEST_PATH)/vlidort_experimental_local_pars.f90 \
   $(TEST_PATH)/vbrdf_sup_inputs_def.f90 \
   $(TEST_PATH)/vbrdf_sup_core_def.f90\
   $(TEST_PATH)/vbrdf_sup_outputs_def.f90\
   $(TEST_PATH)/vbrdf_sup_aux.f90	  \
   $(TEST_PATH)/vbrdf_findpar.f90	  \
   $(TEST_PATH)/vbrdf_sup_kernels_experimental.f90 \
   $(TEST_PATH)/vbrdf_sup_core_routines.f90	  \
   $(TEST_PATH)/vbrdf_sup_core_masters.f90	  \
   $(TEST_PATH)/vbrdf_sup_makerfourier_routines.f90	  \
   $(TEST_PATH)/vbrdf_sup_masters.f90	  \
   $(TEST_PATH)/vbrdf_sup_mod.f90	  \
   $(TEST_PATH)/V2p8p5_VBRDF_exptester_Full.f90

$(OBJ_PATH)/%.o : $(TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files

F90SOURCES_BRDF_EXPTESTER_KERNELS_ONLY := $(notdir $(filter %.f90, $(SOURCES_BRDF_EXPTESTER_KERNELS_ONLY)))
F90OBJECTS_BRDF_EXPTESTER_KERNELS_ONLY := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BRDF_EXPTESTER_KERNELS_ONLY)))

F90SOURCES_BRDF_EXPTESTER_FULL := $(notdir $(filter %.f90, $(SOURCES_BRDF_EXPTESTER_FULL)))
F90OBJECTS_BRDF_EXPTESTER_FULL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BRDF_EXPTESTER_FULL)))

# Define desired targets

All: V2p8p5_VBRDF_exptester_Konly.exe V2p8p5_VBRDF_exptester_Full.exe 

V2p8p5_VBRDF_exptester_Konly.exe: $(F90OBJECTS_BRDF_EXPTESTER_KERNELS_ONLY)
	$(FC) $^ -o $@
V2p8p5_VBRDF_exptester_Full.exe: $(F90OBJECTS_BRDF_EXPTESTER_FULL)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

