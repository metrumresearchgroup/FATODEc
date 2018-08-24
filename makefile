#Configuration flags for linear solvers. Must choose one from the following three options.
#  -DFULL_ALGEBRA            use BLAS and LAPACK library (full algebra)
#  -DSPARSE_UMF              use UMFpack (sparse format) 
#  -DSPARSE_LU               use SuperLU (sparse format)
#  
FATODE_LS_CONFIG  =  -DFULL_ALGEBRA

ARFLAGS = -rs
LIB_FATODE_DIR ?= .

#~~~> intel fortran
#export FC=ifort
#export FFLAGS= -cpp -O2 -nogen-interface -fltconsitency $(FATODE_LS_CONFIG) -warn all
#BLAS=    /opt/ifort_libs/libblas.a
#LAPACK =  /opt/ifort_libs/liblapack.a

#~~~> gfortran (GNU FORTRAN Compiler)
FC = gfortran
FFLAGS = -cpp -O2 -ffree-line-length-none $(FATODE_LS_CONFIG)
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/ERK
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/RK
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/ROS
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/SDIRK
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/ERK_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/RK_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/ROS_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/SDIRK_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/ERK_ADJ
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/RK_ADJ
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/ROS_ADJ
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/SDIRK_ADJ

BLAS=-framework Accelerate
LAPACK=
# BLAS=/opt/lapack/lib/libblas.a
# LAPACK=/opt/lapack/lib/liblapack.a

# FATODE_TESTS := $(subst .cpp,$(EXE),$(shell find test/unit/math/torsten -name *fatode_test.cpp))
# $(FATODE_TESTS) : $(LIB_FATODE)

LIBDIR    = $(FATDIR)/LSS_LIBS/x86_64_Linux
MODEL     = ADJ
FAMILY    = RK_ADJ

INTEGRATOR= $(FAMILY)_f90_Integrator.o

all: $(LIB_FATODE_DIR)/libfatode.a

$(LIB_FATODE_DIR)/src/%.o : $(LIB_FATODE_DIR)/src/%.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -J$(dir $@) -c $< -o $@

FAT_OBJS := $(patsubst %.F90, %.o, $(shell find $(LIB_FATODE_DIR)/src -name *.F90))


$(LIB_FATODE_DIR)/libfatode.a : $(FAT_OBJS)
	$(AR) $(ARFLAGS) $@ $^

LDFLAGS += $(LIB_FATODE_DIR)/libfatode.a

LIB = -lm

FAT_LUPACK = 
FAT_LUWRAP = 
FAT_UMFWRAP =
FAT_LUWRAP =

# default: driver clean

# driver: cbm4_rk_adj_dr.o $(PAR) $(LSSOLVER) $(INTEGRATOR)
# 	$(FC) $(FFLAGS) -o $(APP) $< $(FAT_LUWRAP) $(FAT_LUPACK) $(FAT_UMFWRAP) $(FAT_UMFPACK) $(LSSOLVER) $(INTEGRATOR) $(PAR) $(LAPACK) $(BLAS) $(LIB)


clean:
	@rm $(LIB_FATODE_DIR)/*.a
