#Configuration flags for linear solvers. Must choose one from the following three options.
#  -DFULL_ALGEBRA            use BLAS and LAPACK library (full algebra)
#  -DSPARSE_UMF              use UMFpack (sparse format) 
#  -DSPARSE_LU               use SuperLU (sparse format)
#  
FATODE_LS_CONFIG  =  -DFULL_ALGEBRA

ARFLAGS = -rs
LIB_FATODE_DIR ?= .

FC = gfortran
CXX = clang++

#~~~> intel fortran
#export FC=ifort
#export FFLAGS= -cpp -O2 -nogen-interface -fltconsitency $(FATODE_LS_CONFIG) -warn all
#BLAS=    /opt/ifort_libs/libblas.a
#LAPACK =  /opt/ifort_libs/liblapack.a

FFLAGS = -cpp -O2 -ffree-line-length-none $(FATODE_LS_CONFIG)
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/ERK
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/RK
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/ROS
FFLAGS += -I$(LIB_FATODE_DIR)/src/FWD/SDIRK
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/ERK_TLM
# FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/RK_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/ROS_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/TLM/SDIRK_TLM
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/ERK_ADJ
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/RK_ADJ
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/ROS_ADJ
FFLAGS += -I$(LIB_FATODE_DIR)/src/ADJ/SDIRK_ADJ

debug: FFLAGS += -g -DDEBUG -O0
debug: $(LIB_FATODE_DIR)/libfatode.a

BLAS=-framework Accelerate
LAPACK=
# BLAS=/opt/lapack/lib/libblas.a
# LAPACK=/opt/lapack/lib/liblapack.a

all: $(LIB_FATODE_DIR)/libfatode.a # $(LIB_FATODE_DIR)/libfatode_cc.a

$(LIB_FATODE_DIR)/src/%.o : $(LIB_FATODE_DIR)/src/%.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -J$(dir $@) -c $< -o $@

CXXFLAGS = -I$(LIB_FATODE_DIR)/include -std=c++14
$(LIB_FATODE_DIR)/src/fatode_cc.o : $(LIB_FATODE_DIR)/src/fatode_cc.cpp
	$(COMPILE.cc) -o $@ $^

FAT_FWD_OBJS := $(patsubst %.F90, %.o, $(shell find $(LIB_FATODE_DIR)/src/FWD -name *.F90))
FAT_TLM_OBJS := $(patsubst %.F90, %.o, $(shell find $(LIB_FATODE_DIR)/src/TLM -name *.F90))
FAT_ADJ_OBJS := $(patsubst %.F90, %.o, $(shell find $(LIB_FATODE_DIR)/src/ADJ -name *.F90))
FAT_WRAPPER  := $(LIB_FATODE_DIR)/src/integrate_fatode.o
$(FAT_WRAPPER) : $(FAT_FWD_OBJS) $(FAT_TLM_OBJS)

FAT_OBJS := $(FAT_WRAPPER) $(FAT_FWD_OBJS) $(FAT_TLM_OBJS) # $(FAT_ADJ_OBJS)
FAT_OBJS += $(LIB_FATODE_DIR)/src/fatode_cc.o
$(LIB_FATODE_DIR)/libfatode.a : $(FAT_OBJS)
	$(AR) $(ARFLAGS) $@ $^

LDFLAGS += $(LIB_FATODE_DIR)/libfatode.a

LIB = -lm

FAT_LUPACK = 
FAT_LUWRAP = 
FAT_UMFWRAP =
FAT_LUWRAP =

# default: driver clean

clean:
	find src -name "*.o" -exec rm {} \;
	find src -name "*.mod" -exec rm {} \;
	@rm $(LIB_FATODE_DIR)/*.a
