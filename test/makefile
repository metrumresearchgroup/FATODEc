CXX = clang++

GTEST ?= gtest_1.7.0
GTEST_MAIN = $(GTEST)/src/gtest_main.cc
GTEST_CXXFLAGS = -isystem $(GTEST)/include -isystem $(GTEST) -O$O -DGTEST_USE_OWN_TR1_TUPLE

LIB_FATODE_DIR := ..
CXXFLAGS = -I$(LIB_FATODE_DIR)/include -I$(GTEST)/include
LDFLAGS += -L/usr/local/Cellar/gcc/7.3.0_1/lib/gcc/7 -lquadmath -lgfortran
LDFLAGS += -L/usr/local/Cellar//gcc/7.3.0_1/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.3.0 -lgcc
LDFLAGS += $(LIB_FATODE_DIR)/libfatode.a
LDFLAGS += -framework accelerate # for macOS to link BLAS/LAPACK
LDFLAGS += $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(LIB_FATODE_DIR)/libfatode.a

$(GTEST)/src/gtest-all.o: CXXFLAGS = $(GTEST_CXXFLAGS)

% : CXXFLAGS += $(GTEST_CXXFLAGS) -std=c++14
% : %.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(LIB_FATODE_DIR)/libfatode.a
	$(LINK.cpp) -o $@ $^
