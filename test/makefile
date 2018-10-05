CXX = clang++

GTEST ?= gtest_1.7.0
GTEST_MAIN = $(GTEST)/src/gtest_main.cc
GTEST_CXXFLAGS = -isystem $(GTEST)/include -isystem $(GTEST) -O$O -DGTEST_USE_OWN_TR1_TUPLE

LIB_FATODE_DIR := ..
CXXFLAGS = -I$(LIB_FATODE_DIR)/include -I$(GTEST)/include
LDFLAGS += -L/usr/local/Cellar/gcc/7.3.0_1/lib/gcc/7 -lquadmath -lgfortran
LDFLAGS += -L/usr/local/Cellar//gcc/7.3.0_1/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.3.0 -lgcc
LDFLAGS += -framework accelerate # for macOS to link BLAS/LAPACK

gtest : $(GTEST)/src/gtest-all.o
all : sho_test roberts_test species_test
sho_test : CXXFLAGS += $(GTEST_CXXFLAGS) -std=c++14
roberts_test : CXXFLAGS += $(GTEST_CXXFLAGS) -std=c++14
species_test : CXXFLAGS += $(GTEST_CXXFLAGS) -std=c++14
debug: CXXLAGS += -g -DDEBUG -O0
debug : sho_test roberts_test species_test

$(GTEST)/src/gtest-all.o: CXXFLAGS = $(GTEST_CXXFLAGS) -g

%.o : %.cpp $(LIB_FATODE_DIR)/include/fatodec.hpp
	$(COMPILE.cc) -o $@ $<
sho_test : sho_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(LIB_FATODE_DIR)/libfatode.a
	$(LINK.cpp) -o $@ $^
roberts_test : roberts_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(LIB_FATODE_DIR)/libfatode.a
	$(LINK.cpp) -o $@ $^
species_test : species_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(LIB_FATODE_DIR)/libfatode.a
	$(LINK.cpp) -o $@ $^
clean :
	@rm *.o
	@rm *_test
