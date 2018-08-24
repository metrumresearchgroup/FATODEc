all : sho

.PHONY : clean

CC = clang

LIB_FATODE_DIR := ../..
sho: LDFLAGS += -L/usr/local/Cellar/gcc/7.3.0_1/lib/gcc/7 -lquadmath -lgfortran
sho: LDFLAGS += -L/usr/local/Cellar//gcc/7.3.0_1/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.3.0 -lgcc
sho: LDFLAGS += $(LIB_FATODE_DIR)/libfatode.a
sho: LDFLAGS += -framework accelerate # for macOS to link BLAS/LAPACK
sho: sho.o $(LIB_FATODE_DIR)/libfatode.a
	$(LINK.c) $^ -o $@

sho.o: CFLAGS = -I../../include

clean:
	rm *.o sho
