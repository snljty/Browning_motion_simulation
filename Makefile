# Makefile for C source using fftw3

FILENAME = Browning_motion_simulation

FC = gfortran
FLINKER = $(FC)
INCPATH = C:\\fftw-3.3.5-dll64
LIBPATH = $(INCPATH)
LIB = fftw3-3

.PHONY: all $(FILENAME) clean clean_tmp

all: $(FILENAME)

$(FILENAME): $(FILENAME).exe

$(FILENAME).exe: $(FILENAME).o
	@echo Linking $@ ...
	$(FLINKER) -o $@ $^ -L $(LIBPATH) -l $(LIB)

%.o: %.f03
	@echo Compiling $@ ...
	$(FC) -o $@ -c $< -I $(INCPATH)

clean: clean_tmp
	-del $(FILENAME).exe 1> NUL 2> NUL

clean_tmp:
	-del $(FILENAME).o 1> NUL 2> NUL

