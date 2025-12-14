#--------RTsynth makefile---------
#--------Thomas G. Bisbas-------

FC = gfortran
HDF5_DIR = #Here your hdf5 directory
HDF5_INCLUDE = $(HDF5_DIR)/include
HDF5_LIB = $(HDF5_DIR)/lib
LIBS     = -lhdf5_fortran -lhdf5
FFLAGS = -g -fbacktrace -fcheck=all
HDF5_CFLAGS =
HDF5_LDLIBS =

HDF5 = 0
VELOCITY = 1
OPENMP = 1
DUST = 1
BACKGROUND = 1
CRATTENUATION = 0
HI21CM = 0
MULTILINES = 0
MULTIDIRECTION = 0
CDONLY = 0

ifeq ($(HDF5),1)
   OPT += -DHDF5
   HDF5_CFLAGS = -I$(HDF5_INCLUDE)
   HDF5_LDLIBS = -L$(HDF5_LIB) $(LIBS)
endif
ifeq ($(OPENMP),1)
   OPT += -fopenmp -DOPENMP
endif
ifeq ($(VELOCITY),1)
   OPT += -DVELOCITY
endif
ifeq ($(DUST),1)
   OPT += -DDUST
endif
ifeq ($(BACKGROUND),1)
   OPT += -DBACKGROUND
endif
ifeq ($(CRATTENUATION),1)
   OPT += -DCRATTENUATION
endif
ifeq ($(HI21CM),1)
   OPT += -DHI21CM
endif
ifeq ($(MULTILINES),1)
   OPT += -DMULTILINES
endif
ifeq ($(MULTIDIRECTION),1)
   OPT += -DMULTIDIRECTION
endif
ifeq ($(CDONLY),1)
   OPT += -DCDONLY
endif


SRC = modules.F90 read_command_line.F90 readcoolants.F90 read_input.F90 \
      initializations.F90 writeoutputs.F90 excitation.F90 \
      opticaldepth.F90 columndensities.F90 radiativetransfer.F90 

ifeq ($(HDF5),1)
   SRC += m_readh5.F90 readparams.F90
else
   SRC += readparams.F90 readpdrfiles_los.F90 readspopfile_los.F90
endif

ifeq ($(HI21CM),1)
  SRC += hi21cm.F90
endif

SRC += RTsynth.F90

EXEC = RTsynth

all: $(EXEC)
	mv $(EXEC) ../

$(EXEC): $(SRC)
	$(FC) $(FFLAGS) $(OPT) $(HDF5_CFLAGS) -o $(EXEC) $(SRC) $(HDF5_LDLIBS)

clean:
	rm -f $(EXEC) *.o *.mod

