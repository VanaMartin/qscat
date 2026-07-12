.PHONY: all clean exec coupled

# implementations

IMPLEMENTATION ?= intel
LONGLONG_SUPPORT ?= yes
OPB_ROOT=/opt/OpenBLAS
CONFIGURATION ?= debug

# * * * Output directory structure * * * * * * * * * * * *

CCDIR     = linux_dep/coulcc
OBJDIR    = build/obj
OUTDIR    = build
SRCDIR    = source
TESTDIR   = build/tests

MODELS = N2 NO F2
SMAT = TW TWD TF
APPROX = LCP NRM
APPDIRS = bin CS WF

SMATDIR = $(foreach mod, $(MODELS), $(foreach mat, $(SMAT),  output/$(mod)/$(mat)))
WFDIR = $(foreach mod, $(MODELS), output/$(mod)/WF)
APPDIR = $(foreach mod, $(MODELS), $(foreach app, $(APPROX), $(foreach dirc, $(APPDIRS), output/$(mod)/$(app)/$(dirc))))
TSTDIR = $(foreach dr, T1 T2, $(foreach mat, $(SMAT),  output/tests/$(dr)/$(mat)))
DIRSTRUCT = $(SMATDIR) $(APPDIR) $(WFDIR) output/STATS $(TSTDIR)

# Subdirectories of src
SOURCES = $(SRCDIR) $(SRCDIR)/Arrays $(SRCDIR)/FemDvrEcs $(SRCDIR)/FemDvrEcs2d $(SRCDIR)/Model2d $(SRCDIR)/ModelLCP $(SRCDIR)/python
INCLUDE = -Iinclude/ -Ilinux_dep/ -Isource/ -IlibXcuda/libXcuda/include
# Flatten subdirectories
vpath %.cpp $(SOURCES)

# * * * * * * * Compiler options * * * * * * * * * * * * *

ifeq ($(CONFIGURATION),release)
    CONFFLAG = -DNDEBUG -O3 -g -vec-report
else
    CONFFLAG = -DDEBUG -O0 -g -debug all
endif
# BLAS
ifeq ($(IMPLEMENTATION),intel)
    CC   = icc
    CPP  = icpc
    FC   = ifort
    LINK = icpc
    PCP  = icpc
    BLAS = intel
    WRDEP = -wd981 -wd1418 -wd383 -wd1572
    FFLAGS = -i8 -openmp -I$(MKLROOT)/include -traceback -fPIC -g
    CFLIBS = -L/opt/intel/fc/9.0/lib/ -lifcore
    #ifeq($(LONGLONG_SUPPORT),yes)
        CFLAGS = $(CONFFLAG) -Wall $(WRDEP) $(INCLUDE) -DMKL_ILP64 -openmp -I$(MKLROOT)/include -traceback -g -DINTEL_MKL -fPIC
        MASTER_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm $(CFLIBS) -liomp5
    #else
    #    CFLAGS = $(CONFFLAG) -Wall $(WRDEP) $(INCLUDE) -DMKL_LP64 -openmp -I$(MKLROOT)/include -traceback -g -DINTEL_MKL -fPIC
    #    MASTER_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm $(CFLIBS) -liomp5
    #endif
    MLIBS = -lpthread -liomp5
else
    CC   = gcc
    CPP  = g++
    FC   = gfortran
    LINK = g++
    PCP  = g++
    BLAS = openblas
    CFLAGS = $(CONFFLAG) -Wall -Wno-unknown-pragmas -Wno-reorder $(INCLUDE) -openmp -I$(OPB_ROOT)/include -g -D OPENBLAS_LIB_ -fPIC
    FFLAGS = -openmp -I$(OPB_ROOT)/include -fPIC -g -lgomp
    CFLIBS = -L$(OPB_ROOT)/lib -lgfortran
    MASTER_LIBS += -L$(OPB_ROOT)/lib -lopenblas -lm -lgfortran -lgomp
endif

########################################
# * * * * * * * Objects  * * * * * * * #
########################################
# * * * * common * * * *
COMOBJ = common bessel coulomb Storage Object
ARROBJ = Buffer $(BLAS) Vector Matrix RowCompressedMatrix Arrays
FEMOBJ = FemDvrFunctions DvrGrid FemDvrEcsGrid GridVector Projector KineticEnergy OperatorDiagonal OperatorFull OperatorRowCompressed DiscreteStates FemDvrEcs # poly_grid_vector
F2DOBJ = FemDvrEcsGrid2d GridVector2d OperatorRowCompressed2d OperatorFull2d EquidistantProjector2d FemDvrEcs2d # DoubleGridVector2d DoubleOperator2d MaskGrid2d MaskVector2d MaskOperator2d
# * * * * models * * * *
M2DOBJ = input Potentials2d TestFunction2d DiracTestFunction2d FluxTestFunction2d MultiTestFunction2d TimeDependentModel2d CoupledModel2d #compressed_testfunction_2d coupled_model_2d
LCPOBJ = SMatrix ModelLCP
NRMOBJ = # module_NRM
# * * * * rules  * * * *
OBJ  = $(COMOBJ) $(ARROBJ) $(FEMOBJ) $(F2DOBJ) $(LCPOBJ) $(M2DOBJ) $(NRMOBJ)


BUILDOBJ = $(foreach obj, $(OBJ), $(OBJDIR)/$(obj).o)
#DPSOBJ = $(foreach obj, $(OBJ), $(OBJDIR)/$(obj).deps)
SRC    = $(foreach file, $(OBJ), $(file).cpp)
SRCOBJ = $(foreach obj, $(BUILDOBJ), $(OBJDIR)/$(obj))

# * * * * * * * Python Wrapper definitions * * * * * * * *
PYLIBDIR  = build/python
BOOSTROOT = /usr/local/src/boost_1_55_0
PYTHON_LIB=$(shell python3-config --libs)
PYTHON_INC=$(shell python3-config --includes) -Iinclude
PYFLAGS = -shared $(PYTHON_INC) -I$(BOOST_ROOT) -fPIC -I$(BOOSTROOT) -openmp
PYLIBS  = -L$(BOOSTROOT)/libs -L$(BOOSTROOT)/stage/lib -lboost_python3 $(PYTHON_LIB) $(MASTER_LIBS) -lcudart -lcublas -lcusparse libXcuda/build/libXcuda.so
PYFILES = pyArrays
PYOBJ   = $(foreach obj, $(PYFILES), $(PYLIBDIR)/$(obj).o)
PYSOBJ  = $(foreach obj, $(PYFILES), $(PYLIBDIR)/$(obj).so)
PYPREFIX ?= ~/py3
# * * * * * * * Tests  * * * * * * * * * * * * * * * * * *

TESTSRCDIR = $(SRCDIR)/tests
SRCT = tests/t_arrays.cpp tests/t_fem_dvr_ecs.cpp tests/t_fem_dvr_ecs_2D.cpp
TESTS   = $(TESTDIR)/t_arrays $(TESTDIR)/t_fem_dvr_ecs $(TESTDIR)/t_fem_dvr_ecs_2D
TESTOBJ = $(foreach obj, $(TESTS), $(obj).o)

LIBOBJ = $(OBJDIR)/coulcc.o  $(OBJDIR)/sec_stream.o

# include dependencies

#### BUILDS ######################
EXEC=$(OUTDIR)/eMoScat
COUPLED=$(OUTDIR)/simpleCoupling
COUP1D=$(OUTDIR)/coupled_1D
SCAT3D=$(OUTDIR)/scat3d

all: $(OUTDIR) $(EXEC) $(DIRSTRUCT)

install: pylib_install

pylib_install: pylib
	mkdir -p $(PYPREFIX)
	cp -r python/qscat $(PYPREFIX)
	cp $(PYSOBJ) $(PYPREFIX)/qscat/pyArrays.so

tests: $(TESTDIR) $(LIBOBJ) $(BUILDOBJ) $(TESTOBJ) $(TESTS)
	./build/tests/t_arrays

doc:
	doxygen documentation/config

pylib: $(PYLIBDIR) $(PYSOBJ)

coupled: $(OUTDIR) $(COUPLED)

1dcp: $(OUTDIR) $(COUP1D)

scat3d: $(OUTDIR) $(SCAT3D)
	echo scattering

time_independent: $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/time_independent_model.o
	$(CPP) $(CFLAGS) $^ $(MASTER_LIBS) -o $@

clean:
	rm -rf $(EXEC) $(OBJDIR)/eMoScat.o $(COUPLED) $(OBJDIR)/simpleCoupling.o $(SCAT3D) $(OBJDIR)/scat3d.o $(BUILDOBJ) $(LIBOBJ)
	rm -rf $(TESTS) $(PYOBJ) $(PYSOBJ)
	rm -rf $(subst .o,.deps,$(BUILDOBJ))
	rm -f  $(OUTDIR)/*_test

# deletes all libraries
purge:
	(cd libXcuda && make clean)

output/%:
	mkdir -p $@

$(OUTDIR):
	mkdir $@

$(OBJDIR):
	mkdir $@

$(PYLIBDIR):
	mkdir $@

$(TESTDIR):
	mkdir -p $@

lXcuda:
	(cd libXcuda && make)
	cp libXcuda/build/libXcuda.so build/libXcuda.so

$(EXEC): $(OUTDIR) $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/eMoScat.o
	$(CPP) $(CFLAGS) $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/eMoScat.o $(MASTER_LIBS) -o $(EXEC)

# TODO
$(COUPLED): $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/simpleCoupling.o
	$(CPP) $(CFLAGS) $^ $(MASTER_LIBS) -o $@

$(COUP1D): $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/coupled_1D.o
	$(CPP) $(CFLAGS) $^ $(MASTER_LIBS) -o $@

$(SCAT3D): $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/scat3d.o
	$(CPP) $(CFLAGS) $^ $(MASTER_LIBS) -o $@

$(PYLIBDIR)/%.so: $(LIBOBJ) $(BUILDOBJ) $(PYLIBDIR)/%.o
	$(PCP) $(PYFLAGS) $^ -o $@ $(PYLIBS)

$(TESTS): $(TESTOBJ)

$(TESTDIR)/%: $(TESTDIR)/%.o
	$(CPP) $(CFLAGS) $< $(BUILDOBJ) $(LIBOBJ) $(MASTER_LIBS) -o $@

### SPECIAL: GENERAL TESTS
$(OUTDIR)/%_test: $(OBJDIR)/%_test.o $(BUILDOBJ) $(LIBOBJ)
	$(CPP) $(CFLAGS) $^ $(MASTER_LIBS) -o $@

#$(OUTDIR)/%_test
.PRECIOUS: $(OUTDIR)/%_test $(OBJDIR)/%_test.o
%_test: $(OUTDIR)/%_test
	@echo

# dependencies - no need to make target BUILDOBJ
-include $(subst .o,.deps,$(BUILDOBJ))
-include $(subst .o,.deps,$(TESTOBJ))
-include $(OBJDIR)/simpleCoupling.deps $(OBJDIR)/coupled_1d.deps

-include $(OBJDIR)/test.deps
#cocg: $(OUTDIR) $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/test.o
#	$(CPP) $(CFLAGS) $(LIBOBJ) $(BUILDOBJ) $(OBJDIR)/test.o $(MASTER_LIBS) -o $(OUTDIR)/test

########## Rules

$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

$(OBJDIR)/coulcc.o: linux_dep/coulcc/coulcc.f
	$(FC) $(FFLAGS) -c linux_dep/coulcc/coulcc.f -o $@

$(OBJDIR)/sec_stream.o: linux_dep/sec_stream.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

$(PYLIBDIR)/%.o: %.cpp
	$(PCP) -I$(BOOSTROOT) $(PYTHON_INC) $(CFLAGS) -fPIC -c $< -o $@

# makedepend does not work with vpath so we maintain dependencies manually:
#$(OBJDIR)/%.deps : %.c
#	@echo "make c deps: $@"
#	@echo -n "$@ $(OBJDIR)/" > $@+ && $(CC) -MM $(CFLAGS) $< 2>/dev/null >> $@+ && mv $@+ $@

$(OBJDIR)/%.deps : %.cpp
	@echo "make cpp deps: $@"
	@echo -n "$@ $(OBJDIR)/" > $@+ && $(CPP) -MM $(CFLAGS) $< 2>/dev/null >> $@+ && mv $@+ $@

$(TESTDIR)/%.o: $(TESTSRCDIR)/%.cpp
	$(CPP)  $(CFLAGS) -c $< -o $@

$(TESTDIR)/%.deps : $(TESTSRCDIR)/%.cpp
	@echo "make cpp deps: $@"
	@echo -n "$@ $(OBJDIR)/" > $@+ && $(CPP) -MM $(CFLAGS) $< 2>/dev/null >> $@+ && mv $@+ $@
