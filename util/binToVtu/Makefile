# Default values if they are not set
VTKINCLUDE ?= /usr/include/vtk-6.2
VTKLIBPATH ?=
MYCXX ?= g++
MYCFLAGS ?= -O2 -fopenmp -DOMP
EXENAME ?= ConvertBinToVtu

LIB = -lm -lvtkIOMPIParallel-6.2 -lvtkIOXML-6.2 -lvtkFiltersSources-6.2 -lvtkFiltersCore-6.2 \
		-lvtkCommonExecutionModel-6.2 -lvtkCommonDataModel-6.2 -lvtkCommonCore-6.2

CXX= $(MYCXX)
CFLAGS += $(MYCFLAGS)
CPPFLAGS += -I$(VTKINCLUDE)
LDFLAGS += $(VTKLDFLAGS) $(LIB)
SRC = ./src/bin_to_vtu.cxx

all: $(EXENAME)

$(EXENAME): $(SRC)
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(SRC) -o $(EXENAME) $(LDFLAGS)

clean:
	rm -f $(EXENAME)
