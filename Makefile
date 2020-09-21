ifeq ($(P56xxLIBS),)
$(error P56xxLIBS environment variable is missing, source <path>/p56xxlib/setup.sh)
endif

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -I$(P56xxLIBS)/inc -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -L$(P56xxLIBS)/lib -lP56xx
GXX	   = g++ $(CXXFLAGS)

SRCS = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
EXES = $(SRCS:%.cpp=%)
dep = $(OBJ:.o=.d)  # one dependency file for each source

all: $(EXES)

$(EXES): $(SRCS)
	$(GXX) $(CXXFLAGS) -MMD -c $@.cpp
	$(GXX) $@.cpp -o$@ $(LDFLAGS)

clean:
	rm -f *.o *.so *.d *~ $(EXES)

cleanall: clean
	rm -f $(dep)
