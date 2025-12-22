CXX = g++
CXXFLAGS = -g -Wall -fPIC -Wno-unused-variable
ROOTFLAGS = `root-config --cflags --glibs --libs` -lTreePlayer -lEG -lMinuit -lMathMore

INCLUDE_SRCS = $(wildcard includes/*.cpp)  # all object files in includes/
INCLUDE_OBJS = $(INCLUDE_SRCS:.cpp=.o)

########################################
## main programs to compile
########################################
TARGETS = NC_analyzer univmake        # MODIFY


########################################
########################################

# includes/*.o
includes/%.o: includes/%.cpp
    $(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

# ********************************
# # main program build --> MODIFY
# ********************************           
NC_analyzer: NC_analyzer.cpp $(INCLUDE_OBJS)
    $(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -o $@ $^
# ********************************

# univmake
univmake: univmake.C $(INCLUDE_OBJS)
    $(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -o $@ $^

# root dictionary
stv_root_dict.o:
    $(RM) stv_root_dict*.*
    rootcling -f stv_root_dict.cc -c LinkDef.h
    $(CXX) $(shell root-config --cflags --libs) -O3 -fPIC -o stv_root_dict.o -c stv_root_dict.cc
    $(RM) stv_root_dict.cc


########################################
########################################
# cleaning

clean:
    rm -f $(wildcard *.o) \
          $(INCLUDE_OBJS) \
          $(wildcard includes/*.so) \
          $(wildcard includes/*.pcm) \
          $(TARGETS) \
          stv_root_dict.*

.PHONY: all clean

