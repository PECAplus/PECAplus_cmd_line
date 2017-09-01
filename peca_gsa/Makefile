all: bin/cps

objs = main data Option print Pre Post
objs := $(addsuffix .o, $(objs))
objs := $(addprefix src/, $(objs))
CXXFLAGS = -Wall -O3 -I. -std=c++0x -I..

bin/cps: $(objs)
	$(CXX) $(objs) -o $@

.PHONY : clean
clean :
	$(RM) bin/cps $(objs)
