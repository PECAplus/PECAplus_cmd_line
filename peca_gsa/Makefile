all: gsa_bin

objs = main data Option print Pre Post
objs := $(addsuffix .o, $(objs))
objs := $(addprefix src/, $(objs))
CXXFLAGS = -Wall -O3 -I. -std=c++0x -I../include_dir

gsa_bin: $(objs)
	$(CXX) $(objs) -o $@

.PHONY : clean
clean :
	$(RM) gsa_bin $(objs)
