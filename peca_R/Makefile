all: peca_r_bin

objs = Pre Mcmc main data Option print Module Est
objs := $(addsuffix .o, $(objs))
objs := $(addprefix src/, $(objs))
CXXFLAGS = -Wall -O3 -I. -std=c++0x -I../include_dir/nlopt/include -I../include_dir

$(objs): nlopt

peca_r_bin: nlopt $(objs)
	$(CXX) $(objs) -o $@ -L../include_dir/nlopt/lib -lnlopt 

nlopt:
	cd ../include_dir/nlopt-2.4.2/;\
	    ./configure --prefix=$(CURDIR)/../include_dir/nlopt;\
	    $(MAKE) install

.PHONY : clean
clean :
	$(RM) peca_r_bin $(objs)


