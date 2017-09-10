all: peca_r_bin

objs = Pre Mcmc main data Option print Module Est
objs := $(addsuffix .o, $(objs))
objs := $(addprefix src/, $(objs))
CXXFLAGS = -Wall -O3 -I. -std=c++0x -I../nlopt_dir/nlopt/include -I../include_dir


peca_r_bin: nlopt $(objs)
	$(CXX) $(objs) -o $@ -L../nlopt_dir/nlopt/lib -lnlopt 

nlopt:
	cd ../nlopt_dir/nlopt-2.4.2/;\
	    ./configure --prefix=$(CURDIR)/../nlopt_dir/nlopt;\
	    $(MAKE) install

.PHONY : clean
clean :
	$(RM) peca_r_bin $(objs)


