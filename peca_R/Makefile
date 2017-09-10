all: bin/peca

objs = Pre Mcmc main data Option print Module Est
objs := $(addsuffix .o, $(objs))
objs := $(addprefix src/, $(objs))
CXXFLAGS = -Wall -O3 -I. -std=c++0x -Inlopt/include -I../include_dir


bin/peca: nlopt $(objs)
	$(CXX) $(objs) -o $@ -Lnlopt/lib -lnlopt 

nlopt:
	cd nlopt-2.4.2/;\
	    ./configure --prefix=$(CURDIR)/nlopt;\
	    $(MAKE) install

.PHONY : clean
clean :
	$(RM) bin/peca $(objs)
