all: peca_core_bin peca_r_bin peca_gsa_bin

peca_core_bin: nlopt 
	cd peca_core;\
	    $(MAKE) 

peca_r_bin: nlopt 
	cd peca_r;\
	    $(MAKE) 

peca_gsa_bin:
	cd peca_gsa;\
	    $(MAKE) 

nlopt:
	cd $(CURDIR)/nlopt_dir/nlopt-2.4.2/;\
	    ./configure --prefix=$(CURDIR)/nlopt_dir/nlopt;\
	    $(MAKE) install

