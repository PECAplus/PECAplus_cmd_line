# PECAplus CLI

Command-Line Interface for PECA+

<img src="https://github.com/PECAplus/PECAplus_cmd_line/wiki/images/PECALogo.png" align="left" width="200" height="150">



**Protein Expression Control Analysis** (**PECA**) is a statistical toolbox to analyze time-series multi-omics dataset where molecules in one -omics platform serve as template for synthesis of the other.  For example, PECA can analyze paired RNA and protein time series data in order to identify change points of the impact translation and protein degradation on changes in protein concentration *given* changing RNA concentrations. 



The modules are presented mostly for the analysis of paired protein-RNA data. However, we note that PECA can be applied to any analogous dataset. For example, one can use PECA can analyze paired DNA and RNA concentration data (where the DNA concentration is typically set to 1) and provide a change point analysis for the impact of transcription and RNA degradation on changes in RNA concentration, given DNA.

The principal method was published [here](http://pubs.acs.org/doi/abs/10.1021/pr400855q) (Teo et al, *J. Prot. Res.* 2014). Additional information on the methods and the modules contained in PECAplus are described in a forthcoming manuscript. 

## Full Documentation

See the [Wiki](https://github.com/PECAplus/PECAplus_cmd_line/wiki) for full documentation.

See [Perseus-PluginPECA](https://github.com/PECAplus/Perseus-PluginPECA) for Perseus Plugin version.

## Installing 

### Linux and OSX (Mac)

Type `make` to install the software. The makefile automatically compiles the program. The executable will be in the `bin` folder. Binaries for 64-bit Windows are included. Due to limitations of GNU Make, the path of PECA is not allowed to contain any whitespace characters. The PECA software consists of the Markov chain Monte Carlo (MCMC) sampler and a Python script to summarize the results from the MCMC samples. The Python script requires Python 3.x and matplotlib.

### Windows

Executables for 64-bit Windows are included here (<- insert link()).

## Bugs and Feedback

For bugs, questions and discussions please use the [GitHub Issues](https://github.com/PECAplus/PECAplus_cmd_line/issues).

## License

Copyright (C) <2017> Guoshou Teo < gt49@nyu.edu >, Christine Vogel < cvogel@nyu.edu >, and Hyungwon Choi < hwchoi@nus.edu.sg >, National University of Singapore.

Licensed under the Apache License, Version 2.0 (the "License");

you may not use this file except in compliance with the License.

You may obtain a copy of the License at

[Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software

distributed under the License is distributed on an "AS IS" BASIS,

WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and

limitations under the License.


