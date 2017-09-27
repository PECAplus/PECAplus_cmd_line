

#include"main.hpp"
#include"Option.hpp"


Option::Option(const map<string,string>& opm) :
    d_opm(opm),d_FDRcut(-1),d_syndeg(true)
{
    istringstream iss(d_opm.at("GO_term_table"));
    iss>>d_GO_term;
    if (not ifstream(GO_term().c_str())) throw runtime_error("GO_term_table?");

    d_CPS="data_R_CPS.txt";

    d_module="Adjacency_list.txt";

    double double0;
    iss.clear(); iss.str(d_opm.at("FDR_cutoff"));
    if (ifstream(CPS().c_str()) and iss>>double0) {
        d_FDRcut=double0;
        if (FDRcut()<0 or FDRcut()>1) throw runtime_error("FDR_cutoff?");
    }


    iss.clear(); iss.str(d_opm.at("BACKGROUND"));
    iss>>d_percent;
    iss>>d_minGene;

    if (SEL().empty() and (not CPS().empty())) {
        d_selbool=false;
        if (not ifstream(CPS().c_str())) throw runtime_error("CPS_table not found");
    } else if ((not SEL().empty()) and CPS().empty()) {
        d_selbool=true;
        if (not ifstream(SEL().c_str())) throw runtime_error("SELECTION not found");
    } else {
        throw runtime_error("SELECTION or CPS_table?");
    }

    d_modulebool=static_cast<bool>(ifstream(module().c_str()));

    iss.clear(); iss.str(d_opm.at("SYNTHESIS"));
    iss>>d_syndeg;
}
