

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"


bool read_param(ifstream& ifs,string& lstr,string& rstr,const map<string,string>& opm)
{
    istringstream liss;
    for (string str0;getline(ifs,str0);) {
        str0=str0.substr(0,str0.find("#"));
        const size_t e=str0.find("=");
        if (e!=string::npos) {
            liss.str(str0.substr(0,e));
            liss>>lstr;
            if (opm.find(lstr)==opm.end()) throw runtime_error("Unknown parameter: "+lstr);
            rstr=str0.substr(e+1);
            return true;
        }
    }
    return false;
}


Option set_param(const string& filepath)
{
    ifstream input_ifs(filepath.c_str());
    if (not input_ifs) throw runtime_error("can't open "+filepath);

    map<string,string> opm;
    opm["GO_term_table"];
    opm["FDR_cutoff"];
    opm["BACKGROUND"];
    opm["SYNTHESIS"]; 

    set<string> opset;
    for (string lstr,rstr;read_param(input_ifs,lstr,rstr,opm);opm[lstr]=rstr) {
        if (opset.find(lstr)!=opset.end()) throw runtime_error("Duplicate \""+lstr+" = \"");
        opset.insert(lstr);
    }

    Option op(opm);
    return op;
}
