

#include"main.hpp"
#include"Module.hpp"
#include"Option.hpp"

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


Option set_param(const string& filepath,Module& mo)
{
    ifstream input_ifs(filepath.c_str());

    map<string,string> opm;
    opm["FILE_X"];
    opm["FILE_M"];
    opm["FILE_Y"];
    opm["N_REP"];
    opm["TIME"];
    opm["N_BURN"]; opm["N_THIN"]; opm["N_SAMPLE"];
    opm["SMOOTHING"];
    opm["MODULE"];
    set<string> opset;
    for (string lstr,rstr;read_param(input_ifs,lstr,rstr,opm);opm[lstr]=rstr) {
        if (opset.find(lstr)!=opset.end()) throw runtime_error("Duplicate \""+lstr+" = \"");
        opset.insert(lstr);
    }

    mo.setModulebool(opm.at("MODULE"));

    if (mo.modulebool()) {
        mo.setModule(opm.at("MODULE"));

    }

    Option op(opm);

    return op;
}
