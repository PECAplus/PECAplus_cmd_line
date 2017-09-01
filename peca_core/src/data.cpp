

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
    opm["FILE_Y"];
    opm["N_REP"];
    opm["TIME"];
    opm["N_BURN"]; opm["N_THIN"]; opm["N_SAMPLE"];
    opm["LEVEL"];
    opm["MIN_CORREL"];
    opm["MODULE"]; opm["MODULE_SIZE"];
    opm["PROTEIN_VARIANCE"]; opm["EXPERIMENTAL_DESIGN"];
    opm["SMOOTHING"];

    set<string> opset;
    for (string lstr,rstr;read_param(input_ifs,lstr,rstr,opm);opm[lstr]=rstr) {
        if (opset.find(lstr)!=opset.end()) throw runtime_error("Duplicate \""+lstr+" = \"");
        opset.insert(lstr);
    }

    mo.setModulebool(opm.find("MODULE")->second);

    if (mo.modulebool()) {
        mo.setModule(opm.find("MODULE")->second);
        //mo.setModule_size(opm.find("MODULE_SIZE")->second);
    }

    Option op(opm);

    return op;
}