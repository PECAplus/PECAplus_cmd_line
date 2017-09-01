

#include"main.hpp"
#include"Option.hpp"


Option::Option(const map<string,string>& opm) :
    d_opm(opm),d_xbool(true),d_Mbool(false),d_log_x(true),d_log_m(true),d_log_h(true),
    d_smooth(false),d_level(1)
{
    istringstream iss(d_opm.find("FILE_X")->second);
    iss>>d_filex;
    if (not ifstream(filex().c_str())) throw runtime_error("x");
    iss>>d_log_x;

    if (not d_opm.find("FILE_M")->second.empty()) {
        d_Mbool=true;
        iss.clear(); iss.str(d_opm.find("FILE_M")->second);
        iss>>d_fileM;
        if (not ifstream(fileM().c_str())) throw runtime_error("M");
        iss>>d_log_m;
    }

    iss.clear(); iss.str(d_opm.find("FILE_Y")->second);
    iss>>d_fileH;
    if (not ifstream(fileH().c_str())) throw runtime_error("H");
    iss>>d_log_h;

    if (filex()==fileH()) d_xbool=false;

    iss.clear(); iss.str(d_opm.find("N_REP")->second);
    iss>>d_nrep;

    double double0;

    if (not d_opm.find("SMOOTHING")->second.empty()) {
        d_smooth=true;
        iss.clear(); iss.str(d_opm.find("SMOOTHING")->second);
        iss>>double0;
        d_smoothing.push_back(double0);
        iss>>double0;
        d_smoothing.push_back(double0);
    }

    iss.clear(); iss.str(d_opm.find("TIME")->second);
    for (;iss>>double0;d_timep.push_back(double0));

    d_timei=d_timep;
    for (unsigned i=0;i<timep().size();i++) d_timei.at(i)=i;//




    iss.clear(); iss.str(d_opm.find("N_BURN")->second);
    iss>>d_nburn;

    iss.clear(); iss.str(d_opm.find("N_THIN")->second);
    iss>>d_nthin;

    iss.clear(); iss.str(d_opm.find("N_SAMPLE")->second);
    iss>>d_ns;

    iss.clear(); iss.str(d_opm.find("LEVEL")->second);
    iss>>d_level;

    if (level()==2) {
        iss.clear(); iss.str(d_opm.find("MIN_CORREL")->second);
        if (not (iss>>d_min_correl) or (min_correl()<-1 or min_correl()>1)) throw runtime_error("MIN_CORREL?");
    }
}