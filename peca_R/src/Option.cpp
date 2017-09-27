

#include"main.hpp"
#include"Option.hpp"


Option::Option(const map<string,string>& opm) :
    d_opm(opm),d_xbool(true),d_Mbool(false),d_log_x(true),d_log_m(true),d_log_h(true),
    d_smooth(false),d_level(1)
{
    istringstream iss(d_opm.at("FILE_Y"));
    iss>>d_fileH;
    if (not ifstream(fileH().c_str())) throw runtime_error("H");
    iss>>d_log_h;

    if (not d_opm.at("FILE_M").empty()) {
        d_Mbool=true;
        iss.clear(); iss.str(d_opm.at("FILE_M"));
        iss>>d_fileM;
        if (not ifstream(fileM().c_str())) throw runtime_error("M");
        iss>>d_log_m;
    }

    iss.clear(); iss.str(d_opm.at("FILE_X"));
    iss>>d_filex;
    if (filex().empty()) d_filex=fileH();
    if (not ifstream(filex().c_str())) throw runtime_error("x");
    iss>>d_log_x;

    if (filex()==fileH()) d_xbool=false;

    iss.clear(); iss.str(d_opm.at("N_REP"));
    iss>>d_nrep;

    double double0;

    if (not d_opm.at("SMOOTHING").empty()) {
        d_smooth=true;
        iss.clear(); iss.str(d_opm.at("SMOOTHING"));
        iss>>double0;
        d_smoothing.push_back(double0);
        iss>>double0;
        d_smoothing.push_back(double0);
    }

    iss.clear(); iss.str(d_opm.at("TIME"));
    for (;iss>>double0;d_timep.push_back(double0));

    d_timei=d_timep;
    for (unsigned i=0;i<timep().size();i++) d_timei.at(i)=i;//




    iss.clear(); iss.str(d_opm.at("N_BURN"));
    iss>>d_nburn;

    iss.clear(); iss.str(d_opm.at("N_THIN"));
    iss>>d_nthin;

    iss.clear(); iss.str(d_opm.at("N_SAMPLE"));
    iss>>d_ns;


}
