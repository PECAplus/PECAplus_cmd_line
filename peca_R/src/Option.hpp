



class Option
{
    const map<string,string> d_opm;
    bool d_xbool,d_Mbool,d_log_x,d_log_m,d_log_h,d_smooth;
    string d_filex,d_fileM,d_fileH,d_protvar,d_rep;
    int d_nrep,d_nburn,d_nthin,d_ns,d_level;
    vector<double> d_timep,d_timei,d_smoothing;
    double d_min_correl;

    public:
    Option(const map<string,string>& opm);

    const string& filex() const { return d_filex; }
    const string& fileM() const { return d_fileM; }
    const string& fileH() const { return d_fileH; }
    const string& protvar() const { return d_protvar; }
    const string& rep() const { return d_rep; }
    const vector<double>& timep() const { return d_timep; }
    const vector<double>& timei() const { return d_timei; }
    const vector<double>& smoothing() const { return d_smoothing; }
    int nt() const { return timep().size(); }
    int nrep() const { return d_nrep; }
    int nburn() const { return d_nburn; }
    int nthin() const { return d_nthin; }
    int ns() const { return d_ns; }
    int nc() const { return nrep()*timep().size(); }
    int level() const { return d_level; }
    bool xbool() const { return d_xbool; }
    bool Mbool() const { return d_Mbool; }
    bool log_x() const { return d_log_x; }
    bool log_h() const { return d_log_h; }
    bool log_m() const { return d_log_m; }
    double min_correl() const { return d_min_correl; }
    bool smooth() const { return d_smooth; }
};
