

class Option
{
    const map<string,string> d_opm;
    string d_CPS,d_GO_term,d_SEL,d_module;
    double d_FDRcut,d_percent;
    int d_minGene;//,d_minSel;
    bool d_selbool,d_modulebool;
    bool d_syndeg;

    public:
    Option(const map<string,string>& opm);
    const string& CPS() const { return d_CPS; }
    const string& GO_term() const { return d_GO_term; }
    const string& module() const { return d_module; }
    const string& SEL() const { return d_SEL; }
    double FDRcut() const { return d_FDRcut; }
    double percent() const { return d_percent; }
    int minGene() const { return d_minGene; }
    bool selbool() const { return d_selbool; }
    bool modulebool() const { return d_modulebool; }
    bool syndeg() const { return d_syndeg; }
};
