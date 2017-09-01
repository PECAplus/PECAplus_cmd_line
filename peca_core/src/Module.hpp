

class Module
{
    bool d_modulebool;
    string d_module;//d_module_type,d_mrf;
    int d_min_size,d_max_size;//,d_module_freq;
    vector<set<int> > d_adj;


    public:
    Module() : d_modulebool(false),
    d_min_size(0),d_max_size(INT_MAX){}//,d_module_freq(0) { }
    void setModulebool(const string& str0);
    bool modulebool() const { return d_modulebool; }
    void setModule(const string& str0);
    const string& module() const { return d_module; }
    const vector<set<int> >& adj() const { return d_adj; }
    void module_data(const vector<string>& pidvec);
    void setModule_size(const string& str0);
    int max_size() const { return d_max_size; }
    int min_size() const { return d_min_size; }
};
