

class Option;


struct Row
{
    string pid; //protein id
    vector<double> R; //rate ratio
    vector<double> cp; //change point prob.
};


struct Row0
{
    string GO;//GO id
    string GOf;//GO function
    vector<string> mem;//member genes
};


struct Row1
{
    string pid; //protein id
    vector<int> sel; //selection
};


class Pre
{
    const Option d_op;
    list<Row> d_cdata;
    list<Row0> d_GOdata;
    list<Row1> d_sdata;
    set<string> d_pidset;
    vector<string> d_pidvec,d_GOid,d_GOf;
    vector<vector<double> > d_R,d_cp;
    vector<vector<int> > d_mem,d_selup,d_seldown,d_selsig;
    vector<vector<bool> > d_adj;
    vector<vector<string> > d_mem_all;
    void read_sel();
    void set_sel();
    void read_cp();
    void set_cp();
    void read_GO();
    void set_GO();
    void read_net();


    public:
    Pre(const Option& op);
    const Option& op() const { return d_op; }
    const list<Row>& cdata() const { return d_cdata;}
    const list<Row1>& sdata() const { return d_sdata;}
    const list<Row0>& GOdata() const { return d_GOdata;}
    const set<string>& pidset() const { return d_pidset;}
    const vector<string>& pidvec() const { return d_pidvec;}
    const vector<string>& GOid() const { return d_GOid;}
    const vector<string>& GOf() const { return d_GOf;}
    const vector<vector<double> >& cp() const { return d_cp;};
    const vector<vector<double> >& R() const { return d_R;};
    const vector<vector<int> >& mem() const { return d_mem;};
    const vector<vector<string> >& mem_all() const { return d_mem_all;};
    const vector<vector<int> >& selup() const { return d_selup; }
    const vector<vector<int> >& seldown() const { return d_seldown; }
    const vector<vector<int> >& selsig() const { return d_selsig; }
    const vector<vector<bool> >& adj() const { return d_adj; }
    size_t ng() const { return d_mem.size(); }
    size_t np() const;//{ return d_cp.size(); }
    size_t nl() const;
};
