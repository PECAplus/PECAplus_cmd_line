

class Option;


struct Row
{
    string pid,qid;
    vector<double> in;
    int nrepsi;//number of obs per gene
    double m_cor;
    Row(const int nc) :  pid(""),qid(""),in(vector<double>(nc)),nrepsi(nc) { }
};


class Pre
{
    list<Row> d_gdata0,d_gdata,d_pdata0,d_pdata;
    vector<string> d_xh,d_yh,d_xid,d_pid;
    vector<vector<string> > d_qid;
    vector<vector<double> > d_xx;
    vector<vector<vector<double> > > d_yy;
    const Option d_op;

    void setMc();
    void read();
    void setAll();
    double Mfn(const vector<double>& m,const double xs,const vector<double>& ll);
    double Kfn(const vector<double>&theta,const double xp,const double xq,const int e);

    public:
    Pre(const Option& op);
    void impute_x();
    void impute_y();
    void norm_data();
    const Option& op() const { return d_op; };
    int nprot() const { return d_gdata.size(); }
    const vector<string>& xh() const { return d_xh; }
    const vector<string>& yh() const { return d_yh; }
    const vector<string>& xid() const { return d_xid; }
    const vector<string>& pid() const { return d_pid; }
    const vector<vector<string> >& qid() const { return d_qid; }
    const vector<vector<double> >& xx() const { return d_xx; }
    const vector<vector<vector<double> > >& yy() const { return d_yy; }
    const list<Row>& gdata() const { return d_gdata; }
    const list<Row>& gdata0() const { return d_gdata0; }
    const list<Row>& pdata() const { return d_pdata; }
    const list<Row>& pdata0() const { return d_pdata0; }
};
