

class Option;


struct Row
{
    string pid,qid;
    vector<double> in;
    int nrepsi;
    double m_cor;
    Row(const int nc) :  pid(""),qid(""),in(vector<double>(nc)),nrepsi(nc),m_cor(-1) { }
};


class Pre
{
    list<Row> d_gdata0,d_gdata,d_mdata0,d_mdata,d_hdata0,d_hdata;
    vector<string> d_xh,d_Mh,d_Hh,d_xid,d_pid;
    vector<vector<string> > d_qid;
    vector<int> d_nobs;
    vector<vector<double> > d_xx;
    vector<vector<vector<double> > > d_yM,d_yH;
    const Option d_op;

    void setrepcor(list<Row>& );
    void setMc();
    void read();
    void setAll();
    double Mfn(const vector<double>& m,const double xs,const vector<double>& ll);
    double Kfn(const vector<double>&theta,const double xp,const double xq,const int e);

    public:
    Pre(const Option& op);
    void impute_x();
    void impute_y(const char c);
    void norm_data();
    const Option& op() const { return d_op; };
    int nprot() const { return d_gdata.size(); }
    const vector<int>& nobs() const { return d_nobs; }
    const vector<string>& xh() const { return d_xh; }
    const vector<string>& Mh() const { return d_Mh; }
    const vector<string>& Hh() const { return d_Hh; }
    const vector<string>& xid() const { return d_xid; }
    const vector<string>& pid() const { return d_pid; }
    const vector<vector<string> >& qid() const { return d_qid; }
    const vector<vector<double> >& xx() const { return d_xx; }
    const vector<vector<vector<double> > >& yM() const { return d_yM; }
    const vector<vector<vector<double> > >& yH() const { return d_yH; }
    const list<Row>& gdata() const { return d_gdata; }
    const list<Row>& gdata0() const { return d_gdata0; }
    const list<Row>& mdata() const { return d_mdata; }
    const list<Row>& mdata0() const { return d_mdata0; }
    const list<Row>& hdata() const { return d_hdata; }
    const list<Row>& hdata0() const { return d_hdata0; }
};
