

class Option;
class Pre;

struct fdrs { double fdr; int t; };

class Post
{
    const Option d_op;
    const Pre d_pr;
    vector<vector<int> > d_up,d_selup,d_seldown,d_selsig;
    vector<vector<double> > d_pvup,d_pvdown,d_pvsig;
    vector<map<double,fdrs,greater<double> > > d_fdrup,d_fdrdown,d_fdrsig;

    void fdr(vector<map<double,fdrs,greater<double> > >& fdr,const int up);
    void selection(vector<vector<int> >& sel,const vector<map<double,fdrs,greater<double> > >& fdr,const int up);
    void test(const vector<vector<int> >& sel,vector<vector<double> >& pval,const int up);


    public:
    Post(const Pre& pr);
    const Option& op() const { return d_op; }
    const Pre& pr() const { return d_pr; }
    const vector<map<double,fdrs,greater<double> > >& fdrup() const { return d_fdrup; }
    const vector<map<double,fdrs,greater<double> > >& fdrdown() const { return d_fdrdown; }
    const vector<map<double,fdrs,greater<double> > >& fdrsig() const { return d_fdrsig; }
    const vector<vector<int> >& selup() const { return d_selup; }
    const vector<vector<int> >& seldown() const { return d_seldown; }
    const vector<vector<int> >& selsig() const { return d_selsig; }
    const vector<vector<double> >& pvup() const { return d_pvup; }
    const vector<vector<double> >& pvdown() const { return d_pvdown; }
    const vector<vector<double> >& pvsig() const { return d_pvsig; }
};
