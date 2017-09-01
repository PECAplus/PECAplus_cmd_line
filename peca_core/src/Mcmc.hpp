

class Option;
class Pre;

inline double logit(const double x) { return log(x/(1-x)); }
inline double inv_logit(const double x) { return 1/(1+exp(-x)); }
inline double pq(const double x) { return x*(1-x); }

class Mcmc
{
    const Option d_op;
    const bool d_prior;
    const vector<vector<double> > d_varphi,d_xx;
    const vector<vector<vector<double> > > d_yy;
    vector<double> d_tau2;
    vector<vector<double> > d_y0,d_RR;
    vector<vector<int> > d_CPS,d_CPSum;
    double d_a_tau,d_b_tau;

    void tau_hyperparams();
    void update_y0(const int i);
    void update_R(const int i);
    void update_CPS(const int i);
    void Construct();

    public:
    Mcmc(const Pre& pr);
    Mcmc(const Pre& pr,const vector<vector<double> >& varphi);

    const Option& op() const { return d_op; }
    const vector<vector<double> >& xx() const { return d_xx; }
    const vector<vector<vector<double> > >& yy() const { return d_yy; }
    double a_tau() const { return d_a_tau; }
    double b_tau() const { return d_b_tau; }
    const vector<double>& tau2() const { return d_tau2; }
    const vector<vector<double> >& y0() const { return d_y0; }
    const vector<vector<double> >& RR() const { return d_RR; }
    const vector<vector<int> >& CPS() const { return d_CPS; }
    const vector<vector<int> >& CPSum() const { return d_CPSum; }
};
