

class Option;
class Pre;


class Mcmc
{
    const Option d_op;
    const bool d_prior;
    const vector<vector<double> > d_varphiR,d_varphiD,d_xx;
    const vector<vector<vector<double> > > d_yM,d_yH;
    const vector<int> d_nobs;
    vector<double> d_tau2M,d_tau2H;
    vector<vector<double> > d_M0,d_H0,d_RR,d_DD;
    vector<vector<int> > d_CPR,d_CPRum,d_CPD,d_CPDum;
    double d_a_tauM,d_b_tauM,d_a_tauH,d_b_tauH,upRD,d_mean_R,d_mean_D,d_sigma2_RD;

    void tau_hyperparams(const char c);
    void gibbs_tau2(const int i,const char c);
    void update_y0(const int i,const char c);
    void update_RD(const int i,const char c);
    void update_R(const int i);
    void update_D(const int i);
    void update_CP(const int i,const char c);
    void Construct();
    void loc_RD();

    const double logit(const double x,const double u) const { return log(x/(u-x)); }
    const double inv_logit(const double x,const double u) const { return u/(1+exp(-x)); }
    const double pq(const double x,const double u) const { return x*(u-x); }
    const double mean_RD(const char c) const;
    bool Ey_fn(const int p, const int j, const vector<double>& y0i,const vector<double>& Ri,const vector<double>& Di, vector<double>& Ey_ij,const Mcmc& mc,const char c);
    double sum_sq(const int p, const int j, const vector<double>& y0i, const vector<double>& Ri, const vector<double>& Di,const Mcmc& mc,const char c);
    double sum_sq(const int p, const vector<double>& y0i, const vector<double>& Ri, const vector<double>& Di,const Mcmc& mc,const char c,const int s);

    public:
    Mcmc(const Pre& pr);
    Mcmc(const Pre& pr,const vector<vector<double> >& varphiR,const vector<vector<double> >& varphiD);

    const Option& op() const { return d_op; }
    const vector<vector<double> >& xx() const { return d_xx; }
    const vector<vector<vector<double> > >& yM() const { return d_yM; }
    const vector<vector<vector<double> > >& yH() const { return d_yH; }
    double a_tauM() const { return d_a_tauM; }
    double b_tauM() const { return d_b_tauM; }
    double a_tauH() const { return d_a_tauH; }
    double b_tauH() const { return d_b_tauH; }
    double sigma2_RD() const { return d_sigma2_RD; };
    const vector<double>& tau2M() const { return d_tau2M; }
    const vector<double>& tau2H() const { return d_tau2H; }
    const vector<vector<double> >& M0() const { return d_M0; }
    const vector<vector<double> >& H0() const { return d_H0; }
    const vector<vector<double> >& RR() const { return d_RR; }
    const vector<vector<double> >& DD() const { return d_DD; }
    const vector<vector<int> >& CPR() const { return d_CPR; }
    const vector<vector<int> >& CPRum() const { return d_CPRum; }
    const vector<vector<int> >& CPD() const { return d_CPD; }
    const vector<vector<int> >& CPDum() const { return d_CPDum; }
};
