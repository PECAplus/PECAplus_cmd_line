

class Mcmc;

class Est
{
    const Mcmc d_mc;
    const Module d_mo;
    vector<double> d_Phi;
    vector<vector<double> > d_varphiR,d_varphiD;
    void get_varphi(const char c);

    public:
    Est(const Mcmc& mc,const Module& mo);

    const Mcmc& mc() const { return d_mc; }
    const vector<double>& Phi() const { return d_Phi; }
    const vector<vector<double> >& varphiR() const { return d_varphiR; }
    const vector<vector<double> >& varphiD() const { return d_varphiD; }
};
