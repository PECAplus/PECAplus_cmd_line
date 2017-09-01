

class Mcmc;

class Est
{
    const Mcmc d_mc;
    const Module d_mo;
    vector<double> d_Phi;
    vector<vector<double> > d_varphi;

    public:
    Est(const Mcmc& mc,const Module& mo);

    const Mcmc& mc() const { return d_mc; }
    const vector<double>& Phi() const { return d_Phi; }
    const vector<vector<double> >& varphi() const { return d_varphi; }
};
