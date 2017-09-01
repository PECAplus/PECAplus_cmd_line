

#include"main.hpp"
#include"Option.hpp"
#include"Mcmc.hpp"
#include"Pre.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/gamma_distribution.hpp>
#include"nlopt.hpp"

boost::mt19937 r;

const double sigma2_y0 = 10000;
const bool bounded=false;


void Mcmc::loc_RD()
{
    double protsum=0,mrnasum=0;
    int ncount=0;
    for (unsigned p=0;p<yH().size();p++) for (unsigned q=0;q<yH().at(p).size();q++) for (int j=0;j<op().nrep();j++) for (int t=0;t<op().nt();t++) {
        protsum+=exp(yH().at(p).at(q).at(j*op().nt()+t));
        if (op().Mbool()) protsum+=exp(yM().at(p).at(q).at(j*op().nt()+t));
        mrnasum+=exp(xx().at(p).at(j*op().nt()+t));
        ncount++;
    }
    d_mean_R=0;
    d_mean_D=0;
    d_sigma2_RD=op().Mbool() ? .25 : 1;
}


const double Mcmc::mean_RD(const char c) const
{
    if (c!='R' and c!='D') throw runtime_error("mean_RD");
    return (c=='R'?d_mean_R:d_mean_D);
}


double quantile(const vector<double>& data,const double f)
{
    vector<double> sorted_data(data);
    sort(sorted_data.begin(),sorted_data.end());
    const double index = f * (data.size() - 1) ;
    const size_t lhs = static_cast<size_t>(index) ;
    const double delta = index - lhs ;
    double result;
    if (data.size() == 0) throw runtime_error("quantile");
    if (lhs == data.size() - 1) {
        result = sorted_data.at(lhs) ;
    } else {
        result = (1 - delta) * sorted_data.at(lhs) + delta * sorted_data.at(lhs + 1) ;
    }
    return result ;
}


bool Mcmc::Ey_fn(const int p, const int j, const vector<double>& y0i,const vector<double>& Ri,const vector<double>& Di, vector<double>& Ey_ij,const Mcmc& mc,const char c)
{
    if (c!='M' and c!='H') throw runtime_error("Ey_fn");
    Ey_ij.at(0) = y0i.at(j);
    if (not op().Mbool()) for (int t=1;t<mc.op().nt();t++) {
        Ey_ij.at(t) = Ey_ij.at(t-1) 
            + (mc.op().timep().at(t)-mc.op().timep().at(t-1))*( log(1.+exp(mc.xx().at(p).at(j*mc.op().nt()+t-1))) *Ri.at(t-1)-Ey_ij.at(t-1)*Di.at(t-1) );
        if (Ey_ij.at(t)<=0) return true;
    } else if (c=='M') for (int t=1;t<mc.op().nt();t++) {
        Ey_ij.at(t) = Ey_ij.at(t-1) 
            + (mc.op().timep().at(t)-mc.op().timep().at(t-1))*( -Ey_ij.at(t-1)*Di.at(t-1) );
        if (Ey_ij.at(t)<=0) return true;
    } else for (int t=1;t<mc.op().nt();t++) {
        Ey_ij.at(t) = Ey_ij.at(t-1) 
            + (mc.op().timep().at(t)-mc.op().timep().at(t-1))*( log(1.+exp(mc.xx().at(p).at(j*mc.op().nt()+t-1))) *Ri.at(t-1));
        if (Ey_ij.at(t)<=0) return true;
    }
    return false;
}


double Mcmc::sum_sq(const int p, const int j, const vector<double>& y0i, const vector<double>& Ri, const vector<double>& Di,const Mcmc& mc,const char c)
{
    if (c!='M' and c!='H') throw runtime_error("sum_sqj");
    double sum2=0;
    vector<double> Ey_ij(mc.op().nt()+0);
    if (Ey_fn(p,j,y0i,Ri,Di,Ey_ij,mc,c)) return NAN;
    const vector<vector<double> >& y_p=(c=='M'?mc.yM():mc.yH()).at(p);
    for (int t=0; t<mc.op().nt(); t++) for (unsigned q=0;q<y_p.size();q++) {
        sum2 += pow(y_p.at(q).at(j*mc.op().nt()+t)-log(Ey_ij.at(t)),2);
    }
    return sum2;
}


double Mcmc::sum_sq(const int p, const vector<double>& y0i, const vector<double>& Ri, const vector<double>& Di,const Mcmc& mc,const char c,const int s=0)
{
    if (c!='M' and c!='H') throw runtime_error("sum_sq");
    double sum2=0;
    const vector<vector<double> >& y_p=(c=='M'?mc.yM():mc.yH()).at(p);
    for (int j=0; j<mc.op().nrep(); j++) {
        bool flag=false;
        for (unsigned q=0;q<y_p.size();q++) if (obs(y_p.at(q).at(j*mc.op().nt()))) {
            flag=true;
            break;
        }
        if (flag) {
            vector<double> Ey_ij(mc.op().nt()+0);
            if (Ey_fn(p,j,y0i,Ri,Di,Ey_ij,mc,c)) return NAN;
            for (int t=s; t<mc.op().nt(); t++) for (unsigned q=0;q<y_p.size();q++) {
                sum2 += pow(y_p.at(q).at(j*mc.op().nt()+t)-log(Ey_ij.at(t)),2);
            }
        }
    }
    return sum2;
}



boost::random::normal_distribution<> ugaussian(0,1);
boost::uniform_real<> uniform(0,1);
void Mcmc::update_y0(const int p,const char c)
{
    if (c!='M' and c!='H') throw runtime_error("update_y0");
    const vector<vector<double> >& y_p=(c=='M'?yM():yH()).at(p);
    const double tau2p=(c=='M'?tau2M():tau2H()).at(p);
    const vector<double>& y0p=(c=='M'?M0():H0()).at(p);
    for (int j=0;j<op().nrep();j++) {
        bool flag=false;
        for (unsigned q=0;q<y_p.size();q++) if (obs(y_p.at(q).at(j*op().nt()))) {
            flag=true;
            break;
        }
        if (flag) {
            vector<double> newy0p(y0p);
            newy0p.at(j) *= exp(ugaussian(r));
            double logr = -sum_sq(p,j,newy0p,RR().at(p),DD().at(p),*this,c)/tau2p/2;
            logr -= -sum_sq(p,j,y0p,RR().at(p),DD().at(p),*this,c)/tau2p/2;
            logr += -pow(log(newy0p.at(j)),2)/2/sigma2_y0;
            logr -= -pow(log(y0p.at(j)),2)/2/sigma2_y0;
            if (logr>0 or uniform(r)<exp(logr)) {
                (c=='M'?d_M0: d_H0).at(p).at(j) = newy0p.at(j);
            }
        }
    }
}


void Mcmc::update_RD(const int p,const char c)
{
    if (c!='R' and c!='D') throw runtime_error("update_RD");
    const vector<int>& CPp=(c=='R'?CPR():CPD()).at(p);
    const vector<double>& RDp=(c=='R'?RR():DD()).at(p);
    const double& upl=upRD;
    vector<double> newRDp;
    int s=-1;
    for (int t=0;t<op().nt()-1;t++) {
        if (/*t==0 or */CPp.at(t)==1) {
            s = t+1;
            newRDp = RDp;
            if (bounded) {
                newRDp.at(t) = inv_logit(logit(RDp.at(t),upl)+ugaussian(r),upl);
            } else {
                newRDp.at(t) *= exp(ugaussian(r));
            }
        } else {
            newRDp.at(t) = newRDp.at(t-1);
        }
        if (/*t!=op().nt()-2 and */CPp.at(t+1)==0) continue;
        double logr = -sum_sq(p,H0().at(p),(c=='R'?newRDp:RR().at(p)),(c=='R'?DD().at(p):newRDp),*this,'H',s)/tau2H().at(p)/2;
        logr -= -sum_sq(p,H0().at(p),RR().at(p),DD().at(p),*this,'H',s)/tau2H().at(p)/2;
        if (op().Mbool()) {
            logr += -sum_sq(p,M0().at(p),(c=='R'?newRDp:RR().at(p)),(c=='R'?DD().at(p):newRDp),*this,'M',s)/tau2M().at(p)/2;
            logr -= -sum_sq(p,M0().at(p),RR().at(p),DD().at(p),*this,'M',s)/tau2M().at(p)/2;
        }
        if (bounded) {
            logr += log(pq(newRDp.at(t),upl))-log(pq(RDp.at(t),upl));
        } else {
            logr += -pow(log(newRDp.at(t))-mean_RD(c),2)/2/sigma2_RD();
            logr -= -pow(log(RDp.at(t))-mean_RD(c),2)/2/sigma2_RD();
        }
        if (logr>0 or uniform(r)<exp(logr)) {
            (c=='R'?d_RR:d_DD).at(p) = newRDp;
        }
    }
}


void Mcmc::update_CP(const int p,const char c)
{
    if (c!='R' and c!='D') throw runtime_error("update_CP");
    const vector<double>& RDp=(c=='R'?RR():DD()).at(p);
    const vector<int>& CPp=(c=='R'?CPR():CPD()).at(p);
    const vector<vector<double> >& d_varphi=(c=='R'?d_varphiR:d_varphiD);
    const double& upl=upRD;
    for (int cp=1;cp<op().nt()-1;cp++) {
        vector<double> newRDp=RDp;
        int l_time_index=-1, r_time_index=-1;
        for (int t=cp-1;l_time_index<0;t--) if (/*t==0 or */CPp.at(t)==1) l_time_index=t;
        for (int t=cp+1;r_time_index<0;t++) if (/*t==op().nt()-1 or */CPp.at(t)==1) r_time_index=t;
        const double l_time = op().timep().at(cp)-op().timep().at(l_time_index);
        const double r_time = op().timep().at(r_time_index)-op().timep().at(cp);
        const double t_time = op().timep().at(r_time_index)-op().timep().at(l_time_index);
        double logr=0;
        if (CPp.at(cp)==0) {
            const double u = uniform(r);
            const double newRleft = (bounded ? inv_logit(logit(RDp.at(cp),upl)+r_time/t_time*logit(u,1),upl) : exp(log(RDp.at(cp))+r_time/t_time*logit(u,1)) );
            const double newRright= (bounded ? inv_logit(logit(RDp.at(cp),upl)-l_time/t_time*logit(u,1),upl) : exp(log(RDp.at(cp))-l_time/t_time*logit(u,1)) );
            for (int t=l_time_index;t<cp;t++) newRDp.at(t) = newRleft;
            for (int t=cp;t<r_time_index;t++) newRDp.at(t) = newRright;
            if (d_prior) logr += d_varphi.at(p).at(cp);
            if (bounded) {
                logr += -log(upl);
                logr += 2*log(newRleft*(upl-newRright)+newRright*(upl-newRleft))-log(pq(RDp.at(cp),upl))-log(upl);
            } else {
                logr -= -log(RDp.at(cp))-pow(log(RDp.at(cp))-mean_RD(c),2)/2/sigma2_RD();
                logr += -log(newRleft)-pow(log(newRleft)-mean_RD(c),2)/2/sigma2_RD();
                logr += -log(newRright)-pow(log(newRright)-mean_RD(c),2)/2/sigma2_RD();
                logr += -log(sigma2_RD())-.5*log(2*M_PI);
                logr += 2*log(newRleft+newRright)-log(RDp.at(cp));
            }
        } else {
            const double mergedR=(bounded ? inv_logit((l_time*logit(RDp.at(cp-1),upl)+r_time*logit(RDp.at(cp),upl))/t_time,upl) : exp((l_time*log(RDp.at(cp-1))+r_time*log(RDp.at(cp)))/t_time) );
            for (int t=l_time_index;t<r_time_index;t++) newRDp.at(t)=mergedR;
            if (d_prior) logr -= d_varphi.at(p).at(cp);
            //if ((not op().Mbool()) and c=='D') logr-=logit(.1,1);
            if (bounded) {
                logr -= -log(upl);
                logr -= 2*log(RDp.at(cp-1)*(upl-RDp.at(cp))+RDp.at(cp)*(upl-RDp.at(cp-1)))-log(pq(mergedR,upl))-log(upl);
            } else {
                logr += -log(mergedR)-pow(log(mergedR)-mean_RD(c),2)/2/sigma2_RD();
                logr -= -log(RDp.at(cp-1))-pow(log(RDp.at(cp-1))-mean_RD(c),2)/2/sigma2_RD();
                logr -= -log(RDp.at(cp))-pow(log(RDp.at(cp))-mean_RD(c),2)/2/sigma2_RD();
                logr -= -log(sigma2_RD())-.5*log(2*M_PI);
                logr -= 2*log(RDp.at(cp-1)+RDp.at(cp))-log(mergedR);
            }
        }
        logr += -sum_sq(p,H0().at(p),(c=='R'?newRDp:RR().at(p)),(c=='R'?DD().at(p):newRDp),*this,'H',l_time_index+1)/tau2H().at(p)/2;
        logr -= -sum_sq(p,H0().at(p),RR().at(p),DD().at(p),*this,'H',l_time_index+1)/tau2H().at(p)/2;
        if (op().Mbool()) {
            logr += -sum_sq(p,M0().at(p),(c=='R'?newRDp:RR().at(p)),(c=='R'?DD().at(p):newRDp),*this,'M',l_time_index+1)/tau2M().at(p)/2;
            logr -= -sum_sq(p,M0().at(p),RR().at(p),DD().at(p),*this,'M',l_time_index+1)/tau2M().at(p)/2;
        }
        if (logr>0 or uniform(r)<exp(logr)) {
            if (c=='R') {
                d_RR.at(p) = newRDp;
                d_CPR.at(p).at(cp)=1-CPp.at(cp);
            } else {
                d_DD.at(p) = newRDp;
                d_CPD.at(p).at(cp)=1-CPp.at(cp);
            }
        }
    }
}


void Mcmc::gibbs_tau2(const int p,const char c)
{
    if (c!='M' and c!='H') throw runtime_error("gibbs_tau2");
    const double sumsq=sum_sq(p,(c=='M'?M0():H0()).at(p),RR().at(p),DD().at(p),*this,c);
    if (obs(sumsq)) {
        (c=='M'?d_tau2M:d_tau2H).at(p) = 1 / boost::random::gamma_distribution<> ((c=='M'?a_tauM():a_tauH())+.5*d_nobs.at(p)*op().nt(),1/((c=='M'?b_tauM():b_tauH())+sumsq/2))(r);
    }
}


void Mcmc::tau_hyperparams(const char c)
{
    if (c!='M' and c!='H') throw runtime_error("tau_hyperparams");
    double& d_a_tau=(c=='M'?d_a_tauM:d_a_tauH);
    double& d_b_tau=(c=='M'?d_b_tauM:d_b_tauH);
    vector<double>& d_tau2=(c=='M'?d_tau2M:d_tau2H);
    const vector<double>& tau2=d_tau2;

    vector<vector<vector<double> > > yy(yH());
    if (op().Mbool()) for (unsigned p=0;p<yH().size();p++) for (unsigned q=0;q<yH().at(p).size();q++) for (int j=0;j<op().nrep();j++) for (int t=0;t<op().nt();t++) {
        yy.at(p).at(q).at(j*op().nt()+t)=log(exp(yM().at(p).at(q).at(j*op().nt()+t))+exp(yH().at(p).at(q).at(j*op().nt()+t)));
    }

    d_tau2.assign(yy.size(),0);
    for (unsigned p=0;p<yy.size();p++) {
        for (unsigned q=0; q<yy.at(p).size(); q++) {
            for (int j=0;j<op().nrep();j++) {
                vector<double> mu(op().nt()-1);
                for (int t=0;t<op().nt()-1;t++) {
                    mu.at(t)=yy.at(p).at(q).at(j*op().nt()+t)+yy.at(p).at(q).at(j*op().nt()+t+1);
                    mu.at(t)/=2;
                }
                for (int t=1;t<op().nt()-1;t++) {
                    d_tau2.at(p)+=.5*pow(yy.at(p).at(q).at(j*op().nt()+t)-mu.at(t-1),2);
                    d_tau2.at(p)+=.5*pow(yy.at(p).at(q).at(j*op().nt()+t)-mu.at(t),2);
                }
                d_tau2.at(p)+=pow(yy.at(p).at(q).at(j*op().nt()+0)-mu.at(0),2);
                const int t=op().nt()-2;
                d_tau2.at(p)+=pow(yy.at(p).at(q).at(j*op().nt()+t)-mu.at(t),2);
            }
        }
        d_tau2.at(p) /= op().nc()*yy.at(p).size();
    }
    const double m=quantile(d_tau2,.1);
    for (unsigned p=0;p<yy.size();p++) d_tau2.at(p) += m;
    ofstream var_ofs("varc"+boost::lexical_cast<string>(c)+".txt");

    //std::transform(tau2.begin(),tau2.end(),d_tau2.begin(),std::bind1st(std::multiplies<double>(),2.));

    copy(tau2.begin(),tau2.end(),ostream_iterator<double>(var_ofs,"\n"));
    const double mom1=accumulate(tau2.begin(),tau2.end(),0.)/tau2.size();
    double mom2=0;
    for (vector<double>::const_iterator it=tau2.begin();it!=tau2.end();it++) mom2 += (*it)*(*it);
    mom2 /= tau2.size();
    d_a_tau=(2*mom2-mom1*mom1)/(mom2-mom1*mom1);//MOM estimates
    d_b_tau=mom1*mom2/(mom2-mom1*mom1);         //MOM estimates
}


Mcmc::Mcmc(const Pre& pr) : d_op(pr.op()),d_prior(false),d_xx(pr.xx()),d_yM(pr.yM()),d_yH(pr.yH()),d_nobs(pr.nobs())
{
    Construct();
}


Mcmc::Mcmc(const Pre& pr,const vector<vector<double> >& varphiR,const vector<vector<double> >& varphiD) : d_op(pr.op()),d_prior(true),d_varphiR(varphiR),d_varphiD(varphiD),d_xx(pr.xx()),d_yM(pr.yM()),d_yH(pr.yH()),d_nobs(pr.nobs())
{
    Construct();
}




void Mcmc::Construct()
{
    d_M0.assign(yM().size(),vector<double>(op().nrep()));
    d_H0.assign(yH().size(),vector<double>(op().nrep()));
    for (unsigned p=0;p<yM().size();p++) for (int j=0;j<op().nrep();j++) {
        double sum=0;
        for (unsigned q=0;q<yM().at(p).size();q++) sum+=yM().at(p).at(q).at(j*op().nt());
        d_M0.at(p).at(j)=exp(sum/yM().at(p).size());
    }
    for (unsigned p=0;p<yH().size();p++) for (int j=0;j<op().nrep();j++) {
        double sum=0;
        for (unsigned q=0;q<yH().at(p).size();q++) sum+=yH().at(p).at(q).at(j*op().nt());
        d_H0.at(p).at(j)=exp(sum/yH().at(p).size());
    }
    d_RR.assign(yH().size(),vector<double> (op().nt(),1e-5));
    d_DD.assign(yH().size(),vector<double> (op().nt(),1e-5));

    loc_RD();


    vector<int> setcp(op().nt(),1);setcp.front()=1;setcp.back()=1;
    d_CPR.assign(yH().size(),setcp);
    d_CPRum.assign(yH().size(),vector<int> (op().nt()-1));
    d_CPD.assign(yH().size(),setcp);
    d_CPDum.assign(yH().size(),vector<int> (op().nt()-1));

    if (op().Mbool()) tau_hyperparams('M');
    tau_hyperparams('H');
    ofstream M0_ofs("s_M0");
    ofstream H0_ofs("s_H0");
    ofstream RR_ofs("s_RR");
    ofstream DD_ofs("s_DD");
    ofstream tau2M_ofs("s_tau2M");
    ofstream tau2H_ofs("s_tau2H");
    if (not tau2M_ofs) throw runtime_error("can't open s_tau2M");
    ofstream CPR_ofs("s_CPR");
    ofstream CPD_ofs("s_CPD");
    ofstream etaM_ofs("s_etaM");
    ofstream etaH_ofs("s_etaH");
    ofstream loglike_ofs("s_loglike");

    upRD=10;

    for (int iter=0,ns=0; ns<op().ns(); iter++) {

        if (iter%100==0) cout<<'\r'<<iter<<flush;
        for (unsigned p=0; p<yH().size(); p++) {
            if (op().Mbool()) update_y0(p,'M');
            update_y0(p,'H');
            update_RD(p,'R');
            update_RD(p,'D');
            update_CP(p,'R');
            update_CP(p,'D');
        }

        if (iter>op().nburn() and iter%op().nthin()==0) {
            ns++;

            for (unsigned p=0;p<yH().size();p++) {
                if (op().Mbool()) tau2M_ofs<<d_tau2M.at(p)<<'\t';
                tau2H_ofs<<d_tau2H.at(p)<<'\t';
                for (int j=0;j<op().nrep();j++) {
                    if (op().Mbool()) M0_ofs<<d_M0.at(p).at(j)<<'\t';
                    H0_ofs<<d_H0.at(p).at(j)<<'\t';
                }
                for (int t=0;t<op().nt()-1;t++) {
                    RR_ofs<<d_RR.at(p).at(t)<<'\t';
                    DD_ofs<<d_DD.at(p).at(t)<<'\t';
                }
                for (int t=1;t<op().nt()-1;t++) {
                    CPR_ofs<<d_CPR.at(p).at(t)<<'\t';
                    CPD_ofs<<d_CPD.at(p).at(t)<<'\t';
                    d_CPRum.at(p).at(t)+=d_CPR.at(p).at(t);
                    d_CPDum.at(p).at(t)+=d_CPD.at(p).at(t);
                }

            }

            if (op().Mbool()) {
                M0_ofs<<'\n';
                tau2M_ofs<<'\n';
            }
            H0_ofs<<'\n';
            tau2H_ofs<<'\n';
            RR_ofs<<'\n';
            DD_ofs<<'\n';
            CPR_ofs<<'\n';
            CPD_ofs<<'\n';

            double loglike = 0;
            for (unsigned p=0;p<yH().size();p++) for (unsigned q=0;q<yH().at(p).size();q++) for (int j=0;j<op().nrep();j++) {
                if (op().Mbool()) {
                    if (obs(yM().at(p).at(q).at(j*op().nt()))) {
                        vector<double> Ey_ij(op().nt());
                        if (Ey_fn(p,j,d_M0.at(p),d_RR.at(p),d_DD.at(p),Ey_ij,*this,'M')) throw runtime_error("fit");
                        for (int t=0;t<op().nt();t++) {
                            etaM_ofs << log(Ey_ij.at(t)) << '\t';
                            loglike += -log(d_tau2M.at(p))-pow(yM().at(p).at(q).at(j*op().nt()+t)-log(Ey_ij.at(t)),2)/d_tau2M.at(p)/2;
                        }
                    } else {
                        for (int t=0;t<op().nt();t++) etaM_ofs << "NA\t";
                    }
                }
                if (obs(yH().at(p).at(q).at(j*op().nt()))) {
                    vector<double> Ey_ij(op().nt());
                    if (Ey_fn(p,j,d_H0.at(p),d_RR.at(p),d_DD.at(p),Ey_ij,*this,'H')) throw runtime_error("fit");
                    for (int t=0;t<op().nt();t++) {
                        etaH_ofs << log(Ey_ij.at(t)) << '\t';
                        loglike += -log(d_tau2H.at(p))-pow(yH().at(p).at(q).at(j*op().nt()+t)-log(Ey_ij.at(t)),2)/d_tau2H.at(p)/2;
                    }
                } else {
                    for (int t=0;t<op().nt();t++) etaH_ofs << "NA\t";
                }
            }
            if (op().Mbool()) etaM_ofs<<'\n';
            etaH_ofs<<'\n';
            loglike_ofs << loglike << '\n';
        }

    }
    cout<<'\n';
}





