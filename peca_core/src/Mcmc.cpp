

#include"main.hpp"
#include"Option.hpp"
#include"Mcmc.hpp"
#include"Pre.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/gamma_distribution.hpp>

boost::mt19937 r;

const double sigma2_y0 = 10000;


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


bool Ey_fn(const int p, const int j, const vector<double>& y0i,const vector<double>& Ri, vector<double>& Ey_ij,const Mcmc& mc)
{
    Ey_ij.at(0) = y0i.at(j);
    for (int t=1;t<mc.op().nt();t++) {
        Ey_ij.at(t) = Ey_ij.at(t-1) 
            + (mc.op().timep().at(t)-mc.op().timep().at(t-1))*(exp(mc.xx().at(p).at(j*mc.op().nt()+t-1))*Ri.at(t-1)-Ey_ij.at(t-1)*(1-Ri.at(t-1)));
        if (Ey_ij.at(t)<=0) return true;
    }
    return false;
}


double sum_sq(const int p, const int j, const vector<double>& y0i, const vector<double>& Ri,const Mcmc& mc)
{
    double sum2=0;
    vector<double> Ey_ij(mc.op().nt());
    if (Ey_fn(p,j,y0i,Ri,Ey_ij,mc)) return NAN;
    for (int t=0; t<mc.op().nt(); t++) for (unsigned q=0;q<mc.yy().at(p).size();q++) {
        sum2 += pow(mc.yy().at(p).at(q).at(j*mc.op().nt()+t)-log(Ey_ij.at(t)),2);
    }
    return sum2;
}


double sum_sq(const int p, const vector<double>& y0i, const vector<double>& Ri,const Mcmc& mc,const int s=0)
{
    double sum2=0;
    for (int j=0; j<mc.op().nrep(); j++) {
        bool flag=false;
        for (unsigned q=0;q<mc.yy().at(p).size();q++) if (obs(mc.yy().at(p).at(q).at(j*mc.op().nt()))) {
            flag=true;
            break;
        }
        if (flag) {
            vector<double> Ey_ij(mc.op().nt());
            if (Ey_fn(p,j,y0i,Ri,Ey_ij,mc)) return NAN;
            for (int t=s; t<mc.op().nt(); t++) for (unsigned q=0;q<mc.yy().at(p).size();q++) {
                sum2 += pow(mc.yy().at(p).at(q).at(j*mc.op().nt()+t)-log(Ey_ij.at(t)),2);
            }
        }
    }
    return sum2;
}


boost::random::normal_distribution<> ugaussian(0,1);
boost::uniform_real<> uniform(0,1);
void Mcmc::update_y0(const int p)
{
    for (int j=0;j<op().nrep();j++) {
        bool flag=false;
        for (unsigned q=0;q<yy().at(p).size();q++) if (obs(yy().at(p).at(q).at(j*op().nt()))) {
            flag=true;
            break;
        }
        if (flag) {
            vector<double> newy0p(y0().at(p));
            newy0p.at(j) *= exp(ugaussian(r));
            double logr = -sum_sq(p,j,newy0p,RR().at(p),*this)/tau2().at(p)/2;
            logr -= -sum_sq(p,j,y0().at(p),RR().at(p),*this)/tau2().at(p)/2;
            logr += -pow(log(newy0p.at(j)),2)/2/sigma2_y0;
            logr -= -pow(log(y0().at(p).at(j)),2)/2/sigma2_y0;
            if (logr>0 or uniform(r)<exp(logr)) {
                d_y0.at(p).at(j) = newy0p.at(j);
            }
        }
    }
}


void Mcmc::update_R(const int p)
{
    vector<double> newRp;
    int s=-1;
    for (int t=0;t<op().nt()-1;t++) {
        if (CPS().at(p).at(t)==1) {
            s = t+1;
            newRp = RR().at(p);
            newRp.at(t) = inv_logit(logit(RR().at(p).at(t))+ugaussian(r));
        } else {
            newRp.at(t) = newRp.at(t-1);
        }
        if (CPS().at(p).at(t+1)==0) continue;
        double logr = -sum_sq(p,y0().at(p),newRp,*this,s)/tau2().at(p)/2;
        logr -= -sum_sq(p,y0().at(p),RR().at(p),*this,s)/tau2().at(p)/2;
        logr += log(pq(newRp.at(t)))-log(pq(RR().at(p).at(t)));
        if (logr>0 or uniform(r)<exp(logr)) {
            d_RR.at(p) = newRp;
        }
    }
}




void Mcmc::update_CPS(const int p)
{
    const vector<double>& Rp=RR().at(p);
    for (int cp=1;cp<op().nt()-1;cp++) {
        vector<double> newRp(Rp);
        int l_time_index=-1, r_time_index=-1;
        for (int t=cp-1;l_time_index<0;t--) if (CPS().at(p).at(t)==1) l_time_index=t;
        for (int t=cp+1;r_time_index<0;t++) if (CPS().at(p).at(t)==1) r_time_index=t;
        const double l_time = op().timep().at(cp)-op().timep().at(l_time_index);
        const double r_time = op().timep().at(r_time_index)-op().timep().at(cp);
        const double t_time = op().timep().at(r_time_index)-op().timep().at(l_time_index);
        double logr=0;
        if (CPS().at(p).at(cp)==0) {
            const double u = uniform(r);
            const double newRleft = inv_logit(logit(Rp.at(cp))+r_time/t_time*logit(u));
            const double newRright = inv_logit(logit(Rp.at(cp))-l_time/t_time*logit(u));
            for (int t=l_time_index;t<cp;t++) newRp.at(t) = newRleft;
            for (int t=cp;t<r_time_index;t++) newRp.at(t) = newRright;
            if (d_prior) logr += d_varphi.at(p).at(cp);
            //if (d_prior) logr += varphi_fn(p,cp,*this,d_mo);
            logr += 2*log(newRleft*(1-newRright)+newRright*(1-newRleft))-log(pq(Rp.at(cp)));
        } else {
            const double mergedR = inv_logit((l_time*logit(Rp.at(cp-1))+r_time*logit(Rp.at(cp)))/t_time);
            for (int t=l_time_index;t<r_time_index;t++) newRp.at(t)=mergedR;
            if (d_prior) logr -= d_varphi.at(p).at(cp);
            //if (d_prior) logr -= varphi_fn(p,cp,*this,d_mo);
            logr -= 2*log(Rp.at(cp-1)*(1-Rp.at(cp))+Rp.at(cp)*(1-Rp.at(cp-1)))-log(pq(mergedR));
        }
        logr += -sum_sq(p,y0().at(p),newRp,*this,l_time_index+1)/tau2().at(p)/2;
        logr -= -sum_sq(p,y0().at(p),Rp,*this,l_time_index+1)/tau2().at(p)/2;
        if (logr>0 or uniform(r)<exp(logr)) {
            d_RR.at(p) = newRp;
            d_CPS.at(p).at(cp)=1-CPS().at(p).at(cp);
        }
    }
}




void Mcmc::tau_hyperparams()
{

        d_tau2.assign(yy().size(),0);
            for (unsigned p=0;p<yy().size();p++) {
                for (unsigned q=0; q<yy().at(p).size(); q++) {
                    for (int j=0;j<op().nrep();j++) {
                        vector<double> mu(op().nt()-1);
                        for (int t=0;t<op().nt()-1;t++) {
                            mu.at(t)=yy().at(p).at(q).at(j*op().nt()+t)+yy().at(p).at(q).at(j*op().nt()+t+1);
                            mu.at(t)/=2;
                        }
                        for (int t=1;t<op().nt()-1;t++) {
                            d_tau2.at(p)+=.5*pow(yy().at(p).at(q).at(j*op().nt()+t)-mu.at(t-1),2);
                            d_tau2.at(p)+=.5*pow(yy().at(p).at(q).at(j*op().nt()+t)-mu.at(t),2);
                        }
                        d_tau2.at(p)+=pow(yy().at(p).at(q).at(j*op().nt()+0)-mu.at(0),2);
                        const int t=op().nt()-2;
                        d_tau2.at(p)+=pow(yy().at(p).at(q).at(j*op().nt()+t)-mu.at(t),2);
                    }
                }
                d_tau2.at(p) /= op().nc()*yy().at(p).size();
            }
            const double m=quantile(d_tau2,.1);
            for (unsigned p=0;p<yy().size();p++) d_tau2.at(p) += m;
        ofstream var_ofs("varc.txt");
        copy(tau2().begin(),tau2().end(),ostream_iterator<double>(var_ofs,"\n"));
        const double mom1=accumulate(tau2().begin(),tau2().end(),0.)/tau2().size();
        double mom2=0;
        for (vector<double>::const_iterator it=tau2().begin();it!=tau2().end();it++) mom2 += (*it)*(*it);
        mom2 /= tau2().size();
        d_a_tau=(2*mom2-mom1*mom1)/(mom2-mom1*mom1);//MOM estimates
        d_b_tau=mom1*mom2/(mom2-mom1*mom1);         //MOM estimates
}


Mcmc::Mcmc(const Pre& pr) : d_op(pr.op()),d_prior(false),d_varphi(vector<vector<double> >()),d_xx(pr.xx()),d_yy(pr.yy())//,d_nobs(pr.nobs())
{
    Construct();
}


Mcmc::Mcmc(const Pre& pr,const vector<vector<double> >& varphi) : d_op(pr.op()),d_prior(true),d_varphi(varphi),d_xx(pr.xx()),d_yy(pr.yy())//,d_nobs(pr.nobs())
{
    Construct();
}




void Mcmc::Construct()
{
    d_y0.assign(yy().size(),vector<double>(op().nrep()));
    for (unsigned p=0;p<yy().size();p++) for (int j=0;j<op().nrep();j++) {
        double sum=0;
        for (unsigned q=0;q<yy().at(p).size();q++) sum+=yy().at(p).at(q).at(j*op().nt());
        d_y0.at(p).at(j)=exp(sum/yy().at(p).size());
    }
    d_RR.assign(yy().size(),vector<double> (op().nt(),.9));
    d_CPS.assign(yy().size(),vector<int> (op().nt(),1));
    d_CPSum.assign(yy().size(),vector<int> (op().nt()));

    tau_hyperparams();
    ofstream y0_ofs("s_y0");
    ofstream RR_ofs("s_RR");
    ofstream tau2_ofs("s_tau2");
    if (not tau2_ofs) throw runtime_error("can't open s_tau2");
    ofstream CPS_ofs("s_CPS");
    ofstream xRyD_ofs("s_xRyD");
    ofstream loglike_ofs("s_loglike");


    for (int iter=0,ns=0; ns<op().ns(); iter++) {

        if (iter%100==0) cout<<'\r'<<iter<<flush;
        for (unsigned p=0; p<yy().size(); p++) {
            update_y0(p);
            update_R(p);
            update_CPS(p);
        }

        if (iter>op().nburn() and iter%op().nthin()==0) {
            ns++;

            for (unsigned p=0;p<yy().size();p++) {
                tau2_ofs<<d_tau2.at(p)<<'\t';
                for (int j=0;j<op().nrep();j++) y0_ofs<<d_y0.at(p).at(j)<<'\t';
                for (int t=0;t<op().nt()-1;t++) {
                    RR_ofs<<d_RR.at(p).at(t)<<'\t';
                }
                for (int t=1;t<op().nt()-1;t++) {
                    CPS_ofs<<d_CPS.at(p).at(t)<<'\t';
                    d_CPSum.at(p).at(t)+=d_CPS.at(p).at(t);
                }

            }

            y0_ofs<<'\n';
            tau2_ofs<<'\n';
            RR_ofs<<'\n';
            CPS_ofs<<'\n';

            double loglike = 0;
            for (unsigned p=0;p<yy().size();p++) for (unsigned q=0;q<yy().at(p).size();q++) for (int j=0;j<op().nrep();j++) {
                if (obs(yy().at(p).at(q).at(j*op().nt()))) {
                    vector<double> Ey_ij(op().nt());
                    if (Ey_fn(p,j,d_y0.at(p),d_RR.at(p),Ey_ij,*this)) throw runtime_error("fit");
                    for (int t=0;t<op().nt();t++) {
                        xRyD_ofs << log(Ey_ij.at(t)) << '\t';
                        //likelihood
                        loglike += -log(d_tau2.at(p))-pow(yy().at(p).at(q).at(j*op().nt()+t)-log(Ey_ij.at(t)),2)/d_tau2.at(p)/2;
                    }
                } else for (int t=0;t<op().nt();t++) xRyD_ofs << "NA\t";
            }
            xRyD_ofs<<'\n';
            loglike_ofs << loglike << '\n';
        }
    }
    cout<<'\n';
}
