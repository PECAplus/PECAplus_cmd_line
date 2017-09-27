

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"
#include"Post.hpp"

#include <boost/math/distributions/hypergeometric.hpp>


void Post::fdr(vector<map<double,fdrs,greater<double> > >& fdr,const int up)
{
    fdr.resize(1?1:d_pr.nl());
    for (unsigned l=0;l<d_pr.nl();l++) {
        map<double,fdrs,greater<double> >& fdr_l(fdr.at(fdr.size()==1?0:l));
        for (unsigned p=0;p<d_pr.np();p++) {
            fdr_l[up==-1 or d_up.at(l).at(p)==up?pr().cp().at(p).at(l):0].t++;
        }
    }
    for (unsigned l=0;l<fdr.size();l++) {
        map<double,fdrs,greater<double> >& fdr_l(fdr.at(l));
        int denom=fdr_l.begin()->second.t;
        double numer=(1-fdr_l.begin()->first)*denom;
        map<double,fdrs>::iterator it=fdr_l.begin();
        for (it++;it!=fdr_l.end();it++) {
            it->second.fdr=numer/denom;
            numer += (1-it->first)*it->second.t;
            denom += it->second.t;
        }
    }
}


void Post::selection(vector<vector<int> >& sel,const vector<map<double,fdrs,greater<double> > >& fdr,const int up)
{
    double minCPS=1;
    for (unsigned l=0;l<d_pr.nl();l++) {
        const map<double,fdrs,greater<double> >& fdr_l(fdr.at(fdr.size()==1?0:l));
        sel.at(l).assign(d_pr.np(),0);
        for (unsigned p=0;p<d_pr.np();p++) {
            if (fdr_l.at(up==-1 or d_up.at(l).at(p)==up?pr().cp().at(p).at(l):0).fdr<d_op.FDRcut()) {
                sel.at(l).at(p)=1;
                if (pr().cp().at(p).at(l)<minCPS) minCPS=pr().cp().at(p).at(l);
            }
        }
    }
}


void Post::test(const vector<vector<int> >& sel,vector<vector<double> >& pval,const int up)
{
    for (unsigned l=0;l<d_pr.nl();l++) {
        const int r=accumulate(sel.at(l).begin(),sel.at(l).end(),0);
        for (unsigned g=0;g<d_pr.ng();g++) {
            int sum=0;
            for (unsigned m=0;m<d_pr.mem().at(g).size();m++) {
                sum+=sel.at(l).at(d_pr.mem().at(g).at(m));
            }
            boost::math::hypergeometric a(r, d_pr.mem().at(g).size(), d_pr.np());
            pval.at(l).at(g)=cdf(complement(a,sum))+pdf(a,sum);
        }
    }
}


Post::Post(const Pre& pr) :
    d_op(pr.op()),d_pr(pr),d_up(pr.nl()),d_selup(pr.nl()),d_seldown(pr.nl()),d_selsig(pr.nl()),d_pvup(pr.nl(),vector<double>(pr.ng())),d_pvdown(pr.nl(),vector<double>(pr.ng())),d_pvsig(pr.nl(),vector<double>(pr.ng()))
{
    if (op().selbool()) {
        d_selup=pr.selup();
    } else {
        for (unsigned l=0;l<d_pr.nl();l++) {
            d_up.at(l).assign(d_pr.np(),0);
            for (unsigned p=0;p<d_pr.np();p++) {
                if (d_pr.R().at(p).at(l)<d_pr.R().at(p).at(l+1)) d_up.at(l).at(p)=1;
            }
        }
        fdr(d_fdrup,1);
        selection(d_selup,d_fdrup,1);
    }
    test(d_selup,d_pvup,1);

    if (op().selbool()) {
        d_seldown=pr.seldown();
    } else {
        fdr(d_fdrdown,0);
        selection(d_seldown,d_fdrdown,0);
    }
    test(d_seldown,d_pvdown,0);

    if (op().selbool()) {
        d_selsig=pr.selsig();
    } else {
        fdr(d_fdrsig,-1);
        selection(d_selsig,d_fdrsig,-1);
    }
    test(d_selsig,d_pvsig,-1);
}
