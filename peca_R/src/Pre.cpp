

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

using boost::lexical_cast;
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
typedef multimap<string,Row*>::const_iterator mmit;
typedef multimap<string,const Row*>::const_iterator cmmit;



class pred_mc
{
    const double d_mc_cut;
    public:
    pred_mc(const double mc_cut) : d_mc_cut(mc_cut) { }
    bool operator() (const Row& tmprow) const { return tmprow.m_cor<d_mc_cut; }
};


class pred_g
{
    const set<string> d_set;
    public:
    pred_g(const set<string>& tmpset) : d_set(tmpset) {}
    bool operator() (const Row& tmprow) const { return d_set.end()==d_set.find(tmprow.pid); }
};


class pred_nr
{
    const int d_nr;
    public:
    pred_nr(const int nr) : d_nr(nr) {}
    bool operator() (const Row& tmpco) const { return tmpco.nrepsi<d_nr; }
};


void Pre::setrepcor(list<Row>& mh_data)
{
    for (list<Row>::iterator it=mh_data.begin();it!=mh_data.end();it++) {
        const vector<double> rep1(it->in.begin()+0*op().nt(),it->in.begin()+1*op().nt());
        const vector<double> rep2(it->in.begin()+1*op().nt(),it->in.begin()+2*op().nt());
        vector<double> center(2);
        vector<double> sd(center.size());
        int ncol0=0;
        for (int t=0;t<op().nt();t++) if (obs(rep1.at(t)) and obs(rep2.at(t))) {
            center.at(0) += rep1.at(t);
            center.at(1) += rep2.at(t);
            ncol0++;
        }
        if (ncol0<3) continue;
        center.at(0) /= ncol0;
        center.at(1) /= ncol0;
        for (int t=0;t<op().nt();t++) if (obs(rep1.at(t)) and obs(rep2.at(t))) {
            sd.at(0) += pow(rep1.at(t)-center.at(0),2);
            sd.at(1) += pow(rep2.at(t)-center.at(1),2);
        }
        sd.at(0) = sqrt(sd.at(0));
        sd.at(1) = sqrt(sd.at(1));
        double cortmp=0;
        for (int t=0;t<op().nt();t++) if (obs(rep1.at(t)) and obs(rep2.at(t))) {
            cortmp += (rep1.at(t)-center.at(0))*(rep2.at(t)-center.at(1));
        }
        if (sd.at(0)==0 or sd.at(1)==0) continue;
        cortmp /= sd.at(0);
        cortmp /= sd.at(1);
        it->m_cor=cortmp;
    }
}


void Pre::setMc()
{
    multimap<string,Row*> pid_mm;
    for (list<Row>::iterator it=d_mdata.begin();it!=d_mdata.end();it++) {
        pid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (mmit end1,itp=pid_mm.begin();itp!=pid_mm.end();itp=end1) {
        end1=pid_mm.upper_bound(itp->first);
        vector<vector<double> > ww;
        for (mmit it=itp;it!=end1;it++) ww.push_back(it->second->in);
        if (ww.size()==1) {
            itp->second->m_cor=1;
            continue;
        }
        mmit it=itp;
        for (unsigned f=0;f<ww.size();f++,it++) {
            vector<double> cor;
            for (unsigned g=0;g<ww.size();g++) if (f != g) {
                vector<double> center(2);
                vector<double> sd(center.size());
                int ncol0=0;
                for (int l=0;l<op().nc();l++) if (obs(ww.at(f).at(l)) and obs(ww.at(g).at(l))) {
                    center.at(0) += ww.at(f).at(l);
                    center.at(1) += ww.at(g).at(l);
                    ncol0++;
                }
                if (ncol0<2) continue;
                center.at(0) /= ncol0;
                center.at(1) /= ncol0;
                for (int l=0;l<op().nc();l++) if (obs(ww.at(f).at(l)) and obs(ww.at(g).at(l))) {
                    sd.at(0) += pow(ww.at(f).at(l)-center.at(0),2);
                    sd.at(1) += pow(ww.at(g).at(l)-center.at(1),2);
                }
                sd.at(0) = sqrt(sd.at(0));
                sd.at(1) = sqrt(sd.at(1));
                double cortmp=0;
                for (int l=0;l<op().nc();l++) if (obs(ww.at(f).at(l)) and obs(ww.at(g).at(l))) {
                    cortmp += (ww.at(f).at(l)-center.at(0))*(ww.at(g).at(l)-center.at(1));
                }
                if (sd.at(0)==0 or sd.at(1)==0) continue;
                cortmp /= sd.at(0);
                cortmp /= sd.at(1);
                cor.push_back(cortmp);
            }
            if (cor.size()==0) {
                it->second->m_cor=-1;
            } else {
                it->second->m_cor=median(cor.begin(),cor.end());
            }
        }
    }
}


void Pre::norm_data()
{
    if (op().xbool()) {
        for (int p=0;p<nprot();p++) {
            vector<double> gvec;
            for (int l=0;l<op().nc();l++) if (obs(xx().at(p).at(l))) {
                gvec.push_back(xx().at(p).at(l));
            }
            const double center_g=accumulate(gvec.begin(),gvec.end(),0.)/gvec.size();
            for (int j=0;j<op().nrep();j++) {//centering by reps
                vector<double> xvec;
                for (int t=0; t<op().nt(); t++) if (obs(xx().at(p).at(j*op().nt()+t))) {
                    xvec.push_back(xx().at(p).at(j*op().nt()+t));
                }
                const double center_x=accumulate(xvec.begin(),xvec.end(),0.)/xvec.size();
                for (int t=0; t<op().nt(); t++) if (obs(xx().at(p).at(j*op().nt()+t))) {
                    d_xx.at(p).at(j*op().nt()+t) += center_g-center_x;
                }
            }
        }
    }
    for (int p=0;p<nprot();p++) for (unsigned q=0;q<yM().at(p).size();q++) {
        vector<double> gvec;
        for (int l=0;l<op().nc();l++) if (obs(yM().at(p).at(q).at(l))) {
            gvec.push_back(yM().at(p).at(q).at(l));
        }
        const double center_g=accumulate(gvec.begin(),gvec.end(),0.)/gvec.size();
        for (int j=0;j<op().nrep();j++) {
            vector<double> yvec;
            for (int t=0; t<op().nt(); t++) if (obs(yM().at(p).at(q).at(j*op().nt()+t))) {
                yvec.push_back(yM().at(p).at(q).at(j*op().nt()+t));
            }
            double center_y=accumulate(yvec.begin(),yvec.end(),0.)/yvec.size();
            for (int t=0; t<op().nt(); t++) if (obs(yM().at(p).at(q).at(j*op().nt()+t))) {
                d_yM.at(p).at(q).at(j*op().nt()+t) += center_g-center_y;
            }
        }
    }
    for (int p=0;p<nprot();p++) for (unsigned q=0;q<yH().at(p).size();q++) {
        vector<double> gvec;
        for (int l=0;l<op().nc();l++) if (obs(yH().at(p).at(q).at(l))) {
            gvec.push_back(yH().at(p).at(q).at(l));
        }
        const double center_g=accumulate(gvec.begin(),gvec.end(),0.)/gvec.size();
        for (int j=0;j<op().nrep();j++) {
            vector<double> yvec;
            for (int t=0; t<op().nt(); t++) if (obs(yH().at(p).at(q).at(j*op().nt()+t))) {
                yvec.push_back(yH().at(p).at(q).at(j*op().nt()+t));
            }
            double center_y=accumulate(yvec.begin(),yvec.end(),0.)/yvec.size();
            for (int t=0; t<op().nt(); t++) if (obs(yH().at(p).at(q).at(j*op().nt()+t))) {
                d_yH.at(p).at(q).at(j*op().nt()+t) += center_g-center_y;
            }
        }
    }

}


Pre::Pre(const Option& op) : d_op(op)
{
    read();
    setAll();
}


void Pre::read()
{
    ifstream M_ifs(op().fileM().c_str());
    ifstream H_ifs(op().fileH().c_str());
    ifstream x_ifs(op().filex().c_str());

    string str0;

    if (op().Mbool()) {
        getline(M_ifs,str0);
        split(d_Mh, str0, boost::is_any_of("\t"));
        if (int(Mh().size())<op().nc()) throw runtime_error("file_M");
    }

    getline(H_ifs,str0);
    split(d_Hh, str0, boost::is_any_of("\t"));
    if (int(Hh().size())<op().nc()) throw runtime_error("file_H");

    getline(x_ifs,str0);
    split(d_xh, str0, boost::is_any_of("\t"));
    for (int i=2;getline(x_ifs,str0) and (not str0.empty());i++) {
        Row tmp_f(op().nc());
        tokenizer tok{str0, boost::char_separator<char>{"\t","",boost::keep_empty_tokens}};
        tokenizer::const_iterator tok_it = tok.begin();
        tmp_f.pid = *tok_it; tok_it++;
        for (int l=0;l<op().nc();l++,tok_it++) {
            if (tok_it==tok.end()) throw runtime_error("FILE_X unequal columns: row "+lexical_cast<string>(i));
            if ((*tok_it)=="NA" or (*tok_it).empty() or lexical_cast<double>(*tok_it)==0) {
                tmp_f.in.at(l)=NAN;
            } else {
                try {
                    if (op().log_x()) tmp_f.in.at(l)=log(lexical_cast<double>(*tok_it));
                    else tmp_f.in.at(l)=(lexical_cast<double>(*tok_it));
                } catch (boost::bad_lexical_cast& e) {
                    cerr<<e.what()<<"FILE_X : row "<<i<<", column "<<l+1+op().level();
                    exit(1);
                }
            }
        }
        d_gdata0.push_back(tmp_f);
    }

    for (int i=2;getline(M_ifs,str0) and (not str0.empty());i++) {
        Row tmp_f(op().nc());
        tokenizer tok{str0, boost::char_separator<char>{"\t","",boost::keep_empty_tokens}};
        tokenizer::const_iterator tok_it = tok.begin();
        tmp_f.pid = *tok_it; tok_it++;
        if (op().level()==1) tmp_f.qid=tmp_f.pid;
        else if (op().level()==2) {tmp_f.qid = *tok_it; tok_it++;}
        for (int l=0;l<op().nc();l++,tok_it++) {
            if (tok_it==tok.end()) throw runtime_error("FILE_M unequal columns: row "+lexical_cast<string>(i)+" "+lexical_cast<string>(l));
            if (*tok_it=="NA" or (*tok_it).empty() or lexical_cast<double>(*tok_it)==0) {
                tmp_f.in.at(l)=NAN;
            } else {
                try {
                    if (op().log_m()) tmp_f.in.at(l)=log(lexical_cast<double>(*tok_it));
                    else tmp_f.in.at(l)=(lexical_cast<double>(*tok_it));
                } catch (boost::bad_lexical_cast& e) {
                    cerr<<e.what()<<"FILE_M : row "<<i<<", column "<<l+1+op().level();
                    exit(1);
                }
            }
        }
        d_mdata0.push_back(tmp_f);
    }

    for (int i=2;getline(H_ifs,str0) and (not str0.empty());i++) {
        Row tmp_f(op().nc());
        tokenizer tok{str0, boost::char_separator<char>{"\t","",boost::keep_empty_tokens}};
        tokenizer::const_iterator tok_it = tok.begin();
        tmp_f.pid = *tok_it; tok_it++;
        if (op().level()==1) tmp_f.qid=tmp_f.pid;
        else if (op().level()==2) {tmp_f.qid = *tok_it; tok_it++;}
        for (int l=0;l<op().nc();l++,tok_it++) {
            if (tok_it==tok.end()) throw runtime_error("FILE_H unequal columns: row "+lexical_cast<string>(i));
            if (*tok_it=="NA" or (*tok_it).empty() or lexical_cast<double>(*tok_it)==0) {
                tmp_f.in.at(l)=NAN;
            } else {
                try {
                    if (op().log_h()) tmp_f.in.at(l)=log(lexical_cast<double>(*tok_it));
                    else tmp_f.in.at(l)=(lexical_cast<double>(*tok_it));
                } catch (boost::bad_lexical_cast& e) {
                    cerr<<e.what()<<"FILE_H : row "<<i<<", column "<<l+1+op().level();
                    exit(1);
                }
            }
        }
        d_hdata0.push_back(tmp_f);
    }
}


void Pre::setAll()
{
    multimap<string,const Row*> gid_mm;
    for (list<Row>::const_iterator it=gdata0().begin();it!=gdata0().end();it++) {
        gid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (cmmit end1,itp=gid_mm.begin();itp!=gid_mm.end();itp=end1) {
        end1=gid_mm.upper_bound(itp->first);
        Row tmp_f(op().nc());
        tmp_f.pid=itp->first;
        vector<int> cobs(op().nc());
        for (cmmit it=itp;it!=end1;it++) {
            for (int i=0;i<op().nc();i++) {
                if (obs(it->second->in.at(i))) {
                    tmp_f.in.at(i)+=it->second->in.at(i);
                    cobs.at(i)++;
                }
            }
        }
        for (int i=0;i<op().nc();i++) {
            if (cobs.at(i)>0) tmp_f.in.at(i) /= cobs.at(i); else tmp_f.in.at(i)=NAN;
        }
        d_gdata.push_back(tmp_f);
    }

    multimap<string,const Row*> pqid_mm;
    for (list<Row>::const_iterator it=mdata0().begin();it!=mdata0().end();it++) {
        pqid_mm.insert(make_pair(it->pid+it->qid,&(*it)));
    }
    for (cmmit end1,itp=pqid_mm.begin();itp!=pqid_mm.end();itp=end1) {
        end1=pqid_mm.upper_bound(itp->first);
        Row tmp_f(op().nc());
        tmp_f.pid=itp->second->pid;
        tmp_f.qid=itp->second->qid;
        vector<int> cobs(op().nc());
        for (cmmit it=itp;it!=end1;it++) {
            for (int i=0;i<op().nc();i++) {
                if (obs(it->second->in.at(i))) {
                    tmp_f.in.at(i)+=it->second->in.at(i);
                    cobs.at(i)++;
                }
            }
        }
        for (int i=0;i<op().nc();i++) {
            if (cobs.at(i)>0) tmp_f.in.at(i) /= cobs.at(i); else tmp_f.in.at(i)=NAN;
        }
        d_mdata.push_back(tmp_f);
    }

    pqid_mm.clear();
    for (list<Row>::const_iterator it=hdata0().begin();it!=hdata0().end();it++) {
        pqid_mm.insert(make_pair(it->pid+it->qid,&(*it)));
    }
    for (cmmit end1,itp=pqid_mm.begin();itp!=pqid_mm.end();itp=end1) {
        end1=pqid_mm.upper_bound(itp->first);
        Row tmp_f(op().nc());
        tmp_f.pid=itp->second->pid;
        tmp_f.qid=itp->second->qid;
        vector<int> cobs(op().nc());
        for (cmmit it=itp;it!=end1;it++) {
            for (int i=0;i<op().nc();i++) {
                if (obs(it->second->in.at(i))) {
                    tmp_f.in.at(i)+=it->second->in.at(i);
                    cobs.at(i)++;
                }
            }
        }
        for (int i=0;i<op().nc();i++) {
            if (cobs.at(i)>0) tmp_f.in.at(i) /= cobs.at(i); else tmp_f.in.at(i)=NAN;
        }
        d_hdata.push_back(tmp_f);
    }

    for (list<Row>::iterator it=d_gdata.begin();it!=d_gdata.end();it++) {
        for (int l=0;l<op().nc();l++) {
            if (not obs(it->in.at(l))) { it->nrepsi--; /*break;*/ }
        }
        //cout<<it->nrepsi<<' ';
    }
    for (list<Row>::iterator it=d_mdata.begin();it!=d_mdata.end();it++) {
        for (int l=0;l<op().nc();l++) {
            if (not obs(it->in.at(l))) { it->nrepsi--; /*break;*/ }
        }
    }
    for (list<Row>::iterator it=d_hdata.begin();it!=d_hdata.end();it++) {
        for (int l=0;l<op().nc();l++) {
            if (not obs(it->in.at(l))) { it->nrepsi--; /*break;*/ }
        }
    }

    cout<<"0gdata "<<gdata().size()<<" mdata "<<mdata().size()<<" hdata "<<hdata().size()<<'\n';
    //d_gdata.remove_if(pred_nr(op().nc()-2));
    //d_mdata.remove_if(pred_nr(op().nc()-2));
    //d_hdata.remove_if(pred_nr(op().nc()-2));

    //if (op().Mbool() and op().nrep()>1 and op().fileH()!="synH.txt") {
    //    //cout<<"thisgdata "<<gdata().size()<<" mdata "<<mdata().size()<<" hdata "<<hdata().size()<<'\n';
    //    setrepcor(d_gdata);
    //    d_gdata.remove_if(pred_mc( 0.2));
    //    setrepcor(d_hdata);
    //    d_hdata.remove_if(pred_mc( 0.2));
    //    setrepcor(d_mdata);
    //    d_mdata.remove_if(pred_mc( 0.2));
    //}

    //setMc();
    //d_mdata.remove_if(pred_mc(op().min_correl()));

    set<string> gset;
    for (list<Row>::const_iterator it=gdata().begin();it!=gdata().end();it++) {
        gset.insert(it->pid);
    }
    d_mdata.remove_if(pred_g(gset));
    d_hdata.remove_if(pred_g(gset));

    if (op().Mbool()) {
        gset.clear();
        for (list<Row>::const_iterator it=mdata().begin();it!=mdata().end();it++) {
            gset.insert(it->pid);
        }
        d_gdata.remove_if(pred_g(gset));
        d_hdata.remove_if(pred_g(gset));
    }

    gset.clear();
    for (list<Row>::const_iterator it=hdata().begin();it!=hdata().end();it++) {
        gset.insert(it->pid);
    }
    d_gdata.remove_if(pred_g(gset));
    d_mdata.remove_if(pred_g(gset));

    cout<<"gdata "<<gdata().size()<<" mdata "<<mdata().size()<<" hdata "<<hdata().size()<<'\n';

    map<string,const Row*> gid_m;
    for (list<Row>::const_iterator it=gdata().begin();it!=gdata().end();it++) {
        gid_m[it->pid]=&(*it);
    }

    multimap<string,const Row*> mid_mm;
    for (list<Row>::const_iterator it=mdata().begin();it!=mdata().end();it++) {
        mid_mm.insert(make_pair(it->pid,&(*it)));
    }

    multimap<string,const Row*> hid_mm;
    for (list<Row>::const_iterator it=hdata().begin();it!=hdata().end();it++) {
        hid_mm.insert(make_pair(it->pid,&(*it)));
    }

    cout<<"nprot = "<<nprot()<<'\n';
    d_nobs.resize(nprot());
    d_xid.resize(nprot());
    d_pid.resize(nprot());
    d_qid.resize(nprot());
    d_xx.resize(nprot(),vector<double>(op().nc()));
    d_yM.resize(nprot());
    d_yH.resize(nprot());
    int p=0;
    multimap<string,const Row*>::const_iterator endm,itm=mid_mm.begin();
    multimap<string,const Row*>::const_iterator endh,ith=hid_mm.begin();
    for (map<string,const Row*>::const_iterator itg=gid_m.begin();itg!=gid_m.end();itg++,itm=endm,ith=endh,p++) {
        d_nobs.at(p)=ith->second->nrepsi;
        d_xid.at(p)=itg->second->pid;
        d_pid.at(p)=ith->second->pid;
        for (int l=0;l<op().nc();l++) d_xx.at(p).at(l)=itg->second->in.at(l);
        endm=mid_mm.upper_bound(itm->first);
        for (cmmit it=itm;it!=endm;it++) {
            d_yM.at(p).push_back(it->second->in);
        }
        endh=hid_mm.upper_bound(ith->first);
        for (cmmit it=ith;it!=endh;it++) {
            d_qid.at(p).push_back(it->second->qid);
            d_yH.at(p).push_back(it->second->in);
        }
    }
    if (not op().xbool()) {
        for (int g=0;g<nprot();g++) for (int l=0;l<op().nc();l++) d_xx.at(g).at(l)=0;
    }

}


#include"Eigen/Dense"


double Pre::Kfn(const vector<double>&theta,const double xp,const double xq,const int e=1)
{
    if (not op().smooth()) return exp(-pow((xp-xq),2)/2);
    const static double ll = op().smoothing().at(1);
    return pow(op().smoothing().at(0),2)*exp(-pow((xp-xq)/ll,2)/2)+e*(xp==xq);
}


double Pre::Mfn(const vector<double>& m,const double xs,const vector<double>& adjdiff)
{
    if (op().timei().front()>xs or xs>op().timei().back()) throw runtime_error("xs range");
    int fxs=-1;
    if (xs==op().timei().back()) return m.back();
    for (int i=op().nt()-2;i>=0;i--) if (op().timei().at(i)<=xs) {fxs=i; break;}
    return m.at(fxs)+(m.at(fxs+1)-m.at(fxs))/adjdiff.at(fxs+1)*(xs-op().timei().at(fxs));
}


typedef multimap<int,int>::const_iterator cmmii;




void Pre::impute_x()
{
    double** data = new double*[nprot()];
    int** mask = new int*[nprot()];
    for (int i = 0; i < nprot(); i++) {
        data[i] = new double[op().nc()];
        mask[i] = new int[op().nc()];
    }

    vector<vector<double> > mvec(nprot(),vector<double>(op().nrep()));
    vector<vector<double> > vvec(nprot(),vector<double>(op().nrep()));
    for (int p=0;p<nprot();p++) for (int j=0;j<op().nrep();j++) {
        int nj=0;
        const vector<double>& xxp=xx().at(p);
        for (int t=0;t<op().nt();t++) if (obs(xxp.at(j*op().nt()+t))) {
            mvec.at(p).at(j)+=xxp.at(j*op().nt()+t);
            nj++;
        }
        mvec.at(p).at(j) /= nj;
        for (int t=0;t<op().nt();t++) if (obs(xxp.at(j*op().nt()+t))) {
            vvec.at(p).at(j) += pow(xxp.at(j*op().nt()+t)-mvec.at(p).at(j),2);
        }
        vvec.at(p).at(j) /= nj-1;
        if (vvec.at(p).at(j)==0) vvec.at(p).at(j)=1;
        for (int t=0;t<op().nt();t++) if (obs(xxp.at(j*op().nt()+t))) {
            d_xx.at(p).at(j*op().nt()+t) -= mvec.at(p).at(j);
            d_xx.at(p).at(j*op().nt()+t) /= sqrt(vvec.at(p).at(j));
        }
    }

    for (int p=0;p<nprot();p++) for (int l=0; l<op().nc(); l++) {
        data[p][l]=0;
        int nq=0;
        if (obs(xx().at(p).at(l))) {
            data[p][l] += xx().at(p).at(l);
            nq++;
        }
        if (nq>0) {
            data[p][l] /= nq;
            mask[p][l] = 1;
        } else {
            data[p][l] = NAN;
            mask[p][l] = 0;
        }
    }


    vector<int> clusterid(nprot());

    multimap<int,int> c_mm;
    for (int p=0;p<nprot();p++) c_mm.insert(make_pair(clusterid[p],p));

    vector<vector<Eigen::VectorXd> > xax,yy_;
    xax.resize(xx().size());
    yy_.resize(xx().size());
    for (int p=0;p<nprot();p++) {
        xax.at(p).resize(op().nrep());
        yy_.at(p).resize(op().nrep());
        for (int j=0;j<op().nrep();j++) {
            vector<double> x,y;
            for (int t=0;t<op().nt();t++) {
                if (obs(xx().at(p).at(j*op().nt()+t))) {
                    x.push_back(op().timei().at(t));
                    y.push_back(xx().at(p).at(j*op().nt()+t));
                }
            }
            xax.at(p).at(j).resize(x.size());
            yy_.at(p).at(j).resize(y.size());
            for (unsigned i=0;i<x.size();i++) {
                xax.at(p).at(j)(i)=x.at(i);
                yy_.at(p).at(j)(i)=y.at(i);
            }
        }
    }

    const int Ns=100;
    const double step=double(op().timei().back()-op().timei().front())/Ns;
    Eigen::VectorXd xs(Ns);
    for (int i=0;i<Ns;i++) xs(i)=op().timei().at(0)+step/2+i*step;

    ofstream ofs1("EfsX.txt");
    ofs1<<"p\tq\tj";
    for (unsigned i=0;i<xs.size();i++) ofs1<<'\t'<<xs(i);
    ofs1<<'\n';


    int ijk=0;
    vector<double> adjdiff(op().timei().size());
    adjacent_difference(op().timei().begin(),op().timei().end(),adjdiff.begin());

    for (cmmii end1,itc=c_mm.begin();itc!=c_mm.end()/*,ijk<8*/;itc=end1,ijk++) {
        end1=c_mm.upper_bound(itc->first);

        vector<vector<double> > m(op().nrep(),vector<double>(op().nt()));

        vector<double> theta(3,1);

        for (cmmii it=itc;it!=end1;it++) {
            const int p=it->second;
            for (int j=0;j<op().nrep();j++) {
                Eigen::VectorXd& y = yy_.at(p).at(j);
                Eigen::VectorXd& x = xax.at(p).at(j);
                const int N = y.size();
                //theta.at(0)=max(1,log10(exp(mvec.at(p).at(j))));
                Eigen::MatrixXd Ky(N,N);
                for (int i=0;i<N;i++) for (int k=0;k<N;k++) Ky(i,k)=Kfn(theta,x(i),x(k));
                Eigen::VectorXd ym(y);
                for (int i=0;i<N;i++) ym(i)-=Mfn(m.at(j),x(i),adjdiff);
                Eigen::MatrixXd invKy=Ky.inverse();
                ofs1<<p<<"\t0\t"<<j;
                for (unsigned i=0;i<xs.size();i++) {
                    Eigen::VectorXd ks(N);
                    for (int k=0;k<N;k++) ks(k)=Kfn(theta,xs(i),x(k),0);
                    const double Efs=ks.transpose()*invKy*ym+Mfn(m.at(j),xs(i),adjdiff);
                    ofs1<<'\t'<<Efs*sqrt(vvec.at(p).at(j))+mvec.at(p).at(j);
                }
                ofs1<<'\n';



                for (int t=0;t<op().nt();t++) {
                    Eigen::VectorXd ks(N);
                    for (int k=0;k<N;k++) ks(k)=Kfn(theta,op().timei().at(t),x(k),0);
                    d_xx.at(p).at(j*op().nt()+t)=
                        ks.transpose()*invKy*ym+Mfn(m.at(j),op().timei().at(t),adjdiff);
                    d_xx.at(p).at(j*op().nt()+t)*=sqrt(vvec.at(p).at(j));
                    d_xx.at(p).at(j*op().nt()+t)+=mvec.at(p).at(j);
                }
            }
        }
    }
}


void Pre::impute_y(const char c)
{
    if (c!='M' and c!='H') throw runtime_error("impute_y");
    double** data = new double*[nprot()];
    int** mask = new int*[nprot()];
    for (int i = 0; i < nprot(); i++) {
        data[i] = new double[op().nc()];
        mask[i] = new int[op().nc()];
    }
    const vector<vector<vector<double> > >& yMH(c=='M'?yM():yH());
    vector<vector<vector<double> > >& d_yMH(c=='M'?d_yM:d_yH);

    vector<vector<vector<double> > > mvec(nprot());
    vector<vector<vector<double> > > vvec(nprot());
    for (int p=0;p<nprot();p++) {
        mvec.at(p).resize(yMH.at(p).size(),vector<double>(op().nrep()));
        vvec.at(p).resize(yMH.at(p).size(),vector<double>(op().nrep()));
        for (unsigned q=0;q<yMH.at(p).size();q++) for (int j=0;j<op().nrep();j++) {
            int nj=0;
            const vector<double>& yMHpq=yMH.at(p).at(q);
            for (int t=0;t<op().nt();t++) if (obs(yMHpq.at(j*op().nt()+t))) {
                mvec.at(p).at(q).at(j)+=yMHpq.at(j*op().nt()+t);
                nj++;
            }
            mvec.at(p).at(q).at(j) /= nj;
            for (int t=0;t<op().nt();t++) if (obs(yMHpq.at(j*op().nt()+t))) {
                vvec.at(p).at(q).at(j) += pow(yMHpq.at(j*op().nt()+t)-mvec.at(p).at(q).at(j),2);
            }
            vvec.at(p).at(q).at(j) /= nj-1;
            if (vvec.at(p).at(q).at(j)==0) vvec.at(p).at(q).at(j)=1;
            for (int t=0;t<op().nt();t++) if (obs(yMHpq.at(j*op().nt()+t))) {
                d_yMH.at(p).at(q).at(j*op().nt()+t) -= mvec.at(p).at(q).at(j);
                d_yMH.at(p).at(q).at(j*op().nt()+t) /= sqrt(vvec.at(p).at(q).at(j));
            }
        }
    }

    for (int p=0;p<nprot();p++) for (int l=0; l<op().nc(); l++) {
        data[p][l]=0;
        int nq=0;
        for (unsigned q=0;q<yMH.at(p).size();q++) if (obs(yMH.at(p).at(q).at(l))) {
            data[p][l] += yMH.at(p).at(q).at(l);
            nq++;
        }
        if (nq>0) {
            data[p][l] /= nq;
            mask[p][l] = 1;
        } else {
            data[p][l] = NAN;
            mask[p][l] = 0;
        }
    }


    vector<int> clusterid(nprot());

    multimap<int,int> c_mm;
    for (int p=0;p<nprot();p++) c_mm.insert(make_pair(clusterid[p],p));


    vector<vector<vector<Eigen::VectorXd> > > xax,yy_;
    xax.resize(yMH.size());
    yy_.resize(yMH.size());
    for (int p=0;p<nprot();p++) {
        xax.at(p).resize(yMH.at(p).size());
        yy_.at(p).resize(yMH.at(p).size());
        for (unsigned q=0;q<yMH.at(p).size();q++) {
            xax.at(p).at(q).resize(op().nrep());
            yy_.at(p).at(q).resize(op().nrep());
            for (int j=0;j<op().nrep();j++) {
                vector<double> x,y;
                for (int t=0;t<op().nt();t++) {
                    if (obs(yMH.at(p).at(q).at(j*op().nt()+t))) {
                        x.push_back(op().timei().at(t));
                        y.push_back(yMH.at(p).at(q).at(j*op().nt()+t));
                    }
                }
                xax.at(p).at(q).at(j).resize(x.size());
                yy_.at(p).at(q).at(j).resize(y.size());
                for (unsigned i=0;i<x.size();i++) {
                    xax.at(p).at(q).at(j)(i)=x.at(i);
                    yy_.at(p).at(q).at(j)(i)=y.at(i);
                }
            }
        }
    }

    const int Ns=100;
    const double step=double(op().timei().back()-op().timei().front())/Ns;
    Eigen::VectorXd xs(Ns);
    for (int i=0;i<Ns;i++) xs(i)=op().timei().at(0)+step/2+i*step;

    ofstream ofs1("Efs"+lexical_cast<string>(c)+".txt");
    ofs1<<"p\tq\tj";
    for (unsigned i=0;i<xs.size();i++) ofs1<<'\t'<<xs(i);
    ofs1<<'\n';


    int ijk=0;
    vector<double> adjdiff(op().timei().size());
    adjacent_difference(op().timei().begin(),op().timei().end(),adjdiff.begin());

    for (cmmii end1,itc=c_mm.begin();itc!=c_mm.end()/*,ijk<8*/;itc=end1,ijk++) {
        end1=c_mm.upper_bound(itc->first);

        vector<vector<double> > m(op().nrep(),vector<double>(op().nt()));

        vector<double> theta(3,1);

        for (cmmii it=itc;it!=end1;it++) {
            const int p=it->second;
            for (unsigned q=0;q<yy_.at(p).size();q++) for (int j=0;j<op().nrep();j++) {
                Eigen::VectorXd& y = yy_.at(p).at(q).at(j);
                Eigen::VectorXd& x = xax.at(p).at(q).at(j);
                const int N = y.size();
                //theta.at(0)=max(1,log10(exp(mvec.at(p).at(q).at(j))));
                Eigen::MatrixXd Ky(N,N);
                for (int i=0;i<N;i++) for (int k=0;k<N;k++) Ky(i,k)=Kfn(theta,x(i),x(k));
                Eigen::VectorXd ym(y);
                for (int i=0;i<N;i++) ym(i)-=Mfn(m.at(j),x(i),adjdiff);
                Eigen::MatrixXd invKy=Ky.inverse();
                ofs1<<p<<'\t'<<q<<'\t'<<j;
                for (unsigned i=0;i<xs.size();i++) {
                    Eigen::VectorXd ks(N);
                    for (int k=0;k<N;k++) ks(k)=Kfn(theta,xs(i),x(k),0);
                    const double Efs=ks.transpose()*invKy*ym+Mfn(m.at(j),xs(i),adjdiff);
                    ofs1<<'\t'<<Efs*sqrt(vvec.at(p).at(q).at(j))+mvec.at(p).at(q).at(j);
                }
                ofs1<<'\n';



                for (int t=0;t<op().nt();t++) {
                    Eigen::VectorXd ks(N);
                    for (int k=0;k<N;k++) ks(k)=Kfn(theta,op().timei().at(t),x(k),0);
                    d_yMH.at(p).at(q).at(j*op().nt()+t)=
                        ks.transpose()*invKy*ym+Mfn(m.at(j),op().timei().at(t),adjdiff);
                    d_yMH.at(p).at(q).at(j*op().nt()+t)*=sqrt(vvec.at(p).at(q).at(j));
                    d_yMH.at(p).at(q).at(j*op().nt()+t)+=mvec.at(p).at(q).at(j);
                }
            }
        }
    }
}
