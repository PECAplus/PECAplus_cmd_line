

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"


typedef map<string,Row*>::const_iterator mit;
typedef map<string,const Row*>::const_iterator cmit;
typedef map<string,const Row1*>::const_iterator cmit1;


double median(const vector<double>& data)
{
    vector<double> sorted_data(data);
    sort(sorted_data.begin(),sorted_data.end());
    const size_t lhs = (data.size() - 1) / 2 ;
    const size_t rhs = data.size() / 2 ;
    if (data.size() == 0) throw runtime_error("median");
    if (lhs == rhs) return sorted_data.at(lhs) ;
    else return (sorted_data.at(lhs) + sorted_data.at(rhs))/2.0;
}


size_t Pre::nl() const
{
    if (d_op.selbool()) return sdata().front().sel.size()/3;
    return d_cp.front().size();
}


size_t Pre::np() const
{
    if (d_op.selbool()) return sdata().size();
    return d_cp.size();
}


void Pre::read_sel()
{
    string str0;
    ifstream sel_ifs(op().SEL().c_str());
    if (not sel_ifs) throw runtime_error("can't open "+op().SEL());
    getline(sel_ifs,str0);
    istringstream iss(str0);
    for (int i=2;getline(sel_ifs,str0) and (not str0.empty());i++) {
        Row1 tmp_p;
        iss.clear(); iss.str(str0+"\t");
        getline(iss,tmp_p.pid,'\t');
        for (string str1;getline(iss,str1,'\t');) tmp_p.sel.push_back(boost::lexical_cast<int>(str1)>1);
        tmp_p.sel.pop_back(); tmp_p.sel.pop_back(); tmp_p.sel.pop_back();
        d_sdata.push_back(tmp_p);
    }
}


void Pre::set_sel()
{
    map<string,const Row1*> pid_m;
    for (list<Row1>::const_iterator it=sdata().begin();it!=sdata().end();it++) {
        pid_m[it->pid]=&(*it);
        d_pidset.insert(it->pid);
    }
    d_pidvec.assign(pidset().begin(),pidset().end()); 
    d_selup.resize(nl());
    d_seldown.resize(nl());
    d_selsig.resize(nl());
    for (cmit1 itp=pid_m.begin();itp!=pid_m.end();itp++) {
        for (unsigned l=0;l<nl();l++) {
            d_selup.at(l).push_back(itp->second->sel.at(l));
            d_seldown.at(l).push_back(itp->second->sel.at(nl()+l));
            d_selsig.at(l).push_back(itp->second->sel.at(2*nl()+l));
        }
    }
}


void Pre::read_cp()
{
    string str0;
    ifstream cps_ifs(op().CPS().c_str());
    if (not cps_ifs) throw runtime_error("can't open "+op().CPS());
    getline(cps_ifs,str0);
    istringstream iss(str0);
    int ci=0;
    for (string str1;getline(iss,str1,'\t') and str1!="R0";ci++);
    const int startr=ci;
    ci++;
    const string signedCP=(op().syndeg()?"signedCPS":"signedCPD");
    for (string str1;getline(iss,str1,'\t') and str1.substr(0,9)!=signedCP;ci++);
    const int start=ci;
    for (string str1;getline(iss,str1,'\t') and str1.substr(0,9)==signedCP;ci++);
    const int end=ci+1;
    for (int i=2;getline(cps_ifs,str0) and (not str0.empty());i++) {
        Row tmp_p;
        iss.clear(); iss.str(str0+"\t");
        getline(iss,tmp_p.pid,'\t');
        for (ci=0;ci<end and getline(iss,str0,'\t');ci++) {
            if (ci>=startr and ci<=startr+end-start) {
                tmp_p.R.push_back(boost::lexical_cast<double>(str0));
            } else if (ci>=start) {
                tmp_p.cp.push_back(abs(boost::lexical_cast<double>(str0)));
            }
        }
        d_cdata.push_back(tmp_p);
    }
}


void Pre::set_cp()
{
    map<string,const Row*> pid_m;
    for (list<Row>::const_iterator it=cdata().begin();it!=cdata().end();it++) {
        pid_m[it->pid]=&(*it);
        d_pidset.insert(it->pid);
    }
    d_pidvec.assign(pidset().begin(),pidset().end()); 
    for (cmit itp=pid_m.begin();itp!=pid_m.end();itp++) {
        d_cp.push_back(itp->second->cp);
        d_R.push_back(itp->second->R);
    }
}


void Pre::read_GO()
{
    string str0;
    ifstream go_ifs(op().GO_term().c_str());
    if (not go_ifs) throw runtime_error("can't open "+op().GO_term());
    getline(go_ifs,str0);
    istringstream iss(str0);
    for (int i=2;getline(go_ifs,str0) and (not str0.empty());i++) {

        Row0 tmp_p;
        iss.clear(); iss.str(str0+"\t");
        getline(iss,tmp_p.GO,'\t');
        getline(iss,str0,'\t');
        istringstream iss0(str0);
        for (;getline(iss0,str0,',');) tmp_p.mem.push_back(str0);
        getline(iss,tmp_p.GOf,'\t');
        d_GOdata.push_back(tmp_p);



    }
}


void Pre::set_GO()
{
    map<string,const Row0*> goid_m;
    int i=0;
    for (list<Row0>::const_iterator it=GOdata().begin();it!=GOdata().end();it++,i++) {
        goid_m[it->GO]=&(*it);
    }
    for (map<string,const Row0*>::const_iterator itp=goid_m.begin();itp!=goid_m.end();itp++) {
        const vector<string>& memf=itp->second->mem;
        set<int> tmpvec;
        for (vector<string>::const_iterator it=memf.begin();it!=memf.end();it++) {
            const unsigned j=find(pidvec().begin(),pidvec().end(),*it)-pidvec().begin();
            if (j!=pidvec().size()) tmpvec.insert(j);
        }
        if (tmpvec.size()>=op().minGene() and tmpvec.size()*100./memf.size()>=op().percent()) {
            d_GOid.push_back(itp->first);
            d_GOf.push_back(itp->second->GOf);
            d_mem.push_back(vector<int>(tmpvec.begin(),tmpvec.end()));
            d_mem_all.push_back(memf);
        }
    }
}


void Pre::read_net()
{
    ifstream module_ifs(op().module().c_str());
    if (not module_ifs) throw runtime_error("can't open "+op().module());
    d_adj.assign(pidvec().size(),vector<bool>(pidvec().size()));
    int i=0;
    for (string str0;getline(module_ifs,str0);i++) {
        istringstream iss(str0);
        istream_iterator<string> ii(iss),eos;
        vector<string> strv(ii,eos);
        for (int k=0;k<strv.size();k++) {
            const int j=boost::lexical_cast<int>(strv.at(k));
            d_adj.at(i).at(j)=1;
            d_adj.at(j).at(i)=1;
        }
    }
}


Pre::Pre(const Option& op) :
    d_op(op)
{
    if (d_op.selbool()) {
        read_sel();
        set_sel();
    } else {
        read_cp();
        set_cp();
    }
    read_GO();
    set_GO();
    if (d_op.modulebool()) {
        read_net();
    }
}
