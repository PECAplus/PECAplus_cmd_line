

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"
#include"Post.hpp"

struct Sigdat
{
    int g;
    vector<double> up,down,sig;
    Sigdat(const int g) : g(g) { }
};

void print_analysis(const Pre& pr,const Post& po)
{
    multimap<double,Sigdat,greater<double> > sig_mm;
    for (unsigned g=0;g<pr.ng();g++) {
        Sigdat tmpsig(g);
        for (unsigned l=0;l<pr.nl();l++) {
            tmpsig.up.push_back(-log10(po.pvup().at(l).at(g)));
            tmpsig.down.push_back(-log10(po.pvdown().at(l).at(g)));
            tmpsig.sig.push_back(-log10(po.pvsig().at(l).at(g)));
        }
        sig_mm.insert(make_pair(*max_element(tmpsig.up.begin(),tmpsig.up.end()),tmpsig));
    }
    ofstream ofs("Goterms.txt");
    if (not ofs) throw runtime_error("can't open Goterms.txt");
    ofs<<"MaxSig(Up)\tMaxSig(Down)\tMax(Both)\tGO_id\tGO_name\tGO_size\tGO_size_background";
    if (po.op().modulebool()) ofs<<"\tGO_EdgeCount";
    ofs<<"\tmembers";
    for (unsigned l=0;l<pr.nl();l++) ofs<<"\tUp("<<l+1<<")";
    for (unsigned l=0;l<pr.nl();l++) ofs<<"\tDown("<<l+1<<")";
    for (unsigned l=0;l<pr.nl();l++) ofs<<"\tSig("<<l+1<<")";
    ofs<<'\n';
    typedef multimap<double,Sigdat>::const_iterator sigit;
    for (sigit it=sig_mm.begin();/*it->first>0 and*/ it!=sig_mm.end();it++) {
        const int g=it->second.g;
        ofs<<*max_element(it->second.up.begin(),it->second.up.end())<<'\t';
        ofs<<*max_element(it->second.down.begin(),it->second.down.end())<<'\t';
        ofs<<*max_element(it->second.sig.begin(),it->second.sig.end())<<'\t';
        //ofs<<pr.GOid().at(g)<<'\t'<<pr.GOf().at(g)<<'\t'<<pr.mem().at(g).size()<<'\t';
        ofs<<pr.GOid().at(g)<<'\t'<<pr.GOf().at(g)<<'\t'<<pr.mem_all().at(g).size()<<'\t'<<pr.mem().at(g).size()<<'\t';
        if (po.op().modulebool()) {
            int edge_count=0;
            for (unsigned m=0;m<pr.mem().at(g).size()-1;m++) for (unsigned m1=m+1;m1<pr.mem().at(g).size();m1++) {
                edge_count+=pr.adj().at(pr.mem().at(g).at(m)).at(pr.mem().at(g).at(m1));
            }
            ofs<<edge_count<<"\t";
        }
        for (unsigned m=0;m<pr.mem().at(g).size();m++) {
            if (m>0) ofs<<' ';
            ofs<<pr.pidvec().at(pr.mem().at(g).at(m));
        }
        for (unsigned l=0;l<pr.nl();l++) ofs<<'\t'<<po.pvup().at(l).at(g);
        for (unsigned l=0;l<pr.nl();l++) ofs<<'\t'<<po.pvdown().at(l).at(g);
        for (unsigned l=0;l<pr.nl();l++) ofs<<'\t'<<po.pvsig().at(l).at(g);
        ofs<<'\n';
    }


}
