

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"
#include"Module.hpp"


void print_data(const Pre& pr,const string& str0)
{
    {
        ofstream ofs((str0+"X.txt").c_str());//output normalized data
        ofs<<pr.xh().at(0);
        for (int l=1; l<pr.op().nc()+1; l++) ofs<<'\t'<<pr.xh().at(l);
        ofs<<'\n';
        for (int g=0; g<pr.nprot(); g++) {
            ofs<<pr.xid().at(g);
            for (int l=0; l<pr.op().nc(); l++) {
                if (obs(pr.xx().at(g).at(l))) ofs<<'\t'<<pr.xx().at(g).at(l);
                else ofs<<"\tNA";
            }
            ofs<<'\n';
        }
    }

    if (pr.op().Mbool()) {
        ofstream ofs((str0+"M.txt").c_str());//output normalized data
        ofs<<pr.Mh().at(0);
        for (int l=1; l<pr.op().nc()+pr.op().level(); l++) ofs<<'\t'<<pr.Mh().at(l);
        ofs<<'\n';
        for (int p=0; p<pr.nprot(); p++) {
            for (unsigned q=0; q<pr.yM().at(p).size(); q++) {
                if (pr.op().level()==1) ofs<<pr.pid().at(p);
                else if (pr.op().level()==2) ofs<<pr.pid().at(p)<<'\t'<<pr.qid().at(p).at(q);
                for (int l=0; l<pr.op().nc(); l++) {
                    if (obs(pr.yM().at(p).at(q).at(l))) ofs<<'\t'<<pr.yM().at(p).at(q).at(l);
                    else ofs<<"\tNA";
                }
                ofs<<'\n';
            }
        }
    }
    {
        ofstream ofs((str0+"H.txt").c_str());//output normalized data
        ofs<<pr.Hh().at(0);
        for (int l=1; l<pr.op().nc()+pr.op().level(); l++) ofs<<'\t'<<pr.Hh().at(l);
        ofs<<'\n';
        for (int p=0; p<pr.nprot(); p++) {
            for (unsigned q=0; q<pr.yH().at(p).size(); q++) {
                if (pr.op().level()==1) ofs<<pr.pid().at(p);
                else if (pr.op().level()==2) ofs<<pr.pid().at(p)<<'\t'<<pr.qid().at(p).at(q);
                for (int l=0; l<pr.op().nc(); l++) {
                    if (obs(pr.yH().at(p).at(q).at(l))) ofs<<'\t'<<pr.yH().at(p).at(q).at(l);
                    else ofs<<"\tNA";
                }
                ofs<<'\n';
            }
        }
    }
}


void print_adj(const Module& mo,const string& str0)
{
    ofstream ofs("Adjacency_list.txt");
    for (unsigned i=0;i<mo.adj().size();i++) {
        copy(mo.adj().at(i).begin(),mo.adj().at(i).end(),ostream_iterator<int> (ofs," "));
        ofs<<'\n';
    }
}


int ncol(istream &istr)
{
    istr.clear();
    int pos=istr.tellg();
    istr.seekg(ios_base::beg);
    int n=0;
    string str;
    getline(istr,str);
    istringstream iss(str);
    while (iss>>str) n++;
    istr.clear();
    istr.seekg(pos);
    return n;
}


void posterior_mean(const int nMCsamples, const vector<string>& param_id)
{
    using boost::lexical_cast;
    for (unsigned j=0; j<param_id.size(); j++) {
        string buf;
        ifstream param_ifs(("s_"+param_id.at(j)).c_str());
        if (not param_ifs) throw runtime_error("can't open s_"+param_id.at(j));
        int n_cols = ncol(param_ifs);
        vector<double> param_mean(n_cols);
        for (int i=0; param_ifs>>buf; i++) {
            if (buf!="NA") param_mean.at(i%n_cols) += lexical_cast<double>(buf)/nMCsamples;
            else if (i<n_cols) param_mean.at(i) = NAN;
        }
        ofstream param_mean_ofs(("mean_"+param_id.at(j)).c_str());
        if (not param_mean_ofs) throw runtime_error("can't open mean_"+param_id.at(j));
        for (int i=0; i<n_cols; i++) {
            if (obs(param_mean.at(i))) param_mean_ofs<<param_mean.at(i)<<'\t';
            else param_mean_ofs<<"NA\t";
        }
    }
}
