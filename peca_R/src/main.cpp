

#include"main.hpp"
#include"Module.hpp"
#include"Option.hpp"
#include"Pre.hpp"
#include"Mcmc.hpp"
#include"Est.hpp"//


Option set_param(const string& filepath,Module& mo);

void print_data(const Pre& pre,const string& str0);
void print_adj(const Module& mo,const string& str0);
void posterior_mean(const int nMCsamples, const vector<string>& param_id);

int main(int /*argc*/, char** argv)
{
    Module mo;
    const Option& op=set_param(string(argv[ 1 ]),mo);
    Pre pr(op);

    print_data(pr,"");
    cout<<"impute... "<<flush;
    pr.impute_x();
    if (op.Mbool()) pr.impute_y('M');
    pr.impute_y('H');
    cout<<"done\n";
    //print_data(pr,"impute");
    pr.norm_data();
    print_data(pr,"norm");

    Mcmc mc(pr);

    if (mo.modulebool()) {
        mo.module_data(pr.pid());
        print_adj(mo,"");
        Est es(mc,mo);
        Mcmc mc1(pr,es.varphiR(),es.varphiD());
    }

    ofstream attr_ofs("attr.txt");
    attr_ofs<<op.nrep()<<'\t'<<op.nt()<<'\t'<<mc.a_tauH()<<'\t'<<mc.b_tauH()<<'\t'<<mc.a_tauM()<<'\t'<<mc.b_tauM();
    attr_ofs<<'\n';
    copy(op.timei().begin(),op.timei().end(),ostream_iterator<double>(attr_ofs," "));

    vector<string> param_id;
    if (op.Mbool()) param_id.push_back("etaM");
    param_id.push_back("etaH");
    posterior_mean(op.ns(),param_id);
}
