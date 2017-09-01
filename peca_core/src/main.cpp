

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
    pr.impute_y();
    cout<<"done\n";
    print_data(pr,"impute");
    pr.norm_data();
    print_data(pr,"norm");

    Mcmc mc(pr);

    if (mo.modulebool()) {
        mo.module_data(pr.pid());
        print_adj(mo,"");
        Est es(mc,mo);
        Mcmc mc1(pr,es.varphi());

        //Mcmc mc1(pr,mo);
        //ofstream attr_ofs("attr.txt");
        //attr_ofs<<op.nrep()<<'\t'<<op.nt()<<'\t'<<mc1.a_tau()<<'\t'<<mc1.b_tau();
        //attr_ofs<<'\n';
    } else {
        //Est es(mc,mo);
        //Mcmc mc1(pr,es.varphi());
    }

    ofstream attr_ofs("attr.txt");
    attr_ofs<<op.nrep()<<'\t'<<op.nt()<<'\t'<<mc.a_tau()<<'\t'<<mc.b_tau();
    attr_ofs<<'\n';

    vector<string> param_id;
    param_id.push_back("xRyD");
    posterior_mean(op.ns(),param_id);
}
