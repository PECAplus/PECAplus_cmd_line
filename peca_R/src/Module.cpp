

#include"main.hpp"
#include"Module.hpp"


void Module::setModulebool(const string& str0)
{
    istringstream iss(str0);
    string str1;
    iss>>str1;
    d_modulebool=bool(ifstream(str1.c_str()));
}


void Module::setModule_freq(const string& str0)
{
    istringstream iss(str0);
    if (not (iss>>d_module_freq)) throw runtime_error("MODULE_FREQ?");
}


void Module::setModule_size(const string& str0)
{
    istringstream iss(str0);
    if (not (iss>>d_min_size)) throw runtime_error("MODULE_SIZE?");
    int int0;
    if (iss>>int0) d_max_size=int0;
}


void Module::setModule(const string& str0)
{
    istringstream iss(str0);
    iss>>d_module;
}





void Module::module_data(const vector<string>& pidvec)
{
    ifstream module_ifs(module().c_str());
    if (not module_ifs) throw runtime_error("can't open "+module());
    {
        d_adj.resize(pidvec.size());
        for (string str0,str1;getline(module_ifs,str0);) {
            istringstream iss(str0);
            if (getline(iss,str0,'\t') and getline(iss,str1,'\t') and str0!=str1) {
                const unsigned i=find(pidvec.begin(),pidvec.end(),str0)-pidvec.begin();
                const unsigned j=find(pidvec.begin(),pidvec.end(),str1)-pidvec.begin();
                if (i!=pidvec.size() and j!=pidvec.size()) {//if both found
                    d_adj.at(j).insert(i);
                    d_adj.at(i).insert(j);
                }
            }
        }
    }


    int sum=0;
    for (unsigned p=0;p<pidvec.size();p++) sum += adj().at(p).size();
    cout<<sum/2.<<" edges\n";
}
