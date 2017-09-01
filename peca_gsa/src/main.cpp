

#include"main.hpp"
#include"Option.hpp"
#include"Pre.hpp"
#include"Post.hpp"


Option set_param(const string& filepath);

void print_analysis(const Pre& pr,const Post& po);


int main(int argc, char** argv)
{
    const Option& op=set_param(string(argv[ 1 ]));

    Pre pr(op);
    Post po(pr);

    print_analysis(pr,po);

    cout<<"Done!\n";

}
