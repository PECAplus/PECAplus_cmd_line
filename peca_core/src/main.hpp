
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<string>
#include<list>
#include<vector>
#include<algorithm>
#include<stdexcept>
#include<cfloat>

#include<iomanip>
#include<set>
#include<map>
#include<numeric>
#include<iterator>
#include<functional>
#include"boost/lexical_cast.hpp"
using namespace std;


inline bool obs(const double x) { return not std::isnan(x); }

template <class InputIterator>
double median(InputIterator first, InputIterator last)
{
    vector<double> sorted_data(first,last);
    sort(sorted_data.begin(),sorted_data.end());
    const size_t lhs = (sorted_data.size() - 1) / 2 ;
    const size_t rhs = sorted_data.size() / 2 ;
    if (first == last) throw runtime_error("median");
    if (lhs == rhs) return sorted_data.at(lhs) ;
    else return (sorted_data.at(lhs) + sorted_data.at(rhs))/2.0;
}
