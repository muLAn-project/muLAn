%module esbltaylor
%{
#include <iostream>
#include <vector>

using namespace std;

extern vector< vector<double> > magnifcalc(vector< vector<double> >);
%}

%include "std_vector.i"
namespace std {
  %template(VecDouble) vector<double>;
  %template(VecVecdouble) vector< vector<double> >;
}

extern std::vector< std::vector<double> > magnifcalc (std::vector< std::vector<double> >);
