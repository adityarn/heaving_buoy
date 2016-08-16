#include<vector>
#include<iostream>
#include<cmath>
#include<complex>
using namespace std;

class xyz
{
public:
  double x,y,z;

};
int main()
{
  vector<xyz> ob;
  ob.resize(2);
  ob[0].x = 0;
  ob[1].x = 1;

  ob.resize(3);
  ob[2].x = 2;
  ob.resize(100);
  for(int i=0; i<3; i++)
    cout<<"\n"<<ob[i].x;

  cout<<"\n";
  return 0;

}
