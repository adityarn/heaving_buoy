#include<iostream>
#include<vector>
#include<complex>
using namespace std;


class Panel
{
public:
  double coord[3],norm[3],area;

  int set_area(double x, double y)
  {
    area = x*y;
    return 1;
  }

  /*  double rij(coordj[3])
  {
    return sqrt(pow(coord[0]-coordj[0],2)+pow(coord[1]-coordj[1],2)+pow(coord[2]-coordj[2],2));
    }*/
};

int main()
{
  vector <Panel> xyz;
  xyz.resize(3);
  xyz[0].set_area(2,2);
  cout<<"\n"<<xyz[0].area<<"\n";
}
