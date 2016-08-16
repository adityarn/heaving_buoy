#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
using namespace std;

int main()
{
  vector < vector < complex<double> > > array;
  int n = 10;
  array.resize(n);
  for(int i=0; i<n; i++)
    array[i].resize(n);
  for(int j=0; j<n; j++)
    for(int k=0; k<n; k++)
      array[j][k] = complex<double> (j,k);
  for(int j=0; j<n; j++)
    {
    for(int k=0; k<n; k++)
      cout<<array[j][k]<<" ";
    cout<<"\n";
    }
  
  array[3][3] = complex<double>(1,10);
  double flau = 3.3;
  cout<<"\n"<<imag(array[3][3])<<"\n";
  ofstream jhakaas;
  jhakaas.open("cmplx.txt");
  jhakaas<<array[3][3] + flau<<" "<<array[4][4];
}
