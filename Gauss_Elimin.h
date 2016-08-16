#include<iostream>
#include<vector>
#include<complex>
using namespace std;


class Gauss_Elimin
{
 public:
  
 Gauss_Elimin(vector<vector < complex<double> > > A, vector<complex<double> > X, vector < complex<double> > B, int Copy_is_true): a(A),x(X),b(B)
  {
    copy_a = 1;
    dim = A.size();
  }

 Gauss_Elimin(vector<vector < complex<double> > > &A, vector<complex<double> > &X, vector <complex<double> > &B): a(A),x(X),b(B)
  {
    dim = A.size();
  }

  int gauss_Elimin()
  {
    complex<double> r;
    dim = a.size();
    cout<<"\ndim = a.size() = "<<dim<<"\n";
    write_ab();
    for(int i=0; i < dim - 1 ;i++)
      {
	if(abs(a[i][i]) <= 1.0e-8)
	  {
	    cout<<"\npivot called at i = "<<i<<"\nand a[i][i] ="<<a[i][i]<<"\n";
	    pivot_a(i);
	  }

	for(int j=i+1; j < dim; j++)
	  {
	    if(abs(a[i][i])!=0)
	      {
		r= a[j][i]/a[i][i];

		for(int k=i; k < dim; k++)
		  a[j][k] = a[j][k] - r * a[i][k];

		b[j] = b[j] - r*b[i];
	      }
	  }
      }

    back_sub();
    return 1;
  }

 protected:
  vector < vector < complex<double> > > &a;
  vector<complex <double> > &x, &b;
  int dim, copy_a;

  void write_ab()
  {
    ofstream amatrix, bmatrix;
    
    amatrix.open("Amatrix");
    bmatrix.open("Bmatrix");
    for(int i=0; i<dim; i++)
      {
	for(int j=0; j<dim; j++)
	  amatrix<<a[i][j]<<" ";
	amatrix<<"\n";
	bmatrix<<b[i]<<"\n";
      }

  }

  void back_sub()
  {
    complex<double> sum;
    int j,k;

    x.resize(dim);

    for(j = dim-1; j >= 0; j--)
      {
	sum = 0.0;
	if(j != dim-1)
	  for(k= dim-1; k>j; k--)
	    sum = sum + a[j][k] * x[k];
	x[j] = (b[j] - sum) / a[j][j];
      }
  }

  int pivot_a(int row1)
  {
    complex<double> a1;
    int i,row2;

    row2=row1+1;

    while(abs(a[row2][row1]) <= 1.0e-8)
      {
	if(row2 <= dim)
	  row2+=1;

	else
	  return 0;
      }
    for(i=0; i<dim; i++)
      {
	a1 = a[row1][i];
	a[row1][i] = a[row2][i];

	a[row2][i] = a1;

	a1 = b[row1];
	b[row1] = b[row2];
	b[row2] = a1;
      }
    return 1;
  }
};
