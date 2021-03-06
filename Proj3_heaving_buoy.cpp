#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<time.h>
#include<complex>
#include <vector>
#include"Gauss_Elimin.h"

using namespace std;

//Class with all panel properties (x,y,z,nx,ny,nz,area)
class Panel
{
public:
  double coord[3],norm[3],area;

  int set_area(double x, double y)
  {
    area = x*y;
    return 1;
  }
  
  double rij(double coordj[3])
  {
    return sqrt(pow(coord[0]-coordj[0],2)+pow(coord[1]-coordj[1],2)+pow(coord[2]-coordj[2],2));
  }
};

class Green
{

public:

  // Constructor to initialise values
  // Time_period is that of incoming wave
  Green(double Time_Period, double Depth, double Draft, double Dia, int Nbody, int Nf, int Nb, int Nff)
  {

    Time_period = Time_Period;

    // n_body is the number of panels on the surface of the body
    // Ensure that Ns is 2 x (number with exact square root)
    n_body = Nbody;

    // n_fs is the number of panels on the free surface
    n_fs = Nf;

    // n_bottom is the number of panels on the bottom of the domain (flat sea bed)
    n_bottom = Nb;

    // n_on_ff is the number of far field panels
    n_ff = Nff;

    n_total = Nbody + Nf + Nb + Nff;

    // a is a 2 dimensional coefficient matrix a[N][N] x phi[N] = b[N]
    a.resize(n_total);

    // phi is the unknown values of phi that we wish to obtain
    phi.resize(n_total);

    b.resize(n_total);

    // xyz.coord[] holds x,y,z coordinates of control points of each panel

    for(int i=0; i<n_total; i++)
      a[i].resize(n_total);

    draft = Draft;
    pi = 4*atan(1.0);
    rad = Dia*0.5;

    sigma = 2*pi/Time_period;

    depth = Depth;
  }


  int Solve()
  {
    clock_t t1;
    t1 = clock();

    iterate_L();
    coord_set_body();
    coord_set_fs();
    coord_set_bb();
    coord_set_ff();
    
    write_coord();
    cout<<"\nAll body coords set!\n";

    //    a_b_compute();
    //    write(t1);
    
    return 1;
  }

protected:
  vector < vector < complex <double> > > a;
  vector < complex <double> > phi, b;
  double pi, rad, draft, domain_rad, sigma, L, depth, Time_period;
  int n_body, n_fs, n_bottom, n_ff, n_total;
  vector <Panel> xyz; 


  int iterate_L()
  {
    double Lo, tL;

    Lo = 1.56 * pow(Time_period,2);
    L = Lo;
    tL = 0;

    while( abs(L - tL) >= 1e-8)
      {
	tL = L;
	L = Lo * tanh(2*pi/tL * depth);
      }

    cout<<"\nLo = "<<Lo<<"\nL = "<<L<<"\n";
    return 1;
  }

  // Function which sets the control points
  // coordinate to x[], y[], z[]
  int coord_set_body()
  {
    double dy, dtheta, dr, denom;
    int ndtheta, ndy, ndr, c;

    c = 0;
    ndy = int(pow(n_body*.5,0.5));
    ndtheta = ndy;
    dtheta = 2*pi/ndtheta;
    dy = draft/ndy;

    for(int i=0; i<ndy; i++)
      {
	for(int j=0; j<ndtheta; j++)
	  {
	    xyz.resize(c+1);
	    xyz[c].coord[0] = rad * cos(0.5*dtheta + j*dtheta);
	    xyz[c].coord[1] = rad * sin(0.5*dtheta + j*dtheta);
	    xyz[c].coord[2] = - (dy * 0.5 + i * dy);
	    xyz[c].set_area(rad * dtheta, dy);
	    denom = sqrt(4*pow(xyz[c].coord[0],2) + 4 * pow(xyz[c].coord[1],2));
	    xyz[c].norm[0] = - 2 * xyz[c].coord[0] / denom;
	    xyz[c].norm[1] = - 2 * xyz[c].coord[1] / denom;
	    xyz[c].norm[2] = 0.0;
	    c++;
	  }
      }
    
    if(c != ndy * ndtheta)
      {
	cout<<"Check body discretization! c .ne. (ndy x ndtheta), c = "<<c<<"\n";
	//	return 0;
      }

    // Bottom (of body) discretization
    ndr = ndy;
    dr = rad / ndr;
    cout<<"\ndr = "<<dr<<"\nndr = "<<ndr<<"\n";
    for(int i=0; i<ndr; i++)
      {
	for(int j=0; j<ndtheta; j++)
	  {
	    xyz.resize(c+1);
	    xyz[c].coord[0] = (0.5*dr + i*dr) * cos(0.5*dtheta + j*dtheta);
	    xyz[c].coord[1] = (0.5*dr + i*dr) * sin(0.5*dtheta + j*dtheta);
	    xyz[c].coord[2] = -draft;
	    xyz[c].set_area((0.5*dr + i*dr)*dtheta,dr);
	    xyz[c].norm[0] = 0.0;
	    xyz[c].norm[1] = 0.0;
	    xyz[c].norm[2] = 1.0;
	    c++;
	  }
      }
    
    return 1;

  }

  int coord_set_fs()
  {
    double dtheta, dr, ndtheta,circum, arc_l, ndr, circum_t = 0;
    int c = n_body;
    
    ndr = sqrt(n_fs);
    dr = (10*L) / ndr;

    for(int i=0; i<int(ndr); i++)
      circum_t += 2*pi* (dr+i*dr);

    arc_l = circum_t/n_fs;
    
    for(int i= 0; i < ndr; i++)
      {
	circum = 2*pi * (i+1)*dr;
	ndtheta = circum/arc_l;
	dtheta = 2*pi/ndtheta;

	for(int j=0; j<ndtheta; j++)
	  {
	    xyz.resize(c+1);
	    xyz[c].coord[0] = (0.5 * dr + i*dr) * cos(0.5*dtheta + j*dtheta);
	    xyz[c].coord[1] = (0.5 * dr + i*dr) * sin(0.5*dtheta + j*dtheta);
	    xyz[c].coord[2] = 0;

	    xyz[c].set_area( (0.5 * dr + i*dr) * dtheta , dr);

	    xyz[c].norm[0] = 0;
	    xyz[c].norm[1] = 0;
	    xyz[c].norm[1] = 1.0;

	    c++;
	  }
      }
    
    n_fs = c - n_body;

    return 1;
  }

  int coord_set_bb()
  {
    double dtheta, dr, ndtheta, ndr, circum, arc_l, circum_t = 0;
    int c = n_body + n_fs;
    
    ndr = sqrt(n_fs);
    dr = (10*L) / ndr;

    for(int i=0; i<int(ndr); i++)
      circum_t += 2*pi* (dr+i*dr);

    arc_l = circum_t/n_fs;
    
    for(int i= 0; i < ndr; i++)
      {
	circum = 2*pi * (i+1)*dr;
	ndtheta = circum/arc_l;
	dtheta = 2*pi/ndtheta;

	for(int j=0; j<ndtheta; j++)
	  {
	    xyz.resize(c+1);
	    xyz[c].coord[0] = (0.5 * dr + i*dr) * cos(0.5*dtheta + j*dtheta);
	    xyz[c].coord[1] = (0.5 * dr + i*dr) * sin(0.5*dtheta + j*dtheta);
	    xyz[c].coord[2] = 0;

	    xyz[c].set_area( (0.5 * dr + i*dr) * dtheta , dr);

	    xyz[c].norm[0] = 0;
	    xyz[c].norm[1] = 0;
	    xyz[c].norm[1] = 1.0;

	    c++;
	  }
      }
    
    n_bottom = c - (n_body + n_fs);

    return 1;
  }


  int coord_set_ff()
  {
    double dtheta, dy, ndtheta, ndy, denom;
    int c = n_body + n_fs + n_bottom;

    ndy = int(pow(n_ff , 0.5));
    ndtheta = ndy;
    dtheta = 2*pi/ndtheta;
    dy = depth/ndy;
    
    for(int i=0; i<ndy; i++)
      {
	for(int j=0; j<ndtheta; j++)
	  {
	    xyz.resize(c+1);
	    xyz[c].coord[0] = 10 * L * cos(0.5*dtheta + j*dtheta);
	    xyz[c].coord[1] = 10 * L * sin(0.5*dtheta + j*dtheta);
	    xyz[c].coord[2] = - (dy * 0.5 + i * dy);

	    xyz[c].set_area(L/10 * dtheta, dy);
	    denom = sqrt(4*pow(xyz[c].coord[0],2) + 4 * pow(xyz[c].coord[1],2));

	    xyz[c].norm[0] = 2 * xyz[c].coord[0] / denom;
	    xyz[c].norm[1] = 2 * xyz[c].coord[1] / denom;
	    xyz[c].norm[2] = 0.0;
	    c++;
	  }
      }
    
    if( (c - n_body - n_fs - n_bottom) != ndy * ndtheta - 1)
      {
	cout<<"Check far field discretization! c .ne. (ndy x ndtheta - 1), c = "<<c<<"\n";
	//	return 0;
      }

    return 1;
  }

  // Function which computes a[] and b[]
  // for each mode and calls gauss_Elimin()
  // which in turn calls back_sub() to compute phi
  // for each mode (from 1 to 6)
  int a_b_compute()
  {
    int i,j,k,f;
    double denom,numer,numer2,denom2;
    double r, bsum;

    cout<<"Starting a_b_compute()\n";

    f=0;

    for(j=0; j < n_total; j++)
      {
	bsum = 0;

	for(k=0; k < n_body; k++)
	  {
	    if(j != k)
	      {
		r = pow(xyz[j].coord[0]- xyz[k].coord[0],2) + pow(xyz[j].coord[1] - xyz[k].coord[1],2) + pow(xyz[j].coord[2] - xyz[k].coord[2],2);
		denom = sqrt(r);
		numer = -sigma * xyz[k].norm[2] * xyz[k].area; 
		
		bsum += numer/denom;
	      }
	  }
	b[j] = complex<double>(0,bsum);

	for(k=0; k < n_body; k++)
	  {

	    if(j!= k)
	      {
		r = pow(xyz[j].coord[0]- xyz[k].coord[0],2) + pow(xyz[j].coord[1] - xyz[k].coord[1],2) + pow(xyz[j].coord[2] - xyz[k].coord[2],2);
		denom = pow(r,1.5);

		numer = ( (xyz[j].coord[0]-xyz[k].coord[0])*xyz[k].norm[0] + (xyz[j].coord[1]-xyz[k].coord[1])*xyz[k].norm[1] + (xyz[j].coord[2]-xyz[k].coord[2])*xyz[k].norm[2] ) * xyz[k].area;

		a[j][k] = complex<double>(numer/denom,0);
	      }

	    if(j == k)
	      a[j][k] = complex<double>(2*pi,0);		// All diagonal terms of a matrix are 2 x pi
	  }
	  

	for(k= n_body; k < n_body + n_fs; k++)
	  {
	    if(j!= k)
	      {
		r = pow(xyz[j].coord[0]- xyz[k].coord[0],2) + pow(xyz[j].coord[1] - xyz[k].coord[1],2) + pow(xyz[j].coord[2] - xyz[k].coord[2],2);
		denom = pow(r,1.5);

		numer = ( (xyz[j].coord[0]-xyz[k].coord[0])*xyz[k].norm[0] + (xyz[j].coord[1]-xyz[k].coord[1])*xyz[k].norm[1] + (xyz[j].coord[2]-xyz[k].coord[2])*xyz[k].norm[2] );

		denom2 = 9.8 * pow(r,0.5);

		numer2 = pow(sigma,2);

		a[j][j] = complex<double>( (numer/denom - numer2/denom2) * xyz[k].area, 0);
	      }

		if(j == k)
		  a[j][j] = complex<double>(2*pi,0);		// All diagonal terms of 'a' matrix are 2 x pi
	  }
	  
	for(k= n_body + n_fs; k < n_body + n_fs + n_bottom; k++)
	  {
	    if(j!= k)
	      {
		r = pow(xyz[j].coord[0]- xyz[k].coord[0],2) + pow(xyz[j].coord[1] - xyz[k].coord[1],2) + pow(xyz[j].coord[2] - xyz[k].coord[2],2);
		denom = pow(r,1.5);

		numer = ( (xyz[j].coord[0]-xyz[k].coord[0])*xyz[k].norm[0] + (xyz[j].coord[1]-xyz[k].coord[1])*xyz[k].norm[1] + (xyz[j].coord[2]-xyz[k].coord[2])*xyz[k].norm[2] );

		a[j][k] = complex<double>( (numer/denom) * xyz[k].area, 0);
	      }

		if(j == k)
		  a[j][j] = complex<double>(2*pi,0);		// All diagonal terms of 'a' matrix are 2 x pi
	  }

	for(k= n_body + n_fs + n_bottom; k< n_total; k++)
	  {
	    if(j!= k)
	      {
		r = pow(xyz[j].coord[0]- xyz[k].coord[0],2) + pow(xyz[j].coord[1] - xyz[k].coord[1],2) + pow(xyz[j].coord[2] - xyz[k].coord[2],2);
		denom = pow(r,1.5);

		numer = ( (xyz[j].coord[0]-xyz[k].coord[0])*xyz[k].norm[0] + (xyz[j].coord[1]-xyz[k].coord[1])*xyz[k].norm[1] + (xyz[j].coord[2]-xyz[k].coord[2])*xyz[k].norm[2] ) * xyz[k].area;

		denom2 = pow(r,.5);

		numer2 = -2*pi/L * xyz[k].area;

		a[j][k] = complex<double>(numer/denom , numer2/denom2);
	      }

		if(j == k)
		  a[j][j] = complex<double>(2*pi, 0);		// All diagonal terms of 'a' matrix are 2 x pi
	  }
	
	cout<<"\nCoefficient matrix computed!\n";
	write_ab();
	Gauss_Elimin Elim(a,phi,b);
	f = Elim.gauss_Elimin();

	if(f != 0)
	  return 1;


      }

    return 0;

  }

  double added_mass_coeff(int cmp)
  {
    double mu=0.0;

    for(int i=0; i < n_body; i++)
      {
	if(cmp == 1)
	  mu = mu + imag(phi[i])*xyz[i].norm[0]*xyz[i].area;
	else
	  mu = mu + real(phi[i])*xyz[i].norm[0]*xyz[i].area;
      }

    return mu;
  }

  void write_coord()
  {
    ofstream coord;
    coord.open("xyz.txt");
    coord<<Time_period<<" "<<L<<" "<<n_total<<" "<<n_body<<" "<<n_fs<<" "<<n_bottom<<" "<<n_ff<<"\n";
    for(int i=0; i<n_total; i++)
      coord<<xyz[i].coord[0]<<" "<<xyz[i].coord[1]<<" "<<xyz[i].coord[2]<<" "<<xyz[i].norm[0]<<" "<<xyz[i].norm[1]<<" "<<xyz[i].norm[2]<<" "<<xyz[i].area<<"\n";
  }

  void write_ab()
  {
    ofstream amatrix, bmatrix;
    
    amatrix.open("amatrix");
    bmatrix.open("bmatrix");
    for(int i=0; i<n_total; i++)
      {
	for(int j=0; j<n_total; j++)
	  amatrix<<a[i][j]<<" ";
	amatrix<<"\n";
	bmatrix<<b[i]<<"\n";
      }

  }

  void write(float t1)
  {
    ofstream phi_real, phi_imag, run_info;
    clock_t t2;

    phi_real.open("phi_real");
    phi_imag.open("phi_imag");
    run_info.open("run_info");

    for(int i=0;i<n_total;i++)
      {
	phi_real<<real(phi[i])<<"\n";
	phi_imag<<imag(phi[i])<<"\n";
      }

    t2 = clock();

    run_info<<"Green's boundary integration carried out for a heaving buoy in inviscid fluid with incident wave and reflected wave\n\n";

    run_info<<"n_total =>"<<n_total<<"\n\n";

    run_info<<"Wave Length = "<<L<<"\n\n";

    run_info<<"Mu (real) = "<<added_mass_coeff(0)<<"\n\n";
    run_info<<"Mu (imag) = "<<added_mass_coeff(1)<<"\n\n";

    run_info<<"\nTotal execution time = "<< float(t2) - t1<<" us";
  }


};

int main()
{

  int f=0;
  //  Green(double Time_period, double Depth, double Draft, double Dia, int Nbody, int Nf, int Nb, int Nff)
  // Nbody = 2 x (perfect square)
  // Nf = Nb = Nff = (perfect square)
  Green A(10.0, 30, 3.0, 1.0, 5000, 1681, 1681, 1681);

  A.Solve();

  return 1;
}
