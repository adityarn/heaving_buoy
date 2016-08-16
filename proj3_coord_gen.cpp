#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<time.h>
#include <vector>

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

    // xyz.coord[] holds x,y,z coordinates of control points of each panel
    xyz.resize(n_total);

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
    
    return 1;
  }

protected:
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
    double dtheta, dr, ndtheta, ndr;
    int c = n_body;
    
    ndr = sqrt(n_fs);
    ndtheta = ndr;

    dtheta = 2*pi/ndtheta;
    dr = (10*L) / ndr;
    
    for(int i= 0; i < ndtheta; i++)
      {
	for(int j=0; j<ndr; j++)
	  {
	    xyz[c].coord[0] = (0.5 * dr + j*dr) * cos(0.5*dtheta + i*dtheta);
	    xyz[c].coord[1] = (0.5 * dr + j*dr) * sin(0.5*dtheta + i*dtheta);
	    xyz[c].coord[2] = 0;

	    xyz[c].set_area( (0.5 * dr + j*dr) * dtheta , dr);

	    xyz[c].norm[0] = 0;
	    xyz[c].norm[1] = 0;
	    xyz[c].norm[1] = 1.0;

	    c++;
	  }
      }
    return 1;
  }

  int coord_set_bb()
  {
    double dtheta, dr, ndtheta, ndr;
    int c = n_body + n_fs;
    
    ndr = sqrt(n_fs);
    ndtheta = ndr;

    dtheta = 2*pi/ndtheta;
    dr = (10*L) / ndr;
    
    for(int i= 0; i < ndtheta; i++)
      {
	for(int j=0; j<ndr; j++)
	  {
	    xyz[c].coord[0] = (0.5 * dr + j*dr) * cos(0.5*dtheta + i*dtheta);
	    xyz[c].coord[1] = (0.5 * dr + j*dr) * sin(0.5*dtheta + i*dtheta);
	    xyz[c].coord[2] = -depth;

	    xyz[c].set_area( (0.5 * dr + j*dr) * dtheta , dr);
	    xyz[c].norm[0] = 0;
	    xyz[c].norm[1] = 0;
	    xyz[c].norm[1] = -1.0;

	    c++;
	  }
      }
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

  void write_coord()
  {
    ofstream coord;
    coord.open("xyz.txt");
    coord<<Time_period<<" "<<L<<" "<<n_total<<" "<<n_body<<" "<<n_fs<<" "<<n_bottom<<" "<<n_ff<<"\n";
    for(int i=0; i<n_total; i++)
      coord<<xyz[i].coord[0]<<" "<<xyz[i].coord[1]<<" "<<xyz[i].coord[2]<<" "<<xyz[i].norm[0]<<" "<<xyz[i].norm[1]<<" "<<xyz[i].norm[2]<<" "<<xyz[i].area<<"\n";
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
