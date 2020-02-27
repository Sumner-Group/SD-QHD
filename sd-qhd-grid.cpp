#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <new>
#include <time.h>
#include <algorithm>
//Calculates trajectory information (velocity, lambda, force, etc) on a cartesian grid
using namespace std;


double basiswfx (double cr, double cj, double x, double xr, int lxj, double y, double yr, int lyj, double z, double zr, int lzj, double hj)
{
  //Generates the Gaussian basis function
  double phi,xdum,ydum,zdum,r2;
  r2 = -hj*((x-xr)*(x-xr)+(y-yr)*(y-yr)+(z-zr)*(z-zr));
  phi =  cr*cj*exp(r2);
  return phi;
}

double dbasis(int lxj, double x, double xr,double hj)
{
  //First derivative of the Gaussian basis function  
  double fdx,dxdum,dphi,fdxdum;
  dxdum = pow((x-xr),(lxj-1));
  fdx = lxj*dxdum;
  fdxdum= pow((x-xr),(lxj+1));
  if (lxj == 0)
  {
    dphi=-2*hj*fdxdum; 
  }
  else
  {
    dphi=fdx-2*hj*fdxdum;
  }
  return dphi;
}

double d2basis (int lxj,double hj,double x, double xr)
{
  //Second derivative of the Gaussian basis function
  double d2xdum,d2dxdum,d2phi;

    d2xdum= pow((x-xr),(lxj-2));
  d2dxdum = pow((x-xr),(lxj));
       
  if (lxj == 0)
  {
    d2phi=-2.0*hj;
  }
   else if (lxj==1)
  {
    d2phi=-2*hj*(lxj+1)*d2dxdum;
  }
   else
  {
	d2phi= lxj*(lxj-1)*d2xdum-2*hj*(lxj+1)*d2dxdum;
  }
    return d2phi;
}

double lambda (double m, double fx,double fy,double fz, double fxx, double fyy, double fxy, double vx, double vy)
{
  //Lagrange multiplier
  double lamd;
  lamd = -m*((fxx*vx*vx)+(2*fxy*vx*vy)+ (fyy*vy*vy))/ (fx*fx + fy*fy);
  return lamd;
}


double position (double x, double v, double dt, double a0)
{
  //Updates the position
  double pos; 
  pos = x + v*dt + 0.5*(a0*dt*dt);
 return pos; 
}


double energy (double m, double vx,double vy)
{
  //Calculates the kinetic energy
 double ene;
 ene =  (0.5*m*((vx)*(vx)+(vy)*(vy)));
 return ene;
}

int main ()

{
  double m, vx,vy, a0x=0, a0y=0, dt, dt0,x, y, e, z=0, fxx, fxy, fyy, fx, fy,fz,v1x,v1y, a1y, a1x,v0x,v0y,x1, y1, xhalf, yhalf, vxhalf, vyhalf, axhalf, ayhalf, k1vx,k1x,k1vy,k1y,k2vx,k2x,k3x,k4vx,k4x,k2vy,k2y,k3vy,k3y,k4vy,k4y,k4, k3vx, * cj,inpsi,cr[1],* xr, * yr, * zr,* hj, xi,yi, xf,yf,d,vxi, vyi,dot,vmag,dist,dist0,xrange,yrange,zrange,randnum1,randnum2,randnum3,randnum4,rw,maxval,wrange,phi,xdum,ydum,zdum,dphix,dphiy,phix,dsphi,phiy,d2phi,phixy;
  double maxx,minx,maxy,miny,maxz,minz,dx,dy,dz,x0,y0,z0,lam;
  int n=0, j,k,* lxj,* lyj,* lzj,* nbasis ,natom,ntotbasis,l,* ncount,ndum,numtraj,n1,ifrev,i,dummy, xcount,ycount,zcount,ngridx,ngridy,ngridz,ifgrid;

  //Read input file
  cin >> m;

  cin >> ngridx;
  cin >> ngridy;
  cin >> ngridz;
  
  numtraj=(ngridx+1)*(ngridy+1)*(ngridz+1);

  cin >> xrange;

  cin >> yrange;

  cin >> zrange;

  cin >> natom;

  xr= new double[natom];
  yr= new double[natom];
  zr= new double[natom];
  nbasis= new int[natom];

  ncount= new int[natom];
  ntotbasis=0;
  for (k=0;k<=natom-1; k++)
    {
      cin >>nbasis[k];
      ntotbasis+=nbasis[k];
      
    }

  //Read in coefficient, angular momentum and exponent of Gaussian basis function
  cj = new double[ntotbasis];
  lxj = new int[ntotbasis];
  lyj = new int[ntotbasis];
  lzj = new int[ntotbasis];
  hj = new double[ntotbasis];
 
  cin >>cr[0]; 
  
  for (k=0;k<ntotbasis; k++)
    {
      cin >>cj[k];
      
    }
  
  
  for (k=0;k<ntotbasis; k++)
    {
      cin >>lxj[k];
      
    }
  
  
  for (k=0;k<ntotbasis; k++)
    {
      cin >>lyj[k];
  
    }
  

  for (k=0;k<ntotbasis; k++)
    {
      cin >>lzj[k];
 
    }
  
  for (k=0;k<=natom-1; k++)
    {
      cin >> xr[k];
      
      cin >> yr[k];
      
      cin >> zr[k];  
    }
  
  for (k=0; k<=ntotbasis-1; k++)
    {
      cin >> hj[k];
      
    }
  maxx = *max_element(xr,xr+natom);
  minx = *min_element(xr,xr+natom);
  maxy = *max_element(yr,yr+natom);
  miny = *min_element(yr,yr+natom);
  maxz = *max_element(zr,zr+natom);
  minz = *min_element(zr,zr+natom);
  x0=minx-xrange;
  y0=miny-yrange;
  z0=minz-zrange;
  dx=(maxx+xrange-x0)/ngridx;
  dy=(maxy+yrange-y0)/ngridy;
  dz=(maxz+zrange-z0)/ngridz;
  xcount=0;
  ycount=0;
  zcount=0;
  //Loop over grid points
  for (n1=0; n1<numtraj; n1++)
    {
      x=x0+dx*(xcount);
      y=y0+dy*(ycount);
      z=z0+dz*(zcount);
      xcount++;
      if(xcount == ngridx+1) 
	{
	  xcount=0;
	  ycount++;
	}
    if(ycount == ngridy+1)
      {
	ycount=0;
	zcount++;
      }
    inpsi=0.0;
    dummy=0;
    fx=0.0;
    fy=0.0;
    fxx=0.0;
    fxy=0.0;
    fyy=0.0;
    for (i=0; i<natom; i++)
      {
     
	for (k=0; k<nbasis[i];k++)
	  {


	    phi = basiswfx(cr[0], cj[k+dummy], x, xr[i], lxj[k+dummy], y, yr[i], lyj[k+dummy], z, zr[i], lzj[k+dummy], hj[k+dummy]);

	    //wavefunction
	    xdum = pow((x-xr[i]),lxj[k+dummy]);
	    ydum = pow((y-yr[i]),lyj[k+dummy]);
	    zdum = pow((z-zr[i]),lzj[k+dummy]);
	    inpsi += phi*xdum*ydum*zdum;

	    //first derivative x
	    dphix = dbasis( lxj[k+dummy], x,  xr[i],hj[k+dummy]);
	    fx += phi*dphix*ydum*zdum;

	    //first derivative y
	    dphiy = dbasis( lyj[k+dummy], y,  yr[i],hj[k+dummy]);
	    fy += phi*dphiy*xdum*zdum;

	    //2nd derivative xx
	    phix = phi*ydum*zdum;
	    d2phi= d2basis(lxj[k+dummy], hj[k+dummy], x, xr[i]);
	    fxx+= phix*dphix*(-2.0*hj[k+dummy]*(x-xr[i]))+ phix*d2phi;

	    //2nd derivative yy
	    phiy = phi*xdum*zdum;
	    d2phi= d2basis(lyj[k+dummy], hj[k+dummy], y, yr[i]);
	    fyy+= phiy*dphiy*(-2.0*hj[k+dummy]*(y-yr[i]))+ phiy*d2phi;

	    //2nd derivative xy
	    phixy = phi*zdum;
	    fxy+= phixy*dphix*dphiy;

	  }
	dummy+=nbasis[i];
      }
    inpsi = 1.0/inpsi;

    //Spin-dependent velocity
    vx= fy*inpsi;
    vy= -fx*inpsi;
 
    if (abs(fx) < 1.0e-10 && abs(fy) < 1.0e-10)
      {
	printf ("stationary point \n");
      }
    //memset might be faster
    std:fill( ncount, ncount+natom, 0);
  
    vmag = sqrt(vx*vx+vy*vy);
    if(vmag >= 68.5)
      {
	printf ("vmag >50percent of c %8.6f %8.6f %8.6f %8.6f \n",vmag,x,y,z);
	ifrev = 1;
      }
    
    //kinetic energy
    e = energy (m, vx,vy);
    
    //Lagrange multiplier
    lam = lambda(m, fx,fy,fz, fxx, fyy, fxy, vx, vy);
    
    //print out information
    printf ("position: %8.6f %8.6f %8.6f velocity: %8.6f %8.6f energy %8.6f wfn %14.12f lambda %8.6f force %8.6f %8.6f\n", x, y, z, vx, vy, e, 1.0/inpsi,lam,lam*fx,lam*fy);
    ncount[ndum]++;
   
    }
  return 0;
}

