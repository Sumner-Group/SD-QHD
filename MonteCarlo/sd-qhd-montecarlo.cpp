#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <new>
#include <time.h>
//Randomly places starting points. Uses acceptance rejection method to get nonuniform distribution of points.
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
  double d2xdum,d2dxdum,d2phi;
  //Second derivative of the Gaussian basis function
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

double distance (double xi, double xf, double yi, double yf,double zi,double zf)
{ 
  //Calculates the distance
  double d,dummy1,dummy2,dummy3,dummy4;
  dummy1= (xf-xi)*(xf-xi);
  dummy2=(yf-yi)*(yf-yi);
  dummy3= (zf-zi)*(zf-zi);
  dummy4=dummy1+dummy2+dummy3;
  d=sqrt(dummy4);
  return d;
}

double dotproduct (double vxi,double vx, double vyi, double vy)
{
  //Uses dot product to calculate cos(theta), theta is angle between vectors
  double dot,dum1,dum2,dum3,dum4,dum5,dum6;
  dum1= (vxi*vx)+(vyi*vy);
  dum2= (vxi*vxi)+(vyi*vyi);
  dum3= sqrt(dum2);
  dum4= (vx*vx)+(vy*vy);
  dum5= sqrt(dum4);
  dum6= dum3*dum5;
  dot= dum1/dum6;
  return dot;
}

double cross (double vxi,double vx, double vyi, double vy)
{
  //Uses cross product to calculate sin(theta), theta is angle between vectors
  double cross,dum1,dum2,dum3,dum4,dum5,dum6;
  dum1= (vxi*vy)-(vyi*vx);
  dum2= (vxi*vxi)+(vyi*vyi);
  dum3= sqrt(dum2);
  dum4= (vx*vx)+(vy*vy);
  dum5= sqrt(dum4);
  dum6= dum3*dum5;
  cross= dum1/dum6;
  return cross;
}

int main ()

{
  double m, vx,vy, a0x=0, a0y=0, dt, dt0,x, y, e, z=0, fxx, fxy, fyy, fx, fy,fz,v1x,v1y, a1y, a1x,v0x,v0y,x1, y1, xhalf, yhalf, vxhalf, vyhalf, axhalf, ayhalf, k1vx,k1x,k1vy,k1y,k2vx,k2x,k3x,k4vx,k4x,k2vy,k2y,k3vy,k3y,k4vy,k4y,k4, k3vx, * cj,inpsi,cr[1],* xr, * yr, * zr,* hj, xi,yi, xf,yf,d,vxi, vyi,dot,vmag,dist,dist0,xrange,yrange,zrange,randnum1,randnum2,randnum3,randnum4,rw,maxval,wrange,phi,xdum,ydum,zdum,dphix,dphiy,phix,dsphi,phiy,d2phi,phixy,* posx, * posy, * posz, rho0, dx, dy,dz,w, rhotot,xave,yave,zave,* orbit,* orbit2, d0,xmin,xmax,ymin,ymax,zmin,zmax,lam;
  int n=0, j,k,* lxj,* lyj,* lzj,* nbasis ,natom,ntotbasis,l,* ncount,ndum,numtraj,n1,ifrev,i,dummy,* orbits;
  //Read input file
  cin >> m;

  cin >> dt;
  dt0 = dt;

  cin >> numtraj;

  cin >> xrange;
  cin >> xmax;
  cin >> xmin;

  cin >> yrange;
  cin >> ymax;
  cin >> ymin;

  cin >> zrange;
  cin >> zmax;
  cin >> zmin;

  cin >> wrange;

  cin >> j;

  cin >> natom;

  xr= new double[natom];
  yr= new double[natom];
  zr= new double[natom];
  nbasis= new int[natom];

  ncount= new int[natom];
  orbit = new double[natom*j];
  orbit2 = new double[natom*j];
  orbits = new int[natom];
  ntotbasis=0;
  for (k=0;k<=natom-1; k++)
    {
      cin >>nbasis[k];
      ntotbasis+=nbasis[k];
      orbits[k]=0;
    }
  //Read in coefficient, angular momentum and width of Gaussian basis function
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
  
  //Set seed for random number generator
  //for reproducability, use srand(1)
  srand(time(NULL));
  //srand (1);
  //Loop over starting points
  for (n1=0; n1<numtraj; n1++)
    {
      
      //Calculate random number between the max and min in the x, y and z-directions
      x= xmin-xrange+double(rand() %1000000)/1000000.0*(xmax+2.0*xrange-xmin);

      y= ymin-yrange+double(rand() %1000000)/1000000.0*(ymax+2.0*yrange-ymin);

      z= zmin-zrange+double(rand() %1000000)/1000000.0*(zmax+2.0*zrange-zmin);
      xi=x;
      yi=y;

      //Calculate random number for w function used in acceptance/rejection method
      randnum4= double(rand() %1000000)/1000000.0;
      rw= double(rand() %1000000)/1000000.0*wrange;
      //Calculate initial conditions
      inpsi=0.0;
      dummy=0;
      fx=0.0;
      fy=0.0;
      fxx=0.0;
      fxy=0.0;
      fyy=0.0;      
      //wavefunction
      for (i=0; i<natom; i++)
	{
	  
	  for (k=0; k<nbasis[i];k++)
	    {

	      phi = basiswfx(cr[0], cj[k+dummy], x, xr[i], lxj[k+dummy], y, yr[i], lyj[k+dummy], z, zr[i], lzj[k+dummy], hj[k+dummy]);
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
      //Compare random number to probability density, rho. Reject starting point if random number greater than rho.
      if (rw > (inpsi*inpsi))
	{
	  n1=n1-1;
	  continue;
	}
      inpsi = 1.0/inpsi;

      vx= fy*inpsi;
      vy= -fx*inpsi;
 
      vxi = vx;
      vyi = vy;
      
      if (abs(fx) < 1.0e-10 && abs(fy) < 1.0e-10)
	{
	  printf ("stationary point \n");
	  n1=n1-1;
	  continue;
	}
      cout << "start traj\n";
      //memset might be faster
    std:fill( ncount, ncount+natom, 0);
      for(i=0;i<natom;i++){
	xdum=x-xr[i];
	ydum=y-yr[i];
	x1=xi-xr[i];
	y1=yi-yr[i];
	dot=dotproduct(xdum,x1,ydum,y1);
	orbit[0+i*j]=dot;
	dot=cross(xdum,x1,ydum,y1);
	orbit2[0+i*j]=dot;
      }
      
      v0x=vx;
      v0y=vy;
      lam=lambda(m, fx,fy,fz, fxx, fyy, fxy, vx, vy);
      a0x=lam*fx/m;
      a0y= lam*fy/m;
      printf ("position: %8.6f %8.6f %8.6f velocity: %8.6f %8.6f energy %8.6f wfn %14.12f lambda %8.6f force %8.6f %8.6f\n", x, y, z, vx, vy, e, 1.0/inpsi,lam,lam*fx,lam*fy);
      //Start trajectory
      for (n=1; n<j; n++)
	
	{ 
  
	  ifrev=1;
	  v0x=vx;
	  v0y=vy;
	  
	  dt = dt0;
	  vmag = sqrt(vx*vx+vy*vy);
	  //typically only occurs if the electron is too close to a node
	  if(vmag >= 68.5)
	    {
	      printf ("vmag >50percent of c %8.6f %8.6f %8.6f %8.6f \n",vmag,x,y,z);
	      ifrev = 1;
	      break;
	    }
	  if (vmag>=0.5)
	    {
	      dt=dt0*0.5/vmag;
	      printf ("dt rescaled by 0.5/vmag. dt= %8.6f Vmag= %8.6f \n", dt,vmag);
	    }
	  
	  //update position
	  x= position (x, vx, dt, a0x);
	  
	  y= position (y, vy, dt, a0y);

	  //calculate cos(theta) and sin(theta) , where theta is angle between current distance vector from atom to electron, to initial distance
	  for(i=0;i<natom;i++){
	    xdum=x-xr[i];
	    ydum=y-yr[i];
	    x1=xi-xr[i];
	    y1=yi-yr[i];
	    dot=dotproduct(xdum,x1,ydum,y1);
	    orbit[n+i*j]=dot;
	    dot=cross(xdum,x1,ydum,y1);
	    orbit2[n+i*j]=dot;
	  }
	  //calculate wavefunction, first and second derivatives based on the new updated position
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
	  //recalculate velocity based on new position
	  vx= fy*inpsi;
	  vy= -fx*inpsi;
	  lam=lambda(m, fx,fy,fz, fxx, fyy, fxy, vx, vy);
	  a0x=lam*fx/m;
	  a0y= lam*fy/m;

	  e = energy (m, vx,vy);
	  if(n%10==0&&j>2){
	    printf ("position: %8.6f %8.6f %8.6f velocity: %8.6f %8.6f energy %8.6f wfn %14.12f lambda %8.6f force %8.6f %8.6f\n", x, y, z, vx, vy, e, 1.0/inpsi,lam,lam*fx,lam*fy);
	  }

	  
	  //Calculate distance from initial position
	  d= distance(xi,x,yi,y,0,0);
 
	  //Calculate dotproduct with initial velocity
	  dot= dotproduct(vxi, vx, vyi, vy);
 
	  dist0=1000;
	  for (l=0; l<natom; l++)
	    {
	      
	      dist=distance(xr[l],x,yr[l],y,zr[l],z);
	      if (dist<dist0)
		{
		  dist0=dist;
		  ndum=l;
		}
	      
	    }
	  
	  //printf ("dist: %8.6f dot: %8.6f \n", d, dot);
	  ncount[ndum]++;
	  //Determine if the electron completed a full orbit. Dot product with initial velocity should be 1 and the distance from initial position should be small
	  if (n>100 &&dot> 0.99&&d<0.05&&d>d0)
	    {
	      printf ("went around at least one revolution \n");
	      printf ("ncount[0]: %8.6f %8.6f %8.6f ncount[1] wfn %8.6f\n", double(ncount[0])/double (n), double(ncount[1])/double (n), double(ncount[2])/double (n),1.0/inpsi);
	      printf ("final position: %8.6f %8.6f %8.6f \n", x,y,z);
	      ifrev=0; 
	      break;
	    }
	  d0=d;
	  
	}
      //If turned on, the line below will only count a trajectory if it goes around a full revolution
      if(ifrev && j>2) n1=n1-1;
 
      //Determine if a trajectory orbits an atom. An trajectory orbits the atom if cos(theta) changes sign twice during trajectory and if sin(theta) changes sign once while cos(theta) is less than zero
      //This can likely be placed inside the previous loop to be made more efficient. The orbits arrays don't nned to be compared over the entire trajectory, just between subsequent steps.
      dot=1.0;
      for(k=0;k<natom;k++){
	ifrev=0;
	for(i=0;i<n;i++){
	  dot=abs(orbit[i+k*j]+1);
	  if((orbit[i+k*j]<0.0&&orbit[i+k*j-1]>0.0)||(orbit[i+k*j]>0.0&&orbit[i+k*j-1]<0.0)||(orbit2[i+k*j]>0.0&&orbit2[i+k*j-1]<0.0&&orbit[i+k*j]<0.0)||(orbit2[i+k*j]<0.0&&orbit2[i+k*j-1]>0.0&orbit[i+k*j]<0)){
	    orbits[k]++;
	  }
	}
      }
      i=0;
      for(k=0;k<natom;k++){
	if(orbits[k]>=3){
	  i++;
	  cout << "Electron orbits atom " << k << " ";
	}
      } 
      if(i==0) {
	cout << "Electron orbits NO atoms ";
      }
      cout << "\n";
      for(k=0;k<natom;k++){
	orbits[k]=0;
      }
    }


 return 0;
}

