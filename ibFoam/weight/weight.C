#include "weight.H"
 

namespace Foam
{

}


Foam::scalar Foam::weight(Foam::vector xcloud, Foam::vector xgrid)
{

	scalar hx=.05;  //TODO read from mesh
    scalar hy=.05;
    scalar hz=.05;
	scalar phi1;
	scalar phi2;
    scalar phi3;

///computing scaled distance between grid point and cloud point 
   scalar r1 = (1/hx)*(xgrid[0]-xcloud[0]);
   scalar r2 = (1/hy)*(xgrid[1]-xcloud[1]);
   scalar r3 = (1/hz)*(xgrid[2]-xcloud[2]);

	 r1=fabs(r1);
	 r2=fabs(r2);
     r3=fabs(r3);
     //Info<<"r1 ="<<r1<<endl;
     //Info<<"r2 ="<<r2<<endl;
     //Info<<"r3 ="<<r3<<endl;

///computing phi1

	if (r1<=1.0)
	{
		phi1=(1.0/8.0)*(3.0-2.0*r1+sqrt(1.0+4.0*r1-4.0*r1*r1));
	} 
	else if(r1>=1.0 && r1<=2.0)
	{
		phi1=(1.0/8.0)*(5.0-2.0*r1-sqrt(-7.0+12.0*r1-4.0*r1*r1));
	} 
	else
	{
		phi1=0.0;
	} 


///computing phi2

	if (r2<=1.0)
	{
		phi2=(1.0/8.0)*(3.0-2.0*r2+sqrt(1.0+4.0*r2-4.0*r2*r2));
    } 
	else if(r2>=1.0 && r2<=2.0)
	{
		phi2=(1.0/8.0)*(5.0-2.0*r2-sqrt(-7.0+12.0*r2-4.0*r2*r2));
	} 
	else
	{
		phi2=0.0;
	}

///computing phi3

	if (r3<=1.0)
	{
		phi3=(1.0/8.0)*(3.0-2.0*r3+sqrt(1.0+4.0*r3-4.0*r3*r3));
    } 
	else if(r3>=1.0 && r3<=2.0)
	{
		phi3=(1.0/8.0)*(5.0-2.0*r3-sqrt(-7.0+12.0*r3-4.0*r3*r3));
	} 
	else
	{
		phi3=0.0;
	}
 


 scalar w=(phi1/hx)*(phi2/hy)*(phi3/hz);

 return w;

}


