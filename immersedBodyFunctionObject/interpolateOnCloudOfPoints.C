/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interpolateOnCloudOfPoints.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "surfaceFields.H"
#include "scalar.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interpolateOnCloudOfPoints, 0);
    addToRunTimeSelectionTable(functionObject, interpolateOnCloudOfPoints, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interpolateOnCloudOfPoints::interpolateOnCloudOfPoints
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
     
)
:

     
    fvMeshFunctionObject(name, runTime, dict),
    patchNames_(dict.lookup("patches")),

pointcloud_
    (
        IOobject
        (
            "cloud_of_points",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    ),

current_pointcloud_
    (
        IOobject
        (
            "currentcloudofpoints",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    ),
 
connectivity_matrix_
   (
        IOobject
        (
            "connectivity",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    ),
volume_(connectivity_matrix_.size(),0.0) 	 
{
    read(dict);

	forAll(volume_,i)
	{
		vector  p0=pointcloud_[connectivity_matrix_[i][0]];
	  	vector  p1=pointcloud_[connectivity_matrix_[i][1]];
	  	vector  p2=pointcloud_[connectivity_matrix_[i][2]];
	  	vector  p3=pointcloud_[connectivity_matrix_[i][3]];

		volume_[i]= (1.0/6.0)*fabs( ((p2-p0)^(p3-p0))& (p1-p0) );
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::interpolateOnCloudOfPoints::~interpolateOnCloudOfPoints()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interpolateOnCloudOfPoints::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


Foam::scalar Foam::functionObjects::interpolateOnCloudOfPoints::weight(Foam::vector xcloud, Foam::vector xgrid)
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




//**********************************************************************************************************************************//




bool Foam::functionObjects::interpolateOnCloudOfPoints::execute()
{

Info<<"begin execute"<<endl;

/*******************************************VELOCITY INTERPOLATION********************/

        scalar hx=.05;  //TODO read from mesh
        scalar hy=.05;
        scalar hz=.05;
        scalar delta=0.005;  //TODO read time step
        //FOR LOOP RUNNING OVER ALL CLOUD POINTS

		forAll(current_pointcloud_,i)
			{

		 	 
			//COMPUTING VELOCITY AT THE PROBE POINT

			   vector probePoint=current_pointcloud_[i];               
			   vector velocity(0.0,0.0,0.0);                                     //initialising probe point velocity

			   const volVectorField& U =mesh_.lookupObject<volVectorField>("U");   
			   labelHashSet cell_neighbors=findNeighbourCells(probePoint);

               

			  	
							forAllConstIter(labelHashSet, cell_neighbors,iter) 
					   	 		{
									  label celli = iter.key();
									  vector cellicentre=mesh_.cellCentres()[celli];
									  scalar W=weight(probePoint,cellicentre);

								 	  //Info<<"U of neighbour cell"<<celli<<"is"<<U[celli]<<endl;
									  //Info<<"Weight of neighbour cell"<<celli<<"is"<<W<<endl;

									   

									  velocity=velocity+(W*U[celli]*hx*hy*hz);

						 		}
				 
				  

                //MOVE THE POINT FOR THE FIRST HALF TIME STEP//

                current_pointcloud_[i]=probePoint+ velocity*(delta/2);     
                //Info<<"current point"<<current_pointcloud_[i]<<endl;
                //Info<<"velocity of the cloud point"<<velocity<<endl;
             }



  
/*********************************************3D FORCE**************************************************/

           vectorField force(pointcloud_.size(),Zero);           ///Initialising Force at each vertex

           forAll(volume_,i)                                     ///For loop runs over each element
				{

                    //Info<<"element"<<i<<endl;

                    //Coordinates of reference tetrahedron:
					vector  s0=pointcloud_[connectivity_matrix_[i][0]];   
	  				vector  s1=pointcloud_[connectivity_matrix_[i][1]];
	  				vector  s2=pointcloud_[connectivity_matrix_[i][2]];
                    vector  s3=pointcloud_[connectivity_matrix_[i][3]];
 

                    //Coordinates of deformed tetrahedron:
					vector  x0=current_pointcloud_[connectivity_matrix_[i][0]];   
	  				vector  x1=current_pointcloud_[connectivity_matrix_[i][1]];
	  				vector  x2=current_pointcloud_[connectivity_matrix_[i][2]];
                    vector  x3=current_pointcloud_[connectivity_matrix_[i][3]];

                    //difference between coordinates of reference nodes:
                    vector  ds1=s1-s0;
                    vector  ds2=s2-s0;
                    vector  ds3=s3-s0;


                    //Deformation gradient:

                    //First row:

				    vector  a0=(x1[0]-x0[0])*(ds2^ds3)+(x2[0]-x0[0])*(ds3^ds1)+(x3[0]-x0[0])*(ds1^ds2);
                    vector  a1=(x1[1]-x0[1])*(ds2^ds3)+(x2[1]-x0[1])*(ds3^ds1)+(x3[1]-x0[1])*(ds1^ds2);
                    vector  a2=(x1[2]-x0[2])*(ds2^ds3)+(x2[2]-x0[2])*(ds3^ds1)+(x3[2]-x0[2])*(ds1^ds2);
                    
                    //Info<<"volume"<<volume_[i]<<endl;

                    a0=(1.0/(6.0*volume_[i]))*a0;
                    a1=(1.0/(6.0*volume_[i]))*a1;
                    a2=(1.0/(6.0*volume_[i]))*a2;

                    //Info<<"def_grad"<<a0<<endl;
                    //Info<<"def_grad"<<a1<<endl;
                    //Info<<"def_grad"<<a2<<endl;

                    //Piola Kirchoff Stress:

                    //For simplicity we set deformationgradient=Piola Kirchoff stress

                    vector P0=a0;
                    vector P1=a1;
                    vector P2=a2;
                    
                    //Derivative of deformation gradient wrt position:


	               vector dela_delx0 = -(1.0/(6.0*volume_[i]))*((s1^s2) + (s2^s3) + (s3^s1));
	               vector dela_delx1 =  (1.0/(6.0*volume_[i]))*(ds2^ds3);
	               vector dela_delx2 =  (1.0/(6.0*volume_[i]))*(ds3^ds1);
	               vector dela_delx3 =  (1.0/(6.0*volume_[i]))*(ds1^ds2);

                    
                   vector f0( P0&dela_delx0, P1&dela_delx0, P2&dela_delx0);
                   vector f1( P0&dela_delx1, P1&dela_delx1, P2&dela_delx1);
                   vector f2( P0&dela_delx2, P1&dela_delx2, P2&dela_delx2);
                   vector f3( P0&dela_delx3, P1&dela_delx3, P2&dela_delx3);

                   f0 = -volume_[i]*f0;
                   f1 = -volume_[i]*f1;
                   f2 = -volume_[i]*f2;
                   f3 = -volume_[i]*f3;

                   //Info<<"f0"<<f0<<endl;
                   //Info<<"f1"<<f1<<endl;
                   //Info<<"f2"<<f2<<endl;
                   //Info<<"f3"<<f3<<endl;

                   force[connectivity_matrix_[i][0]]=force[connectivity_matrix_[i][0]] + f0;
                   force[connectivity_matrix_[i][1]]=force[connectivity_matrix_[i][1]] + f1;
				   force[connectivity_matrix_[i][2]]=force[connectivity_matrix_[i][2]] + f2;
                   force[connectivity_matrix_[i][3]]=force[connectivity_matrix_[i][3]] + f3;

 
				}  

                    
/****************************FORCE SPREADING OPERATION************************************/
 volVectorField& f = const_cast<volVectorField&>(mesh_.lookupObject<volVectorField>("f"));

    //FOR LOOP RUNNING OVER ALL CLOUD POINTS

		forAll(current_pointcloud_,i)
			{

			 
			 

			   vector probePoint=current_pointcloud_[i];               
			   vector probePoint_force=force[i];

			    
			   labelHashSet cell_neighbors=findNeighbourCells(probePoint);

			  	
							forAllConstIter(labelHashSet, cell_neighbors,iter) 
					   	 		{
									  label celli = iter.key();
									  vector cellicentre=mesh_.cellCentres()[celli];
									  scalar W=weight(cellicentre,probePoint);


                                      //Info<<"Weight of neighbour cell"<<celli<<"is"<<W<<endl;
								 	    
                                      f[celli]=f[celli]+W*probePoint_force;
                                      //Info<<"f of neighbour cell"<<celli<<"is"<<f[celli]<<endl;
							
						 		}
 

             }

//Info<<"I am running till this point"<<endl;
//Info<<"grid_force"<<f<<endl;

/***********************************************************************************************/

 
    
return true;
}



//TODO check whether the cloud point is in mesh or not

Foam::labelHashSet Foam::functionObjects::interpolateOnCloudOfPoints::findNeighbourCells(const vector probePoint) const
{
    
     

    // find cell containing this point
    label celli = mesh_.findCell(probePoint);

    // container for neighbours set by dumping the cell containing it
    labelHashSet neighbourCellSet(0);
    neighbourCellSet.set(celli);

     
    // number of layers
    int nLayers = 2;
    for (int n = 0; n < nLayers; n++)
    {
        // make a copy of marked cells
        labelHashSet markedNeighbours = neighbourCellSet;

        // loop over all marked cells
        forAllConstIter(labelHashSet, markedNeighbours,iter)
        {
            celli = iter.key();
             
            // get points of celli
            labelList celliPoints = mesh_.cellPoints()[celli];

            forAll(celliPoints,j)
            {
                // get neighbor cells of j th point
                labelList cellJNeighbours = mesh_.pointCells()[celliPoints[j]];

                // append these cells in neighbourCellSet
                forAll(cellJNeighbours, k)
                {
                    neighbourCellSet.set(cellJNeighbours[k]);
                 
                }
        
            }
        }  

    }
 
    return neighbourCellSet;

}







bool Foam::functionObjects::interpolateOnCloudOfPoints::write()
{
    return true;
}


// ************************************************************************* //
