/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

Application
     ibFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * GEOMETRY* * * * * * * * * * * * * * * * * * * //

        #include "volume.H"                     //volume of solid mesh


        scalar hx=.05;                          //TODO read from mesh
        scalar hy=.05;
        scalar hz=.05;
        scalar delta=0.005;                     //TODO read time step

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

         #include "CourantNo.H"

         #include "deform1.H"              //deform solid        X(n+1/2)       
         #include "Force.H"                //compute solid force F(n+1/2)        
         #include "spread.H"               //spread force        f(n+1/2)  
         #include "predictor.H"            //fluid solver        v(n+1/2) and p(n+1/2)                          
         #include "createPhi.H"            //convective term     v(n+1/2)          
         #include "deform2.H"              //deform solid        X(n+1)              
         #include "corrector.H"            //fluid solver        v(n+1)   and p(n+1) 
   

        runTime.write();

        runTime.printExecutionTime(Info);

		Info<<"fluid solver end"<<endl;
    }

    Info<< "Hey Don't worry I'm running\n" << endl;

    return 0;
}


// ************************************************************************* //
