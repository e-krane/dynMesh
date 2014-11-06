/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "myTurboFvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "ZoneIDs.H"
#include "forces.H"
#include "IFstream.H"
#include <iostream>

//#include "OFstream.H"
//#include "sixDoFRigidBodyMotion.H"
//#include "sixDOFbodies.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myTurboFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, myTurboFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::myTurboFvMesh::calcMovingPoints(scalar rpmOld)
{
    if (debug)
    {
        Info<< "void myTurboFvMesh::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

// Remove to be able to update every time step
//    if (movingPointsPtr_)
//    {
//        FatalErrorIn("void myTurboFvMesh::calcMovingMasks() const")
//            << "point mask already calculated"
//            << abort(FatalError);
//    }

    // Retrieve the cell zone Names
    const wordList cellZoneNames = cellZones().names();

    // Set the points
    movingPointsPtr_ = new vectorField(allPoints().size(), vector::zero);

    vectorField& movingPoints = *movingPointsPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    // Claculate rotational speed from fluid forces
    scalar rpm;
    #include "variableRpm.H"	

    // Write data to file readable by matlab etc.
    fileName myDir = time().path();
    word tracingFileName = "rpmData.dat";
    fileName myFile = myDir/tracingFileName;

    std::ofstream os
    ( 
        myFile.c_str(),
        ios_base::app
    );

    os << time().value() << tab 
       << rpm   << tab 
       << Mflow << tab 
       << Mpow  << tab
       << alpha << tab
       << Mpress  << tab 
       << Mtau <<"\n";
    //-----------------------------------------

    const labelList& cellAddr =
          cellZones()[cellZones().findZoneID("rotor")];
	
     //rpm = readScalar
     //(
     //    dict_.subDict("rpm").lookup("rotor")
     //);
     Info<< "Moving Cell Zone Name: " << "rotor" 
         << " rpm: " << rpm << endl;

     forAll (cellAddr, cellI)
     {
         const cell& curCell = c[cellAddr[cellI]];
         
         forAll (curCell, faceI)
         {
                	
	     const face& curFace = f[curCell[faceI]];
             
             forAll (curFace, pointI)
             {
          	         
                 movingPoints[curFace[pointI]] =
                 vector(0, rpm/60.0*360.0, 0);
             }
          }
     }
return rpm;
}

//Foam::scalar Foam::myTurboFvMesh::rpmCalc(scalar rpmOld)
//{
//    #include "variableRpm.H"
//    return rpm;
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::myTurboFvMesh::myTurboFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "runner",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        )
    ),
    rpmWrite_
    (
        IOdictionary
        (
            IOobject
            (
                "rpm",
		time().constant(),
                //time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            )
        )
    ),
    cs_
    (
        "coordinateSystem",
        dict_.subDict("coordinateSystem")
    ),
    movingPointsPtr_(NULL)
{
    // Make sure the coordinate system does not operate in degrees
    // Bug fix, HJ, 3/Oct/2011
    if (!cs_.inDegrees())
    {
        WarningIn("myTurboFvMesh::myTurboFvMesh(const IOobject& io)")
            << "Mixer coordinate system is set to operate in radians.  "
            << "Changing to rad for correct calculation of angular velocity."
            << nl
            << "To remove this message please add entry" << nl << nl
            << "inDegrees true;" << nl << nl
            << "to the specification of the coordinate system"
            << endl;

        cs_.inDegrees() = true;
    }


    Info<< "Turbomachine Mixer mesh:" << nl
        << "    origin: " << cs().origin() << nl
        << "    axis  : " << cs().axis() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myTurboFvMesh::~myTurboFvMesh()
{
    deleteDemandDrivenData(movingPointsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return moving points mask
Foam::vectorField& Foam::myTurboFvMesh::movingPoints()
{
    scalar rpm;
    scalar rpmOld;
    rpmOld = readScalar( rpmWrite_.lookup("rpm") );
    if (time().value()<0.0018)
    {    
        rpmOld = readScalar( dict_.subDict("variableRpm").lookup("startRpm") );
    }
    rpm = calcMovingPoints(rpmOld);
    rpmWrite_.add("rpm", rpm, true);

    return *movingPointsPtr_;
}


bool Foam::myTurboFvMesh::update()
{
    movePoints
    (
        cs_.globalPosition
        (
            cs_.localPosition(allPoints())
          + movingPoints()*time().deltaT().value()
        )
    );
    
// The mesh is not morphing
    return false;
}


// ************************************************************************* //
