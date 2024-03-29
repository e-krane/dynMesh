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

Class
    myTurboFvMesh

Description
    Simple mixer mesh using an ggi interfaces

    This tool is used to have multiple rotating regions around the same origin
    with different rpms. Creating the cellZones is not implemented
    in this tool.
    The steps to obtain the cellZones and faceZones are:

    1) use regionCellSets utility. With this command you can have different
       cellSets for each region.

    2) run command "setsToZones -noFlipMap".  After this command the
       cellSets are transformed to cellZones.

    4) for each rotating cellZone add an entry in rpm subDict of
       constant/dynamicMeshDict

    5) in parallel you need to create from all rotating coupled-interface
       patches faceZones (if you use a ggi interface between you have already
       created these faceZones), Then you need to specify these faceZones with
       the corresponding rpm in the additional "slider" subDict in
       constant/dynamicMeshDict

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

SourceFiles
    myTurboFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef myTurboFvMesh_H
#define myTurboFvMesh_H

#include "dynamicFvMesh.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class myTurboFvMesh Declaration
\*---------------------------------------------------------------------------*/

class myTurboFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

	//- sixDOFbodies dictionary
	dictionary motionDict_;

	//- rpm dictionary
	dictionary rpmWrite_;

        //- Coordinate system
        cylindricalCS cs_;

        //- Markup field for points.
        mutable vectorField* movingPointsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        myTurboFvMesh(const myTurboFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const myTurboFvMesh&);

        //- Calculate moving Points
        Foam::scalar calcMovingPoints(scalar rpmOld);

	//- Calculate rpm (if variable)
	//Foam::scalar rpmCalc(scalar rpmOld);

        //- Return moving points
        vectorField& movingPoints();


public:

    //- Runtime type information
    TypeName("myTurboFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit myTurboFvMesh(const IOobject& io);


    // Destructor

        virtual ~myTurboFvMesh();


    // Member Functions

        //- Return coordinate system
        const cylindricalCS& cs() const
        {
            return cs_;
        }

        //- Update the mesh for both mesh motion
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
