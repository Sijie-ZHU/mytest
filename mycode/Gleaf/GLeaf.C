/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "GLeaf.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(GLeaf, 0);
    addToRunTimeSelectionTable(functionObject, GLeaf, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::GLeaf::GLeaf
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::GLeaf::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        return true;
    }

    return false;
}


bool Foam::functionObjects::GLeaf::execute()
{
    // Assign and build fields
    const auto& T = lookupObject<volScalarField>("T");
    const auto& G = lookupObject<volScalarField>("G");

    // Calculate the temperature of the leaf (GLeaf) quantity
    Info<< "Calculating the Incident radiation on the leaf (GLeaf)" << endl;

    tmp<volScalarField> tGLeaf
    (
        volScalarField::New
        (
            "GLeaf",
            G.mesh(),
            dimPower/sqr(dimLength)
        )
    );
    volScalarField& GLeaf = tGLeaf.ref();

    const dimensionedScalar StefanBoltzmann
    (
        dimPower/(sqr(dimLength)*pow4(dimTemperature)),
        5.670374e-8
    );

    GLeaf = G - 4 * StefanBoltzmann * pow4(T);    

    word fieldNameGLeaf = "GLeaf";

    return store(fieldNameGLeaf, tGLeaf);
}

bool Foam::functionObjects::GLeaf::write()
{
    return writeObject("GLeaf");
}


// ************************************************************************* //
