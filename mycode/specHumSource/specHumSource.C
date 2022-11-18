/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "specHumSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = e97fe6ced1b0ca40c9289d9da8b1add0a8482aa2
//
// unique function name that can be checked if the correct library version
// has been loaded

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(specHumSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    specHumSource,
    dictionary
);

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::specHumSource::specHumSource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
	Cp0_("Cp0", dimSpecificHeatCapacity, 1003.5),
	L_("L", dimLength, 0.1),
	C_("C", sqrt(dimTime)/dimLength, 130),
	L_v_("L_v", dimensionSet(0,2,-2,0,0,0,0), 2.5e+6)

{
    read(dict);           // <<-- Add
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::specHumSource::~specHumSource()
{
    Info << "Destructor: specHumSource" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::fv::specHumSource::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }
	
    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);

    return true;
}


void Foam::fv::specHumSource::addSup
(
    const volScalarField& rho,
	fvMatrix<scalar>& eqn,
    const label fieldi
)
    {
        
        Info<<"specHumSource::addSup()\n";
        
		const dimensionedScalar timeUnit
    (
        dimTime,
        1
    );
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
		const volScalarField& TLeaf = mesh_.lookupObject<volScalarField>("TLeaf");
		const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");
		const volScalarField& GLeaf = mesh_.lookupObject<volScalarField>("GLeaf");
		const dimensionedScalar Umin(dimVelocity, 0.001);
        const dimensionedScalar Umax(dimVelocity, 1000);
        volScalarField Umag(mag(U));
        Umag.clip(Umin, Umax);
        volScalarField r_a = C_ * sqrt(L_ / Umag);//[s/m]
		volScalarField qPlantSen = (TLeaf - T) * Cp0_ * rho *2 / r_a;// [W/m^2]
		volScalarField qPlantLat = GLeaf - qPlantSen; // [W/m^2]
		
        eqn += qPlantLat * LAD * rho / L_v_ ;
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

