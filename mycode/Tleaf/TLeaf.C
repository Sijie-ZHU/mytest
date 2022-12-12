/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "TLeaf.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TLeaf, 0);
    addToRunTimeSelectionTable(functionObject, TLeaf, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::TLeaf::converged
(
    const volScalarField& phi
) const
{
    return
        max(mag(phi.primitiveField() - phi.prevIter().primitiveField()))
      < tolerance_;
}

Foam::volScalarField Foam::functionObjects::TLeaf::plantFilter
(
    const volScalarField& LAD
) const
{
    volScalarField plantFilter = LAD * dimensionedScalar(dimLength, 1.0);

    forAll(plantFilter, i)
    {
        if(plantFilter[i] != 0)
        {
            plantFilter[i] = 1;
        } 
		else 
		{
            plantFilter[i] = 0;
        }
    }

    return plantFilter;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TLeaf::TLeaf
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    Cp0_("Cp0", dimSpecificHeatCapacity, 1003.5),
    //rho0_("rho0", dimDensity, 1.225),
    //p0_("p0", dimPressure, 101325),
    L_("L", dimLength, 0.1),
    C_("C", sqrt(dimTime)/dimLength, 130),
    r_sMin_("r_sMin", dimTime/dimLength, 150),
    R_a_("R_a", dimGasConstant, 287.042),
    R_v_("R_v", dimGasConstant, 461.524),
    L_v_("L_v", dimensionSet(0,2,-2,0,0,0,0), 2.5e+6),
    //G0_("G0", dimPower/sqr(dimLength), 500),
    tolerance_(1e-4),
    maxIter_(100)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TLeaf::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Cp0_.readIfPresent(dict);
        //rho0_.readIfPresent(dict);
        //p0_.readIfPresent(dict);
        L_.readIfPresent(dict);
        C_.readIfPresent(dict);
        r_sMin_.readIfPresent(dict);
        R_a_.readIfPresent(dict);
        R_v_.readIfPresent(dict);
        L_v_.readIfPresent(dict);
        //G0_.readIfPresent(dict);
        tolerance_ = dict.getOrDefault("tolerance", 1e-4);
        maxIter_ = dict.getOrDefault("maxIter", 100);

        return true;
    }

    return false;
}


bool Foam::functionObjects::TLeaf::execute()
{
    // Assign and build fields
    const auto& T = lookupObject<volScalarField>("T");
    const auto& GLeaf = lookupObject<volScalarField>("GLeaf");
    const auto& LAD = lookupObject<volScalarField>("LAD");
    const auto& p = lookupObject<volScalarField>("p");
	const auto& rho = lookupObject<volScalarField>("rho");
    const auto& w = lookupObject<volScalarField>("specHum");
    const auto& U = lookupObject<volVectorField>("U");

    // Calculate the temperature of the leaf (TLeaf) quantity
    Info<< "Calculating the Temperature of the leaf (TLeaf)" << endl;

    tmp<volScalarField> tTLeaf
    (
        volScalarField::New
        (
            "TLeaf",
            T.mesh(),
            dimTemperature
        )
    );
    volScalarField& TLeaf = tTLeaf.ref();

    // Initial guess
    TLeaf = T;

    label i = 0;

    TLeaf.storePrevIter();

    const dimensionedScalar StefanBoltzmann
    (
        dimPower/(sqr(dimLength)*pow4(dimTemperature)),
        5.670374e-8
    );

    const dimensionedScalar pressureUnit
    (
        dimPressure,
        1
    );

    do
    {
        TLeaf.storePrevIter();
        
        // Radiative flux to the leaf (+ if entering the leaf) (as if all the radiation is absorbed)
        // [W/m^2] = [W/m^2]
        volScalarField qPlantRad = GLeaf;

        // Aerodinamic resistance
        // [s/m] = [s^0.5/m] * [s^0.5]
        const dimensionedScalar Umin(dimVelocity, 0.001);
        const dimensionedScalar Umax(dimVelocity, 1000);
        volScalarField Umag(mag(U));
        Umag.clip(Umin, Umax);
        volScalarField r_a = C_ * sqrt(L_ / Umag);//[s/m]

        // Partial vapour pressure at T and w
        // [Pa]
        const dimensionedScalar A(dimless, 17.27);
        const dimensionedScalar B(dimTemperature, 273.15);
        const dimensionedScalar C(dimTemperature, 237.3);
        const dimensionedScalar E(dimless, 0.61078);
		const dimensionedScalar F(dimPressure/dimTemperature/dimDensity, 461.5);
        
        // [Pa]/[kg/(m s^2)]
		
		//Tetens equation [kPa]
		//volScalarField pVSatAir = 1000 * pressureUnit * E * exp((A * (T - B)) / ((T - B) + C));
        volScalarField pVAir = w * R_v_ * T;

        // Saturation vapour pressure at TLeaf (Tetens equation)
        // [Pa]/[kg/(m s^2)]
        volScalarField pVSatLeaf = 1000 * pressureUnit * E * exp((A * (TLeaf - B)) / ((TLeaf - B) + C));

        // Stomatal resistance [r_s = r_sMin * f1(G0) * f2(D)]
        // [s/m] = [s/m] * [-] * [-]
        // const dimensionedScalar a_1(dimPower/sqr(dimLength), 169); // [W/m^2]
        // const dimensionedScalar a_2(dimPower/sqr(dimLength), 18); // [W/m^2]
        // const dimensionedScalar a_3(dimless/sqr(dimPressure), 0.005); // [1/kPa^2]
        // const dimensionedScalar D0(dimPressure, 1.2); // [kPa]
        
        // volScalarField D = (pVSatAir - pVAir) / 1000; // [kPa]

        // volScalarField r_s = r_sMin_ * (a_1 + G0_) / (a_2 + G0_) * (1 + a_3 * sqr(D - D0));

        volScalarField r_s = r_a * 0 + r_sMin_;

        // e) Vapour mass flux from the leaf
        // [kg/(s m^2)] = ([J/(kg K)] * /[kg/(m s^2)]) / ([J/(kg K)] * [s/m])
        volScalarField g_vLeaf = R_a_ * rho * (pVSatLeaf - pVAir) / ( p * R_v_ * (r_a + r_s) );

        //    latent flux from the leaf
        // [W/m^2] = [J/kg] * [kg/(s m^2)]
        volScalarField qPlantLat = L_v_ * g_vLeaf;

        // f) Leaf energy balance
        // [K] = [K] + [W/m^2] / ([J/(kg K)] * [s/m])
        TLeaf = T + (qPlantRad - qPlantLat) / (2 * Cp0_ * rho / r_a);

        // reset TLeaf to T where there is not any plant
        TLeaf = (1 - plantFilter(LAD)) * T + plantFilter(LAD) * TLeaf;

    } while (!converged(TLeaf) && i++ < maxIter_);

    if (i < maxIter_) {
        Info
            << "The leaf temperature converged within " << i
            << " iterations" << endl;
    } else if (i == maxIter_) {
        WarningInFunction
            << "The leaf temperature did not converge within " << i
            << " iterations" << nl;
    }

    

    word fieldNameTLeaf = "TLeaf";

    return store(fieldNameTLeaf, tTLeaf);
}

bool Foam::functionObjects::TLeaf::write()
{
    return writeObject("TLeaf");
}


// ************************************************************************* //
