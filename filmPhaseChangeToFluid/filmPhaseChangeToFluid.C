/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "filmPhaseChangeToFluid.H"
#include "fluidPhaseChangeToFilm.H"
#include "mappedFvPatchBaseBase.H"
#include "multicomponentFluid.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"
#include "ops.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(filmPhaseChangeToFluid, 0);
        
        addToRunTimeSelectionTable
        (
            fvModel,
            filmPhaseChangeToFluid,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::filmPhaseChangeToFluid::filmPhaseChangeToFluid
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel
    (
        sourceName, 
	    modelType, 
	    mesh, 
	    dict
    ),
    /*
    outputProperties_
    (
        IOobject
        (
            sourceName + "OutputProperties",
            this->mesh().time().name(),
            "uniform"/sourceName,
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    */
    film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName)),
    curTimeIndex_(-1),
    dmTotal_(0),
    dmTotalOld_(0),
    dmTotalCumulative_(0),
    mDot_
    (
        IOobject
        (
            this->name() + ":mDot",
            film_.db().time().name(),
            film_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimArea/dimTime, 0)
    ),
    heLatDot_
    (
        IOobject
        (
            this->name() + ":heLatDot",
            film_.db().time().name(),
            film_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimArea/dimTime, 0)
    ),
    vaporisationMode_
    (
        IOobject
        (
            this->name() + ":vaporisationMode",
            film_.db().time().name(),
            film_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    filmMassOld_(0),
    nbrPatchi_(film_.surfacePatchMap().nbrFvPatch().index()),
    overrideLRef_(dict.lookupOrDefault("overrideLRef", false)),
    LRef_(dict.lookupOrDefault<scalar>("LRef", 0)),
    Lwet_(0),
    deltaMin_(dict.lookupOrDefault<scalar>("deltaMin", 1e-13)),
    dm_
    (
        volScalarField::Internal::New
        (
            "massPhaseChangeRate",
            mesh,
            dimensionedScalar(dimMass, 0)
        )
    ),
    heLat_
    (
        volScalarField::Internal::New
        (
            "heLat",
            mesh,
            dimensionedScalar(dimEnergy/dimMass, 0)
        )
    ),
    debug_(dict.lookupOrDefault("debug", false))
{
    const word& liqName = dict.lookup("activeLiquid");
    activeLiquid_ = liquidProperties::New(liqName);
    
    Info<< "    Active liquid specie: \t" 
        << activeLiquid_->name()
        << endl;
    
    if (overrideLRef_)
    {
        Info<< "    Using user-defined LRef: \t" << LRef_ << endl;
    }
    else
    {
        if (LRef_ != 0)
        {
            Info<< "    Note! Dictionary entry 'LRef' has no effect unless"
                << " the entry 'overrideLRef' is set to true" << endl;
        }
        Info<< "    Using film wetted length as LRef" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::filmPhaseChangeToFluid::addSupFields() const
{
    return wordList
    {
        film_.alpha.name(),
        film_.thermo.he().name()
    };
}

const Foam::word Foam::fv::filmPhaseChangeToFluid::liquidName() const
{
    return activeLiquid_->name();
}

Foam::scalar Foam::fv::filmPhaseChangeToFluid::calculateSh
(
    const scalar Re,
    const scalar Sc
) const
{
    // Correlations from OpenFOAM-10 code

    if (Re < 5.0e+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


void Foam::fv::filmPhaseChangeToFluid::correctHeLat()
{
    // Film properties
    const labelList& faceCells = film_.surfacePatch().faceCells();
    const volScalarField& pc = film_.thermo.p();
    const volScalarField& Tc = film_.thermo.T();

    // Fluid properties
    const solvers::multicomponentFluid& mcf
    (
        film_.surfacePatchMap().nbrMesh()
                               .lookupObject<solvers::multicomponentFluid>
        (
            solver::typeName
        )
    );

    const label nbrPatchi = film_.surfacePatchMap().nbrFvPatch().index();
    
    const scalarField pInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.p().boundaryField()[nbrPatchi]
        )
    );
    const scalarField TInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.T().boundaryField()[nbrPatchi]
        )
    );


    scalar heLatMax = 0;
    scalar heLatMin = 0;

    Switch writeHes = true;

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        
        // Latent heat approach
        //heLat_[celli] = activeLiquid_->hl(p[celli], T[celli]);
        //

        // Enthalpy difference approach
        const label cid = mcf.thermo.species()[activeLiquid_->name()];
        const scalar hf = activeLiquid_->ha(pInf[celli], TInf[celli]);
        const scalar hc = mcf.thermo.hai(cid, pInf[celli], TInf[celli]);
        heLat_[celli] = hc - hf;
        //
        
        if (heLat_[celli] < heLatMin) { heLatMin = heLat_[celli]; }
        if (heLat_[celli] > heLatMax) { heLatMax = heLat_[celli]; }

        if (writeHes)
        {
            Info<<activeLiquid_->name() << nl
                <<"Film ha 300K 1 bar: "<< activeLiquid_->ha(100000, 300) << nl
                <<"Vapour ha 300K 1 bar: "<< mcf.thermo.hai(cid, 100000, 300) << endl;
            
            writeHes = false;
        }
    }

    // DEBUG:

    reduce(heLatMax, maxOp<scalar>());
    reduce(heLatMin, minOp<scalar>());

    Info<< "Latent heat field: (max, min):" << nl
        << "(" << heLatMax << ", " << heLatMin << ")"
        << endl;
}

void Foam::fv::filmPhaseChangeToFluid::correctLwet()
{
    // Film properties
    const labelList& faceCells = film_.surfacePatch().faceCells();
    const scalarField& V = mesh().V();
    const scalarField& delta = film_.delta();
    const scalarField& alpha = film_.alpha();
    const scalarField& magSf = film_.magSf;
    
    scalar wetV = 0;
    scalar wetA = 0;
    
    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        
        if (delta[celli] < deltaMin_)
        {
            continue;
        }
        wetV += V[celli]*alpha[celli];
        wetA += magSf[celli];
    }
    
    // Sum data from all processors and scatter
    reduce(wetV, sumOp<scalar>());
    reduce(wetA, sumOp<scalar>());
    Pstream::scatter(wetV);
    Pstream::scatter(wetA);
    
    // Calculate
    if (wetA > 0)
    {
        Lwet_ = wetV / wetA;
    }
}

void Foam::fv::filmPhaseChangeToFluid::solvePhaseChange()
{
    const scalar deltaT = mesh().time().deltaTValue();

    // Film properties
    
    const labelList& faceCells = film_.surfacePatch().faceCells();
    const scalarField& delta = film_.delta();
    const scalarField& p = film_.thermo.p();
    const scalarField& T = film_.thermo.T();
    const scalarField& rho = film_.thermo.rho();

    // Coupled region (multicomponentFluid) properties

    const solvers::multicomponentFluid& mcf
    (
        film_.surfacePatchMap().nbrMesh()
                               .lookupObject<solvers::multicomponentFluid>
        (
            solver::typeName
        )
    );

    const label nbrPatchi = film_.surfacePatchMap().nbrFvPatch().index();
    
    // Map and get fluid region quantitites at interface

    const scalarField pInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.p().boundaryField()[nbrPatchi]
        )
    );
    const vectorField UInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.U.boundaryField()[nbrPatchi]
        )
    );
    const scalarField rhoInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.rho()().boundaryField()[nbrPatchi]
        )
    );
    const scalarField muInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.mu().boundaryField()[nbrPatchi]
        )
    );
    const scalarField WInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.W()().boundaryField()[nbrPatchi]
        )
    );
    const scalarField YInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.Y(activeLiquid_->name()).boundaryField()[nbrPatchi]
        )
    );

    // Reset fields and recalculate
    dm_      = Zero;
    mDot_    = Zero;
    heLatDot_ = Zero;
    vaporisationMode_ = Zero;

    dmTotal_ = 0;

    // Mass transfer coefficient counters
    scalar hmMax = 0;
    scalar hmMin = 0;
    scalar hmAvg = 0;

    scalar BmMax = 0;
    scalar BmMin = 0;
    scalar BmAvg = 0;

    label cellCounter = 0;
    label evapCellCounter = 0;
    label boilCellCounter = 0;

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];

        // ------ Based on OpenFOAM-10 implementation ----- // 

        // Cell surface area [m^2]
        const scalar Ac = film_.magSf[celli];

        // Available evaporating mass in cell [kg]
        const scalar limMassc = max(Ac*rho[celli]*(delta[celli] - deltaMin_),
                                    Zero);

        if (delta[celli] > deltaMin_)
        {
            // Cell pressure [Pa]
            const scalar pc = p[celli];

            // Far away pressure [Pa]
            const scalar pInfc = pInf[celli];
 
            // Cell boiling temperature [K]
            const scalar Tb = activeLiquid_->pvInvert(pc);           

            // Cell temperature [K]
            // Tb coefficient (1.1) from OpenFOAM-10 source code
            // (standardPhaseChange.C)
	        // Lower bound 220 K for stability
            const scalar Tc = min(max(220, Tb*1.1), max(220, T[celli]));

            // Cell saturation pressure [Pa]
            const scalar pSat = activeLiquid_->pv(pc, Tc);

            // Boiling
            if (pSat >= 0.999*pc)
            {
                // boiling,     limited by heat transfer

                const scalar Cp = activeLiquid_->Cp(pc, Tc);
                const scalar TCorr = max(Zero, Tc - Tb);
                const scalar qCorr = limMassc*Cp*TCorr;

                dm_[celli] += qCorr/(heLat_[celli] + rootVSmall);
                
                mDot_[celli] += (qCorr/(heLat_[celli] + rootVSmall))
                                     /
                                (Ac*deltaT);
                
                vaporisationMode_[celli] = 2;
                
                boilCellCounter++;
            }
            
            // Evaporation
            else
            {
                // evaporation, limited by mass transfer 
                
                // Reference length [m], by default Lwet 
                // (corrected on each iteration)
                const scalar LRef
                (
                    overrideLRef_ ? LRef_
                                  : Lwet_
                );

                // Molecular mass of vapour [kg/kmol]
                const scalar Wvap = activeLiquid_->W();

                // Average molecular mass of bulk fluid [kg/kmol]
                const scalar Wbulk = WInf[celli];
                
                // Bulk fluid density [kg/m^3]
                const scalar rhoInfc = rhoInf[celli];

                // Bulk fluid dynamic viscosity [kg/(m*s)]
                const scalar muInfc = muInf[celli];

                // Film and the bulk fluid velocity difference at 
                // the surface [m/s]
                const vector deltaUc = (UInf[celli] - film_.U[celli]);

                // Reynolds number [-]
                const scalar Re = rhoInfc*mag(deltaUc)*LRef/muInfc;

                // Vapour mass fraction at the surface
                //   Assuming (at surface):
                //   - Fully saturated with vapour
                //   - Vapour temperature and bulk fluid temperature are equal
                //   - Vapour and bulk fluid are ideal gases
                const scalar Ys = (Wvap*pSat)/(Wvap*pSat + Wbulk*(pc-pSat));

                // Vapour mass fraction far from surface [-]
                const scalar YInfc = YInf[celli];
                
                // Vapour diffusivity [m^2/s]
                const scalar Dab = activeLiquid_->D(pInfc, Tc);
                
                // Schmidt number [-]
                const scalar Sc = muInfc/(rhoInfc*(Dab + rootVSmall));

                // Sherwood number [-]
                const scalar Sh = calculateSh(Re, Sc);

                // Mass transfer coefficient [m/s]
                const scalar hm = Sh*Dab/(LRef + rootVSmall);

                // Mass transfer driving force [-]
                scalar Bm = (Ys - YInfc)/(1.0 - Ys);

                // no condensation (for now)
                if (true)
                {
                    Bm = max(Bm, rootVSmall);
                }

                // Corrected mass transfer coefficient [m/s] accounting
                // for blowing and suction effects
                //const scalar gmCorr = hm*(log(1 + Bm)/Bm);
                //const scalar gmCorr = hm; //disabled

                // Calculate evaporation [kg]
                dm_[celli] += deltaT*Ac*rhoInfc*hm*Bm;
                
                mDot_[celli] += rhoInfc*hm*Bm;

                mDot_[celli] > 0 ? vaporisationMode_[celli] = 1
                                 : vaporisationMode_[celli] = 0;

                evapCellCounter++;

                // log heat transfer coefficient
                hmMax = max(hm, hmMax);
                hmMin = min(hm, hmMin);
                hmAvg += hm;

                // log mass transfer driving force
                BmMax = max(Bm, BmMax);
                BmMin = min(Bm, BmMin);
                BmAvg += Bm;
            }
        }

        // Bound phase change to mass limit
        dm_[celli] = min(limMassc, max(dm_[celli], Zero));
        
        mDot_[celli] = min(limMassc/(Ac*deltaT), 
                           max(mDot_[celli], Zero)
                          );

        heLatDot_[celli] = heLat_[celli]*mDot_[celli];
        
        //dmTotal_ += dm_[celli];
        dmTotal_ += mDot_[celli]*Ac*deltaT;
        
        cellCounter++;
    }

    // -------------------------- PRINTING FUNCTIONALITY ------------------------//

    //dmTotal_ = gSum(dm_);

    reduce(dmTotal_, sumOp<scalar>());
    dmTotalCumulative_ += dmTotal_;

    // Calculate sums across all processors
    reduce(cellCounter, sumOp<label>());
    reduce(boilCellCounter, sumOp<label>());
    reduce(evapCellCounter, sumOp<label>());
    reduce(hmAvg, sumOp<scalar>());
    reduce(BmAvg, sumOp<scalar>());
    
    reduce(hmMax, maxOp<scalar>());
    reduce(BmMax, maxOp<scalar>());
    reduce(hmMin, minOp<scalar>());
    reduce(BmMin, minOp<scalar>());

    // Scatter sum totals to all processors
    Pstream::scatter(evapCellCounter);
    Pstream::scatter(boilCellCounter);
    Pstream::scatter(cellCounter);

    Pstream::scatter(hmAvg);
    Pstream::scatter(BmAvg);
    Pstream::scatter(hmMax);
    Pstream::scatter(hmMin);
    Pstream::scatter(BmMax);
    Pstream::scatter(BmMin);

    if (evapCellCounter > 0) 
    {    
        hmAvg /= evapCellCounter;
        BmAvg /= evapCellCounter;
    }

    // debug printing
    if (debug)
    {
        // C++ casts are used here, because foam casts were not working
        const scalar evapCellsFrac
        (
            evapCellCounter > 0 ? static_cast<double>(evapCellCounter) 
                                  / static_cast<double>(cellCounter)
                                : 0.0
        );

        const scalar boilCellsFrac
        (
            boilCellCounter > 0 ? static_cast<double>(boilCellCounter) 
                                  / static_cast<double>(cellCounter)
                                : 0.0
        );

        tmp<volScalarField::Internal> tdhe
        (
            volScalarField::Internal::New
            (
                "dhe",
                dm_*heLat_
            )
        );

        scalar dheMax = gMax(tdhe());
        scalar dheMin = gMin(tdhe());
        scalar dheAvg = gAverage(tdhe());
        scalar dheSum = gSum(tdhe());

        // Print data
        Info<<endl;
        Info<< name() << ": " << endl
            << "  LRef: \t = "
            << (overrideLRef_ ? LRef_ : Lwet_)
            << endl

            << "  mass change:" << endl
            << "    this timestep: \t= " << dmTotal_ << endl
            << "    cumulative:    \t= " << dmTotalCumulative_ << endl
            << endl
            << "  modes of evaporation:" << endl
            << "    cells:       " << cellCounter << endl
            
            << "    evaporating: " << evapCellCounter << " (" 
            << evapCellsFrac * 100 << " \% of surface cells)" << endl
            
            << "    boiling:     " << boilCellCounter << " (" 
            << boilCellsFrac * 100 << " \% of surface cells)" << endl
            << endl

            << "  energy sources: " << endl
            << "    (min, avg, max): " << endl
            << "    (" << dheMin << ", " << dheAvg << ", " << dheMax << ")" << endl
            << "    sum: " << dheSum << endl
            << endl

            << "  mass transfer coefficient hm" << endl
            << "    (min, avg, max): " << endl
            << "    (" << hmMin << ", " << hmAvg << ", " << hmMax << ")" << endl
            << endl
            << "  mass transfer driving force Bm" << endl
            << "    (min, avg, max): " << endl
            << "    (" << BmMin << ", " << BmAvg << ", " << BmMax << ")" << endl
            <<endl;
    }

    // standard printing
    else
    {
        // Print data
        Info<<endl;
        Info<< name() << ": " << endl
            << "  mass change:" << endl
            << "    this timestep: \t= " << dmTotal_ << endl
            << "    cumulative:    \t= " << dmTotalCumulative_ << endl
            << endl

            << "  mass transfer coefficient hm" << endl
            << "    (min, avg, max): " << endl
            << "    (" << hmMin << ", " << hmAvg << ", " << hmMax << ")" << endl
            << endl
            << "  mass transfer driving force Bm" << endl
            << "    (min, avg, max): " << endl
            << "    (" << BmMin << ", " << BmAvg << ", " << BmMax << ")" << endl
            <<endl;
    }

    /*
    if (mesh().db().time().writeTime())
    {
        // Write to output properties
        if (outputProperties_.found(name()))
        {
            //TODO, write to output properties
        }
    }
    */

    //---------------------------------------------------------------------------//
}

void Foam::fv::filmPhaseChangeToFluid::printFilmMass() const
{
    const scalarField& delta = film_.delta();
    const scalarField& magSf = film_.magSf;
    const volScalarField& rho = film_.thermo.rho();

    scalar filmMass = gSum(delta*magSf*rho);

    scalar rhoMax = gMax(rho.primitiveField());
    scalar rhoMin = gMin(rho.primitiveField());
    scalar rhoAvg = gAverage(rho.primitiveField());

    scalar filmMassDelta = filmMassOld_ - filmMass;
    scalar filmMassDeltaRatio
        = Foam::mag(dmTotalOld_) > 0 ? dmTotalOld_ / filmMassDelta 
                                    : 0;

    Info<< "  film mass: \t= " << filmMass << endl
        << "      delta: \t= " << filmMassDelta << endl
        << "  ratio (evaporation / film mass change): \t= " << filmMassDeltaRatio << endl
        << endl
        << "  rho:" << endl
        << "    (min, avg, max)" << endl
        << "    (" << rhoMin << ", " << rhoAvg << ", " << rhoMax << ")" << endl
        << endl;
    
    filmMassOld_ = filmMass;
}

void Foam::fv::filmPhaseChangeToFluid::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();

    if (!overrideLRef_)
    {
        correctLwet();
    }

    correctHeLat();
    solvePhaseChange();
    printFilmMass();

    dmTotalOld_ = dmTotal_;
}

const Foam::fv::fluidPhaseChangeToFilm& Foam::fv::filmPhaseChangeToFluid::FluidToFilm
(
    const Foam::fvModels& fvModels
) const
{
    const fluidPhaseChangeToFilm* fluidToFilmPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<fluidPhaseChangeToFilm>(fvModels[i]))
        {
            const fluidPhaseChangeToFilm& fluidToFilm
            (
                refCast<const fluidPhaseChangeToFilm>(fvModels[i])
            );

            const labelList fluidToFilmPatchIndices
            (
                fluidToFilm.filmPatchIndices()
            );

            forAll(fluidToFilmPatchIndices, j)
            {
                if
                (
                    fluidToFilmPatchIndices[j]
                        == 
                    film_.surfacePatchMap().nbrFvPatch().index()
                )
                {
                    fluidToFilmPtr = &fluidToFilm;
                }
            }
        }
    }

    if (!fluidToFilmPtr)
    {
        FatalErrorInFunction
            << "Cannot find fluidPhaseChangeToFilm fvModel for this film "
               "in neighbour region '" << film_.surfacePatchMap().nbrMesh().name()
            << "'" << exit(FatalError);
    }

    return *fluidToFilmPtr;
}

void Foam::fv::filmPhaseChangeToFluid::addSup
(
    const volScalarField& rho,
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn
) const
{
    if (true)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&alpha == &film_.alpha && &eqn.psi() == &film_.alpha)
    {
        const scalar deltaT = mesh().time().deltaTValue();

        // Base source term [kg/m2/s]
        volScalarField::Internal tSu
        (
            volScalarField::Internal::New
            (          
                "SuRho", 
                mDot_.clone()
            )
        );
        
        // Modify source term into the correct form
        tSu *= film_.magSf;
        tSu /= mesh().V(); // [kg/m3/s]

        // Source term [kg/m3/s]
        // Explicit source term (because eqn is the alpha-eqn)
        eqn += fvm::Su(-tSu, eqn.psi());

        /*
        // Explicit transfer in/out of the film
        volScalarField::Internal tSu
        (
            volScalarField::Internal::New
            (          
                "SuRho", 
                dm_.clone()
            )
        );
        tSu /= dimensionedScalar(dimTime, deltaT);
        tSu /= mesh().V();

        // Source term [kg/m3/s]
        eqn += fvm::Su(-tSu, eqn.psi()); // explicit, since eqn is alpha eq.
        */
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << alpha.name() << " is not implemented"
            << exit(FatalError);
    }
}

void Foam::fv::filmPhaseChangeToFluid::addSup
(
    const volScalarField& rho,  
    const volScalarField& alpha,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    if (true)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }
    
    if (&he == &film_.thermo.he() && &eqn.psi() == &he)
    {
        const scalar deltaT = mesh().time().deltaTValue();

        // Mass source term [kg/m2/s]
        volScalarField::Internal tSuRho
        (
            volScalarField::Internal::New
            (          
                "SuRho", 
                mDot_.clone()
            )
        );

        // Latent heat of vaporisation source term [J/m2/s]
        volScalarField::Internal tSuHe
        (
            volScalarField::Internal::New
            (          
                "SuHeLat", 
                heLatDot_.clone()
            )
        );

        // Modify source terms into the correct form
        tSuRho *= film_.magSf;
        tSuRho /= mesh().V(); // [kg/m3/s]
        
        tSuHe *= film_.magSf;
        tSuHe /= mesh().V();  // [J/m3/s]
        
        // Add energy contained in the leaving mass to the source term
        tSuHe += tSuRho*he;

        // Source term [J/m3/s]
        // Explicit source term
        // Latent heat of vaporisation is supplied fully by the film
        eqn += fvm::Su(-tSuHe, eqn.psi());  

        /*
        // Multiply dm with heLat
        volScalarField::Internal tSu
        (
            volScalarField::Internal::New
            (          
                "SuHe", 
                dm_.clone()
            )
        );
        tSu /= dimensionedScalar(dimTime, deltaT);
        tSu /= mesh().V();

        // Transfer out of the film
        
        // Phase change energy is supplied by the film
        // Source term [J/m3/s]
        eqn += fvm::Su(-tSu*(he + heLat_), eqn.psi());
        
        // Phase change energy is supplied by the bulk fluid
        //eqn += fvm::Su(-tSu*he, eqn.psi());
        */

    }   
    else
    {
        FatalErrorInFunction
            << "Support for field " << he.name() << " is not implemented"
            << exit(FatalError);
    }
}

/*
Foam::tmp<Foam::volScalarField> Foam::fv::filmPhaseChangeToFluid::rhoTrans() const
{
    //tmp<volScalarField> trhoTrans =
    //    volScalarField::New
    //    (
    //        "trhoTrans",
    //        mesh(),
    //        dimensionedScalar(dimMass/dimVolume, Zero)
    //    );

    //trhoTrans.ref().internalFieldRef() += (dm_/mesh().V());

    //return trhoTrans;

    return rhoTransV()/mesh().V();
}
*/

Foam::tmp<Foam::volScalarField> Foam::fv::filmPhaseChangeToFluid::rhoTransV() const
{
    const dimensionedScalar deltaT
    (
        dimTime, 
        mesh().time().deltaTValue()
    );

    tmp<volScalarField> trhoTrans =
        volScalarField::New
        (
            "trhoTrans",
            mesh(),
            dimensionedScalar(dimMass, Zero)
        );

    //trhoTrans.ref().internalFieldRef() += dm_;
    trhoTrans.ref().internalFieldRef() += mDot_*film_.magSf*deltaT;

    return trhoTrans;
}

/*
Foam::tmp<Foam::volScalarField> Foam::fv::filmPhaseChangeToFluid::hTrans() const
{
    
    //tmp<volScalarField> thTrans =
    //    volScalarField::New
    //    (
    //        "thTrans",
    //        mesh(),
    //        dimensionedScalar(dimEnergy/dimVolume, Zero)
    //    );

    // Phase change energy is supplied from the film
    // (zero source)

    // Phase change energy is supplied from the bulk fluid
    //thTrans.ref().internalFieldRef() -= (dm_*heLat_/mesh().V());
    
    //return thTrans;

    return hTransV()/mesh().V();
}
*/

Foam::tmp<Foam::volScalarField> Foam::fv::filmPhaseChangeToFluid::hTransV() const
{
    const dimensionedScalar deltaT
    (            
        dimTime, 
        mesh().time().deltaTValue()
    );

    tmp<volScalarField> thTrans =
        volScalarField::New
        (
            "thTrans",
            mesh(),
            dimensionedScalar(dimEnergy, Zero)
        );

    // Do this stuff inside the phase change calculations? ----------------------
    
    const solvers::multicomponentFluid& mcf
    (
        film_.surfacePatchMap().nbrMesh()
                               .lookupObject<solvers::multicomponentFluid>
        (
            solver::typeName
        )
    );

    const label nbrPatchi = film_.surfacePatchMap().nbrFvPatch().index();

    const scalarField pInf
    (
        film_.surfacePatchMap().fromNeighbour
        (
            mcf.thermo.p().boundaryField()[nbrPatchi]
        )
    );

    const label speciei = mcf.thermo.specieIndex(mcf.thermo.Y(activeLiquid_->name()));

    // Mass change [kg]
    const scalarField dm = mDot_*film_.magSf*deltaT;
    // Vapour sensible enthalpy [J/kg]
    const scalarField hs = mcf.thermo.hsi(speciei, pInf, film_.thermo.T());

    // --------------------------------------------------------------------------

    // Energy source [J]
    // Sensible enthalpy of the mass transferred into bulk fluid
    thTrans.ref().primitiveFieldRef() += dm*hs;

    // Phase change energy is supplied from the film
    // (zero source)

    // Phase change energy is supplied from the bulk fluid  
    //thTrans.ref().internalFieldRef() -= (dm_*heLat_);
    
    return thTrans;
}

void Foam::fv::filmPhaseChangeToFluid::topoChange(const polyTopoChangeMap&)
{
    dm_.setSize(mesh().nCells());
    heLat_.setSize(mesh().nCells());
    mDot_.setSize(mesh().nCells());
    heLatDot_.setSize(mesh().nCells());
    vaporisationMode_.setSize(mesh().nCells());
}

void Foam::fv::filmPhaseChangeToFluid::mapMesh(const polyMeshMap&)
{
    dm_.setSize(mesh().nCells());
    heLat_.setSize(mesh().nCells());
    mDot_.setSize(mesh().nCells());
    heLatDot_.setSize(mesh().nCells());
    vaporisationMode_.setSize(mesh().nCells());
}

void Foam::fv::filmPhaseChangeToFluid::distribute(const polyDistributionMap&)
{
    dm_.setSize(mesh().nCells());
    heLat_.setSize(mesh().nCells());
    mDot_.setSize(mesh().nCells());
    heLatDot_.setSize(mesh().nCells());
    vaporisationMode_.setSize(mesh().nCells());
}

bool Foam::fv::filmPhaseChangeToFluid::movePoints()
{
    return true;
}


// ************************************************************************* //

