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

#include "fluidPhaseChangeToFilm.H"
#include "filmPhaseChangeToFluid.H"
#include "mappedFvPatchBaseBase.H"
#include "multicomponentFluid.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(fluidPhaseChangeToFilm, 0);
        
        addToRunTimeSelectionTable
        (
            fvModel,
            fluidPhaseChangeToFilm,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fluidPhaseChangeToFilm::fluidPhaseChangeToFilm
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
    coupled_(dict.lookupOrDefault("coupled", true)),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    mcf_(mesh.lookupObject<solvers::multicomponentFluid>(solver::typeName)),
    curTimeIndex_(-1),
    FilmToFluidsCachedList_(),
    FilmToFluidsCached_(false)
{
    // Get all patches
    if (coupled_)
    {
        if (dict.found("patches") && dict.found("patch"))
        {
            Info<<name()<<": Both 'patches' and 'patch' defined. "
                <<"Ignoring the 'patch' entry."
                <<endl;
        }
        
        if (dict.found("patches"))
        {
            labelHashSet patchiSet = mesh.boundaryMesh().patchSet
            (
                dict.lookup<wordReList>("patches")
            );

            if (patchiSet.size() > 0)
            {
                forAll(mesh.boundary(), i)
                {
                    if (patchiSet[i])
                    {
                        filmPatchiList_.append((mesh.boundary()[i]).index());
                    }
                }
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Unable to find any patches with the regular expression: " 
                    << dict.lookup<wordReList>("patches")
                    << exit(FatalIOError);
            }
        }

        else if (dict.found("patch"))
        {
            const word patchName(dict.lookup("patch"));
            const label patchIndex = mesh.boundaryMesh().findIndex(patchName);
            
            if (patchIndex >= 0)
            {
                filmPatchiList_.append(patchIndex);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Unable to find patch " << patchName 
                    << exit(FatalIOError);
            }
        }

        else
        {
            FatalIOErrorInFunction(dict)
                << "Neither 'patch' or 'patches' specified"
                << exit(FatalIOError);
        }

        // Check that all patches are coupled
        Info<<name()<<endl;
        Info<<name()<<": Coupled patches: "<<endl;
        Info<<name()<<": (patch)\t|    (neighbour region)"<<endl;

        forAll(filmPatchiList_, i)
        {
            const fvPatch& filmPatch
            (
                mesh.boundary()[filmPatchiList_[i]]
            );

            Info<<name()<<": "<<filmPatch.name()<<"\t";

            if (!isA<mappedFvPatchBaseBase>(filmPatch))
            {
                FatalErrorInFunction
                    << "Patch '" << filmPatch.name()
                    << "' is not coupled! (Cannot be cast into"
                    << " 'mappedFvPatchBaseBase')"
                    << exit(FatalError);
            }

            Info<<"|    "
                <<refCast<const mappedFvPatchBaseBase>(filmPatch)
                    .nbrMesh().name()
                <<endl;
        }
        Info<<name()<<endl;
    }
    else
    {
        Info<<name()<<" 'coupled' set to false,"
            <<" transfer from film to bulk fluid will not be solved!"
            <<endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::fluidPhaseChangeToFilm::addSupFields() const
{
    wordList fieldNames(0);

    if (!coupled_)
    {
        return fieldNames;
    }

    fieldNames.append(mcf_.thermo.he().name());
    fieldNames.append(rhoName_);
    fieldNames.append(liquidNames());

    return fieldNames;
}

void Foam::fv::fluidPhaseChangeToFilm::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();

    cacheFilmToFluids();
}

void Foam::fv::fluidPhaseChangeToFilm::cacheFilmToFluids()
{
    if (FilmToFluidsCached_ || !coupled_)
    {
        return;
    }

    Info<<name()<<": Caching filmPhaseChangeToFluid models "
        <<"from coupled neighbour regions"
        <<endl;

    // Determine all coupled models
    
    List<const filmPhaseChangeToFluid*> ptrList(filmPatchiList_.size());

    label n_models = 0;

    forAll(filmPatchiList_, patchesi)
    {
        // initialise with nullptr
        ptrList[patchesi] = nullptr;
        
        const filmPhaseChangeToFluid* filmToFluidPtr = nullptr;

        const fvPatch& filmPatch
        (
            mesh().boundary()[filmPatchiList_[patchesi]]
        );

        const mappedFvPatchBaseBase& filmPatchMap
        (
            refCast<const mappedFvPatchBaseBase>(filmPatch)
        );

        const fvModels& filmFvModels
        (
            fvModels::New(filmPatchMap.nbrMesh())
        );

        forAll(filmFvModels, modeli)
        {
            if (isType<filmPhaseChangeToFluid>(filmFvModels[modeli]))
            {
                const filmPhaseChangeToFluid& filmToFluid
                (
                    refCast<const filmPhaseChangeToFluid>(filmFvModels[modeli])
                );

                if
                (
                    filmToFluid.nbrPatchIndex()
                        ==
                    filmPatch.index()
                )
                {
                    filmToFluidPtr = &filmToFluid;
                }
            }
        }
        
        if (!filmToFluidPtr)
        {
            if (requireNbrModel_)
            {
                FatalErrorInFunction
                    << "Cannot find filmPhaseChangeToFluid fvModel "
                    << "in the film region " 
                    << filmPatchMap.nbrMesh().name()
                    << " via patch "
                    << filmPatch.name()
                    << exit(FatalError);
            }
            else
            {
                Info<<name()<<": Did not find filmPhaseChangeToFluid "
                    <<"fvModel in the film region "
                    << filmPatchMap.nbrMesh().name()
                    << " via patch "
                    << filmPatch.name()
                    << ". Ignoring. (set 'requireNbrModel' to 'true' "
                    << "to force error)"
                    <<endl;
            }
            continue;
        }

        ptrList[patchesi] = filmToFluidPtr;
        n_models++;
    }

    FilmToFluidsCachedList_.clear();
    FilmToFluidsCachedList_.transfer(ptrList);
    FilmToFluidsCached_ = true;

    Info<<name()<<": Cached "<<n_models<<" models."<<endl;
}

const Foam::List<const Foam::fv::filmPhaseChangeToFluid*>&
  Foam::fv::fluidPhaseChangeToFilm::FilmToFluids() const
{
    return FilmToFluidsCachedList_;
}

const Foam::wordList Foam::fv::fluidPhaseChangeToFilm::liquidNames() const
{
    wordList names(0);

    const List<const filmPhaseChangeToFluid*>& ptrList = FilmToFluids();

    forAll(ptrList, i)
    {
        if (ptrList[i])
        {
            names.append(ptrList[i]->liquidName());
        }
    }

    return names;
}

Foam::tmp<Foam::volScalarField> Foam::fv::fluidPhaseChangeToFilm::mapFilmField
(
    const filmPhaseChangeToFluid& filmToFluid,
    const volScalarField& nbrField
) const
{
    // Determine the correct mapped patch

    label coupledPatchi = -1;

    forAll(filmPatchiList_, patchesi)
    {
        const fvPatch& filmPatch
        (
            mesh().boundary()[filmPatchiList_[patchesi]]
        );

        if 
        (
            filmToFluid.nbrPatchIndex()
                ==
            filmPatch.index()
        )
        {
            coupledPatchi = filmPatchiList_[patchesi];
            break;
        }
    }

    if (coupledPatchi < 0)
    {
        FatalErrorInFunction
            << "Could not determine the patch mapping to the "
            << "given filmPhaseChangeToFluid model! (No patch "
            << "index match found!)"
            << exit(FatalError);
    }

    const fvPatch& coupledPatch
    (
        mesh().boundary()[coupledPatchi]
    );

    const mappedFvPatchBaseBase& coupledPatchMap
    (
        refCast<const mappedFvPatchBaseBase>(coupledPatch)
    );

    // Create new empty field (shape of this mesh)
    tmp<volScalarField> tField =
        volScalarField::New
        (
            "tField",
            mesh(),
            dimensionedScalar(nbrField.dimensions(), Zero)
        );

    // Get neighbour field
    tmp<scalarField> mappedField =
        coupledPatchMap.fromNeighbour
        (
            nbrField
        );

    const labelList& faceCells = coupledPatch.faceCells();

    // Map into this mesh
    forAll (faceCells, facei)
    {
        const label celli = faceCells[facei];
        tField.ref().primitiveFieldRef()[celli] += mappedField.ref()[facei];
    }

    return tField;
}

void Foam::fv::fluidPhaseChangeToFilm::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (true)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
        //Info<< eqn.psi().name() << endl;
        //Info<< mcf_.thermo.rho()().name() << endl;
    }

    if (rho.name() == rhoName_)
    {
        // Explicit transfer into fluid
        volScalarField::Internal tSu
        (
            volScalarField::Internal::New
            (
                "SuRho",
                mesh(),
                dimensionedScalar(dimMass, Zero)
            )
        );

        const List<const filmPhaseChangeToFluid*>& filmToFluids = FilmToFluids();

        // Add all boundary sources into the temporary source field
        forAll(filmToFluids, modeli)
        {
            if (filmToFluids[modeli])
            {
                tSu += mapFilmField
                (
                    *filmToFluids[modeli],
                    filmToFluids[modeli]->rhoTransV().ref()
                ).ref().internalField();
            }
        }

        tSu /= mesh().time().deltaT();

        //eqn += fvm::Su(tSu, eqn.psi());
        eqn.source() -= tSu;

        return;
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << eqn.psi().name() << " is not implemented"
            << exit(FatalError);
    }
}

void Foam::fv::fluidPhaseChangeToFilm::addSup
(
    const volScalarField& rho,
    const volScalarField& heOrYi,
    fvMatrix<scalar>& eqn
) const
{
    if (true)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }
   
    //- Energy source
    if (&heOrYi == &mcf_.thermo.he() && &eqn.psi() == &mcf_.thermo.he())
    {
        // Explicit transfer into fluid
        volScalarField::Internal tSu
        (
            volScalarField::Internal::New
            (
                "SuHe",
                mesh(),
                dimensionedScalar(dimEnergy, Zero)
            )
        );

        const List<const filmPhaseChangeToFluid*>& filmToFluids = FilmToFluids();

        // Add all boundary sources into the temporary source field
        forAll(filmToFluids, modeli)
        {
            if (filmToFluids[modeli])
            {
                tSu += mapFilmField
                (
                    *filmToFluids[modeli],
                    filmToFluids[modeli]->hTransV().ref()
                ).ref().internalField();
            }
        }

        tSu /= mesh().time().deltaT();

        // Latent heat (film phase change)
        eqn.source() -= tSu;

        //eqn += fvm::Su(tSu, eqn.psi());

        return;
    }

    //- Mass fraction source
    else
    {
        // Iterate all coupled models and find the models with the 
        // corresponding liquid specie
        
        bool oneFound = false;

        const List<const filmPhaseChangeToFluid*>& filmToFluids = FilmToFluids();

        forAll(filmToFluids, modeli)
        {
            if (!filmToFluids[modeli])
            {
                continue;
            }

            const filmPhaseChangeToFluid& filmToFluid
            (
                *FilmToFluids()[modeli]
            );

            if 
            (
                heOrYi.name() == filmToFluid.liquidName()
                   && 
                eqn.psi().name() == filmToFluid.liquidName()
            )
            {
                // Explicit transfer into fluid
                volScalarField::Internal tSu
                (
                    volScalarField::Internal::New
                    (
                        "SuYi",
                        mapFilmField
                        (
                            filmToFluid,
                            filmToFluid.rhoTransV().ref()
                        ).ref().internalField()
                    )
                );

                tSu /= mesh().time().deltaT();

                //eqn += fvm::Su(tSu, eqn.psi());
                eqn.source() -= tSu;

                oneFound = true;

                Info<<"    from film region '"
                    <<filmToFluid.regionName()<<"'"
                    <<endl;
                continue;
            }
        }

        // Return to avoid the error below
        if (oneFound)
        {
            return;
        }
    }
    
    {
        FatalErrorInFunction
            << "Support for field " << heOrYi.name() << " is not implemented"
            << exit(FatalError);
    }
}

void Foam::fv::fluidPhaseChangeToFilm::topoChange(const polyTopoChangeMap&)
{
    FilmToFluidsCached_ = false;
}

void Foam::fv::fluidPhaseChangeToFilm::mapMesh(const polyMeshMap&)
{
    FilmToFluidsCached_ = false;
}

void Foam::fv::fluidPhaseChangeToFilm::distribute(const polyDistributionMap&)
{
    FilmToFluidsCached_ = false;
}

bool Foam::fv::fluidPhaseChangeToFilm::movePoints()
{
    return true;
}


// ************************************************************************* //

