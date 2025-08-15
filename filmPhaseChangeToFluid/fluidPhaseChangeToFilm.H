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
 
Class
    Foam::fv::fluidPhaseChangeToFilm
 
Description
    film <---> multicomponentFluid multi-region phase change model.

    Phase change is handled by the film region counterpart model; this
    model only reads the phase change fields from it, and applies the
    corresponding source terms into the continuity, mass-fraction and
    energy equations.

    Supports multiple coupled film regions.

    Supports one active specie (determined by the film region counterpart 
    model) per coupled film region. In case of multiple films, does support 
    films with different active species (e.g., film1 = IC8H18 and film2 = CH3OH)

    Requires a multi-region setup, where this region is coupled with one or
    more film regions (with the film solver being either 'isothermalFilm'
    or its derivative). The interface patches between the regions must be 
    mapped patches.
 
Usage
    Example usage:

    \verbatim
    fluidPhaseChangeToFilm
    {
        type    fluidPhaseChangeToFilm;
 
                                       //  | required | default value |
        patches ( ... ... );           //  |    yes   |               |

        coupled         true;          //  |    no    |     true      |
        requireNbrModel true;          //  |    no    |     true      |
        rho             rho;           //  |    no    |     rho       |
    }
    \endverbatim
 
SourceFiles
    fluidPhaseChangeToFilm.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef fluidPhaseChangeToFilm_H
#define fluidPhaseChangeToFilm_H
 
#include "fvModel.H"
#include "multicomponentFluid.H"
#include "mappedFvPatchBaseBase.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
namespace fv
{
 
//- Forward declarations
class filmPhaseChangeToFluid;

/*---------------------------------------------------------------------------*\
                  Class fluidPhaseChangeToFilm Declaration
\*---------------------------------------------------------------------------*/
 
class fluidPhaseChangeToFilm
:
    public fvModel
{
private:
    
    // Private Data

        //- Switch to enable/disable transfer between film and fluid
        const Switch coupled_;

        //- Switch to force FatalError if a filmToFluid model cannot be found
        //  in the nbrRegion of any coupled patch
        const Switch requireNbrModel_;

        //- Name of rho field
        const word rhoName_;

        //- The fluid model
        const solvers::multicomponentFluid& mcf_;

        //- Regular expression of coupled patch names
        wordReList filmPatchNamesList_;

        //- List of the coupled patch indices in this region
        labelList filmPatchiList_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;

        //- Cached coupled phase change models. Use FilmToFluids() for access,
        //  instead of directly accessing this variable.
        //
        //  Note: It is assumed that the destruction of the pointers contained
        //  within is handled by the solver ( this is also why smart memory, 
        //  i.e., PtrList<T> was not used here, because it caused segfaults 
        //  with double-destruction )
        List<const filmPhaseChangeToFluid*> FilmToFluidsCachedList_;
        
        //- Is the FilmToFluidsCachedList_ up to date?
        mutable Switch FilmToFluidsCached_;

    // Private Member Functions
        
        //- Get and cache all coupled phase change models
        void cacheFilmToFluids();

        //- Return all coupled phase change models
        const List<const filmPhaseChangeToFluid*>& FilmToFluids() const;

        //- Return the names of the active liquids in the coupled phase change models
        const wordList liquidNames() const;
 
        //- Map a film field 'nbrField' from the film mesh into this mesh
        tmp<volScalarField> mapFilmField
        (
            const filmPhaseChangeToFluid& filmToFluid,
            const volScalarField& nbrField
        ) const;

public:

    //- Runtime type information
    TypeName("fluidPhaseChangeToFilm");

    // Constructors

        fluidPhaseChangeToFilm
	    (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
	    );

    // Member Functions
        
        // Checks
            
            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;
 
            //- Return the indices of coupled patches in this region
            const labelList filmPatchIndices() const
            {
                return filmPatchiList_;
            }

        // Correct

            //- Cache coupled models, if necessary
            virtual void correct();

	    // Sources

            //- Add source term to continuity equation
	        virtual void addSup
	        (
                const volScalarField& rho,
	            fvMatrix<scalar>& eqn
            ) const;

            //- Add source term to energy/mass fraction equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& heOrYi,
                fvMatrix<scalar>& eqn
            ) const;

        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluidPhaseChangeToFilm&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

