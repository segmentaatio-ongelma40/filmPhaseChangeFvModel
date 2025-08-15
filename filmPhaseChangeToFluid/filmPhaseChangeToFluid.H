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
    Foam::fv::filmPhaseChangeToFluid
 
Description
    film <---> multicomponentFluid multi-region phase change model.

    The model is based on OpenFOAM-10 'standardPhaseChange' model.

    Features mass transfer driven convective evaporation and heat transfer
    driven boiling.

    Supports only one active specie (a limitation imposed by the film
    thermodynamic model, which does not support multiple species).

    Requires a multi-region setup, where this (film) region is coupled with
    a neighbouring fluid region (with fluid solver being either 
    'multicomponentFluid' or its derivative). The interface patch
    between the regions must be a mapped patch.

    Reference length LRef is automatically determined to be the film wetted
    length (Lwet = Vwet / Awet). LRef can be manually overriden by setting
    'overrideLRef' to true and providing a value for LRef.

Usage
    Example usage:

    \verbatim
    filmPhaseChangeToFluid
    {
        type          filmPhaseChangeToFluid;
 
                                            // | required | default value |
        activeLiquid  IC8H18;               // |   yes    |               |
        overrideLRef  false;                // |   no     |     false     |
        LRef          0;                    // |   no     |  film length  |
        deltaMin      1e-13;                // |   no     |     1e-13     |
        debug         false;                // |   no     |     false     |
    }
    \endverbatim
 
SourceFiles
    filmPhaseChangeToFluid.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef filmPhaseChangeToFluid_H
#define filmPhaseChangeToFluid_H
 
#include "fvModel.H"
#include "isothermalFilm.H"
#include "liquid.H"
//#include "timeIOdictionary"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
namespace fv
{

//- Forward declarations
class fluidPhaseChangeToFilm;
 
/*---------------------------------------------------------------------------*\
                  Class filmPhaseChangeToFluid Declaration
\*---------------------------------------------------------------------------*/
 
class filmPhaseChangeToFluid
:
    public fvModel
{
private:
    
    // Private Data

        //- The film model
        const solvers::isothermalFilm& film_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;

        //- Mass phase change counters

            //- Mass phase change counter, this time step [kg]
            scalar dmTotal_;

            //- Mass phase change counter, previous time step [kg]
            scalar dmTotalOld_;

            //- Cumulative mass phase change counter [kg]
            scalar dmTotalCumulative_;

        //- Phase change mass flux [kg/m2/s] OUT of the film
        volScalarField::Internal mDot_;

        //- Total phase change latent heat flux [J/m2/s] OUT of the film
        volScalarField::Internal heLatDot_;

        //- Vaporisation mode field [-]
        //  - 2: boiling
        //  - 1: evaporation
        //  - 0: none
        volScalarField::Internal vaporisationMode_;

        //- Film mass, previous time step [kg]
        mutable scalar filmMassOld_;

        //- The liquid, of which the film is composed of
        autoPtr<liquidProperties> activeLiquid_;

        //- Name of the coupled film surface patch
        word filmPatchName_;

        //- Index of the coupled patch in the neighbour region
        const label nbrPatchi_;

        //- Use user-defined LRef instead of Lwet?
        const Switch overrideLRef_;

        //- Characteristic length [m]
        const scalar LRef_;

        //- Film wetted length (volume / wetted area) [m]
        scalar Lwet_;

        //- Minimum delta [m], above which the vaporisation model activates
        //  (default 1e-13)
        const scalar deltaMin_;

        //- Mass phase change (out of the film) field (extensive) [kg]
        volScalarField::Internal dm_;

        //- Latent heat of vaporisation field (extensive) [J/kg],
        //  positive sign = vaporisation drains energy
        volScalarField::Internal heLat_;

        //- Enable additional debug information?
        const Switch debug_;

    // Private Member Functions

        //- Calculate and return Sherwood number
        scalar calculateSh
        (
            const scalar Re, 
            const scalar Sc
        ) const;

        //- Correct the hLat field
        void correctHeLat();

        //- Recalculate the film wetted length
        void correctLwet();

        //- Solve phase change equations and determine vaporisation rate
        void solvePhaseChange();

        //- Calculates and prints film mass
        void printFilmMass() const;

        //- Reference to the coupled model
        const fluidPhaseChangeToFilm& FluidToFilm(const fvModels&) const;

        // TODO: remove
        //- Return a transfer rate from coupledModel
        template<class Type, class TransferRateFunc>
        inline tmp<VolInternalField<Type>> fluidToFilmTransferRate
        (
            TransferRateFunc transferRateFunc,
            const dimensionSet& dimProp
        ) const;
        // TODO

protected:
    
    //- Protected data

        //- Dictionary of output properties (TODO)
        //timeIOdictionary outputProperties_;

public:

    //- Runtime type information
    TypeName("filmPhaseChangeToFluid");

    // Constructors

        filmPhaseChangeToFluid
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

            //- Return the name of the active liquid specie
            const word liquidName() const;

            //- Return the index of the coupled patch in the neighbour region
            label nbrPatchIndex() const
            {
                return nbrPatchi_;
            }

            //- Return the name of the mesh region
            const word regionName() const
            {
                return mesh().name();
            }

        // Correct

            //- Solve phase change rate and update the sources
            virtual void correct();

	    // Sources

            //- Add source term to continuity equation
	        virtual void addSup
	        (
	            const volScalarField& rho,
                const volScalarField& alpha,
	            fvMatrix<scalar>& eqn
            ) const;

            //- Add source term to energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& alpha,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- TODO: momentum source terms?

        // Transfer to fluid
            
            //- Mass transfer (intensive) [kg/m3]
            //tmp<volScalarField> rhoTrans() const;

            //- Mass transfer (extensive) [kg]
            tmp<volScalarField> rhoTransV() const;

            //- Enthalpy transfer (intensive) [J/m3]
            //tmp<volScalarField> hTrans() const;

            //- Enthalpy transfer (extensive) [J]
            tmp<volScalarField> hTransV() const;

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
        void operator=(const filmPhaseChangeToFluid&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

