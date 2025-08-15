Film-multicomponentFluid standard phase change fvModel for OpenFOAM-13
===============

A set of two coupled finite-volume models (fvModels) for OpenFOAM-13,
used for modelling liquid wall film phase-change in multi-region cases.

Port of the *standardPhaseChange* model from OpenFOAM-10.

I created this model originally for my master's thesis, and I'm sharing it here in case somebody 
finds it useful. The code still needs some clean-up, and has several limitations/inaccuracies, as
listed in [Section 5](#5-known-limitations).  

## 1. Technical description

Contains two fvModels, *filmPhaseChangeToFluid* and *fluidPhaseChangeToFilm*, used for modelling liquid wall
film phase change in multi-region cases (examples: the *rivuletBox* and *hotBoxes* tutorial cases).
<ins>Requires a multi-region case, in which the gaseous phase is modelled using *multicomponentFluid* solver,
and the liquid wall film phase is modelled using the *film* solver.
</ins>

In principle, the model functions similar to the [*filmVoFTransfer*](https://cpp.openfoam.org/v13/classFoam_1_1fv_1_1filmVoFTransfer.html) or [*filmCloudTransfer*](https://cpp.openfoam.org/v13/classFoam_1_1fv_1_1filmCloudTransfer.html) models, on which
this model is based upon. The phase change code is based on the OpenFOAM-10 wall film phase change model,
[*standardPhaseChange*](https://cpp.openfoam.org/v10/classFoam_1_1regionModels_1_1surfaceFilmModels_1_1standardPhaseChange.html).

The model considers two modes of wall film phase change: evaporation and boiling.
Evaporation is modelled using convective mass transfer coefficient 
correlations. Boiling is modelled using estimations of available heat
flux. Latent heat of vaporisation is calculated using the enthalpy difference method 
($\Delta h_{lat}=h_a(g)-h_a(l)$)

The *filmPhaseChangeToFluid* fvModel is applied to the film region(s). It calculates the phase change rate
and applies the mass and energy source terms. <ins>Only single-specie films are supported</ins>, limited
to the pre-defined species in OpenFOAM (e.g., CH3OH, IC8H18, H2O). <ins>Requires the film surface to be a 
conformal mapped patch, coupled to the gaseous phase.</ins>

The *filmPhaseChangeToFilm* fvModel is applied to the multicomponentFluid region. It applies the appropriate
mass and energy source terms (equal but opposite to the film source terms). 
Theoretically, it should support an arbitrary number of coupled film regions and different species,
though it has only been extensively tested with a singular coupled film region.



## 2. Installation and compilation

### 2.1 Linux, OpenFOAM compiled from source

1. Ensure that the packages required for compiling OpenFOAM are installed. 
See [OpenFOAM Repo: 1. Software for Compilation](https://openfoam.org/download/source/software-for-compilation/)
for package list.
2. Clone the repository to a directory (e.g., to /home/\<username\>/OpenFOAM/\<username\>-13/dev/).
```
git clone https://github.com/segmentaatio-ongelma40/filmPhaseChangeFvModel.git
cd filmPhaseChangeFvModel/filmPhaseChangeToFluid
```
3. Source the OpenFOAM environment (example using bash)
```
source /path/to/OpenFOAM/etc/bashrc
```
4. Run wmake to compile the model
```
wmake
```

The resulting library, **libfilmPhaseChangeToFluid.so** will be located at the 
directory given by the environment variable **$FOAM_USER_LIBBIN**.

### 2.2 Linux, OpenFOAM installed from repository (e.g., using apt)

**TODO**
I do not have personal experience about this, but I assume that
the instructions of [Section 2.1](#21-linux-compiled-openfoam)
should apply here, if the necessary packages for compilation are installed.

### 2.3 Linux, high-performance computing clusters

Depending on the cluster, packages necessary for compilation may be available as loadable modules.
In that case, load them and follow the instructions in 
[Section 2.1](#21-linux-compiled-openfoam) or [Section 2.2](#22-linux-installed-openfoam-eg-from-apt-repositories).


### 2.4 Windows (WSL)

**TODO**
if using WSL (Windows Subsystem for Linux) to run OpenFOAM, 
[Section 2.1](#21-linux-compiled-openfoam) and [Section 2.2](#22-linux-installed-openfoam-eg-from-apt-repositories)
should apply here.

## 3. Usage

The models are declared in the fvModels file (located at
*/constant/\<regionname\>/fvModels*). 
See the test case for examples.
The library (libfilmPhaseChangeToFluid.so) needs to be included via the 'libs' argument.

Example declarations:
```
filmPhaseChangeToFluid
{
    type    filmPhaseChangeToFluid;
    libs    ("libfilmPhaseChangeToFluid.so");

                                        // | required | default value |
    activeLiquid  IC8H18;               // |   yes    |               |
    overrideLRef  false;                // |   no     |     false     |
    LRef          0;                    // |   no     |  film length  |
    deltaMin      1e-13;                // |   no     |     1e-13     |
    debug         false;                // |   no     |     false     |
}
```

```
fluidPhaseChangeToFilm
{
    type    fluidPhaseChangeToFilm;
    libs    ("libfilmPhaseChangeToFluid.so");

                                   //  | required | default value |
    patches ( ... ... );           //  |    yes   |               |

    coupled         true;          //  |    no    |     true      |
    requireNbrModel true;          //  |    no    |     true      |
}
```


Available options:

**filmPhaseChangeToFluid**

Option | required | default value | description
------ | -------- | ------------- | -----------
activeLiquid | yes |    | The liquid species the film consists of (e.g., IC8H18, CH3OH, H2O). <ins>Currently only supports single-specie liquids</ins>.
overrideLRef | no | false | Switch. If set true, states that an user-defined *L* in evaporation *Re* and *Sh* equations should be used, defined in the 'LRef' entry. If set false, *L* will be calculated from $L = V_{film} / A_{film}$ on each time-step.
LRef | no |             | If 'overrideLRef' is set true, use this value for *L*. 
deltaMin | no | 1e-13   | Minimum film thickness $\delta$, beyond which the phase change model activates. Cells with film thickness less than this are not considered wetted, and phase change is not calculated for them.
debug | no | false | Switch for printing increased debug information.


**fluidPhaseChangeToFilm**

Option | required | default value | description
------ | -------- | ------------- | -----------
patches | yes |    | List of patch name(s), corresponding to the name(s) of the coupled film surface patch(es).
coupled | no  | true | Switch. If set false, no source terms are applied to this region.
requireNbrModel | no | true | Switch. If set true, forces an error if no *filmPhaseChangeToFluid* model is found in one or more of the coupled film regions defined by the 'patches' entry. If set false, the model will log warning and ignore the region.

### Log information
The model logs:
* mass change (phase change), this time-step and cumulative
* mass transfer coefficient $h_m$ and spalding mass transfer number $B_m$(labelled as mass transfer driving force).
* film mass, delta, and ratio. Delta refers to change in film mass between this and previous time-step. Ratio refers to ratio of calculated phase change to film mass change. The ratio may deviate from 1.0 as a result of droplet absorption, film mass out-flow, or film separation (if modelled), leading to additional film mass gain/loss outside of phase change. This was originally logged for debug purposes - may be deprecated in the future.


## 4. Test cases

**TODO**

### Test Case I: Heated plate injection (Jüngst et al., 2021)

![aaa](/assets/testCase1/combined.png)
*Figure: Mesh, spray, and films visualisations*

![aaa](/assets/testCase1/resolution_comparison.png)
*Figure: Mesh resolution (base cell size) comparisons*

* 100 bar injection of iso-octane against a 300 K heated plate.
* Initial mass over-predicted possibly due to spray parameter choices.
* Agreeable evaporation rate results at medium grid resolution. Under-prediction observed
at finer resolutions.
* A relatively heavy case; should run fine on 32-128 CPUs.
* Jüngst N., Frapolli N., Wright Y.M., Boulouchos K., Kaiser S.A., 2021.
"Experimental and numerical investigation of evaporating fuel films in 
combustion". Applications in Energy and Combustion Science. 7. [doi.org/10.1016/j.jaecs.2021.100033](https://doi.org/10.1016/j.jaecs.2021.100033)

**TODO: Include other test cases** 

## 5. Known limitations

* The model does not consider momentum source terms. This is arguably not
accurate at high rates of phase change.
* The evaporation mode is known to be grid-sensitive. At high resolutions,
the evaporation rate is under-predicted, most likely due to decrease
in film-gas $U$ difference in the first wall cell layer.
Likewise, the evaporation rate is over-predicted at low resolutions.
This would imply that the accuracy is dependent on $y^+$.
* The boiling mode does not consider high temperature boiling regimes,
nor the associated reduction in boiling rate due to the Leidenfrost effect.
* The boiling mode has not yet been validated against an experimental 
study.

## 6. Future plans

* Clean-up the code and improve log output. Remove/hide debug stuff.
* Modify the code to utilise submodels for phase change, allowing for different models to be more robustly implemented.
* Implement uniform time dictionary for keeping track of phase change statistics (similar to
the *liquidEvaporation* model).
* Implement cell volume correction for source term transfer, and transfer the source terms in 
intensive form, instead of extensive form.
* Implement momentum source terms.
* Remove/deprecate overrideLRef option.
* Investigate improvements to the grid-sensitivity of the evaporation mode.
* Investigate improvements to the boiling mode.
* Include a validation case for the boiling mode and test it.


