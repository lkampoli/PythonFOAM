/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "PykOmegaSSTGEP.H"
#include "fvOptions.H"
#include "bound.H"

#include "wallDist.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
//#include "evalGeno.H"
template<class BasicMomentumTransportModel>
tmp<volScalarField>
PykOmegaSSTGEP<BasicMomentumTransportModel>::PykOmegaSSTGEP::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
PykOmegaSSTGEP<BasicMomentumTransportModel>::PykOmegaSSTGEP::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
PykOmegaSSTGEP<BasicMomentumTransportModel>::PykOmegaSSTGEP::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
PykOmegaSSTGEP<BasicMomentumTransportModel>::PykOmegaSSTGEP::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicMomentumTransportModel>
void PykOmegaSSTGEP<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    ////BasicMomentumTransportModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void PykOmegaSSTGEP<BasicMomentumTransportModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> PykOmegaSSTGEP<BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> PykOmegaSSTGEP<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> PykOmegaSSTGEP<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> PykOmegaSSTGEP<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> PykOmegaSSTGEP<BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
PykOmegaSSTGEP<BasicMomentumTransportModel>::PykOmegaSSTGEP
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
	
	y_(wallDist::New(this->mesh_).y()),
	
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    // This may not be necessary but maybe usefull at some point
    // for example for plotting purposes ...
    cluster_
    (
        IOobject
        (
            IOobject::groupName("cluster", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ, // NO_READ = 0, MUST_READ = 1, MUST_READ_IF_MODIFIED = 3, READ_IF_PRESENT
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("cluster",dimensionSet(0,0,0,0,0,0,0), 0)
    ),

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --start
////////////////////////////////////////////////////////////////////////    

     nonlinearStress_
    (
        IOobject
        (
            IOobject::groupName("nonlinearStress", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	//0.0*symm(fvc::grad(this->U_)/this->omega_)
        this->mesh_,
        dimensionedSymmTensor
        (
            "nonlinearStress",
            sqr(dimVelocity),
            Zero
        )
    ),
//    bij_
//    (
//        IOobject
//        (
//            IOobject::groupName("bij", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//      //0.0*symm(fvc::grad(this->U_)/this->omega_)
//      this->mesh_
//    ),
     Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
	//0.0*symm(fvc::grad(this->U_))/this->omega_
        this->mesh_,
        dimensionedSymmTensor
        (
            "Rij",
            sqr(dimVelocity),
            Zero
        )
    ),
     Rall_
    (
        IOobject
        (
            IOobject::groupName("Rall", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	//((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + 2*k_*aijx_
	//0.0*( ((2.0/3.0)*I)*this->k_ - this->nut_*twoSymm(fvc::grad(this->U_)) )
        this->mesh_,
        dimensionedSymmTensor
        (
            "Rall",
            sqr(dimVelocity),
            Zero
        )
    )
//    T1
//    (
//        IOobject
//        (
//            "T1",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T2
//    (
//        IOobject
//        (
//            "T2",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T3
//    (
//        IOobject
//        (
//            "T3",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T4
//    (
//        IOobject
//        (
//            "T4",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T5
//    (
//        IOobject
//        (
//            "T5",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T6
//    (
//        IOobject
//        (
//            "T6",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T7
//    (
//        IOobject
//        (
//            "T7",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T8
//    (
//        IOobject
//        (
//            "T8",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T9
//    (
//        IOobject
//        (
//            "T9",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    T10
//    (
//        IOobject
//        (
//            "T10",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    I1
//    (
//      IOobject
//      (
//            IOobject::groupName("I1", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//    ),
//    I2
//    (
//      IOobject
//      (
//            IOobject::groupName("I2", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//    ),
//    I3
//    (
//      IOobject
//      (
//            IOobject::groupName("I3", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//    ),
//    I4
//    (
//      IOobject
//      (
//            IOobject::groupName("I4", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//    ),
//    I5
//    (
//      IOobject
//      (
//            IOobject::groupName("I5", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//    ),
//    Eta1
//    (
//        IOobject
//        (
//            "Eta1",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//      //0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    Eta2
//    (
//        IOobject
//        (
//            "Eta2",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//        //0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    Eta3
//    (
//        IOobject
//        (
//            "Eta3",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//        //0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    Eta4
//    (
//        IOobject
//        (
//            "Eta4",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//        //0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    Eta5
//    (
//        IOobject
//        (
//            "Eta5",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//        //0.0*symm(fvc::grad(this->U_))/omega_
//    ),
//    Sijt 
//    (
//        IOobject
//        (
//            "Sijt",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//        //0.0*symm(fvc::grad(this->U_))/omega_
//        //0.5*symm(gradU + gradU.T())
//    ),
//    Oijt 
//    (
//        IOobject
//        (
//            "Oijt",
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//	//-0.5*(gradU - gradU.T())
//    )


//    xswtch_
//    (
//        IOobject
//        (
//            IOobject::groupName("xswtch", alphaRhoPhi.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        this->mesh_
//    )

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --end
////////////////////////////////////////////////////////////////////////    

{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    // This is where the Python interface is instantiated //
    //#include "PythonCreate.H"
    
    Info<< "initialize python" << endl;
    Py_Initialize();
    Info<< "python initialize successful" << endl;

    Info<< "import sys" << endl;
    PyRun_SimpleString("import sys");

    Info<< "sys.path.append" << endl;
    PyRun_SimpleString("sys.path.append(\".\")");

    // initialize numpy array library
    // init_numpy();
    //Info<< "numpy initialize successful" << endl;
    // if (PyErr_Occurred()) {
    //     Info<< "Failed to import numpy Python module(s)." << endl;
    //     return NULL; // Or some suitable return value to indicate failure.
    // }
    // assert(PyArray_API);
    import_array1();

    pName = PyUnicode_DecodeFSDefault("python_module"); // Python filename

    pModule = PyImport_Import(pName);
    // Py_DECREF(pName);
    
    //if (!pModule)
    //{
    //    FatalErrorInFunction
    //    << "Errors loading python_module (missing imports?)" << nl
    //    << exit(FatalError);
    //}

    // The function "clustering" should be present in the "python_module.py"
    // at run-time, otherwise you will get an error. At compile-time is ok
    // not to have it.
    Info<< " get clustering python function" << endl;
    clustering = PyObject_GetAttrString(pModule, "clustering");

    Info<< " define args for clustering" << endl;
    clustering_args = PyTuple_New(2); // Two arrays sent to Python

    // Get process id
    rank = Pstream::myProcNo();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --start
////////////////////////////////////////////////////////////////////////

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
PykOmegaSSTGEP<BasicMomentumTransportModel>::sigma() const
{
    tmp<volSymmTensorField> tR
    (
        eddyViscosity<RASModel<BasicMomentumTransportModel>>::sigma()
    );
    tR.ref() += nonlinearStress_;
    return tR;
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
PykOmegaSSTGEP<BasicMomentumTransportModel>::devTau() const
{
    tmp<volSymmTensorField> tdevTau
    (
        eddyViscosity<RASModel<BasicMomentumTransportModel>>::devTau()
    );
    tdevTau.ref() += this->rho_*nonlinearStress_;
    return tdevTau;
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
PykOmegaSSTGEP<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    return
    (
        fvc::div(this->rho_*nonlinearStress_) +
        eddyViscosity<RASModel<BasicMomentumTransportModel>>::divDevTau(U)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
PykOmegaSSTGEP<BasicMomentumTransportModel>::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (
        fvc::div(rho*nonlinearStress_) +
        eddyViscosity<RASModel<BasicMomentumTransportModel>>::divDevTau(rho, U)
    );
}

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --end
////////////////////////////////////////////////////////////////////////

template<class BasicMomentumTransportModel>
bool PykOmegaSSTGEP<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void PykOmegaSSTGEP<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();
	
    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );
    dimensionedScalar Nutsmall(
        "0", dimensionSet(0,2,-1,0,0,0,0), 1e-20
    );
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));

    // Moved down because they depend on variables computed before ...
    //volScalarField::Internal GbyNu((dev(twoSymm(tgradU()()))) && tgradU()());
    //volScalarField::Internal GbyNuaijx(((dev(twoSymm(tgradU()()))) - (nonlinearStress_/nut())) && tgradU()());
    //volScalarField::Internal G(this->GName(), nut()*GbyNuaijx);
    
////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --start
////////////////////////////////////////////////////////////////////////

    volScalarField S = sqrt(S2);
	
    //volScalarField tau = 1./max( S/a1_ + omegaMin_,omega_ + omegaMin_);
    volScalarField tau =  1./max( S/0.31 + this->omegaMin_,(omega_ + this->omegaMin_));
    volScalarField tau2 = sqr(tau);
	
    volSymmTensorField sij(dev(symm(tgradU()))); 
    volTensorField omegaij((skew(tgradU())));

    volSymmTensorField Sijt = sij * tau;
    volTensorField     Oijt = omegaij * tau;

    volScalarField I1 = tr(Sijt & Sijt);
    volScalarField I2 = tr(Oijt & Oijt);
    volScalarField I3 = tr(Sijt & (Sijt & Sijt));
    volScalarField I4 = tr(Oijt & (Oijt & Sijt));
    volScalarField I5 = tr(Oijt & (Oijt & (Sijt & Sijt)));

    volSymmTensorField T1 = Sijt;
    volSymmTensorField T2 = symm((Sijt & Oijt) - (Oijt & Sijt));
    volSymmTensorField T3 = symm( Sijt & Sijt) - scalar(1.0/3.0)*I*I1;
    volSymmTensorField T4 = symm( Oijt & Oijt) - scalar(1.0/3.0)*I*I2;
    volSymmTensorField T5 = symm((Oijt & (Sijt & Sijt)) - ((Sijt & Sijt) & Oijt));
    volSymmTensorField T6 = symm(((Oijt & Oijt) & Sijt) + ((Sijt & Sijt) & Oijt) - (2.0/3.0) * I * tr(Sijt & (Oijt & Oijt)));
    volSymmTensorField T7 = symm((Oijt & (Sijt & (Oijt & Oijt))) - (Oijt & (Oijt & (Sijt & Oijt))));
    volSymmTensorField T8 = symm((Sijt & (Oijt & (Sijt & Sijt))) - (Sijt & (Sijt & (Oijt & Sijt))));
    volSymmTensorField T9 = symm((Oijt & (Oijt & (Sijt & Sijt))) + (Sijt & (Sijt & (Oijt & Oijt))) - (2.0/3.0) * I * tr(Sijt & (Sijt & (Oijt & Oijt))));
    volSymmTensorField T10= symm((Oijt & (Sijt & (Sijt & (Oijt & Oijt)))) -(Oijt & (Oijt & (Sijt & (Sijt & Oijt)))));
    
    // imported from kOmegaSSTGEP_Rijhat_added.C
    // This is just for compatibility, with nonLinearModel.H file
    word model (
      this->coeffDict_.lookup ("model")
    );

    word aijx_term (
      this->coeffDict_.lookup("aijx_term")
    );

    word Rterm_corr (
      this->coeffDict_.lookup("Rterm_corr")
    );
    /////////////////////////////////////////////

    // library of GEP models
    if (model == "linear") {
	nonlinearStress_ =  2*k_ *(T1 - T1);
    }
    #include "nonLinearModel.H"
    else {
	FatalError << "unknown model " << model << nl << exit(FatalError);
    }


////////////////////////////////////////////////////////////////////////
//  added by Lorenzo --start
//  compute all derived variables to be passed to python
//  and then call it and get the cluster id back to be used
//  to select the GEP model in each cluster
////////////////////////////////////////////////////////////////////////

  if (this->runTime_.outputTime())
  {
    Info<< "Time = " << this->runTime_.timeName() << nl << endl;
    Info<< "Calculating fields for ML" << nl << endl;

    dimensionedScalar kSMALL("0",dimLength*dimLength/dimTime/dimTime, 1e-10);
    dimensionedScalar rSMALL ("0", dimensionSet(0,2,-2,0,0,0,0),1e-10);
    dimensionedScalar USMALL ("0", dimensionSet(0,1,-1,0,0,0,0),1e-10);
    dimensionedScalar qSMALL ("0", dimensionSet(0,0,-1,0,0,0,0),1e-10);
    dimensionedScalar NewSMALL("0", dimensionSet(0,0,0,0,0,0,0), 1e-10);
    dimensionedScalar constantInRe("0", dimensionSet(0,0,0,0,0,0,0), 2);
    dimensionedScalar osmall ("0",dimensionSet(0,0,-1,0,0,0,0),1e-10);
    dimensionedScalar SMALL_CONVDER( "small", dimensionSet(0,1,-2,0,0,0,0), scalar(1E-30) );

    //volTensorField     gradU = fvc::grad(U);

    //volSymmTensorField S_ = symm(fvc::grad(U));
    //volTensorField     R_ = -skew(fvc::grad(U)); 
    
    // grad(U) produces the transpose of the Jacobian, R is defined based on the Jacobian
    //gradU_ = fvc::grad(U);
    //gradU_ = gradU_.T(); // Output the Jacobian, a more common form of the velocity gradient tensor
    //gradk_ = fvc::grad(k_);
    //DUDt_  = U & gradU_.T();
    //gradp_ = fvc::grad(p_);


//    // Q1
//    volScalarField Q1_org =  0.5*(Foam::sqr(tr(gradU)) - tr(((gradU) & (gradU))));
//    volScalarField Q1 (
//    IOobject (
//        "Q1",
//        this->runTime_.timeName(),
//        this->mesh_,
//        IOobject::NO_READ
//        ),
//        (Q1_org-Foam::min(Q1_org))/(Foam::max(Q1_org)-Foam::min(Q1_org)),
//        "zeroGradient"
//    );
//    Info << "--> Q1_org Min:" << Foam::min(Q1_org).value() << " Max: " << Foam::max(Q1_org).value() << endl;
//    Info << "--> Q1     Min:" << Foam::min(Q1).value() << " Max: " << Foam::max(Q1).value() << endl;
//
//
//    // Q2
//    volScalarField Q2 (
//    IOobject (
//        "Q2",
//        this->runTime_.timeName(),
//        this->mesh_,
//        IOobject::NO_READ
//        ),
//        k_/ (0.5* (U&U) + k_ + kSMALL),
//        "zeroGradient"
//    );
//    Info << "--> Q2    Min:" << Foam::min(Q2).value() << " Max :" << Foam::max(Q2).value() << endl;
//
//
//
//    // nu
//    Info<< "    Reading transport Properties" <<endl;
//    IOdictionary transportProperties
//    (
//    IOobject
//            (
//                "transportProperties",
//                this->runTime_.constant(),
//                this->mesh_,
//                IOobject::MUST_READ,
//                IOobject::NO_WRITE,
//                false
//            )
//    );
//    dimensionedScalar nu (transportProperties.lookup("nu"));
//    Info << "--> nu  :" << nu.value() << endl;
//
//    //    // Q3
////    // wall-distance based Reynolds number
////    volScalarField Q3 (
////    IOobject (
////      "Q3",
////        this->runTime_.timeName(),
////        this->mesh_,
////        IOobject::NO_READ
////        ),
////        mag(y_ * Foam::sqrt(k_+kSMALL) / nut())/10000,
////        "zeroGradient"
////    );
////    Info << "--> Q3    Min:" << Foam::min(Q3).value() << " Max :" << Foam::max(Q3).value() << endl;
//
//
//    // Q4
//    volScalarField Q4_org = U & gradp_;
//    volScalarField Q4 (
//    IOobject (
//        "Q4",
//        this->runTime_.timeName(),
//        this->mesh_,
//        IOobject::NO_READ
//        ),
//        (Q4_org-Foam::min(Q4_org))/(Foam::max(Q4_org)-Foam::min(Q4_org)),
//        "zeroGradient"
//     );
//     Info << "--> Q4_org Min:" << Foam::min(Q4_org).value() << " Max: " << Foam::max(Q4_org).value() << endl;
//     Info << "--> Min Q4:" << Foam::min(Q4).value() << "Max :" << Foam::max(Q4).value() << endl;
//
//     // Q5
//     volScalarField Q5_org = mag(fvc::curl(U));
//     volScalarField Q5 (
//     IOobject (
//        "Q5",
//        this->runTime_.timeName(),
//        this->mesh_,
//        IOobject::NO_READ
//        ),
//        (Q5_org-Foam::min(Q5_org))/(Foam::max(Q5_org)-Foam::min(Q5_org)),
//        "zeroGradient"
//     );
//     Info << "--> Q5_org Min:" << Foam::min(Q5_org).value() << " Max: " << Foam::max(Q5_org).value() << endl;
//     Info << "--> Min Q5:" << Foam::min(Q5).value() <<  "Max :" << Foam::max(Q5).value() << endl;
//
//
//          // Q6
////     volScalarField Q6 (
////            IOobject (
////                "Q6",
////                this->runTime_.timeName(),
////                this->mesh_,
////                IOobject::NO_READ
////            ),
////            mut / (rho * nu*100 +mut),
////            "zeroGradient"
////   );
////   Info << "--> Q6    Min:" << Foam::min(Q6).value() << " Max :" << Foam::max(Q6).value() << endl;
//
//
//
//     // Q7
//     volScalarField Q7_org = Foam::sqrt(gradp_ & gradp_);
//     volScalarField Q7 (
//     IOobject (
//        "Q7",
//        this->runTime_.timeName(),
//        this->mesh_,
//        IOobject::NO_READ
//        ),
//        (Q7_org-Foam::min(Q7_org))/(Foam::max(Q7_org)-Foam::min(Q7_org)),
//        "zeroGradient"
//     );
//     Info << "--> Q7_org Min:" << Foam::min(Q7_org).value() << " Max: " << Foam::max(Q7_org).value() << endl;
//     Info << "--> Min/Max Q7:" << Foam::min(Q7).value() << Foam::max(Q7).value() << endl;
//
//
//     // Q8
//     volScalarField Q8_org = mag((U * U) && gradU);
//     volScalarField Q8 (
//            IOobject (
//                "Q8",
//                this->runTime_.timeName(),
//                this->mesh_,
//                IOobject::NO_READ
//            ),
//            (Q8_org-Foam::min(Q8_org))/(Foam::max(Q8_org)-Foam::min(Q8_org)),
//            "zeroGradient"
//     );
//     Info << "--> Q8_org Min:" << Foam::min(Q8_org).value() << " Max: " << Foam::max(Q8_org).value() << endl;
//     Info << "--> Min/Max Q8:" << Foam::min(Q8).value() <<  Foam::max(Q8).value() << endl;
//
//
//      // Q9
//     volScalarField Q9_org = (fvc::grad(k_) && U);
//     volScalarField Q9 (
//            IOobject (
//                "Q9",
//                this->runTime_.timeName(),
//                this->mesh_,
//                IOobject::NO_READ
//            ),
//            (Q9_org-Foam::min(Q9_org))/(Foam::max(Q9_org)-Foam::min(Q9_org)),
//           "zeroGradient"
//     );
//     Info << "--> Q9_org Min:" << Foam::min(Q9_org).value() << " Max: " << Foam::max(Q9_org).value() << endl;
//     Info << "--> Min/Max Q9:" << Foam::min(Q9).value() << Foam::max(Q9).value() << endl;
//
//
//     //     volScalarField GQ10 (
////            IOobject (
////                "Q10",
////                this->runTime_.timeName(),
////                this->mesh_,
////                IOobject::NO_READ
////            ),
////            Foam::sqrt(R && R)/(k + kSMALL + Foam::sqrt(R && R)),
////            "zeroGradient"
////     );
////     Info << "--> Q10    Min:" << Foam::min(Q10).value() << " Max :" << Foam::max(Q10).value() << endl;

    // Cell centroid coordinates
    const volVectorField& C = this->mesh_.C();
    //volScalarField x_ = pos(this->mesh_.C().component(vector::X));
    //volScalarField y_ = pos(this->mesh_.C().component(vector::Y));
    //volScalarField z_ = pos(this->mesh_.C().component(vector::Z));

    // Define an array of doubles to pass to python module
    int num_cells = this->mesh_.cells().size();
    int num_points = this->mesh_.nPoints();
    Info << " num_cells = " << num_cells << endl;
    Info << " num_points = " << num_points << endl;
    
    // the second index should change according to the number of fields 
    // to be passed to python
    //double input_vals[num_cells][41];
    double input_vals[num_cells][2];
    //double input_vals[num_cells][5];

    // ClusterId flag is used to group all points in one initial cluster
    // which by default has value is set to 0
    int ClusterId[num_cells];
    forAll(k_.internalField(), id)
    {
        ClusterId[id] = 0;
    }
    forAll(k_.boundaryField(), id)
    {
        ClusterId[id] = 0;
    }
    Info << " ClusterId initialized " << endl;

    forAll(k_.internalField(), id) // for boundary field use u_.boundaryField()
    {
        input_vals[id][0] = 0;
        input_vals[id][1] = 0.0;
    }

    // #pragma omp parallel for
    forAll(k_.internalField(), id) // for boundary field use u_.boundaryField()
    {

        input_vals[id][0] = ClusterId[id]; // cluster id label
	input_vals[id][1] = k_[id]; 	   // TKE
//        input_vals[id][1] = C[id].x();	   // x
//        input_vals[id][2] = C[id].y();     // y
//        input_vals[id][3] = C[id].z();     // z
//	input_vals[id][4] = k_[id]; 	   // TKE
//
//        input_vals[id][4] = U[id].x();
//        input_vals[id][5] = U[id].y();
//        input_vals[id][6] = U[id].z();
//
//        input_vals[id][7] = k_[id];
//        input_vals[id][8] = omega_[id];
//
//        input_vals[id][9]  = gradU_[id].xx();
//        input_vals[id][10] = gradU_[id].xy();
//        input_vals[id][11] = gradU_[id].xz();
//        input_vals[id][12] = gradU_[id].yx();
//        input_vals[id][13] = gradU_[id].yy();
//        input_vals[id][14] = gradU_[id].yz();
//        input_vals[id][15] = gradU_[id].zx();
//        input_vals[id][16] = gradU_[id].zy();
//        input_vals[id][17] = gradU_[id].zz();
//
//        input_vals[id][18] = S_[id].xx();
//        input_vals[id][19] = S_[id].xy();
//        input_vals[id][20] = S_[id].xz();
//        input_vals[id][21] = S_[id].yy();
//        input_vals[id][22] = S_[id].yz();
//        input_vals[id][23] = S_[id].zz();
//
//        input_vals[id][24] = R_[id].xx();
//        input_vals[id][25] = R_[id].xy();
//        input_vals[id][26] = R_[id].xz();
//        input_vals[id][27] = R_[id].yx();
//        input_vals[id][28] = R_[id].yy();
//        input_vals[id][29] = R_[id].yz();
//        input_vals[id][30] = R_[id].zx();
//        input_vals[id][31] = R_[id].zy();
//        input_vals[id][32] = R_[id].zz();
//
//        input_vals[id][33] = gradk_[id].x();
//        input_vals[id][34] = gradk_[id].y();
//        input_vals[id][35] = gradk_[id].z();
//
//        input_vals[id][36] = DUDt_[id].x();
//        input_vals[id][37] = DUDt_[id].y();
//        input_vals[id][38] = DUDt_[id].z();
//
//        input_vals[id][39] = Q1[id];
//        input_vals[id][40] = Q8[id];

	  Info << input_vals[id][0] << " " << input_vals[id][1] << endl;
    }
    Info << " input_vals initialized " << endl;

    // Numpy array dimensions
    //npy_intp dim[] = {num_cells, 41};
    npy_intp dim[] = {num_cells, 2};
    //npy_intp dim[] = {num_cells, 5};

    // create a new array using 'buffer'
    array_2d = PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals[0]);
    PyTuple_SetItem(clustering_args, 0, array_2d);
    Info << " create a new array using 'buffer' " << endl;

    // To pass rank to Python interpreter
    rank_val = PyLong_FromLong(rank);
    PyTuple_SetItem(clustering_args, 1, rank_val);
    Info << " rank defined " << endl;

    // Casting to PyArrayObject
    pValue = (PyArrayObject*)PyObject_CallObject(clustering, clustering_args);
    Info << " Casting to PyArrayObject successful " << endl;

    // if you were defining a new array_2d with PyObject *array_2d = ..; here
    // PyArray_ENABLEFLAGS((PyArrayObject*)array_2d, NPY_ARRAY_OWNDATA);

    // Retrieve the number of clusters ////////////////////////////////
    int clust = 0;
    clust = *(int*) PyArray_GETPTR2(pValue, 0, 0); // row 0, column 0
    std::cout<<" How many clusters? " << clust << std::endl;
    ///////////////////////////////////////////////////////////////////

    // Get data back from the internal field
    forAll(k_.internalField(), id) // for boundary field use u_.boundaryField()
    {
        //double* current = (double*) PyArray_GETPTR2(pValue, id, 0);    // row id, column 0

        ClusterId[id] = *(int*) PyArray_GETPTR2(pValue, id, 1); // row id, column 1
        std::cout<<" Internal Cell " << id << " belongs to cluster " << ClusterId[id] << std::endl;

        //ux_[id] = *(double*) PyArray_GETPTR2(pValue, id,             0); // row id, column 0
        //uy_[id] = *(double*) PyArray_GETPTR2(pValue, id+num_cells,   1); // row id, column 1
        //uz_[id] = *(double*) PyArray_GETPTR2(pValue, id+2*num_cells, 2); // row id, column 2

        //ClusterId[id] = *(int*) PyArray_GETPTR2(pValue, id+2*num_cells, 2); // row id, column 2

    }


    // Get data back from the boundary field
//    forAll(k_.boundaryField(), id)
//    {
//        ClusterId[id] = *(int*) PyArray_GETPTR2(pValue, id, 0);
//        std::cout<<" Boundary Cell " << id << " belongs to cluster " << ClusterId[id] << std::endl;
//    }


    // https://www.cfd-online.com/Forums/openfoam-programming-development/95274-writing-output-simple-data-into-ascii-file.html
    // write coordinate points of clusters
    fileName outputFile0("cluster0.txt");
    fileName outputFile1("cluster1.txt");
    fileName outputFile2("cluster2.txt");
    fileName outputFile3("cluster3.txt");
    fileName outputFile4("cluster4.txt");
    fileName outputFile5("outsiders.txt");

    OFstream os0(outputFile0);
    OFstream os1(outputFile1);
    OFstream os2(outputFile2);
    OFstream os3(outputFile3);
    OFstream os4(outputFile4);
    OFstream os5(outputFile5);
    
    // use clusterId to assign GEP models
    forAll(nonlinearStress_,celli)
   {
      // so, here I hardcode the max (5) number of clusters possibly available
      // this part is not ok, since the number of cluster should be also
      // retrieved from python but as the number of clusters is also
      // usually user defined in the clustering algorithm at the moment
      // we suppose to have this info ... TODO: make number cluster retrival automatic
      Info << " ClusterId  " << ClusterId[celli] << endl;
      if (ClusterId[celli] == 0)
      {
              Info<< "\n For ClusterId 0 apply model 0 \n" << endl;
	      nonlinearStress_[celli] = 2 * k_[celli] * (0.018245*I3[celli]*I5[celli]*T1[celli] + I4[celli]*T2[celli] + I5[celli]*T5[celli]);
              os0 << C[celli].x() << "   " << C[celli].y(); 
	      os0 << endl;
      }
      else if (ClusterId[celli] == 1)
      {
              Info<< "\n For ClusterId 1 apply model 1 \n" << endl;
	      nonlinearStress_[celli] = 2 * k_[celli] * (-1.0*I2[celli]*I4[celli]*T3[celli] + I4[celli]*T4[celli] - 0.097*I5[celli]*T1[celli]*(0.03585*I5[celli] - 0.0225));
              os1 << C[celli].x() << "   " << C[celli].y(); 
	      os1 << endl;
      }
      else if (ClusterId[celli] == 2)
      {
              Info<< "\n For ClusterId 2 apply model 2 \n" << endl;
	      nonlinearStress_[celli] = 2 * k_[celli] * (-1.0*I1[celli]*I3[celli]*I4[celli]*T2[celli]
                                                          - 1.0*I3[celli]*I5[celli]*T5[celli]
                                                          - 0.01335*T3[celli] - 0.01455*T4[celli]);
              os2 << C[celli].x() << "   " << C[celli].y(); 
	      os2 << endl;
      }
      else if (ClusterId[celli] == 3)
      {
              Info<< "\n For ClusterId 3 apply model 3 \n" << endl;
	      nonlinearStress_[celli] = 2 * k_[celli] * (I4[celli]*T1[celli]*(0.01453815*I3[celli] - 0.0010099275));
              os3 << C[celli].x() << "   " << C[celli].y(); 
	      os3 << endl;
      }
      else if (ClusterId[celli] == 4)
      {
              Info<< "\n For ClusterId 4 apply model 4 \n" << endl;
	      nonlinearStress_[celli] = 2 * k_[celli] * (-I1[celli]*I4[celli]*T1[celli]*(0.15*I2[celli]*I2[celli]*I2[celli] + I3[celli])
                                + 2.0*I3[celli]*I4[celli]*T4[celli] + I4[celli]*T2[celli]*(I1[celli] - 0.097)
                                + I4[celli]*T3[celli] + I5[celli]*T5[celli]);
              os4 << C[celli].x() << "   " << C[celli].y(); 
	      os4 << endl;
      }
      else
      {
              Info<< "\n More clusters than GEP models ... applying baseline model to them! \n" << endl;
	      nonlinearStress_ =  2*k_ *(T1 - T1);
              os5 << C[celli].x() << "   " << C[celli].y(); 
	      os5 << endl;
      }

      // to visualize clusters on plot
      cluster_[celli] = ClusterId[celli];
      //int y = (int)x;

   }

//   fileName outputFile(this->runTime_.path()/this->runTime_.timeName());
//   OFstream os(outputFile);
//   os << "This is the first line in the file.\n" << endl;
//   os << "scalarField area (" << areaField.size() << ";)" << endl;

  }

////////////////////////////////////////////////////////////////////////
//  added by Lorenzo --end
////////////////////////////////////////////////////////////////////////

    // dimensional anisotropy tensor term
    nonlinearStress_.correctBoundaryConditions();

    // non-dimensional anisotropy tensor term
    //bij_ = nonlinearStress_ / (2 * k_);

    Rall_= ((2.0/3.0)*I)*k_ - this->nut_*dev(twoSymm(tgradU())) + nonlinearStress_;
    Rall_.correctBoundaryConditions();
 
    volScalarField Rterm(Rij_ && symm(tgradU()));	
    //tgradU.clear();

    dimensionedScalar small_val
    (
       "small_val",
       dimensionSet(0, 2, -1, 0, 0, 0, 0),
       1e-25
    );
    
    volScalarField::Internal GbyNu((dev(twoSymm(tgradU()()))) && tgradU()());
    volScalarField::Internal GbyNuaijx(((dev(twoSymm(tgradU()()))) - (nonlinearStress_/nut())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNuaijx);

    // additional invariants involving anisotropy tensor a_ij
    volScalarField Eta1 = tr(nonlinearStress_ & nonlinearStress_);
    volScalarField Eta2 = tr(nonlinearStress_ & Sijt);
    volScalarField Eta3 = tr(nonlinearStress_ & Sijt & Oijt);
    volScalarField Eta4 = tr(nonlinearStress_ & Sijt & Sijt);
    volScalarField Eta5 = tr(nonlinearStress_ & Oijt & Oijt);

////////////////////////////////////////////////////////////////////////
//  added by Yuan Fang --end
////////////////////////////////////////////////////////////////////////

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNuaijx,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          + alpha()*rho()*gamma*Rterm/(nut()+small_val)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

     	// Turbulent kinetic energy equation
     	tmp<fvScalarMatrix> kEqn
     	(
     	    fvm::ddt(alpha, rho, k_)
     	  + fvm::div(alphaRhoPhi, k_)
     	  - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     	 ==
     	    alpha()*rho()*Pk(G)
     	  + alpha()*rho()*Rterm
     	  - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
     	  - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_)
     	  + kSource()
     	  + fvOptions(alpha, rho, k_)
     	);

     	kEqn.ref().relax();
     	fvOptions.constrain(kEqn.ref());
     	solve(kEqn);
     	fvOptions.correct(k_);
     	bound(k_, this->kMin_);

     	correctNut(S2, F23);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace RASModels
// ************************************************************************* //
