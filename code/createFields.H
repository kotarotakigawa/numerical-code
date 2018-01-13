Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
const volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
const volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*
Info<< "Reading/calculating face flux field Uphi\n" << endl;

const surfaceScalarField Uphi
(
    IOobject
    (
        "Uphi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf()
);
*/
Info<< "Reading field p\n" << endl;
const volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field rho\n" << endl;
const volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Creating MW\n" << endl;
const volScalarField PMmix
(
    IOobject
    (
        "meanMolecularWeight",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
/*
Info<< "Reading face flux field phi\n" << endl;

const surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(U*rho) & mesh.Sf()
);
*/
Info<< "Reading mu\n" << endl;
const volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

