/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    Y18Na1Foam0523s

Description
    A solver for vapor-phase flame synthesis using Gauss-Radau CQMOM.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "particleConstants.H"
#include "particleDatas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "createCQMOM.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
/*
    for(label iNode=0; iNode<nNode; ++iNode)
    {
        forAll(abscissav[iNode], cid)
        {
            if(abscissav[iNode][cid] <= v_mol || abscissaa[iNode][cid] <= a_mol)
            {
                shape[iNode][cid] = 1.;
            }
            else
            {
                shape[iNode][cid] = 3*(1-Foam::log(abscissaa[iNode][cid]/a_mol)/Foam::log(abscissav[iNode][cid]/v_mol));
            }
        }
        shape[iNode].correctBoundaryConditions();
    }
*/
    for(label iNode=0; iNode<nNode; ++iNode)
    {
        forAll(mesh.cells(), cid)
        {

            if(abscissav[iNode/Nf][cid] <= v_mol || abscissaa[iNode][cid] <= a_mol)
            {
                shape[iNode][cid] = 1.;
            }
            else
            {
                shape[iNode][cid] = 3*(1-Foam::log(abscissaa[iNode][cid]/a_mol)/(Foam::log(abscissav[iNode/Nf][cid]/v_mol)+SMALL));
            }
        }
        shape[iNode].correctBoundaryConditions();
    }


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
		#include "nodeEqn.H"

		#include "particleCellLoop.H"
/*
        for(label iNode=0; iNode<nNode; ++iNode)
        {
            abscissav[iNode].correctBoundaryConditions();
            abscissaa[iNode].correctBoundaryConditions();
            abscissaf[iNode].correctBoundaryConditions();
            weightv[iNode].correctBoundaryConditions();
            weightf[iNode].correctBoundaryConditions();

            forAll(abscissav[iNode], cid)
            {
                if(abscissav[iNode][cid] <= v_mol || abscissaa[iNode][cid] <= a_mol)
                {
                    shape[iNode][cid] = 1.;
                }
                else
                {
                    shape[iNode][cid] = 3*(1-Foam::log(abscissaa[iNode][cid]/a_mol)/Foam::log(abscissav[iNode][cid]/v_mol));
                }
            }

            shape[iNode].correctBoundaryConditions();
        }
        */
        for(label iNode=0; iNode<Nv; ++iNode)
        {
            abscissav[iNode].correctBoundaryConditions();
            weightv[iNode].correctBoundaryConditions();
        }

        for(label iNode=0; iNode<nNode; ++iNode)
        {
            abscissaa[iNode].correctBoundaryConditions();
            abscissaf[iNode].correctBoundaryConditions();
            weightf[iNode].correctBoundaryConditions();

            forAll(mesh.cells(), cid)
            {
                if(abscissav[iNode/Nf][cid] <= v_mol || abscissaa[iNode][cid] <= a_mol)
                {
                    shape[iNode][cid] = 1.;
                }
                else
                {
                    shape[iNode][cid] = 3*(1-Foam::log(abscissaa[iNode][cid]/a_mol)/(Foam::log(abscissav[iNode/Nf][cid]/v_mol)+SMALL));
                }
            }
            shape[iNode].correctBoundaryConditions();
        }
        Info << MOMv[0][1000] << endl;
        Info << abscissav[0][1000] << endl;
        Info << weightv[0][1000] << endl;
        Info << weightf[0][1000] << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
