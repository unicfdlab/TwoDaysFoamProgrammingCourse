/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    SimplexNelderMead

Description
    Simplex algorithmn for minimisation

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "SimplexNelderMead.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::SimplexNelderMead<Func>::SimplexNelderMead
(
    Func& f,
    const scalarField& p0,
    const scalarField& lambda
)
:
    f_(f),
    points_(f_.nArgs() + 1)
{
    for(label rI = 0; rI < f_.nArgs(); rI++)
    {
        // c equals zero in all apart from the rI-th corner
        // and is offset by p0 over all coordinates
        scalarField c(f_.nArgs(), 0.0);

        c[rI] = lambda[rI];
        c += p0;

        Tuple2<scalar,bool> val = f_(c);
        if (val.second())
        {
            points_.set(rI, new simplexCorner(c, val.first()));
        }
        else
        {
            FatalErrorIn
            (
                "SimplexNelderMead<Func>::SimplexNelderMead\n"
                "(\n"
                "    const Func& f,\n"
                "    const scalarField& p0,\n"
                "    const scalarField& lambda\n"
                ")"
            )   << "failed during initialisation"
                << abort(FatalError);
        }
    }

    Tuple2<scalar, bool> val = f_(p0);

    // Initialise if evaluation was successful
    if (val.second())
    {
        points_.set(f_.nArgs(), new simplexCorner(p0, val.first()));
    }
    else
    {
        FatalErrorIn
        (
            "SimplexNelderMead<Func>::SimplexNelderMead"
            "(const Func& f, const scalarField& p0, const scalarField& lambda)"
        )   << "failed during initialisation"
            << abort(FatalError);
    }


//     for(label rI=0; rI<f_.nArgs()+1; rI++)
//     {
//         Info << "coor = " << points_[rI].coord() << endl;
//     }

    updateOrdering();

//     for(label rI=0; rI<f_.nArgs()+1; rI++)
//     {
//         Info << "coor = " << points_[rI].coord()
//             << " v = " << points_[rI].value() << endl;
//     }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Func>
autoPtr<typename Foam::SimplexNelderMead<Func>::simplexCorner>
Foam::SimplexNelderMead<Func>::newFromWorst
(
    const scalar coeff
)
{
    // Creates a new point by moving a simplex corner scaled by coeff
    // (negative value represents mirroring by the middle point of
    // the "other" corner points)

    scalarField mp(f_.nArgs(), 0.0);

    for(label pI = 0; pI < f_.nArgs(); pI++)
    {
        mp += points_[pI].coord();
    }

    mp /= f_.nArgs();

    scalarField xc = mp - coeff*(mp - points_[f_.nArgs()].coord());

    Tuple2<scalar,bool> val = f_(xc);
    if (val.second())
    {
        return autoPtr<simplexCorner>(new simplexCorner(xc, val.first()));
    }
    else
    {
        return autoPtr<simplexCorner>(NULL);
    }
}


template<class Func>
void Foam::SimplexNelderMead<Func>::contractByBest()
{
    // Contract the simplex in respect to best valued corner. That
    // is, all corners besides the best corner are moved.

    // Starting from 1 - 0 is the best
    for(label pI = 1; pI < f_.nArgs()+1; pI++)
    {
        scalarField newCoord =
            (points_[pI].coord() + minCoord())/2.0;

        Tuple2<scalar,bool> val = f_(newCoord);
        if ( val.second() )
        {
            points_[pI].reset(newCoord, val.first());
        }
        else
        {
            FatalErrorIn
            (
                    "SimplexNelderMead<Func>::contractByBest"
            )   << "cannot contract because of illegal values"
            << abort(FatalError);
        }
    }
}


template<class Func>
scalarField Foam::SimplexNelderMead<Func>::center() const
{
    scalarField center(f_.nArgs(), 0.0);

    for(label pI = 0; pI < f_.nArgs() + 1; pI++)
    {
        center += points_[pI].coord();
    }
    center /= f_.nArgs() + 1;

    return center;
}


template<class Func>
void Foam::SimplexNelderMead<Func>::swapPoints(label i, label j)
{
    simplexCorner tmp = points_[i];
    points_[i] = points_[j];
    points_[j] = tmp;
}


template<class Func>
scalar Foam::SimplexNelderMead<Func>::size()
{
    scalarField c = center();
    scalar ss = 0.0;

    for (label i = 0; i < f_.nArgs()+1; i++)
    {
        ss += sqrt(sumSqr(c - points_[i].coord()));
    }
    ss /= f_.nArgs()+1;

    return ss;
}


template<class Func>
void Foam::SimplexNelderMead<Func>::updateOrdering()
{
    // Sort the points partially by value, so that the smallest value
    // comes first; seconds largest at nArgs() - 1; largest at nArgs()

    // Make sure the first point is smaller than the last
    if (min() > points_[f_.nArgs()].value())
    {
        swapPoints(0,f_.nArgs());
    }

    // Make sure the last is larger than the second last
    if(points_[f_.nArgs() - 1].value() > points_[f_.nArgs()].value())
    {
        swapPoints(f_.nArgs()-1, f_.nArgs());
    }

    for (label i = 1; i < f_.nArgs() - 1; i++)
    {
        if (points_[i].value() < min())
        {
            swapPoints(i,0);
        }

        if (points_[i].value() > points_[f_.nArgs() - 1].value())
        {
            if(points_[i].value() > points_[f_.nArgs()].value())
            {
                swapPoints(f_.nArgs(),f_.nArgs() - 1);
                swapPoints(i,f_.nArgs());
            }
            else
            {
                swapPoints(i,f_.nArgs() - 1);
            }
        }
    }
}


template<class Func>
void Foam::SimplexNelderMead<Func>::iterate()
{
    // Simplex iteration tries to minimize function f value

    // Reflect the highest value

    autoPtr<simplexCorner> np = newFromWorst(-1.0);

    if ( np.valid() && np->value() < min())
    {
        // reflected point becomes lowest point, try expansion

        autoPtr<simplexCorner> np2 = newFromWorst(-2.0);

        if (np2.valid() && np2->value() < min())
        {
            //Info << "reflect 1" << endl;
            points_.set(f_.nArgs(), np2);
        }
        else
        {
            //Info << "reflect 2" << endl;
            points_.set(f_.nArgs(), np);
        }
    }

    // reflection does not improve things enough

    else if (!np.valid() || np->value() > points_[f_.nArgs()-1].value())
    {

        if (np.valid() && np->value() <= points_[f_.nArgs()].value())
        {
            //Info << "reflect 3" << endl;
            points_.set(f_.nArgs(), np);
        }

        // try one dimensional contraction

        autoPtr<simplexCorner> np2 = newFromWorst(0.5);

        if (np2.valid() && np2->value() <= points_[f_.nArgs()].value())
        {
            //Info << "contract" << endl;
            points_.set(f_.nArgs(), np2);
        }
        else
        {
            //contract the whole simplex in respect to the best point

            //Info << "contract by best" << endl;
            contractByBest();
        }
    }
    else
    {
        // trial point is better than second highest point. Replace
        // highest point by it

        //Info << "replace" << endl;
        points_.set(f_.nArgs(), np);
    }

    updateOrdering();

//     for(label rI=0; rI<f_.nArgs()+1; rI++)
//     {
//         Info << "coor = " << points_[rI].coord()
//             << " v = " << points_[rI].value() << endl;
//     }
}


// ************************************************************************* //
