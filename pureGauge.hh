#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
// custom
#include "auxiliary.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//
// NAVIGATON
// generalized row-major format is used for indexing the sites and links instead of a 5-dimensional array
//
// index function for easier understanding and moving around in the flattened array
// taking into account different possible lattice extents and different possible matrix sizes too
// (fixed to 4 spacetime dimensions)
//
// get index from coordinates
int GetIndex(vector<int> const &coordinates, vector<int> const &dims)
{
    // check if mu > 4
    if ((coordinates[4] > 3) || (coordinates[4] < 0))
    {
        cout << "ERROR in GetIndex: mu > 3 or mu < 0" << endl;
        exit(-1);
    }
    // using periodic boundary conditions
    vector<int> redefined_coordinates = coordinates;
    for (int i = 0; i < 4; i++)
    {
        // addig an extra term because C++ is stupid...
        redefined_coordinates[i] = (coordinates[i] + dims[i]) % dims[i];
    }
    // flattened index
    int flattenedIndex = redefined_coordinates[0];
    // strides ~ bigger steps
    int stride = 1;
    for (int i = 1; i < 5; i++)
    {
        stride *= dims[i - 1];
        flattenedIndex += redefined_coordinates[i] * stride;
    }

    // return flattened index
    return flattenedIndex;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// get coordinates and spacetime direction from index ~ (x, y, z, t, mu)
vector<int> GetCoordinates(int const &index, vector<int> const &dims)
{
    // array of coordinates and spacetime direction
    vector<int> coordinates = {0, 0, 0, 0, 0};
    // offset
    int offset = index;
    for (int i = 0; i < 5; i++)
    {
        coordinates[i] = offset % dims[i];
        offset /= dims[i];
    }

    // return coordinates and spacetime direction
    return coordinates;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//
// MONTE CARLO
//

// upper staple (i.e. with given nu) at given site with given mu
// generate arguments of links
vector<Eigen::MatrixXi> StapleUpper_args(vector<int> const &coordinates, int const &nu)
{
    // given direction mu
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do nothing
    Eigen::RowVectorXi coordinates_and_ifDagger(6);
    coordinates_and_ifDagger << coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], 1;

    // used link coordinate combinations in the rectangular staples
    Eigen::MatrixXi initArguments(coordinates_and_ifDagger.colwise().replicate(3));

    // used link coordinate combinations in the staple
    vector<Eigen::MatrixXi> arguments(1, initArguments);

    // arguments
    // U_nu(n + mu)
    arguments[0](0, mu) += 1;
    arguments[0](0, 4) = nu;
    // U_mu(n + nu)^dagger
    arguments[0](1, nu) += 1;
    arguments[0](1, 5) = 0;
    // U_nu(n)^dagger
    arguments[0](2, 4) = nu;
    arguments[0](2, 5) = 0;

    // return arguments
    return arguments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// upper staple (i.e. with given nu) at given site with given mu
// generate arguments of links
vector<Eigen::MatrixXi> StapleLower_args(vector<int> const &coordinates, int const &nu)
{
    // given direction mu
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do nothing
    Eigen::RowVectorXi coordinates_and_ifDagger(6);
    coordinates_and_ifDagger << coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], 1;

    // used link coordinate combinations in the rectangular staples
    Eigen::MatrixXi initArguments(coordinates_and_ifDagger.colwise().replicate(3));

    // used link coordinate combinations in the staple
    vector<Eigen::MatrixXi> arguments(1, initArguments);

    // arguments
    // U_nu(n + mu - nu)^dagger
    arguments[0](0, mu) += 1;
    arguments[0](0, nu) -= 1;
    arguments[0](0, 4) = nu;
    arguments[0](0, 5) = 0;
    // U_mu(n - nu)^dagger
    arguments[0](1, nu) -= 1;
    arguments[0](1, 5) = 0;
    // U_nu(n - nu)
    arguments[0](2, nu) -= 1;
    arguments[0](2, 4) = nu;

    // return arguments
    return arguments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// sum of general staples from sets of link arguments
Eigen::Matrix3cd StapleSum(vector<Eigen::Matrix3cd> const &links, vector<Eigen::MatrixXi> const &arguments, vector<int> const &dims)
{
    // initialize sum of staples
    Eigen::Matrix3cd stapleSum = Eigen::Matrix3cd::Zero();
    // auxiliary matrix
    Eigen::Matrix3cd tmpMat = ID3;
    // number of staples in the sum
    int NumStaples = static_cast<int>(arguments.size());
    // number of links in a staple
    int NumLinks = arguments[0].rows();
    // loop over staples
    for (int staple = 0; staple < NumStaples; staple++)
    {
        // construct staple
        for (int link = 0; link < NumLinks; link++)
        {
            // adjoint or not?
            if (arguments[staple](link, 5) == 0)
                tmpMat *= links[GetIndex(eigen_to_std(arguments[staple].row(link)), dims)].adjoint();
            else
                tmpMat *= links[GetIndex(eigen_to_std(arguments[staple].row(link)), dims)];
        }
        // add up staples
        stapleSum += tmpMat;
        // reset auxiliary matrix
        tmpMat = ID3;
    }
    // return sum of staples
    return stapleSum;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// generating SU(2) matrix
Eigen::Matrix2cd SU2(double const &eps, function<double()> const &randFunc)
{
    // generating four random numbers from [-1/2, 1/2]
    vector<double> r(4, 0);
    for (int i = 0; i < 4; i++)
        r[i] = randFunc() - 0.5;

    // length of 3-vector
    double norm = sqrt(sq(r[0]) + sq(r[1]) + sq(r[2]));

    // generate coefficients
    vector<double> x(4, 0);
    for (int i = 0; i < 3; i++)
        x[i] = eps * r[i] / norm;
    x[3] = ((r[3] > 0.) - (r[3] < 0)) * sqrt(1. - sq(eps));

    // parameterizing SU(2) matrices
    Eigen::Matrix2cd X2 = x[3] * ID2;
    for (int i = 0; i < 3; i++)
        X2 += I * Pauli[i] * x[i];

    // return the SU(2) matrix
    return X2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// candidate link generation through updating matrix X
Eigen::Matrix3cd UpdatingMatrix(double const &eps, function<double()> randFunc)
{
    // generating three SU(2) matrices
    vector<Eigen::Matrix2cd> SU2Matrices(3, Eigen::Matrix2cd::Zero());
    for (int i = 0; i < 3; i++)
        SU2Matrices[i] = SU2(eps, randFunc);

    // generating three 3x3 matrices ~ SU(2) subgropups method
    Eigen::Matrix3cd R{{SU2Matrices[0](0, 0), SU2Matrices[0](0, 1), complex(0., 0.)},
                       {SU2Matrices[0](1, 0), SU2Matrices[0](1, 1), complex(0., 0.)},
                       {complex(0., 0.), complex(0., 0.), complex(1., 0.)}};

    Eigen::Matrix3cd S{{SU2Matrices[1](0, 0), complex(0., 0.), SU2Matrices[1](0, 1)},
                       {complex(0., 0.), complex(1., 0.), complex(0., 0.)},
                       {SU2Matrices[1](1, 0), complex(0., 0.), SU2Matrices[1](1, 1)}};

    Eigen::Matrix3cd T{{complex(1., 0.), complex(0., 0.), complex(0., 0.)},
                       {complex(0., 0.), SU2Matrices[2](0, 0), SU2Matrices[2](0, 1)},
                       {complex(0., 0.), SU2Matrices[2](1, 0), SU2Matrices[2](1, 1)}};

    // return the X matrix
    if (randFunc() > 0.5)
        return R * S * T;
    else
        return (R * S * T).adjoint();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// local change in Wilson action (single site)
double DeltaGaugeAction(double const &beta, vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, Eigen::Matrix3cd const &U, Eigen::Matrix3cd const &X)
{
    // Lorentz index of the given link
    int mu = coordinates[4];

    // sum of staples
    // Eigen::Matrix3cd SumOfStaples = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd A = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd B = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd AUpper = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd ALower = Eigen::Matrix3cd::Zero();

    for (int nu = 0; nu < 4; nu++)
    {
        if (nu == mu)
            continue;
        else
        {
            // A += StapleSum(links, coordinates, dims, nu);
            AUpper += StapleSum(links, StapleUpper_args(coordinates, nu), dims);
            ALower += StapleSum(links, StapleLower_args(coordinates, nu), dims);
        }
    }
    // sum of lower and upper plaquette staples
    A = AUpper + ALower;

    // differences upon proposal in the plaquette and the improved rectangular terms
    Eigen::Matrix3cd diff = (X - ID3) * U * A;

    // return local change in action
    return -beta / 3. * diff.real().trace();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// acceptance rate for the Metropolis algorithm
double Rate(double const &deltaAction)
{
    if (deltaAction <= 0)
        return 1.;
    else
        return exp(-deltaAction);
}

#pragma warning(pop)