#pragma once

#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
#include <random>

// custom
#include "auxiliary.hh"
#include "pureGauge.hh"
#include "observables.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// main function
int main(int, char **argv)
{
    // lattice sizes
    const int Ns = atoi(argv[1]);
    const int Nt = atoi(argv[2]);
    // number of links for 4 spacetime dimensions
    const int NSite = cb(Ns) * Nt;
    const int NLink = NSite * 4;

    // coupling
    const double beta = atof(argv[3]);
    // spread parameter in the Metropolis algorithm
    // probably can be set to ~0.25
    const double eps = atof(argv[4]);
    // number of MC sweeps
    const int T = atoi(argv[5]);
    // number of discarded configurations between measurements
    const int tau = atoi(argv[6]);

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------

    // random number generation
    std::random_device rd{};
    std::mt19937 gen(rd());
    // uniform random from [0, 1)
    std::uniform_real_distribution<double> uniformDistr(0., 1.);
    auto RandomNumber = [&]()
    {
        return uniformDistr(gen);
    };
    // uniform random from [0, Ns)
    std::uniform_int_distribution uniformNsDistr(0, Ns - 1);
    auto RandomNsNumber = [&]()
    {
        return uniformNsDistr(gen);
    };
    // uniform random from [0, Nt)
    std::uniform_int_distribution uniformNtDistr(0, Nt - 1);
    auto RandomNtNumber = [&]()
    {
        return uniformNtDistr(gen);
    };
    // unform random from [0, 4)
    std::uniform_int_distribution uniformDirDistr(0, 3);
    auto RandomDirNumber = [&]()
    {
        return uniformDirDistr(gen);
    };

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------

    // dimensions
    vector<int> dims = {Ns, Ns, Ns, Nt, 4};
    // initial coordinates
    vector<int> coordinates = {0, 0, 0, 0, 0};

    // cold start
    vector<Eigen::Matrix3cd> links(NLink, ID3);

    // RUN simulation
    // SWEEPS
    for (int t = 0; t < T; t++)
    {
        // single sweep via Metropolis steps over the whole lattice
        for (int step = 0; step < NLink; step++)
        {
            // Metropolis step
            //
            // choosing a site and a direction randomly
            for (int coord = 0; coord < 3; coord++)
                coordinates[coord] = RandomNsNumber();
            coordinates[3] = RandomNtNumber();
            coordinates[4] = RandomDirNumber();

            // updating matrix X
            Eigen::Matrix3cd X = UpdatingMatrix(eps, RandomNumber);

            // current link U
            Eigen::Matrix3cd U = links[GetIndex(coordinates, dims)];

            // local change in the Wilson action DeltaS
            double deltaGaugeAction = DeltaGaugeAction(beta, links, coordinates, dims, U, X);
            // rate
            double rate = Rate(deltaGaugeAction);
            // random number from [0, 1)
            double r = RandomNumber();

            // decide if new link is accepted or not
            if (r < rate)
                links[GetIndex(coordinates, dims)] = X * U;
        }

        // MEASUREMENTS
        if ((t % tau) == 0)
        {
            // TASK
        }
    }
}

#pragma warning(pop)