// OFFICE HOURS: 3-4 THURSDAY
//
//  Monte Carlo Ising Model Simulator
//  
//
//  Created by Jack Brown on 10/04/17.
//
//
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include <ctime>
#include <input.h>

class lattice {
    public: std::vector< std::vector<int>> spins;
    public: int energy = 0;
    public: double magnetization = 0;
    public: double m_sqrd = 0;
    public: double m_avg = 0;
    public: double variance = 0;

    void init (int x, int y) {
        spins.resize(y);
        for(int i =0; i<x; ++i) {
            spins[i].resize(x);
            std::fill (spins[i].begin(), spins[i].end(), 1);
        }
    }

    void init_config (int x, int y, int run) {
        spins.resize(y);
        int configs = std::pow(2, x*y);
        int config_number = run % configs;

        std::vector < int > config_vector;
        config_vector.resize(configs);
        std::fill(config_vector.begin(), config_vector.end(), 1);
        for(int change =0; change<=config_number; ++change) {
            for (int &i : config_vector) {
                if (i == 1) {
                    i = -1;
                }
                else {
                    i = 1;
                    break;
                }
            }
        }
        for(int j = 0; j<y; ++j){
            spins[j].resize(x);
            for(int i = 0; i<x; ++i){
                spins[j][i] = config_vector[j*x+i];
            }
        }
    }

    void eval_energy(int x, int y) {
        energy = 0;
        for(int i = 0; i < x - 1; i++) {
            for(int j = 0; j < x - 1; j++) {
                energy = energy + spins[i][j]*spins[i][j+1] + spins[i][j+1]*spins[i+1][j+1]
                    + spins[i+1][j+1]*spins[i+1][j] + spins[i+1][j]*spins[i][j];
            }
        }
    }

    void eval_magnetization(int x, int y){
        magnetization = 0;
        for(int i = 1; i < x - 1; i++) {
            for(int j = 1; j < x - 1; j++) {
                magnetization += spins[i][j];
            }
        }
        m_sqrd = pow((magnetization / (x*y)), 2.0);
        m_avg = magnetization/(x*y);
    }


};

int main() {
    std::ofstream txtoutput;
    txtoutput.open ("snapshots_spins_up.txt");
    int num_runs = 1;
    int step = 0;


    for(int run=0; run < num_runs; run++) {
        lattice grid;
        int size = 10;
        grid.init(size, size);

        std::vector<double> m_sqrd_avg_arr;
        double beta = 0.25;
        int num_steps = 500;
        srand(run);

        for (int i = 0; i < num_steps; i++) {
            int rnd_col = rand() % size;
            int rnd_row = rand() % size;

            grid.eval_energy(size, size);
            int E_before = grid.energy;
            if (grid.spins[rnd_row][rnd_col] == -1) {
                grid.spins[rnd_row][rnd_col] = 1;
                grid.eval_energy(size, size);
                if (exp(-1 * beta * (grid.energy - E_before)) <
                    static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) {
                    grid.spins[rnd_row][rnd_col] = -1;
//                    std::cout << "Reject" << std::endl;
                    //Reject
                }
//                std::cout << "Accept" << std::endl;
                //Accept
                step = step+1;
            } else {
                grid.spins[rnd_row][rnd_col] = -1;
                grid.eval_energy(size, size);
                if (exp(-1 * beta * (grid.energy - E_before)) <
                    static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) {
                    grid.spins[rnd_row][rnd_col] = 1;
//                    std::cout << "Reject" << std::endl;
                    //Reject
                }
//                std::cout << "Accept" << std::endl;
                //Accept

                step = step+1;
            }
            grid.eval_energy(size, size);
            grid.eval_magnetization(size, size);
            m_sqrd_avg_arr.push_back(grid.m_sqrd);
            // txtoutput << ',' << exp(-1*beta*grid.energy);
            for(const auto& row : grid.spins) {
                std::copy(row.cbegin(), row.cend(), std::ostream_iterator<int>(txtoutput, " "));
                txtoutput << '\n';
            }
            txtoutput << ',';
        }
        // if(run != num_runs-1) {txtoutput << std::endl;}

        // avg_magnetization/lattice_area
        double m_sqrd_avg = std::accumulate(m_sqrd_avg_arr.begin(), m_sqrd_avg_arr.end(), 0.000)/m_sqrd_avg_arr.size();

    }
    txtoutput.close();
    return 0;

}

