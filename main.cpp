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

class lattice {
    public: std::vector< std::vector<int>> spins;
    public: int energy = 0;
    public: double magnetization = 0;
    public: double m_sqrd = 0;
    public: double m_avg = 0;
    public: double variance = 0;

    void init (int x, int y) {
        spins.resize(x);
        for(int i =0; i<x; ++i) {
            spins[i].resize(y);
            std::fill (spins[i].begin(), spins[i].end(), 1);
        }
    }
    void eval_energy(int x, int y) {
        energy = 0;
        for(int i = 1; i < x - 1; i++) {
            for(int j = 1; j < x - 1; j++) {
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
    lattice grid;
    int size = 10;
    grid.init(size, size);
    std::vector <double> m_sqrd_avg_arr;
    int step = 0;
    double beta = 0.25;
    int num_steps = 4000;
    srand(time(NULL));
    for (int i = 0; i < num_steps; i++) {
            int rnd_col = rand() % size;
            int rnd_row = rand() % size;

            grid.eval_energy(size, size);
            int E_before = grid.energy;
            if(grid.spins[rnd_row][rnd_col] == -1) {
                grid.spins[rnd_row][rnd_col] = 1;
                grid.eval_energy(size, size);
                if(exp(-1*beta*(grid.energy - E_before)) < static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) {
                    grid.spins[rnd_row][rnd_col] = -1;
                    // std::cout << "\nReject\n";
                }
                //else {std::cout << "\nAccept\n";}
                step++;
            }
            else {
                grid.spins[rnd_row][rnd_col] = -1;
                grid.eval_energy(size, size);
                if(exp(-1*beta*(grid.energy - E_before)) < static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) {
                    grid.spins[rnd_row][rnd_col] = 1;
                    //std::cout << "\nReject";
                }
                //else {std::cout << "\nAccept";}

                step++;
            }
//        for (int j = 0; j < size; ++j) {
//            for (int i = 0; i < size; ++i) {
//                std::cout << '\t' << grid.spins[j][i];
//            }
//            std::cout << '\n';
//        }
        grid.eval_energy(size, size);
        grid.eval_magnetization(size, size);
        m_sqrd_avg_arr.push_back(grid.m_sqrd);
        //std::cout << '\n' << grid.magnetization << '\n';
    }
    // avg_magnetization/lattice_area
    double m_sqrd_avg = std::accumulate(m_sqrd_avg_arr.begin(), m_sqrd_avg_arr.end(), 0.000)/m_sqrd_avg_arr.size();
    std::cout << m_sqrd_avg_arr.size() << std::endl;
//std::pow(std::accumulate(m_sqrd_avg_arr.begin(), m_sqrd_avg_arr.end(), 0) /m_sqrd_avg_arr.size(), 2) / std::pow(size, 4);
    // avg_magnetization
//    float m_avg_sqrd = (std::accumulate(m_vals.begin(), m_vals.end(), 0));
    std::cout << m_sqrd_avg;
    return 0;

}

