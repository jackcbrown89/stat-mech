//  Monte Carlo Ising Model Simulator
//  
//
//  Created by Jack Brown on 10/04/17.
//
//
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>
#include <ctime>
#include "input.h"
#include "Communication.h"

class lattice {
    public: std::vector< std::vector<int> > spins;
    public: int energy = 0;
    public: double magnetization = 0;
    public: double m_sqrd = 0;
    public: double m_avg = 0;
    public: double variance = 0;
    public: double J = 1.0;
    public: double beta = 0.25;
    public: int const_count = 0;
    public: double T = 7.25e+23;

    void init (int x, int y) {
        spins.resize(y);
        for(int i =0; i<x; ++i) {
            spins[i].resize(x);
            std::fill (spins[i].begin(), spins[i].end(), 1);
        }
    }

    void init_existing (vector < vector <int> > a) {
        spins = a;
    }

    void init_config (int x, int y, unsigned long long run) {
        spins.resize(y);
        double configs = pow(2, static_cast<float> (x*y));
        unsigned long long config_number = run % static_cast<unsigned long long> (configs); // identify config type
        cout << run << endl;
        std::vector < int > config_vector;
        config_vector.resize(x*y);
        std::fill(config_vector.begin(), config_vector.end(), 1);
        for(int change = 0; change<=config_number; ++change) {
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

    void init_rand_config (int x, int y) {
        srand(13);
        spins.resize(y);
        for(int i = 0; i < y; i++){
            spins[i].resize(x);
            fill(spins[i].begin(), spins[i].end(), 1);
        }
        int r_changes = rand() % (x*y);
        std::cout << r_changes << endl;
        for(int j=0; j< r_changes; j++) {
            spins[rand() % y][rand() % x] = -1;
        }
    }

    void eval_energy() {
        energy = 0;
        for(int i = 0; i < spins.size() - 1; i++) {
            for(int j = 0; j < spins.size() - 1; j++) {
                energy = energy + spins[i][j]*spins[i][j+1] + spins[i][j+1]*spins[i+1][j+1]
                    + spins[i+1][j+1]*spins[i+1][j] + spins[i+1][j]*spins[i][j];
            }
        }
    }

    int energy_parallel(int i, int j) {
        int energy_p = 0;
        int j_edge_inc = 0;
        int i_edge_inc = 0;
        if(i == spins.size()-1){i_edge_inc = -1;}else if(i == 0){i_edge_inc=1;}
        if(j == spins.size()-1){j_edge_inc = -1;}else if(j == 0){j_edge_inc=1;}
        bool j_edge = (j == spins.size()-1 || j == 0);
        bool i_edge = (i == spins.size()-1 || i == 0);

        if(!i_edge && !j_edge) {
            energy_p = spins[i][j] * spins[i][j + 1] + spins[i][j] * spins[i + 1][j]
                     + spins[i][j] * spins[i - 1][j] + spins[i][j] * spins[i][j - 1];
        }
        if(i_edge && !j_edge) {
            energy_p = spins[i][j] * spins[i][j + 1] + spins[i][j] * spins[i + i_edge_inc][j]
                       + spins[i][j] * spins[i][j - 1];
        }
        if(!i_edge && j_edge) {
            energy_p = spins[i][j] * spins[i + 1][j] + spins[i][j] * spins[i][j + j_edge_inc]
                       + spins[i][j] * spins[i - 1][j];
        }
        if(i_edge && j_edge) {
            energy_p = spins[i][j] * spins[i][j + j_edge_inc] + spins[i][j] * spins[i + i_edge_inc][j];
        }
        return energy_p;
    }

    void eval_magnetization(){
        magnetization = 0;
        for(int i = 0; i < spins.size() - 1; i++) {
            for(int j = 0; j < spins.size() - 1; j++) {
                magnetization += spins[i][j];
            }
        }
        m_sqrd = pow((magnetization / (spins.size()*spins.size())), 2.0);
        m_avg = magnetization/(spins.size()*spins.size());
    }

    int transition_prob_flip (int spin, int i, int j) {
        int E_before = energy_parallel(i, j);
        if (spin == -1) {
            spins[i][j] = 1;
            if (exp(-1 * J * beta * (energy_parallel(i, j) - E_before)) <
                static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) {
                const_count = 0;
                return -1;
            }
            {  const_count++; return 1;}
        }
        else {
            spins[i][j] = -1;
            if (exp(-1 * J * beta * (energy_parallel(i, j) - E_before)) <
                static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) {
                const_count = 0;
                return 1;
            }
            {  const_count++; return -1; }
        }
    }

    int transition_prob_flip_reverse (int spin, int i, int j) {
        int E_before = energy_parallel(i, j);
        if (spin == -1) {
            spins[i][j] = 1;
            if (energy_parallel(i, j) - E_before > 1) {
                const_count++;
                return -1;
            }
            const_count=0;
            return 1;
        }
        else {
            spins[i][j] = -1;
            if (energy_parallel(i, j) - E_before > 1) {
                const_count++;
                return 1;
            }
            const_count=0;
            return -1;
        }
    }

    vector < vector <int> > course_grain(vector < vector <int> > a){
        vector < vector <int> > cg_grid;
        int new_size = a.size()/3;
        cg_grid.resize(new_size);
        for(int i =0; i<new_size; ++i) {
            cg_grid[i].resize(new_size);
            std::fill (cg_grid[i].begin(), cg_grid[i].end(), 1);
        }
        for(int i = 0; i < a.size(); i=i+3){
            for(int j = 0; j < a.size(); j=j+3) {
                int vote = 0;
                vote += std::accumulate(spins[i].begin()+j, spins[i].begin()+j+3, 0)+
                        std::accumulate(spins[i+1].begin()+j, spins[i+1].begin()+j+3, 0)+
                        std::accumulate(spins[i+2].begin()+j, spins[i+2].begin()+j+3, 0);
                if(vote >= 1) {vote=1;} else {vote=-1;}
                cg_grid[floor(i/3)][floor(j/3)] = vote;
            }
        }
        return cg_grid;
    };
};

void print2DArray(vector < vector <int> > a) {
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < a.size(); ++j) {
            cout << a[i][j] << "\t";
        }
        cout << endl;
    }
    cout << "\n\n\n";
}

void txtout_grid(std::ofstream &txtout, vector<vector <int> > a){
    for(const auto& row : a) {
        std::copy(row.cbegin(), row.cend(), std::ostream_iterator<int>(txtout, " "));
        txtout << '\n';
    }
    txtout << ',';
}

void txtout_value(std::ofstream &txtout, double val) {
        txtout << val <<'\n';
}

vector <int> identify_config(vector<vector<int> > a) {
    vector <int> binary;
    for(int i = 0; i < a.size(); i++) {
        for(int j = 0; j < a.size(); j++) {
            if (a[i][j] == -1) { binary.push_back(0); }
            else{ binary.push_back(1); }
        }
    }
    return binary;
}

int main(int argc, char* argv[]) {
    srand(time(NULL));
    std::ifstream f;
    f.open("../input.txt");
    InputClass input;
    input.Read(f);

    double beta=input.toDouble(input.GetVariable("beta"));
    int Lx=input.toInteger(input.GetVariable("Lx"));
    int Ly=input.toInteger(input.GetVariable("Ly"));
    int size = 20;
    int num_runs = 5;
    double K = 1.381e-23;
    vector <double> betas;
    vector <double> coupling_const;
    coupling_const.push_back(1);
    betas.push_back(0.1);
    for(int b = 0; b < 400; b++){
        betas.push_back(static_cast<double>((rand() % 1000))/ static_cast<double>(1000));
    }
    int beta_steps = 100;
    lattice grid;
    lattice cg_single;
    lattice cg_double;
    std::ofstream txtoutput;
    txtoutput.open ("sim_ann_bs100_512_log.txt");

    // run = 26265 with size 4, run=432 with size 3
    grid.init(size, size);
    vector < int > a_spins;
    vector < int > b_spins;

    for(int i = 0; i<size*size; ++i) {
        if((static_cast <int> (floor(i/size)+1.0*(i % (size)))) % 2 == 0) { a_spins.push_back(i); }
        else { b_spins.push_back(i); }
    }
    std::vector<double> m_sqrd_avg_arr;
    std::vector<double> m_sqrd_avg_arr_cg_single;
    std::vector<double> m_sqrd_avg_arr_cg_double;
    grid.beta = betas[0];
    grid.J = pow(-1, rand() % 2)*coupling_const[0];
    grid.init_rand_config(size, size);
    print2DArray(grid.spins);
    double dT = 250e20;
    double T_init = grid.T;
    //for (int beta_iter = 0; beta_iter < beta_steps; beta_iter++) {
      while(exp(-1*grid.J*grid.beta)/exp(-1*coupling_const[0]*betas[0]) < 1.1 || grid.const_count < 20) {
          double diff = exp(-1*grid.J*grid.beta)/exp(-1*coupling_const[0]*betas[0]);
          //grid.beta = cur_beta;
        //auto cur_J = coupling_const[beta_iter];
        //std::cout << cur_beta << endl;
        if(grid.const_count >= 20) {
            std::cout << diff << "hi" << endl;
            grid.T = grid.T - dT;
            grid.beta = 1/(K*(T_init/log(grid.T)));
            //grid.J = grid.J*-1;
            grid.eval_energy();
            //std::cout << grid.beta << endl;
            std::cout << grid.energy << endl;
        }

        //std::cout << grid.beta << endl;
        //grid.J = cur_J;
        for (int run = 0; run < num_runs; run++) {
            grid.eval_energy();
            txtout_value(txtoutput, grid.energy);

            #pragma omp parallel for
            for (int s : a_spins) {
                grid.eval_energy();
                auto i = static_cast<int>(floor(s / size));
                int j = s % size;
                int spin = grid.spins[i][j];

                grid.spins[i][j] = grid.transition_prob_flip(spin, i, j);

                //cg_single.init_existing(grid.course_grain(grid.spins));
                //cg_double.init_existing(grid.course_grain(grid.course_grain(grid.spins)));

                //grid.eval_magnetization();
                //cg_single.eval_magnetization();
                //txtout_value(txtoutput, abs(cg_single.magnetization));
                //auto bin_repr = identify_config(cg_single.spins);
                //for (const auto &b : bin_repr) txtoutput << b;
                //txtoutput << endl;
                //cg_double.eval_magnetization();

                //m_sqrd_avg_arr.push_back(grid.m_sqrd);
                //m_sqrd_avg_arr_cg_single.push_back(cg_single.m_sqrd);
                //m_sqrd_avg_arr_cg_double.push_back(cg_double.m_sqrd);

//                if(grid.const_count > 16 && grid.beta == 0) {// && cur_beta == betas[betas.size()-1]) {
//                    std::cout << "\nLocal Minimum:\t" << endl;
//                    print2DArray(grid.spins);
//                    return 0;
//                }

            }
            #pragma omp barrier

            #pragma omp parallel for
            for (int s : b_spins) {
                grid.eval_energy();
                auto i = static_cast<int>(floor(s / size));
                int j = s % size;
                int spin = grid.spins[i][j];

                grid.spins[i][j] = grid.transition_prob_flip(spin, i, j);

                //cg_single.init_existing(grid.course_grain(grid.spins));
                //cg_double.init_existing(grid.course_grain(grid.course_grain(grid.spins)));

                //grid.eval_magnetization();
                //cg_single.eval_magnetization();
                //txtout_value(txtoutput, abs(cg_single.magnetization));
                //auto bin_repr = identify_config(cg_single.spins);
                //for (const auto &b : bin_repr) txtoutput << b;
                //txtoutput << endl;
                //cg_double.eval_magnetization();

                //m_sqrd_avg_arr.push_back(grid.m_sqrd);
                //m_sqrd_avg_arr_cg_single.push_back(cg_single.m_sqrd);
                //m_sqrd_avg_arr_cg_double.push_back(cg_double.m_sqrd);

//                if(grid.const_count > 16 && grid.beta == 0) {// && cur_beta == betas[betas.size()-1]) {
//                    std::cout << "\nLocal Minimums:\t" << endl;
//                    print2DArray(grid.spins);
//                    return 0;
//                }
            }

            ///double m_sqrd_avg =
              ///      std::accumulate(m_sqrd_avg_arr.begin(), m_sqrd_avg_arr.end(), 0.000) / m_sqrd_avg_arr.size();
            /*double m_sqrd_avg_cg_single =
                    std::accumulate(m_sqrd_avg_arr_cg_single.begin(),
                                    m_sqrd_avg_arr_cg_single.end(), 0.000) / m_sqrd_avg_arr_cg_single.size();
            double m_sqrd_avg_cg_double =
                    std::accumulate(m_sqrd_avg_arr_cg_double.begin(),
                                    m_sqrd_avg_arr_cg_double.end(), 0.000) / m_sqrd_avg_arr_cg_double.size();
            */

            //std::cout << m_sqrd_avg << std::endl;
            //std::cout << m_sqrd_avg_cg_single << std::endl;
            //std::cout << m_sqrd_avg_cg_double << std::endl;

            ///txtout_value(txtoutput, m_sqrd_avg);
            ///txtout_value(txtoutput, grid.beta);
            /*txtout_value(txtoutput, m_sqrd_avg_cg_single);
            txtout_value(txtoutput, m_sqrd_avg_cg_double);
             */
        }

    }
    txtoutput.close();
    std::cout << "\nLocal Minimums:\t" << endl;
    print2DArray(grid.spins);
    return 0;
}
