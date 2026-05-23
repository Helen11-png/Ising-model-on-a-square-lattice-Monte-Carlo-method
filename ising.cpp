#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
using namespace std;

const int L = 128;
const int N = L * L;
const double J = 1.0;
const double kB = 1.0;

const int thermal_steps = 500;
const int measure_steps = 500;

const double T_min = 1.0;
const double T_max = 4.0;
const int T_steps = 31;

vector<vector<int>> spin(L, vector<int>(L));
double M = 0;
double E = 0;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0.0, 1.0);
uniform_int_distribution<> dis_site(0, L - 1);

void save_lattice(const string& filename = "lattice.txt") {
    ofstream out(filename);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            out << spin[i][j] << (j == L - 1 ? '\n' : ' ');
        }
    }
}

double deltaE(int x, int y) {
    int s = spin[(x+1)%L][y] + spin[(x-1+L)%L][y]
          + spin[x][(y+1)%L] + spin[x][(y-1+L)%L];
    return 2.0 * J * spin[x][y] * s;
}

void init_ferro() {
    E = 0;
    M = N;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            spin[i][j] = 1;
            if (i + 1 < L) E += -J;
            if (j + 1 < L) E += -J;
        }
    for (int i = 0; i < L; i++) {
        E += -J; // граница i, L-1 -> 0
        E += -J; // граница L-1, i -> 0
    }
}

void metropolis_step(double T) {
    int x = dis_site(gen);
    int y = dis_site(gen);
    double dE = deltaE(x, y);
    if (dE <= 0 || dis(gen) < exp(-dE / (kB * T))) {
        M -= 2 * spin[x][y];
        E += dE;
        spin[x][y] = -spin[x][y];
    }
}

void run_monte_carlo(double T, double& avg_M, double& avg_M2,
                     double& avg_E, double& avg_E2) {
    for (int step = 0; step < thermal_steps * N; step++)
        metropolis_step(T);

    double sum_M = 0, sum_M2 = 0, sum_E = 0, sum_E2 = 0;
    int cnt = measure_steps * N;

    for (int step = 0; step < cnt; step++) {
        metropolis_step(T);
        double m = fabs(M);
        sum_M  += m;
        sum_M2 += M * M;
        sum_E  += E;
        sum_E2 += E * E;
    }

    avg_M  = sum_M  / cnt / N;
    avg_M2 = sum_M2 / cnt / (double)(N * N);
    avg_E  = sum_E  / cnt / N;
    avg_E2 = sum_E2 / cnt / (double)(N * N);
}

int main() {
    // Результаты пишем в results.txt (читается HTML)
    ofstream out("results.txt");
    out << "T M C Chi\n";

    for (int t_idx = 0; t_idx < T_steps; t_idx++) {
        double T = T_min + (T_max - T_min) * t_idx / (T_steps - 1);

        init_ferro();
        double M_avg, M2_avg, E_avg, E2_avg;
        run_monte_carlo(T, M_avg, M2_avg, E_avg, E2_avg);

        double C   = N * (E2_avg - E_avg * E_avg) / (kB * T * T);
        double Chi = (M2_avg - M_avg * M_avg) / (kB * T);

        out << T << " " << M_avg << " " << C << " " << Chi << "\n";
        cout << "T = " << T
             << "  M = " << M_avg
             << "  C = " << C
             << "  Chi = " << Chi << "\n";
    }
    out.close();

    // Сохраняем решётку при последней температуре
    save_lattice("lattice.txt");
    cout << "Готово. Файлы: results.txt, lattice.txt\n";
    return 0;
}