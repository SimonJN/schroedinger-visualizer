#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

struct NumericalFunction {
    std::vector<double> x;
    std::vector<double> y;
};

void writeCSV(NumericalFunction nf, std::string path) {
    std::ofstream output_file(path);

    for (size_t i = 0; i < nf.x.size(); i++)
    {
        output_file << nf.x[i] << "," << nf.y[i] << std::endl;
    }

    output_file.close();
}

double func(double t, double y) {
    return std::pow(std::sin(t), 2)*y;
}

NumericalFunction rungekutta(double y0, double h, uint max_n) {
    NumericalFunction constructed_function;

    double y_n = y0;
    double t_n = 0;

    constructed_function.x.push_back(t_n);
    constructed_function.y.push_back(y_n);

    std::cout << t_n << "\t" << y_n << std::endl;

    for (size_t n = 1; n < max_n; n++)
    {
        double k1 = func(t_n, y_n);
        double k2 = func(t_n + h/2.0, y_n + h*k1/2.0);
        double k3 = func(t_n + h/2.0, y_n + h*k2/2.0);
        double k4 = func(t_n + h, y_n + h*k3);

        y_n = y_n + h/6.0*(k1 + 2*k2 + 2*k3 + k4);
        t_n = t_n + h;

        constructed_function.x.push_back(t_n);
        constructed_function.y.push_back(y_n);

        std::cout << t_n << "\t" << y_n << std::endl;
    }
    
    return constructed_function;
}

int main() {
    NumericalFunction result = rungekutta(2, 0.5, 11);
    writeCSV(result, "numerical.csv");
    return 0;
}