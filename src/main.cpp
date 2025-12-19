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
    double x_n = 0;

    constructed_function.x.push_back(x_n);
    constructed_function.y.push_back(y_n);

    std::cout << x_n << "\t" << y_n << std::endl;

    for (size_t n = 1; n < max_n; n++)
    {
        double k1 = func(x_n, y_n);
        double k2 = func(x_n + h/2.0, y_n + h*k1/2.0);
        double k3 = func(x_n + h/2.0, y_n + h*k2/2.0);
        double k4 = func(x_n + h, y_n + h*k3);

        y_n = y_n + h/6.0*(k1 + 2*k2 + 2*k3 + k4);
        x_n = x_n + h;

        constructed_function.x.push_back(x_n);
        constructed_function.y.push_back(y_n);

        std::cout << x_n << "\t" << y_n << std::endl;
    }
    
    return constructed_function;
}
double param = 1;

double funcSecond(double x, double y0, double y1) {
    double hbar = 1.054571817e-34; // Js
    double m_e = 9.1093837139e-31; // kg

    return -param*y0; // kg*J/(Js)^2
}

NumericalFunction rungekuttanystrom() {
    NumericalFunction constructed_function;

    double hbar = 1.054571817e-34; // Js
    double m_e = 9.1093837139e-31; // kg

    double h = 0.00001e-3; // m

    double y0_n = 0; // unitless
    double y1_n = 1; // 1/m

    double x_n = 0; // m

    constructed_function.x.push_back(x_n);
    constructed_function.y.push_back(y0_n);

    //std::cout << x_n << "\t" << y0_n << std::endl;

    for (size_t n = 1; n < 100001; n++)
    {
        double k1 = funcSecond(x_n, y0_n, y1_n);

        double ydk1 = y1_n + k1*h/2;
        double yk1 = y0_n + h/2*((y1_n+ydk1)/2);
        double k2 = funcSecond(x_n + h/2, yk1, ydk1);

        double ydk2 = y1_n + k2*h/2;
        double yk2 = y0_n + h/2*((y1_n + ydk2)/2);
        double k3 = funcSecond(x_n + h/2, yk2, ydk2);

        double ydk3 = y1_n + k3*h;
        double yk3 = y0_n + h*((y1_n + ydk3)/2);
        double k4 = funcSecond(x_n + h, yk3, ydk3);

        double ydk4 = y1_n + k4*h;
        
        y1_n = y1_n + h/6*(k1 + 2*k2 + 2*k3 + k4);
        y0_n = y0_n + h/6*(ydk1 + 2*ydk2 + 2*ydk3 + ydk4);
        x_n = x_n + h;

        constructed_function.x.push_back(x_n);
        constructed_function.y.push_back(y0_n);

        //std::cout << x_n << "\t" << y0_n << std::endl;
    }
    
    return constructed_function;
}

NumericalFunction normalizeWavefunction(NumericalFunction input) {
    NumericalFunction output;
    double L = 0.0001e-3*1001;
    double total_integral = 0;
    for (size_t i = 1; i < input.x.size(); i++)
    {
        double previous = std::pow(input.y[i-1], 2);
        double current = std::pow(input.y[i], 2);
        total_integral += (input.x[i]-input.x[i-1])/L*(previous+current)/2;
    }

    double correction_factor = std::sqrt(total_integral);
    for (size_t i = 0; i < input.x.size(); i++)
    {
        output.x.push_back(input.x[i]);
        output.y.push_back(input.y[i]/correction_factor);
    }
    //std::cout << output.y[input.x.size()-1] << std::endl;
    return output;
}

int numRoots(NumericalFunction input) {
    int zeros = 0;
    for (size_t i = 1; i < input.x.size(); i++)
    {
        if (input.y[i]*input.y[i-1] < 0)
        {
            zeros += 1;
        }
    }
    return zeros;    
}

int main() {
    NumericalFunction params;

    for (size_t i = 1; i < 10; i++)
    {
        NumericalFunction result;
        result = rungekuttanystrom();
        result = normalizeWavefunction(result);

        params.x.push_back(param);
        params.y.push_back(result.y[result.y.size() - 1]);

        NumericalFunction as;
        as.x.push_back(param);
        as.y.push_back(result.y[result.y.size() - 1]);

        NumericalFunction bs;

        while (bs.x.size() < 1)
        {
            int roots = numRoots(result);
            if (roots >= i)
            {
                bs.x.push_back(param);
                bs.y.push_back(result.y[result.y.size() - 1]);
            } else
            {
                as.x.push_back(param);
                as.y.push_back(result.y[result.y.size() - 1]);

                param = param * 2;
                result = rungekuttanystrom();
                result = normalizeWavefunction(result);
            }
        }
        //std::cout << "Found new bounds: " << as.x[as.x.size() - 1] << "," << bs.x[bs.x.size() - 1] << std::endl;
        // Ensure only one root
        int roots = -1;
        while (roots != i)
        {
            double middle = (bs.x[bs.x.size() - 1] + as.x[as.x.size() - 1]) / 2;
            param = middle;

            result = rungekuttanystrom();
            result = normalizeWavefunction(result); 

            params.x.push_back(param);
            params.y.push_back(result.y[result.y.size() - 1]);

            roots = numRoots(result);
            if (roots >= i)
            {
                bs.x.push_back(param);
                bs.y.push_back(result.y[result.y.size() - 1]);
            } else if (roots < i)
            {
                as.x.push_back(param);
                as.y.push_back(result.y[result.y.size() - 1]);
            }
        }
        
        // Binary search
        while (std::abs(result.y[result.y.size() - 1]) > 0.00001)
        {
            double middle = (bs.x[bs.x.size() - 1] + as.x[as.x.size() - 1]) / 2;
            param = middle;

            result = rungekuttanystrom();
            result = normalizeWavefunction(result); 

            params.x.push_back(param);
            params.y.push_back(result.y[result.y.size() - 1]);
            
            if (result.y[result.y.size() - 1]*as.y[0] < 0)
            {
                bs.x.push_back(param);
                bs.y.push_back(result.y[result.y.size() - 1]);
            } else
            {
                as.x.push_back(param);
                as.y.push_back(result.y[result.y.size() - 1]);
            }
        }
        
        // result = rungekuttanystrom();
        // result = normalizeWavefunction(result);
        writeCSV(result, "numerical.csv");
        std::cout << "Eigenstate: " << i << ", Eigenvalue: " << param << std::endl;
        param += 10000;
    }

    writeCSV(params, "params.csv");

    return 0;
}