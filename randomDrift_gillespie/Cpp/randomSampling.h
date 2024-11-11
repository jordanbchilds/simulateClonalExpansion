//
//  randomSampling.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 15/05/2024.
//

#ifndef randomSampling_h
#define randomSampling_h

#include <vector>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<double> rStdNorm(0.0, 1.0);

namespace constants {
    constexpr double INV_SQRT_2PI = 0.3989422804014327;
    constexpr double LOG_INV_SQRT_2PI = -0.9189385332046727;
    constexpr double SECONDS_IN_YEAR = 3600.0*24.0*365.0;
    constexpr double SECONDS_IN_CENTURY = SECONDS_IN_YEAR * 100.0;
} // namespace constants

unsigned int choose(const unsigned int n, unsigned int k) {
    if (k>n) return 0;
    if (k*2>n) k = n-k;
    if (k==0) return 1;
    if (k==1) return n;
    int result = n;
    for(int i=2; i<=(n-k); ++i){
        result *= (n-i+1);
        result /= i;
    }
    return result;
} // function choose

inline double runif01(){
    return (double) rand() / (double) RAND_MAX ;
} // function runif01

double runif(const double& lower=0.0, const double& upper=1.0){
    return runif01()*(upper-lower) + lower;
} // funciton runif


double rexp(const double& lambda) {
    return -1.0*log(1.0-runif01())/lambda;
} // function rexp

long unsigned int rdiscrete (const std::vector<double>& weights) {
    double norm = 0.0;
    for (const double& w : weights) {
        norm += w;
    }
    
    double u = runif01() * norm;
    long unsigned int iter = 0;
    for (const double& w : weights) {
        u -= w;
        if (u < 0.0) return iter;
        iter++;
    }
    return iter;
} // function rdiscrete(weights)

long unsigned int rdiscrete (const std::vector<double>& weights, const double& norm) {
    
    double u = runif01() * norm;
    long unsigned int iter = 0;
    for (const double& w : weights) {
        u -= w;
        if (u < 0.0) return iter;
        iter++;
    }
    return iter;
} // function rdiscrete(weights, normalising constant)

double binomialPMF(const unsigned int& k, const unsigned int& n, const double& p, const bool log=false) {
    if (log) {
        return std::log(choose(n, k)) + k*std::log(p) * (n-k)*std::log(1.0 - p);
    } else {
        return (double)choose(n, k) * std::pow(p, k) * std::pow(1.0 - p, n - k);;
    }
} // binomial PMF

double normalPDF(const double x, const double mean, const double stddev, const bool logLikelihood=false) {
    if (std::isinf(stddev)) {
        if (logLikelihood) {
            return -INFINITY;
        } else {
            return 0;
        }
    }
    const double tmp = (x - mean) / stddev;
    double ll;
    if (logLikelihood) {
        ll = constants::LOG_INV_SQRT_2PI - std::log(stddev) - 0.5*tmp*tmp;
    } else {
        ll = constants::INV_SQRT_2PI / stddev * std::exp( -tmp*tmp/2.0 );
    }
    return ll;
} // function normalPDF

double dnorm(const double x, const double m, const double sd) noexcept {
    if (sd<=0) { return 0; }
    if (std::isinf(sd) | std::isinf(m))  { return 0; }
    const double std_diff = (x - m) / sd;
    return 1/(sd*std::sqrt(2*M_PI)) * std::exp( - std_diff*std_diff / 2.0 );
}

double gammaPDF(const double x, const double shape, const double rate, const bool log=false) {
    if (log) {
        const double norm = shape * std::log(rate) - std::log(std::tgamma(shape));
        return norm + (shape-1.0)*std::log(x) -1.0*rate*x;
    } else {
        const double norm = std::pow(rate, shape) / std::tgamma(shape);
        return norm * std::pow(x, shape - 1) * std::exp(-1.0*rate*x);
    }
} // function gammaPDF

std::vector<std::vector<double>> choleskyDecomposition(const std::vector<std::vector<double>>& matrix, const bool upper=false) {
    const std::size_t nrow = matrix.size();
    for (std::size_t i=0; i<nrow; ++i) {
        if (matrix[i].size()!=nrow) throw std::domain_error("poop.");
    }
    
    std::vector<std::vector<double>> result = matrix;
    
    for (size_t i = 0; i < nrow; ++i) {
        for (size_t k = 0; k < i; ++k) {
            result[i][k] = 0.0;
            result[k][i] = 0.0;
            
            double value = matrix[i][k];
            for (std::size_t j = 0; j < k; ++j) {
                value -= result[i][j] * result[k][j];
            }
            result[i][k] = value/result[k][k];
        }
        double value = matrix[i][i];
        for (std::size_t j = 0; j < i; ++j) {
            value -= result[i][j] * result[i][j];
        }
        result[i][i] = std::sqrt(value);
    }
    return result;
} // function cholesky decomposition

std::vector<double> rMultivariateNormal_cholesky(const std::vector<double>& mu, const std::vector<std::vector<double>>& lowerTriangular) {
    
    const std::size_t dim = mu.size();
    std::vector<double> rand(dim, 0.0);
    for (std::size_t i=0; i<dim; ++i) {
        rand[i] = rStdNorm(gen);
    }
    std::vector<double> randomDraw(mu.begin(), mu.end());
    for (std::size_t j=0; j<dim; ++j) {
        for (std::size_t k=0; k<dim; ++k) {
            randomDraw[j] += rand[k]*lowerTriangular[j][k];
        }
    }
    return randomDraw;
} // function rMultivariateNormal_cholesky

std::vector<double> rMultivariateNormal(const std::vector<double>& mu, const std::vector<std::vector<double>>& Sigma) {
    const std::vector<std::vector<double>> lowerTriangular = choleskyDecomposition(Sigma);
    return rMultivariateNormal_cholesky(mu, lowerTriangular);
} // function rMultivariateNormal

#endif /* randomSampling_h */
