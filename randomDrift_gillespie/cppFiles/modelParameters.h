//
//  ModelParameters.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 30/05/2024.
//

#ifndef ModelParameters_h
#define ModelParameters_h

#include <iostream>
#include <fstream>
#include "randomSampling.h"

namespace trueParameters {
    const double genRate = 6.00E-8;
    const double mutationProb = 0.02;
    const double defThreshold = 0.85;
    const double modelError = 5.0;
    const double modelPrecision = 1 / (modelError * modelError);
    const unsigned int initialPopulation = 1000;
    const unsigned int nTrials = 100;
}

namespace priorParameters {
    const double logGenRate_lowerBound = std::log(1.0E-12);
    const double logGenRate_upperBound = std::log(1.0E-3);
    const double logGenRate_priorMean = std::log(6.0E-6);
    const double logGenRate_priorSD = 10.0;
    const double logGenRate_proposalSD = 1.0;
    
    const double logMutationProb_lowerBound = -INFINITY;
    const double logMutationProb_upperBound = std::log(0.5);
    const double logMutationProb_priorMean = std::log(trueParameters::mutationProb);
    const double logMutationProb_priorSD = 10.0;
    const double logMutationProb_proposalSD = 0.5;
    
    const double logDefThreshold_lowerBound = std::log(0.4);
    const double logDefThreshold_upperBound = std::log(1.0);
    const double logDefThreshold_priorMean = std::log(trueParameters::defThreshold);
    const double logDefThreshold_priorSD = 1.0;
    const double logDefThreshold_proposalSD = 0.05;
    
    const double modelPrecision_mean = trueParameters::modelPrecision;
    const double modelPrecision_var = 5.0;
    const double modelPrecision_shape = modelPrecision_mean*modelPrecision_mean / modelPrecision_var;
    const double modelPrecision_rate = modelPrecision_mean / modelPrecision_var;
    const double modelPrecision_scale = 1.0 / modelPrecision_rate;
    std::gamma_distribution<double> rPriorGamma(modelPrecision_shape, modelPrecision_scale);
}

class ModelParameters {
public:
    ModelParameters () noexcept {
        _genRate = trueParameters::genRate;
        _mutationProb = trueParameters::mutationProb;
        _defThreshold = trueParameters::defThreshold;
        _modelPrecision = trueParameters::modelPrecision;
        updateLogs();
    }
    
    ModelParameters (const double genRate, const double mutationProb, const double defThreshold, const double modelPrecision) noexcept {
        _genRate = genRate;
        _mutationProb = mutationProb;
        _defThreshold = defThreshold;
        _modelPrecision = modelPrecision;
        updateLogs();
    }

    ModelParameters (const double* ptr) noexcept {
        _genRate = *ptr;
        _mutationProb = *(ptr+1);
        _defThreshold = *(ptr+2);
        _modelPrecision = *(ptr+3);
        updateLogs();
    }

    ModelParameters (const std::vector<double>::const_iterator itr) noexcept {
        _genRate = *itr;
        _mutationProb = *(itr+1);
        _defThreshold = *(itr+2);
        _modelPrecision = *(itr+3);
        updateLogs();
    }
    
    void priorDraw () noexcept ;
    double priorSKMdensity () const noexcept ;
    void priorDraw (std::fstream& file, unsigned int nOut=1000, const std::string& sep="\t") noexcept ;
    bool checkLogSupport () const noexcept ;
    void updateLogs(const bool& inverse=false) noexcept ;
    void shiftLogSKMparameters (const std::vector<double>& delta) ;
    void shiftLogSKMparameters (const double* const delta) noexcept;
    void shiftLogSKMparameters (const double deltaRepRate, const double deltaMutProb, const double deltaDefThresh) noexcept;

    // PRINT
    void printToFile (std::ofstream& file, const std::string& sep="\t") const noexcept ;
    void print () const noexcept {
        std::cout << _logGenRate << ", " << _logMutationProb << ", " << _logDefThreshold << ", " << _modelPrecision << "\n";
    }
    
    // SETTERs
    void setGenRate (const double r) {_genRate = r; _logGenRate = std::log(r); }
    void setMutationProb (const double p) { _mutationProb = p; _logMutationProb = std::log(p); }
    void setDefThreshold (const double p) { _defThreshold = p; _logDefThreshold = std::log(p); }
    void setModelPrecision (const double sd) { _modelPrecision = sd; } 
    void setInitialPopulation (const unsigned int x) { _initialPopulation = x; }
    void setSKMparameters (const ModelParameters& other) noexcept ;
    
    void setLogGenRate( const double logr) { _logGenRate = logr; _genRate = std::exp(logr);}
    void setLogMutationProb (const double logp) { _logMutationProb = logp; _mutationProb = std::exp(logp); }
    void setLogDefThreshold (const double logp) { _logDefThreshold = logp; _defThreshold = std::exp(logp); }
    
    // GETTERs
    double getGenRate () const noexcept { return _genRate; }
    double getMutationProb () const noexcept { return _mutationProb; }
    double getDefThreshold () const noexcept { return _defThreshold; }
    double getModelPrecision () const noexcept { return _modelPrecision; }
    double getModelError () const noexcept { return 1 / std::sqrt(_modelPrecision); }
    
    double getLogGenRate () const noexcept { return _logGenRate; }
    double getLogMutationProb () const noexcept { return _logMutationProb; }
    double getLogDefThreshold () const noexcept{ return _logDefThreshold; }
    
    unsigned int getInitialPopulation () const noexcept { return _initialPopulation; }
private:
    double _genRate;
    double _mutationProb;
    double _defThreshold;
    double _modelPrecision;
    unsigned int _initialPopulation = trueParameters::initialPopulation;
    
    double _logGenRate;
    double _logMutationProb;
    double _logDefThreshold;
};

double ModelParameters::priorSKMdensity () const noexcept {
    // don't forget the Jacobian terms!!!!
    // the priors are truncated. However do not need to calculate the normalising constant as they are cancelled out in the acceptance probability ratio
    const double rateDensity = normalPDF(_logGenRate, priorParameters::logGenRate_priorMean, priorParameters::logGenRate_priorSD, true) - _logGenRate ;
    const double mutDensity = normalPDF(_logMutationProb, priorParameters::logMutationProb_priorMean, priorParameters::logMutationProb_priorSD, true) - _logMutationProb;
    const double threshDensity = normalPDF(_logDefThreshold, priorParameters::logDefThreshold_priorMean, priorParameters::logDefThreshold_priorSD, true) - _logDefThreshold;
    return rateDensity + mutDensity + threshDensity;
}

void ModelParameters::shiftLogSKMparameters (const std::vector<double>& delta) {
    if (delta.size() != 3) {
        throw std::runtime_error("Dimension mismatch. delta is not of length three.");
    }
    _logGenRate += delta[0];
    _logMutationProb += delta[1];
    _logDefThreshold += delta[2];
};
void ModelParameters::shiftLogSKMparameters (const double* const delta) noexcept {
    _logGenRate += *delta;
    _logMutationProb += *(delta + 1);
    _logDefThreshold += *(delta + 2);
    updateLogs(true);
}
void ModelParameters::shiftLogSKMparameters (const double deltaRepRate, const double deltaMutProb, const double deltaDefThresh) noexcept {
    _logGenRate += deltaRepRate;
    _logMutationProb += deltaMutProb;
    _logDefThreshold += deltaDefThresh;
    updateLogs(true);
}

void ModelParameters::priorDraw () noexcept {
    bool outsideSupport = true;
    _modelPrecision = priorParameters::rPriorGamma(gen);
    while (outsideSupport) {
        _logGenRate = rStdNorm(gen) * priorParameters::logGenRate_priorSD + priorParameters::logGenRate_priorMean;
        _logMutationProb = rStdNorm(gen) * priorParameters::logMutationProb_priorSD + priorParameters::logMutationProb_priorMean;
        _logDefThreshold = rStdNorm(gen) * priorParameters::logDefThreshold_priorSD + priorParameters::logDefThreshold_priorMean;
        outsideSupport = !checkLogSupport();
    }
    updateLogs(true);
}

void ModelParameters::priorDraw (std::fstream& file, unsigned int nOut, const std::string& sep) noexcept {
    for (unsigned int i=0; i<nOut; ++i) {
        priorDraw();
        file << _logGenRate << sep << _logMutationProb << sep << _logDefThreshold << sep << _modelPrecision << "\n";
    }
}

bool ModelParameters::checkLogSupport () const noexcept {
    return _logGenRate > priorParameters::logGenRate_lowerBound && _logGenRate < priorParameters::logGenRate_upperBound && 
        _logMutationProb > priorParameters::logMutationProb_lowerBound && _logMutationProb < priorParameters::logMutationProb_upperBound && 
        _logDefThreshold > priorParameters::logDefThreshold_lowerBound && _logDefThreshold < priorParameters::logDefThreshold_upperBound;
}

void ModelParameters::updateLogs (const bool& inverse) noexcept {
    if (inverse) {
        _genRate = std::exp(_logGenRate);
        _mutationProb = std::exp(_logMutationProb);
        _defThreshold = std::exp(_logDefThreshold);
    } else {
        _logGenRate = std::log(_genRate);
        _logMutationProb = std::log(_mutationProb);
        _logDefThreshold = std::log(_defThreshold);
    }
}

void ModelParameters::setSKMparameters (const ModelParameters& other) noexcept  {
    // UPDATE SKM PARAMS
    _genRate = other.getGenRate();
    _mutationProb = other.getMutationProb();
    _defThreshold = other.getDefThreshold();
    updateLogs();
};

// PRINT
void ModelParameters::printToFile (std::ofstream& file, const std::string& sep) const noexcept {
    file << _logGenRate << sep << _logMutationProb << sep << _logDefThreshold << sep << _modelPrecision << "\n";
}

#endif /* ModelParameters_h */
