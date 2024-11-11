//
//  MetropolisHastings.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 24/05/2024.
//

#ifndef MetropolisHastings_h
#define MetropolisHastings_h

#include <cmath>
#include <array>

#include <iostream>
#include <fstream>

#include "randomSampling.h"
#include "dataClass.h"
#include "ModelParameters_MH.h"
#include "simulate.h"

// --- --- --- 
// --- PRIOR PARAMETERS
// --- --- ---
namespace priorParameters {
    const double logGenRate_lowerBound = std::log(1.0E-10);
    const double logGenRate_upperBound = std::log(5.0E-7);
    const double logGenRate_priorMean = std::log(trueParameters::genRate);
    const double logGenRate_priorSD = 10.0;
    const double logGenRate_proposalSD = 1.0;
    
    const double logMutationProb_lowerBound = -INFINITY;
    const double logMutationProb_upperBound = std::log(0.5);
    const double logMutationProb_priorMean = std::log(trueParameters::mutationProb);
    const double logMutationProb_priorSD = 10.0;
    const double logMutationProb_proposalSD = 1.0;
    
    const double logDefThreshold_lowerBound = std::log(0.5);
    const double logDefThreshold_upperBound = std::log(1.0);
    const double logDefThreshold_priorMean = std::log(trueParameters::defThreshold);
    const double logDefThreshold_priorSD = 1.0;
    const double logDefThreshold_proposalSD = 0.1;

    const double logModelError_lowerBound = -INFINITY;
    const double logModelError_upperBound = INFINITY;
    const double logModelError_priorMean = std::log(trueParameters::modelError);
    const double logModelError_priorSD = 1.0;
    const double logModelError_proposalSD = 0.1;
}

// --- --- ---
// --- MODEL PARAMETERS CLASS
// --- --- ---

class ModelParameters {
public:
    ModelParameters () noexcept {
        _genRate = trueParameters::genRate;
        _mutationProb = trueParameters::mutationProb;
        _defThreshold = trueParameters::defThreshold;
        _modelError = trueParameters::modelError;
        updateLogs();
    }
    
    ModelParameters (const double genRate, const double mutationProb, const double defThreshold, const double modelError, const bool logScale = false) noexcept {
        if (logScale) {
            _logGenRate = genRate;
            _logMutationProb = mutationProb;
            _logDefThreshold = defThreshold;
            _logModelError = modelError;
            updateLogs(true);
        } else  {
            _genRate = genRate;
            _mutationProb = mutationProb;
            _defThreshold = defThreshold;
            _modelError = modelError;
            updateLogs();
        }
    }

    ModelParameters (const double* ptr, const bool logScale = false) noexcept {
        if (logScale) {
            _logGenRate = *ptr;
            _logMutationProb = *(ptr+1);
            _logDefThreshold = *(ptr+2);
            _logModelError = *(ptr+3);
            updateLogs(true);
        } else  {
            _genRate = *ptr;
            _mutationProb = *(ptr+1);
            _defThreshold = *(ptr+2);
            _modelError = *(ptr+3);
            updateLogs();
        }
    }

    ModelParameters (std::vector<double>::iterator itr, const bool logScale = false) noexcept {
        if (logScale) {
            _logGenRate = *itr;
            _logMutationProb = *(itr+1);
            _logDefThreshold = *(itr+2);
            _logModelError = *(itr+3);
            updateLogs(true);
        } else  {
            _genRate = *itr;
            _mutationProb = *(itr+1);
            _defThreshold = *(itr+2);
            _modelError = *(itr+3);
            updateLogs();
        }
    }

    ModelParameters (std::vector<double>& theta, const bool logScale = false) noexcept {
        if (logScale) {
            _logGenRate = theta[0];
            _logMutationProb = theta[1];
            _logDefThreshold = theta[2];
            _logModelError = theta[3];
            updateLogs(true);
        } else  {
            _genRate = theta[0];
            _mutationProb = theta[1];
            _defThreshold = theta[2];
            _modelError = theta[3];
            updateLogs();
        }
    }
    
    double priorLogDensity () const noexcept ;
    void priorDraw () noexcept ;
    void priorDraw (std::fstream& file, unsigned int nOut=1000, const std::string& sep="\t") noexcept ;
    bool checkLogSupport () const noexcept ;
    void updateLogs(const bool inverse=false) noexcept ;

    void shiftLogParameters (const std::vector<double>& delta) noexcept;

    // PRINT
    void printToFile (std::ofstream& file, const std::string& sep="\t") const noexcept ;
    void print () const noexcept {
        std::cout << _logGenRate << ", " << _logMutationProb << ", " << _logDefThreshold << ", " << _logModelError << "\n";
    }
    
    // SETTERs
    void setGenRate (const double r) {_genRate = r; _logGenRate = std::log(r); }
    void setMutationProb (const double p) { _mutationProb = p; _logMutationProb = std::log(p); }
    void setDefThreshold (const double p) { _defThreshold = p; _logDefThreshold = std::log(p); }
    void setModelError (const double sd) { _modelError = sd; _logModelError = std::log(sd); }
    void setInitialPopulation (const unsigned int x) { _initialPopulation = x; }
    
    void setLogGenRate( const double logr) { _logGenRate = logr; _genRate = std::exp(logr);}
    void setLogMutationProb (const double logp) { _logMutationProb = logp; _mutationProb = std::exp(logp); }
    void setLogDefThreshold (const double logp) { _logDefThreshold = logp; _defThreshold = std::exp(logp); }
    void setLogModelError (const double logsd) { _logModelError = logsd; _modelError = std::exp(logsd);}
    
    // GETTERs
    double getGenRate () const noexcept { return _genRate; }
    double getMutationProb () const noexcept { return _mutationProb; }
    double getDefThreshold () const noexcept { return _defThreshold; }
    double getModelError () const noexcept { return _modelError; }
    
    double getLogGenRate () const noexcept { return _logGenRate; }
    double getLogMutationProb () const noexcept { return _logMutationProb; }
    double getLogDefThreshold () const noexcept{ return _logDefThreshold; }
    double getLogModelError () const noexcept { return _logModelError; }
    
    unsigned int getInitialPopulation () const noexcept { return _initialPopulation; }
private:
    double _genRate;
    double _mutationProb;
    double _defThreshold;
    double _modelError;
    unsigned int _initialPopulation = trueParameters::initialPopulation;
    
    double _logGenRate;
    double _logMutationProb;
    double _logDefThreshold;
    double _logModelError;
};

// --- ---
// --- MEMBER FUNCTIONS
// --- ---

// --- 
// PUBLIC MEMBER FUNCTIONS
// ---

double ModelParameters::priorLogDensity () const noexcept {
    /*
        Returns the log density up to proportionality (the truncation constants are ignored). The priors used are normal for all parameters, on the log-scale. 
    */
    double rateDensity = normalPDF(_logGenRate, priorParameters::logGenRate_priorMean, priorParameters::logGenRate_priorSD, true) - _logGenRate ;
    double mutDensity = normalPDF(_logMutationProb, priorParameters::logMutationProb_priorMean, priorParameters::logMutationProb_priorSD, true) - _logMutationProb;
    double threshDensity = normalPDF(_logDefThreshold, priorParameters::logDefThreshold_priorMean, priorParameters::logDefThreshold_priorSD, true) - _logDefThreshold;
    double errorDensity = normalPDF(_logModelError, priorParameters::logModelError_priorMean, priorParameters::logModelError_priorSD, true) - _logModelError;
    return rateDensity + mutDensity + threshDensity + errorDensity;
} // ModelParameters::priorLogDensity


void ModelParameters::priorDraw () noexcept {
    /*
        Alters the object's parameter values to be random draws from their prior distributions.
    */
    bool outsideSupport = true;
    while (outsideSupport) {
        _logGenRate = rStdNorm(gen) * priorParameters::logGenRate_priorSD + priorParameters::logGenRate_priorMean;
        _logMutationProb = rStdNorm(gen) * priorParameters::logMutationProb_priorSD + priorParameters::logMutationProb_priorMean;
        _logDefThreshold = rStdNorm(gen) * priorParameters::logDefThreshold_priorSD + priorParameters::logDefThreshold_priorMean;
        _logModelError = rStdNorm(gen) * priorParameters::logModelError_priorSD + priorParameters::logModelError_priorMean;
        outsideSupport = !checkLogSupport();
    }
    updateLogs(true);
}

void ModelParameters::priorDraw (std::fstream& file, unsigned int nOut, const std::string& sep) noexcept {
    /*
        Prints `nOut` prior draws to a file.
    */
    for (unsigned int i=0; i<nOut; ++i) {
        priorDraw();
        file << _logGenRate << sep << _logMutationProb << sep << _logDefThreshold << sep << _logModelError << "\n";
    }
} // ModelParameters::priorDraw 


bool ModelParameters::checkLogSupport () const noexcept {
    /*
        Returns a boolean indicating if the current values are within the bounds defined in the priorParameter namespace. 
    */
    return _logGenRate > priorParameters::logGenRate_lowerBound && _logGenRate < priorParameters::logGenRate_upperBound && 
        _logMutationProb > priorParameters::logMutationProb_lowerBound && _logMutationProb < priorParameters::logMutationProb_upperBound && 
        _logDefThreshold > priorParameters::logDefThreshold_lowerBound && _logDefThreshold < priorParameters::logDefThreshold_upperBound;
} // ModelParameters::checkLogSupport

void ModelParameters::shiftLogParameters (const std::vector<double>& delta) {
    /*
        Jitters the current parameter values on the log-scale. The vector is checked to be the correct size and alters in the order; replication rate, mutation probability, deficient threshold, and model error.
    */
    if (delta.size() != 4) {
        throw std::runtime_error("Dimension mismatch. `delta` is not of length four.");
    }
    _logGenRate += delta[0];
    _logMutationProb += delta[1];
    _logDefThreshold += delta[2];
    _logModelError += delta[3];
} // ModelParameters::shiftLogParameters

void ModelParameters::updateLogs (const bool inverse) noexcept {
    /*
        Updates the parameter values to be in agreement with the log or raw values of the current object. 
    */
    if (inverse) {
        _genRate = std::exp(_logGenRate);
        _mutationProb = std::exp(_logMutationProb);
        _defThreshold = std::exp(_logDefThreshold);
        _modelError = std::exp(_logModelError);
    } else {
        _logGenRate = std::log(_genRate);
        _logMutationProb = std::log(_mutationProb);
        _logDefThreshold = std::log(_defThreshold);
        _logModelError = std::log(_modelError);
    }
} // ModelParameters::updateLogs 

void ModelParameters::printToFile (std::ofstream& file, const std::string& sep) const noexcept {
    /*
        Print model parameters on the log-scale to a `file`.
    */
    file << _logGenRate << sep << _logMutationProb << sep << _logDefThreshold << sep << _logModelError << "\n";
} // ModelParameters::printToFile


class MetropolisHastings {
private:
    std::vector<ProportionDeficientData> _proportionData;
    std::vector<ProportionDeficientData> _prevSimulatedData;
    std::vector<ProportionDeficientData> _simulatedDataStar;

    unsigned int _nSims;
    unsigned int _nPat;

    ModelParameters _thetaStar;
    double _prevDataLogLikelihood;
    double _prevThetaLogLikelihood;

    std::vector<ModelParameters> _posterior;
    unsigned int _totalProposed;
    unsigned int _totalAccepted;
    std::vector<std::vector<double>> _proposalCovariance;
    std::vector<std::vector<double>> _proposalCovariance_cholesky;
    
    void _genericSetUp ();
    void _simulateData (const ModelParameters& theta) noexcept;
    double _logLikelihood (const std::vector<ProportionDeficientData>& simData, const double modelError) const;
    void _proposalDistribution (ModelParameters& thetaStar, const std::vector<ModelParameters>::const_iterator prevTheta) const noexcept ;

public:
    ~ MetropolisHastings () {}
    MetropolisHastings (const unsigned int nSims, const std::vector<ProportionDeficientData>& proportionData) : _nSims(nSims), _proportionData(proportionData) {
        _genericSetUp();
    }
    
    void inferParameters (const unsigned int nPost) noexcept;
    void modelPrediction (std::ofstream& file, const unsigned int nRep=1, const std::string& sep="\t") const noexcept;
    void proportionPrediction (std::ofstream& file, const unsigned int nRep=1, const std::string& sep="\t") const noexcept ;
    
    // SETTERs
    void setProposalCovarianceMatrix (const std::vector<std::vector<double>>& matrix) {
        _proposalCovariance = matrix;
        _proposalCovariance_cholesky = choleskyDecomposition(matrix);
    }
    
    // GETTERs
    double getAcceptanceRatio () const { return (double)_totalAccepted / (double)_totalProposed; }

    // PRINTERS
    void printPosteriorToFile (std::ofstream& file, const std::string& sep="\t") const ;
};

// --- ---
// MEMBER FUNCTIONS
// --- ---


// --- PRIVATE MEMBER FUNCTIONS

void MetropolisHastings::_genericSetUp () {
    _nPat = static_cast<unsigned int>(_proportionData.size());

    for (auto iter = _proportionData.begin(); iter!=_proportionData.end(); ++iter) {
        _simulatedDataStar.emplace_back(1, iter->getTimeVector(), _nSims);
        _prevSimulatedData.emplace_back(1, iter->getTimeVector(), _nSims);
        
        iter->setNsim(_nSims);
        iter->calculateDataTransform();
        iter->checkDataIsFinite();
    }
    _proposalCovariance = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0} };
    _proposalCovariance_cholesky = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0} };
    _totalAccepted = 0;
    _totalProposed = 0;
} // MetropolisHastings::_genericSetUp

void MetropolisHastings::_proposalDistribution (ModelParameters& thetaStar, const std::vector<ModelParameters>::const_iterator prevTheta) const noexcept {
    /*
        Propose a new theta value (for all parameters) based upon the previously accepted theta value
    */
    bool withinSupport = false;
    while (!withinSupport) {
        const std::vector<double> deltaParam = rMultivariateNormal_cholesky({0,0,0,0}, _proposalCovariance_cholesky);
        thetaStar.shiftLogParameters(deltaParam);
        withinSupport = thetaStar.checkLogSupport();
        thetaStar.updateLogs(true);
    }
} // MetropolisHastings::_proposalDistribution


void MetropolisHastings::_simulateData (const ModelParameters& theta) noexcept {
    for (auto iter=_simulatedDataStar.begin(); iter!=_simulatedDataStar.end(); ++iter) {
        iter->zeroCountData();
        // for (auto iter=_proportionData.begin(); iter!=_proportionData.end(); ++iter)
        for (std::size_t j=0; j<iter->getNsim(); ++j) {
            stochasticFibre mySystem = initialiseFibre(theta);
            std::size_t tt = 0;
            // while (tt<_timeDataSteps[i].size()) {
            while (tt<iter->getNobs()) {
                // simulate over next time step
                gillespieSimulation(mySystem, iter->getTimeJump(tt));
                // Check mutation passed the pathogenic threshold
                if (mySystem.getMutationLoad() > theta.getDefThreshold()) {
                    iter->increaseCount(tt, 0);
                }
                tt++;
            }
        }
        iter->updateProportionData();
        iter->checkDataIsFinite();
    }
} // MetropolisHastings::_simulateData

double MetropolisHastings::_logLikelihood (const std::vector<ProportionDeficientData>& simData, const double modelError) const {
    double logLikelihood = 0.0;
    for (std::size_t i=0; i<_proportionData.size(); ++i) {
        for (unsigned int k=0; k<_proportionData[i].getNobs(); ++k) {
            const double simulatedTransformedProportion = simData[i].getTransformedValue(k, 0);
            for (unsigned int j=0; j<_proportionData[i].getNrep(k); ++j) {
                logLikelihood += normalPDF(_proportionData[i].getTransformedValue(k,j), simulatedTransformedProportion, modelError, true);
            }
        }
    }
    return logLikelihood;
} // MetropolisHastings::_loglikelihood

// --- --- 
// --- PUBLIC MEMBER FUNCTIONS
// --- ---

void MetropolisHastings::inferParameters (const unsigned int nPost) noexcept {
    _posterior.resize(_totalProposed + nPost);
    
    if (_totalProposed == 0) {
        ModelParameters theta0; // initialise at true parameters
        // theta.priorDraw(); 
        *_posterior.begin() = theta0; 
        _prevThetaLogLikelihood = theta0.priorLogDensity(); 
        _simulateData(theta0);
        _prevSimulatedData = _simulatedDataStar;
        _prevDataLogLikelihood = _logLikelihood(_prevSimulatedData, theta0.getModelError());
    }

    for (auto iter=(_posterior.begin() + _totalProposed + 1); iter<_posterior.end(); ++iter) {
        
        ModelParameters thetaStar;
        _proposalDistribution(thetaStar, iter-1);
        _simulateData(thetaStar);

        const double thetaStarLogLikelihood = thetaStar.priorLogDensity();
        const double simulatedDataLogLikelihood = _logLikelihood(_simulatedDataStar, thetaStar.getModelError());
        // const double prevDataLoglikelihood = _logLikelihood(*postIter);
        
        // prior proposal
        const double logAcceptanceFraction = simulatedDataLogLikelihood - _prevDataLogLikelihood + thetaStarLogLikelihood - _prevThetaLogLikelihood;
        const double logAcceptanceProb = std::min(0.0, logAcceptanceFraction);
        const double logRand = std::log(runif01());
        
        _totalProposed++;
        if (logRand < logAcceptanceProb) {
            *(iter) = thetaStar;
            _prevThetaLogLikelihood = thetaStarLogLikelihood;
            _prevDataLogLikelihood = simulatedDataLogLikelihood;
            _prevSimulatedData = _simulatedDataStar;
            _totalAccepted++;
        } else {
            *(iter) = *(iter-1);
        }
    }
} // ModelParameters::inferParameters


void MetropolisHastings::modelPrediction (std::ofstream& file, const unsigned int nRep, const std::string& sep) const noexcept {

    /*
        Using parameter values in the `_posterior` vector simulate the system to make find a posterior predictive. The system is automatically simulated up to 100 years and the state is saved at 100 equi-distant time-intervals up to Tmax.
    */
    double tMax = 100.0;
    // update the system at 100 equi-distant time points between 0.0 and tMax
    const double deltaT = tMax / 100.0;
    
    unsigned int simID = 0;
    
    for (const ModelParameters& theta : _posterior) {
        for (unsigned int i=0; i<nRep; ++i) {
            const std::string initLine = std::to_string(simID) + "." +  std::to_string(i) + sep;
            
            stochasticFibre mySystem = initialiseFibre(theta);
            mySystem.printToFile(file, sep, initLine);
            
            while (mySystem.getSystemTime()<tMax) {
                gillespieSimulation(mySystem, deltaT);
                mySystem.printToFile(file, sep, initLine);
            }
            simID++;
        }
    }
}

void MetropolisHastings::printPosteriorToFile (std::ofstream& file, const std::string& sep="\t") const {
    file << "log_genRate" << sep << "log_mutationProb" << sep << "log_defThreshold" << sep << "log_modelError" << "\n";
    
    for (auto it=_posterior.begin(); it!=_posterior.end(); ++it) {
        file << it->getLogGenRate() << sep << it->getLogMutationProb() << sep << it->getLogDefThreshold() << sep << it->getLogModelError() << "\n";
    }
}

void MetropolisHastings::proportionPrediction (std::ofstream& file, const unsigned int nRep, const std::string& sep) const noexcept {
    // make posterior predictions based upon parameter values
    // find the maximum _timeData value, use this as the simulation maximum
    // change this at a later date to be able to decide tMax and deltaT's
    double tMax = 100.0;
    const unsigned int nPred = 100;
    // update the system at 100 equi-distant time points between 0.0 and tMax
    const double deltaT = tMax / (double)nPred;
    
    unsigned int simID = 0;
    unsigned int thetaID = 0;
    for (const ModelParameters& theta : _posterior) {
        
        for (unsigned int i=0; i<nRep; ++i) {
            std::vector<double> defCount(nPred);
            
            for (unsigned int j=0; j<_nSims; ++j) {
                stochasticFibre mySystem = initialiseFibre(theta);
                unsigned int tt = 0;
                while (tt<nPred) {
                    // simulate over next time step
                    gillespieSimulation(mySystem, deltaT);
                    // Check mutation passed the pathogenic threshold
                    if (mySystem.getMutationLoad() > theta.getDefThreshold()) {
                        defCount[tt]++;
                        /*
                        // The fibre is "deficient" indefinitely, so add one to deficient count for all future time points
                        for (std::size_t l=tt; l<nPred; ++l) {
                            defCount[l]++;
                        }
                        // increase tt to stop the simulation and start again
                        tt = nPred;
                        */
                    }
                    tt++;
                }
            }
            const std::string predID_str = std::to_string(thetaID) + "." + std::to_string(i);
            for (std::size_t i=0; i<nPred; ++i) {
                file << predID_str << sep << (double)i*(double)deltaT << sep << (double)defCount[i] / (double)nPred << "\n";
            }
        }
        thetaID++;
    }
}


void proportionPrediction (std::ofstream& file, const std::string& thetaFile, const unsigned int nRep, const unsigned int nSim, const std::string& sep, const double deltaT, const double tMax) {
    
    std::vector<std::vector<double>> posterior = readMatrixFromFile<double>(thetaFile);
    const unsigned int nPred = tMax / deltaT;
    unsigned int simID = 0;

    for (std::size_t i=0; i<posterior.size(); ++i) {
        ModelParameters theta(&posterior[i][0], true);
        for (unsigned int j=0; j<nRep; ++j) {
            std::vector<unsigned int> defCount(nPred);

            for (unsigned int k=0; k<nSim; ++k) {
                stochasticFibre mySystem = initialiseFibre(theta);
                unsigned int tt = 0;
                while (mySystem.getSystemTime() <= tMax) {
                    gillespieSimulation(mySystem, deltaT);
                    if (mySystem.getMutationLoad() > theta.getDefThreshold()) {
                        defCount[tt]++;
                        /*
                        // The fibre is "deficient" indefinitely, so add one to deficient count for all future time points
                        for (std::size_t l=tt; l<nPred; ++l) {
                            defCount[l]++;
                        }
                        // increase tt to stop the simulation and start again
                        tt = nPred;
                        */
                    }
                    tt++;
                }
            }

            const std::string predID_str = std::to_string(i) + "." + std::to_string(j);
            for (std::size_t ii=0; ii<nPred; ++ii) {
                file << predID_str << sep << (double)ii*deltaT << sep << (double)defCount[ii] / nSim << "\n";
            }
        }
    }
}
#endif /* MetropolisHastings_h */
