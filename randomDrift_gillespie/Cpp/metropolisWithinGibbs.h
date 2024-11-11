//
//  metropolisHastings.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 24/05/2024.
//

#ifndef metropolisHastings_h
#define metropolisHastings_h

#include <cmath>
#include <vector>

#include "dataClass.h"
#include "modelParameters.h"
#include "simulate.h"

class MetropolisGibbs {
private:
    std::vector<ProportionDeficientData> _proportionData;
    std::vector<ProportionDeficientData> _prevSimulatedData;
    std::vector<ProportionDeficientData> _simulatedDataStar;

    unsigned int _nPat;
    unsigned int _nSims;
    double _halfTotalObs;
    double _inv_2nSims;

    ModelParameters _thetaStar;
    double _prevDataLogLikelihood = 1e99;
    double _prevThetaLogLikelihood = 1e99;

    std::vector<ModelParameters> _posterior;
    unsigned int _totalProposed;
    unsigned int _totalAccepted;
    std::vector<std::vector<double>> _proposalCovariance;
    std::vector<std::vector<double>> _proposalCovariance_cholesky;
    
    void _genericSetUp ();
    void _simulateData (const ModelParameters& theta) noexcept;
    double _logLikelihood(const std::vector<ProportionDeficientData>& simData, const double modelError) const;
    void _proposalDistribution (ModelParameters& thetaStar, const std::vector<ModelParameters>::const_iterator prevTheta) const noexcept ;
    void _updateModelParameters (const std::vector<ModelParameters>::iterator postIter) ;
    void _updateModelPrecision (const std::vector<ModelParameters>::iterator postIter) ;

public:
    ~ MetropolisGibbs () {}
    MetropolisGibbs (const unsigned int nSims) : _nSims(nSims) {
        _genericSetUp();
    }
    MetropolisGibbs (const unsigned int nSims, const std::vector<ProportionDeficientData> data) : _nSims(nSims), _proportionData(data) {
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
    
    void printPosteriorToFile (std::ofstream& file, const std::string& sep="\t") const {
        file << "log_genRate" << sep << "log_mutationProb" << sep << "log_defThreshold" << sep << "modelPrecision" << "\n";
        
        for (auto it=_posterior.begin(); it!=_posterior.end(); ++it) {
            file << it->getLogGenRate() << sep << it->getLogMutationProb() << sep << it->getLogDefThreshold() << sep << it->getModelPrecision() << "\n";
        }
    }
};

// --- --- --- --- 
// meptropolisWithinGibbs MEMBER FUNCTIONS 
// --- --- --- ---

// --- ---
// --- PRIVATE MEMBER FUNCTIONS
// --- --- 
void MetropolisGibbs::_genericSetUp () {
    _nPat = static_cast<unsigned int>(_proportionData.size());

    for (auto iter = _proportionData.begin(); iter!=_proportionData.end(); ++iter) {
        _halfTotalObs += iter -> getTotalDataSize();
        _simulatedDataStar.emplace_back(1, iter->getTimeVector(), _nSims);
        _prevSimulatedData.emplace_back(1, iter->getTimeVector(), _nSims);
        
        iter->setNsim(_nSims);
        iter->calculateDataTransform();
        iter->checkDataIsFinite();
    }
    _halfTotalObs /= 2.0;
    _inv_2nSims = 0.5 / (double)_nSims;

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
} // MetropolisGibss:_genericSetUp

void MetropolisGibbs::_proposalDistribution (ModelParameters& thetaStar, const std::vector<ModelParameters>::const_iterator prevTheta) const noexcept {
    /*
        Propose new values for the stochastic model parameters.
        Joint random walk proposal.
    */
    bool withinSupport = false;

    while (!withinSupport) {

        const std::vector<double> deltaParam = rMultivariateNormal_cholesky({0,0,0}, _proposalCovariance_cholesky);
        thetaStar.shiftLogSKMparameters(deltaParam[0], deltaParam[1], deltaParam[2]);
        withinSupport = thetaStar.checkLogSupport();
        thetaStar.updateLogs(true);
    }
} // MetropolisGibbs::_proposalDistribution


void MetropolisGibbs::_simulateData (const ModelParameters& theta) noexcept {
    for (auto iter=_simulatedDataStar.begin(); iter!=_simulatedDataStar.end(); ++iter) {
        iter->zeroCountData();
        // for (auto iter=_proportionData.begin(); iter!=_proportionData.end(); ++iter)
        for (std::size_t j=0; j<iter->getNsim(); ++j) {
            stochasticFibre mySystem = initialiseFibre(theta);
            std::size_t tt = 0;
            // while (tt<_timeDataSteps[i].size()) {
            while (tt<iter->getNobs()) {
                // simulate over next time step
                tauLeapSimulation(mySystem, iter->getTimeJump(tt), 100);
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
} // MetropolisGibbs::_simulateData

double MetropolisGibbs::_logLikelihood (const std::vector<ProportionDeficientData>& simData, const double modelError) const {
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
} // MetropolisGibbs::_loglikelihood

void MetropolisGibbs::_updateModelPrecision (const std::vector<ModelParameters>::iterator postIter) {
    /*
        Gibbs update the model error parameter.
    */
    double squaredError = 0.0;
    for (std::size_t i = 0; i<_nPat; ++i) {
        for (unsigned int k=0; k<_proportionData[i].getNobs(); ++k) {
            const double trueTransformedValue = _prevSimulatedData[i].getTransformedValue(k, 0);
            for (unsigned int j=0; j<_proportionData[i].getNrep(k); ++j) {
                const double error =  (_proportionData[i].getTransformedValue(k, j) -  trueTransformedValue);
                squaredError += error*error;
            }
        }
    }
    if (!isfinite(squaredError)) {
        throw std::runtime_error("Squared error is non-finite.");
    }
    const double shape = priorParameters::modelPrecision_shape + _halfTotalObs ;
    const double scale = 1.0 / (priorParameters::modelPrecision_rate + 0.5*squaredError);
    std::gamma_distribution<double> rgamma (shape, scale);
    const double newModelPrecision = rgamma(gen);
    postIter->setModelPrecision(newModelPrecision);
} // MetropolisGibbs::_updateModelError

void MetropolisGibbs::_updateModelParameters (const std::vector<ModelParameters>::iterator postIter) {
    /*
        Update stochastic model parameters using a metropolis step. 
    */
   // Propose new parameters
    ModelParameters thetaStar(*postIter);

    _proposalDistribution(thetaStar, postIter-1);
    
    // Simulate true proportion of deficiency
    _simulateData(thetaStar);

    // calculate prior log-likelihood of proposed values
    const double thetaStarLogLikelihood = thetaStar.priorSKMdensity();
    
    // calculate log-likelihood of data given current simulated true state
    const double dataStarLogLikelihood = _logLikelihood(_simulatedDataStar, thetaStar.getModelError() );

    const double logAcceptanceFraction = dataStarLogLikelihood - _prevDataLogLikelihood + thetaStarLogLikelihood - _prevThetaLogLikelihood;
    const double logAcceptanceProb = std::min(0.0, logAcceptanceFraction);
    const double logRand = std::log(runif01());

    if (logRand < logAcceptanceProb) {
        postIter->setSKMparameters(thetaStar);
        _prevThetaLogLikelihood = thetaStarLogLikelihood;
        _totalAccepted++;
        _prevSimulatedData = _simulatedDataStar;
        _prevDataLogLikelihood = dataStarLogLikelihood;
    } else {
        postIter -> setSKMparameters( *(postIter-1) );
    }
} // MetropolisGibbs::_updateModelParameters

//
// --- PUBLIC MEMBER FUNCTIONS
//

void MetropolisGibbs::inferParameters (const unsigned int nPost) noexcept {
    _posterior.resize(_totalProposed + nPost);
    if (_totalProposed == 0) {
        ModelParameters theta0; // initialise at true parameters
        theta0.priorDraw(); 
        *_posterior.begin() = theta0;
        _prevThetaLogLikelihood = theta0.priorSKMdensity(); 
        _simulateData(theta0);
        _prevSimulatedData = _simulatedDataStar;
        _prevDataLogLikelihood = _logLikelihood(_prevSimulatedData, theta0.getModelError());
    }

    for (auto iter=(_posterior.begin()+_totalProposed+1); iter!=_posterior.end(); ++iter) {
        _updateModelPrecision(iter);
        _updateModelParameters(iter);
        _totalProposed++;
    }
} // Metropolis::inferParameters


void MetropolisGibbs::modelPrediction (std::ofstream& file, const unsigned int nRep, const std::string& sep) const noexcept {
// make posterior predictions based upon parameter values
    // find the maximum _timeData value, use this as the simulation maximum
    // change this at a later date to be able to decide tMax and deltaT's
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
} // MetropolisGibbs::modelPrediction

void MetropolisGibbs::proportionPrediction (std::ofstream& file, const unsigned int nRep, const std::string& sep) const noexcept {
    // make posterior predictions based upon parameter values
    // find the maximum _timeData value, use this as the simulation maximum
    // change this at a later date to be able to decide tMax and deltaT's
    const double tMax = 100.0;
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
} // Metropolis::proporitonPrediction

#endif /* metropolisHastings_h */
