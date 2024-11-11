//
//  simulate.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 25/07/2024.
//

#ifndef simulate_h
#define simulate_h

#include <cmath>
#include <vector>
#include <fstream>
#include <random>

#include "controller.h"
#include "randomSampling.h"
#include "modelParameters.h"

class stochasticFibre {
private:
    const std::vector<std::vector<int>> _PRE = {{1,0}, {0,1}, {1,0}, {0,1}, {1,0}};
    const std::vector<std::vector<int>> _POST = {{0,0}, {0,0}, {2,0}, {0,2}, {1,1}};
    const std::vector<std::vector<int>> _STOI = {{-1,0}, {0,-1}, {1,0}, {0,1}, {0,1}};
    const unsigned int _reactionNumber = 5;
    const unsigned int _speciesNumber = 2;
    
    std::vector<double> _reactionRates;
    std::vector<double> _reactionHazards;
    std::vector<unsigned int> _populationVector;
    double _mutationProb;
    double _notMutationProb;
    unsigned int _systemPopulation;
    double _systemTime;
    double _nextEventTime;
    double _systemHazard;
    replicationRate_controller _controller;
    unsigned int _targetPopulation;
    double _outputTimeUnit = 60.0*60.0*24.0*365.0;
    double _baseReactionRate;
    
public:
    stochasticFibre () {}
    ~stochasticFibre () {}

    stochasticFibre (const ModelParameters& theta) {

        unsigned int initTotalPopulation = theta.getInitialPopulation();
        const double genRate = theta.getGenRate();
        const double mutationProb = theta.getMutationProb();
        
        replicationRate_controller myController(theta.getGenRate(), initTotalPopulation, 0.1*initTotalPopulation, -0.1*initTotalPopulation);
        myController.setDegradationRate(theta.getGenRate());

        _populationVector = {initTotalPopulation, 0};
        _reactionRates = {genRate, genRate, genRate*(1-mutationProb), genRate, genRate*mutationProb};
        _controller = myController;
        _systemTime = 0.0;
        _nextEventTime = 0.0;
        _systemPopulation = 0.0;
        _reactionHazards = std::vector<double>(_reactionNumber, 0.0);
        _targetPopulation = myController.getTarget();
        _baseReactionRate = _reactionRates[0];
        
        for (auto it=_populationVector.begin(); it!=_populationVector.end(); ++it) {
            _systemPopulation += *it;
        }
        calculateHazards();
        _mutationProb = _reactionRates[4] / (_reactionRates[2] + _reactionRates[4]);
        _notMutationProb = 1.0 - _mutationProb;
    }
    
    stochasticFibre (std::vector<unsigned int> initPopulationVector,
                    const std::vector<double>& reactionRates, 
                    const replicationRate_controller& controller) : _populationVector(initPopulationVector), _reactionRates(reactionRates), _controller(controller) {
        
        _systemTime = 0.0;
        _nextEventTime = 0.0;
        _systemPopulation = 0.0;
        _reactionHazards = std::vector<double>(_reactionNumber, 0.0);
        _targetPopulation = controller.getTarget();
        _baseReactionRate = _reactionRates[0];
        
        for (auto it=_populationVector.begin(); it!=_populationVector.end(); ++it) {
            _systemPopulation += *it;
        }
        calculateHazards();
        _mutationProb = reactionRates[4] / (reactionRates[2] + reactionRates[4]);
        _notMutationProb = 1.0 - _mutationProb;
    }
    
    void calculateHazards () {
        _systemHazard = 0.0;
        for (long unsigned int i=0; i<_reactionNumber; i++) {
            _reactionHazards[i] = _reactionRates[i];
            for (long unsigned int j=0; j<_speciesNumber; j++) {
                _reactionHazards[i] *= choose(_populationVector[j], _PRE[i][j]);
            }
            _systemHazard += _reactionHazards[i];
        }
    }
    
    void updateRates () noexcept {
        double newRate = _controller.calc_replicationRate(_systemPopulation);
        _reactionRates[2] = newRate*_notMutationProb;
        _reactionRates[3] = newRate;
        _reactionRates[4] = newRate*_mutationProb;
    }
    
    std::size_t simulateReactionType () const noexcept {
        double u = runif01() * _systemHazard;
        auto iter = 0;
        while (iter<_reactionNumber) {
            u -= *(_reactionHazards.begin() + iter);
            if (u < 0.0) return iter;
            iter++;
        }
        return iter;
    }
    
    void simulateReaction (const std::size_t reactionInd, const int nRep=1) noexcept {
        for (std::size_t j=0; j<_speciesNumber; j++) {
            const unsigned int populationChange =_STOI[reactionInd][j] * nRep;
            _populationVector[j] += populationChange;
            _systemPopulation += populationChange;
        }
    }
    void updateNextEventTime (const double timeStep) noexcept {
        _systemTime = _nextEventTime;
        _nextEventTime += timeStep;
    }
    
    // GETTERs
    double getSystemTime () const noexcept { return _systemTime; }
    double getNextEventTime () const noexcept { return _nextEventTime; }
    std::vector<unsigned int> getPopulations () const noexcept { return _populationVector; }
    unsigned int getSystemPopulation () const noexcept { return _systemPopulation; }
    double getMutationLoad () const noexcept { return (double)_populationVector[1] / (double)_systemPopulation; }
    double getHazard () const noexcept {return _systemHazard;}
    double getHazard (const std::size_t i) const noexcept { return _reactionHazards[i]; }
    long unsigned int getReactionNumber () const noexcept { return _reactionNumber; }
    long unsigned int getSpeciesNumber () const noexcept { return _speciesNumber; }
    const replicationRate_controller& getController () const noexcept { return _controller; }
    unsigned int getTargetPopulation () const noexcept { return _targetPopulation; }
    double getBaseReactionRate () const noexcept { return _baseReactionRate; }
    // SETTERs
    void setOutputTimeUnit (const double timeUnit) noexcept { _outputTimeUnit = timeUnit; }
    void setController (replicationRate_controller controller) noexcept { _controller = controller; }
    unsigned int getWild () const noexcept { return _populationVector[0]; }
    unsigned int getMutant () const noexcept { return _populationVector[1]; }
    void setNextEventTime (const double tt) noexcept { _nextEventTime = tt; }
    void setSystemTime (const double tt) noexcept { _systemTime = tt; }
    
    // PRINTERs
    void printSummary (const std::string& sep=", ") const {
        std::cout << _systemTime / _outputTimeUnit ;
        for (auto it=_populationVector.begin(); it!=_populationVector.end(); ++it) {
            std::cout << sep << *it;
        }
        std::cout << "\n";
    }
    
    void printToFile (std::ofstream& file, const std::string& sep="\t", const std::string& initLine="") const {
        if (initLine!="") {
            file << initLine << sep << _systemTime / _outputTimeUnit;
        } else {
            file << _systemTime / _outputTimeUnit;
        }
        
        for (auto it=_populationVector.begin(); it!=_populationVector.end(); ++it) {
            file << sep << *it;
        }
        file << "\n";
    }
};

stochasticFibre initialiseFibre (const ModelParameters& theta) {
    
    unsigned int initTotalPopulation = theta.getInitialPopulation();
    std::vector<unsigned int> initPopulationVector({initTotalPopulation, 0});
    const double genRate = theta.getGenRate();
    const double mutationProb = theta.getMutationProb();
    const std::vector<double> reactionRates({genRate, genRate, genRate*(1-mutationProb), genRate, genRate*mutationProb});
    
    replicationRate_controller myController(theta.getGenRate(), initTotalPopulation, 0.1*initTotalPopulation, -0.1*initTotalPopulation);
    myController.setDegradationRate(theta.getGenRate());

    stochasticFibre mySystem(initPopulationVector, reactionRates, myController);
    
    return mySystem;
}


void gillespieSimulation (stochasticFibre& sys, const double timeStep) {
        double Tmax = sys.getSystemTime() + timeStep;
        sys.calculateHazards();
        
        if (sys.getNextEventTime() <= sys.getSystemTime()) {
            const double timeJump = rexp(sys.getHazard());
            sys.setNextEventTime(sys.getSystemTime() + timeJump);
        }
    
        while (sys.getNextEventTime() < Tmax) {
            sys.updateRates();
            const std::size_t reactInd = sys.simulateReactionType();
            sys.simulateReaction(reactInd);
            sys.calculateHazards();
            const double timeJump = rexp(sys.getHazard());
            sys.updateNextEventTime(timeJump);
            if (sys.getHazard() < 1e-10) sys.setSystemTime(1e99);
        }
        sys.setSystemTime(Tmax);
}

inline double calculateLeap (const stochasticFibre& sys, const double increaseFactor) {
    /*
     Define a tau (timeStep over which to simulate) based on the _systemhazard.
     In exact simulation the expected inter-event-time is 1 / systemHazard. Increasing the expected time by a factor of k, we'd expect k reactions to take place
     */
    return increaseFactor / (2.0 * sys.getBaseReactionRate() * sys.getSystemPopulation());
} // function calculateLeap


void tauLeapSimulation (stochasticFibre& sys, const double timeStep, const double eventTimeIncreaseFactor=100.0) {
        /*
         simulate the system over a time length of `timeStep`.
         Ends the simulation at _systemTime + timeStep, saving the time of the next event to continue the simulation.
         */
        const double Tmax = sys.getSystemTime() + timeStep;
        sys.calculateHazards();
        double timeJump = calculateLeap(sys, eventTimeIncreaseFactor);
        
        // if _nextEventTime has not been set:
        if (sys.getNextEventTime() <= sys.getSystemTime()) {
            sys.setNextEventTime(sys.getSystemTime() + timeJump);
        }
        // if _systemTime is between two tau leaps
        // --- move _systemTime back to the previous tau leap time point
        // --- to be able to simulate over that time step.
        // --- adding _tau to the existing _systemTime would be leaving a window
        // --- of un simulated time. It is possible to simulate over that time window but
        // --- it would be a different time-step to the usual _tau. However, it would be
        // --- smaller than _tau so probably fine. This current method is fine for a
        // --- constant _tau (or if the _tau is saved).
        if (sys.getSystemTime() + timeJump > sys.getNextEventTime()) {
            sys.setSystemTime( sys.getNextEventTime() - timeJump);
        }
        
        while (sys.getNextEventTime() <= Tmax) {
            sys.updateRates();
            for (std::size_t i=0; i<sys.getReactionNumber(); ++i) {
                std::poisson_distribution<unsigned int> pd(timeJump*sys.getHazard(i));
                const unsigned int reactRep = pd(gen);
                sys.simulateReaction(i, reactRep);
            }
            sys.calculateHazards();
            timeJump = calculateLeap(sys, eventTimeIncreaseFactor);
            sys.updateNextEventTime(timeJump);
            if (sys.getHazard() < 1e-10) sys.setSystemTime( 1e99);
        }
        sys.setSystemTime(Tmax);
}

#endif /* simulate_h */
