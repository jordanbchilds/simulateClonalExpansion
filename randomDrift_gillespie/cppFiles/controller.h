//
//  controller.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 15/05/2024.
//

#ifndef controller_h
#define controller_h

class LinearController {
private:
    double _target;
    double _baseRate;
    double _controlParameter;
public:
    LinearController () : _target(0), _controlParameter(0), _baseRate(0) {}
    ~ LinearController () {}
    LinearController (const double target, const double baseRate, const double controlParam) : _target(target), _baseRate(baseRate), _controlParameter(controlParam) {}

    double calculateReplicationRate (const double currentState) {
        return _baseRate + _controlParameter*(_target - currentState);
    }
};

class replicationRate_controller {
private:
    double _target;
    double _upperBound;
    double _lowerBound;
    double _genRate = trueParameters::genRate;
    double _inflationFactor = 20.0;
    double _defaultProb = 0.5;
    
public:
    replicationRate_controller () : _target(0), _upperBound(0), _lowerBound(0) {}
    ~ replicationRate_controller ()  {}
    
    replicationRate_controller(const double genRate, const double target, const double upr, const double lwr) : _genRate(genRate), _target(target), _upperBound(upr), _lowerBound(lwr) {}
    
    double calc_replicationRate (const int systemPopulation) {
        int error =  _target - systemPopulation;
        
        if (error <= _lowerBound) return _genRate / _inflationFactor;
        if (error >= _upperBound) return _genRate * _inflationFactor;
        
        double probRep = _defaultProb;
        if (error > 0) {
            probRep = (0.5 * error) / (_upperBound) + _defaultProb;
        } else if ( error < 0 ) {
            probRep = (-0.5 * error) / (_lowerBound) + _defaultProb;
        }
        return (probRep * _genRate) / ( (1.0 - probRep) );
    }
    
    // SETTERs
    constexpr void setTarget (const int target) noexcept { _target = target; }
    constexpr void setLowerbound (const int lwr) noexcept { _lowerBound = lwr; }
    constexpr void setUpperbound (const int upr) noexcept { _upperBound = upr; }
    constexpr void setDegradationRate (const double k) noexcept { _genRate = k; }
    constexpr void setDefaultProb (double rr) noexcept { _defaultProb = rr; }
    
    // GETTERs
    double getTarget () const noexcept { return _target; }
    double getLowerBound () const noexcept { return _lowerBound; }
    double getUpperBound () const noexcept { return _upperBound; }
};



#endif /* controller_h */
