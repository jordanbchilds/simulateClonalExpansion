#ifndef dataClass_h
#define dataClass_h

#include <iostream>
#include <vector>

#include "randomSampling.h"

class ProportionDeficientData {
private:
    unsigned int _totalDataSize; // total data size, to be saved on the heap.
    unsigned int _nObs; // Number of time points in dataset for observations
    unsigned int _nSim; // number of simulations used to simulate data
    const double _dataTransformDelta = 0.5 / 100000.0; // data transform parameter

    double* _time; // time data
    double* _timeJump; // time jump data

    unsigned int* _nRep; // number observations at each time point 
    unsigned int* _obsInitIndex; // index for observations 

    unsigned int* _countData; // count data for simulations i.e. number of simulations which have passed the threshold out of the _nSim
    double* _proportionData; // proportion of deficiency i.e. _countData / _nSim for simulations or real observations 
    double* _transformedData; // transformed proportion data 

    std::string _ID;
public:
    // ----------------------
    // ---  CONSTRUCTORS  ---
    // ----------------------
    ProportionDeficientData () { 
        _totalDataSize = 0;
        _nObs = 0;
        _nSim = 0;
        _time = nullptr;
        _timeJump = nullptr;
        _nRep = nullptr;
        _obsInitIndex = nullptr;
         
        _countData = nullptr;
        _proportionData = nullptr; 
        _transformedData = nullptr;
        _ID = "";
    }
    // ----------------------------------------
    // ---  INITIALISE FOR SIMULATION DATA  ---
    // ----------------------------------------
    ProportionDeficientData (const unsigned int nRep, const double time, const unsigned int nSim) {
        /*
            Initialise for a simulations with a single time point. 
        */
        _nSim = nSim;
        _nObs = 1;
        _nRep = new unsigned int[_nObs];
        for (unsigned int i=0; i<_nObs; ++i) {
            *(_nRep + i) = nRep;
        }
        _totalDataSize = nRep ;

        _obsInitIndex = new unsigned int[_nObs];
        *(_obsInitIndex + 0) = 0;
        for (unsigned int i=1; i<_nObs; ++i) {
            *(_obsInitIndex + i) = nRep*(i-1);
        }

        _time = new double[_nObs];
        _timeJump = new double[_nObs];
        *_time = time;
        *_timeJump = time;
        _countData = new unsigned int[_totalDataSize];
        _proportionData = new double[_totalDataSize];
        _transformedData = new double[_totalDataSize];
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_countData + i) = -1;
            *(_proportionData + i) = -99.0;
            *(_transformedData + i) = -99.0;
        }
    }

    ProportionDeficientData (const unsigned int nRep, const std::vector<double> time, const unsigned int nSim) {
        _nSim = nSim;
        _nObs = time.size();
        _nRep = new unsigned int[_nObs];
        for (unsigned int i=0; i<_nObs; ++i) {
            *(_nRep + i) = nRep;
        }
        _totalDataSize = nRep * _nObs;

        _obsInitIndex = new unsigned int[_nObs];
        *(_obsInitIndex + 0) = 0;
        for (unsigned int i=1; i<_nObs; ++i) {
            *(_obsInitIndex + i) = nRep*(i-1);
        }

        _time = new double[_nObs];
        _timeJump = new double[_nObs];
        *_time = time[0];
        *_timeJump = time[0];
        for (std::size_t i=1; i<_nObs; ++i) {
            *(_time + 0) = time[i];
            *(_timeJump + 0) = time[i] - time[i-1];
        }

        _countData = new unsigned int[_totalDataSize];
        _proportionData = new double[_totalDataSize];
        _transformedData = new double[_totalDataSize];
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_countData + i) = -1;
            *(_proportionData + i) = -99.0;
            *(_transformedData + i) = -99.0;
        }
    }

    // --- INITIALISE FOR OBSERVED DATA
    ProportionDeficientData (const double time, const std::vector<double>& data) {
        _totalDataSize = data.size();
        _nObs = 1;
        _nSim = 0;
        
        _time = new double[_nObs];
        *(_time)= time;
        _timeJump = new double[_nObs];
        *(_timeJump) = *_time;

         _nRep = new unsigned int[_nObs];
        *(_nRep) = _totalDataSize;

        _obsInitIndex = new unsigned int[_nObs];
        *(_obsInitIndex) = 0;
        
        _countData = new unsigned int[_totalDataSize];
        _proportionData = new double[_totalDataSize];
        _transformedData = new double[_totalDataSize];
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_countData + i) = -1;
            *(_proportionData + i) = data[i];
            *(_transformedData + i) = -99.0;
        }
    }

    // DESTRUCTORS
    ~ProportionDeficientData () {
        delete[] _nRep;
        delete[] _obsInitIndex;
        delete[] _time;
        delete[] _timeJump;
        delete[] _countData;
        delete[] _proportionData;
        delete[] _transformedData;
        
    }
    
    // COPY OPERATOR
    ProportionDeficientData (const ProportionDeficientData& other) {
        _totalDataSize = other.getTotalDataSize();
        if (_totalDataSize>1000) {
            throw std::runtime_error("Data size suspiciously large.");
        }
        _nObs = other.getNobs();
        if (_nObs>1) {
            throw std::runtime_error("No. of observations suspiciously large.");
        }
        _nSim = other.getNsim();
        if (_nSim>100) {
            throw std::runtime_error("No. of simulations suspiciously large.");
        }
         _time = new double[_nObs];
         _timeJump = new double[_nObs];
         _nRep = new unsigned int[_nObs];
         _obsInitIndex = new unsigned int[_nObs];

        
        _countData = new unsigned int[_totalDataSize];
        _proportionData = new double[_totalDataSize];
        _transformedData = new double[_totalDataSize];
        
        *(_obsInitIndex) = 0;
        unsigned int cumRep = 0;
        std::size_t it = 0;
        for (std::size_t i=0; i<_nObs; ++i) {
            if (i>0) {
                cumRep += other.getNrep(i-1);
                *(_obsInitIndex + i) = cumRep;
            }
            
            *(_nRep + i) = other.getNrep(i);
            *(_time + i) = other.getTime(i);
            *(_timeJump + i) = other.getTimeJump(i);
           
            for (std::size_t j=0; j<other.getNrep(i); ++j) {
                *(_countData + it) = other.getCountValue(i,j);
                *(_proportionData + it) = other.getValue(i, j);
                *(_transformedData + it) = other.getTransformedValue(i, j);
                it++;
            }
            
        }
        _ID = other.getID();
    }

    // COPY ASSIGNMENT OPERATOR
    ProportionDeficientData operator=(const ProportionDeficientData& other) {
        if (this!=&other) {
            // check data sizes and delete and allocate as necessary
            if (_totalDataSize != other.getTotalDataSize()) {
                _totalDataSize = other.getTotalDataSize();
                delete[] _countData;
                delete[] _proportionData;
                delete[] _transformedData;

                _countData = new unsigned int[_totalDataSize];
                _proportionData = new double[_totalDataSize];
                _transformedData = new double[_totalDataSize];
            }
            if (_nObs!=other.getNobs()) {
                _nObs = other.getNobs();
                delete[] _nRep;
                delete[] _obsInitIndex;
                delete[] _time;
                delete[] _timeJump;
                _nRep = new unsigned int[_nObs];
                _time = new double[_nObs];
                _timeJump = new double[_nObs];
                _obsInitIndex = new unsigned int[_nObs];
                
            }
            _nSim = other.getNsim();
            // --- fill memory blocks
            *(_obsInitIndex + 0) = 0;
            unsigned int cumRep = 0;
            std::size_t it = 0;
            for (std::size_t i=0; i<_nObs; ++i) {
                if (i>0) {
                    cumRep += other.getNrep(i-1);
                    *(_obsInitIndex + i) = cumRep;
                }
                *(_nRep + i) = other.getNrep(i);
                *(_time + i) = other.getTime(i);
                *(_timeJump + i) = other.getTimeJump(i);

                for (std::size_t j=0; j<other.getNrep(i); ++j) {
                    *(_countData + it) = other.getCountValue(i, j);
                    *(_proportionData + it) = other.getValue(i, j);
                    *(_transformedData + it) = other.getTransformedValue(i, j);
                    it++;
                }
            }
            _ID = other.getID();
        }
        return *this;
    }

    void printData () const noexcept ;

    void zeroCountData () noexcept {
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_countData + i) = 0;
        } 
    }

    void resetData () noexcept {
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_countData + i) = -1;
            *(_proportionData + i) = -99.0;
            *(_transformedData + i) = -99.0;
        }
    }
    void increaseCount (const unsigned int timeIndex, const unsigned int obsIndex) noexcept {
       *(_countData + *(_obsInitIndex + timeIndex) + obsIndex) += 1;
    }
    void updateProportionData () noexcept {
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_proportionData + i) = *(_countData + i) / (double)_nSim;
            *(_transformedData + i) = _dataTransformFunction(_proportionData + i);
        }
    }
    void calculateDataTransform () {
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            *(_transformedData + i) = _dataTransformFunction(_proportionData + i);
            if (!isfinite(*(_transformedData + i))) {
                throw std::runtime_error("Transformation non-finite.");
            }
        }
    }
    void checkDataIsFinite () const {
        for (unsigned int i=0; i<_totalDataSize; ++i) {
            if (!isfinite(*(_proportionData + i))) {
                throw std::runtime_error("Proportion data is non-finite.");
            }
        }

        for (unsigned int i=0; i<_totalDataSize; ++i) {
            if (!isfinite(*(_transformedData + i))) {
                throw std::runtime_error("Transformed proportion data is infinite.");
            }
        }
    }
    // GETTERs
    unsigned int getDataSize () const noexcept { return _totalDataSize; }
    unsigned int getNobs () const noexcept { return _nObs; }

    unsigned int getNrep (const unsigned int i=0) const noexcept {
        _checkTimeIndex(i);
        return *(_nRep + i); 
    }

    unsigned int getNsim () const noexcept { return _nSim; }
    unsigned int getTotalDataSize () const noexcept { return _totalDataSize; }

    
    double getTime (const unsigned int timeIndex) const { 
        _checkTimeIndex(timeIndex);
        return *(_time + timeIndex); 
    }
    std::vector<double> getTimeVector () const noexcept { return std::vector<double>(_time, _time+_nObs); }
    double getTimeJump (const unsigned int timeIndex) const {
        _checkTimeIndex(timeIndex);
        return *(_timeJump + timeIndex);
    }
    unsigned int getCountValue(const unsigned int timeIndex, const unsigned int obsIndex) const {
        _checkDataIndex(timeIndex, obsIndex);
        return *(_countData + _dataIndexCalculator(timeIndex, obsIndex));  
    }
    double getValue (const unsigned int timeIndex, const unsigned int obsIndex) const { 
        _checkDataIndex(timeIndex, obsIndex);
        return *(_proportionData + _dataIndexCalculator(timeIndex, obsIndex)); 
    }
    
    double getTransformedValue (const unsigned int timeIndex, const unsigned int obsIndex) const {
        _checkDataIndex(timeIndex, obsIndex);
        return *(_transformedData + _dataIndexCalculator(timeIndex, obsIndex));
    }
    
    std::string getID () const noexcept { return _ID; }

    // SETTERS
    void setID (const std::string& id) noexcept { _ID = id; }
    void setNsim (const unsigned int nSim) noexcept { _nSim = nSim; }

private:
    void _checkTimeIndex (const unsigned int timeIndex) const {
        if (timeIndex > _nObs) {
            throw std::runtime_error("Index out of bounds.");
        }
    }
    
    void _checkDataIndex (const unsigned int timeIndex, const unsigned int obsIndex) const {
        if (obsIndex > *(_nRep + timeIndex)) {
            throw std::runtime_error("Index out of bounds.");
        }
    }
    unsigned int _dataIndexCalculator (const unsigned int timeIndex, const unsigned int obsIndex) const noexcept {
        return *(_obsInitIndex + timeIndex) + obsIndex;
    }
    
    double _dataTransformFunction (const double x) {
        return std::log( (x + _dataTransformDelta) / (1.0 - x + _dataTransformDelta) );
    }
    double _dataTransformFunction (const double* ptr) {
        return std::log( (*ptr + _dataTransformDelta) / (1.0 - *ptr + _dataTransformDelta) );
    }
};

// --- --- ---
// PATIENT DATA MEMBER FUNCTIONS
// --- --- ---

// --- PUBLIC MEMBER FUNCTIONS
void ProportionDeficientData::printData () const noexcept {
    std::cout << "nObs: " << _nObs << "\n";
    std::cout << "nSim: " << _nSim << "\n";
    std::cout << "Time: " << *(_time + 0) / constants::SECONDS_IN_YEAR;
    for (unsigned int i=1; i<_nObs; ++i) {
         std::cout << ", " << *(_time + i) / constants::SECONDS_IN_YEAR;
    }
    std::cout << "\n";

    std::cout << "Time Jumps: " << *(_timeJump + 0) / constants::SECONDS_IN_YEAR;
    for (unsigned int i=1; i<_nObs; ++i) {
         std::cout << ", " << *(_timeJump + i) / constants::SECONDS_IN_YEAR;
    }
    std::cout << "\n";

    std::cout << "Repitions: " << *(_nRep + 0);
    for (unsigned int i=1; i<_nObs; ++i) {
        std::cout << ", " << *(_nRep + i) ;
    }
    std::cout << "\n";

    std::cout << "Obsered data:\n";
    for (unsigned int i=0; i<_nObs; ++i) {
        std::cout << "Time: " << getTime(i) / constants::SECONDS_IN_YEAR << "\n"; 
        std::cout << getValue(i, 0) ;
        for (unsigned int j=1; j< *(_nRep+i); ++j) {
            std::cout << ", " << getValue(i, j);
        }
    }

    std::cout << "\n\nCount data:\n" ;
    for (unsigned int i=0; i<_nObs; ++i) {
        std::cout << "Time: " << getTime(i) / constants::SECONDS_IN_YEAR << "\n"; 
        std::cout << getCountValue(i, 0) ;
        for (unsigned int j=1; j< *(_nRep+i); ++j) {
            std::cout << ", " << getCountValue(i, j);
        }
    }
    
    std::cout << "\n\nTranformed data:\n" ;
    for (unsigned int i=0; i<_nObs; ++i) {  
        std::cout << "Time: " << getTime(i) / constants::SECONDS_IN_YEAR << "\n";
        std::cout << getTransformedValue(i, 0) ;
        for (unsigned int j=1; j< *(_nRep+i); ++j) {
            std::cout << ", " << getTransformedValue(i, j);
        }
    }
    std::cout << "\n\n";
} // ProportionDeficientData::printData

#endif /* dataClass_h */
