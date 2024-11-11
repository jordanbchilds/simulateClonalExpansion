
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <chrono>

#include "readFiles.h"
#include "randomSampling.h"
#include "modelParameters.h"
#include "simulate.h"

using namespace std;

int main () {
    
    srand(static_cast<unsigned int>(time(NULL)));

    const string timeDataFp = "../Data/synData_RGD/timeDataMatrix.txt";
    vector<double> myTimeData;
    if (filesystem::exists(timeDataFp)) {
        // Read time data
        myTimeData = readVectorFromFile<double>(timeDataFp);
        cout << "Time data read from: " << timeDataFp << "\n";
    } else {
        for (int i=0; i<6; ++i) {
            myTimeData.push_back( (2 + i) * 10 * constants::SECONDS_IN_YEAR );
        }
        cout << "Unable to read time data file. Default times used: 20, 30, 40, 50, 60, 70 years.\n";
    }

    const unsigned int nTimes = myTimeData.size();

    const unsigned int nSim = 100;
    const unsigned int nRep = 100;

    cout << "Number of repitions for each time point: " << nRep << ".\n";
    cout << "Number of simulations per time point: " << nSim << ".\n";

    const string sep = "\t";
    const string outputFile = "./Output/synData_RGD/rawData.txt";
    cout << "Saving output to: " << outputFile << "\n";
    // ofstream outFile(outputFile);

    // outFile << "simID" << sep << "systemTime" << sep << "nTrials" << sep << "defCount" << "\n";

    const ModelParameters trueTheta; 

    for (int i=0; i<nTimes; ++i) {
        for (int j=0; j<nRep; ++j) {
            unsigned int defCount = 0;
            for (int k=0; k<nSim; ++k) {
                stochasticFibre mySystem = initialiseFibre(trueTheta);
                gillespieSimulation(mySystem, myTimeData[i]);

                if (mySystem.getMutationLoad() > trueTheta.getDefThreshold()) {
                    defCount++ ;
                }
            }
            cout << i+1 << "-" << j+1 << sep << myTimeData[i] / constants::SECONDS_IN_YEAR << sep << nSim << sep << defCount << "\n";
            // outFile << i+1 << "-" << j+1 << sep << myTimeData[i] / constants::SECONDS_IN_YEAR << sep << nSim << sep << defCount << "\n";
        }
    }
    // outFile.close();
}
