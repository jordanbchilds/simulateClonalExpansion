
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <chrono>

#include "readFiles.h"
#include "metropolisWithinGibbs.h"

using namespace std;

int main () {
    srand(static_cast<unsigned int>(time(NULL)));
    
    const unsigned int nSim = 100;
    const unsigned int nOut = 500;
    
    const vector<vector<double>> proportionData = readMatrixFromFile<double>("Data/synData_SoS/proportionDataMatrix.txt");
    const vector<double> timeData = readVectorFromFile<double>("Data/synData_SoS/timeData.txt");

    std::vector<ProportionDeficientData> myDataVector;
    for (std::size_t i=0; i<proportionData.size(); ++i) {
        myDataVector.emplace_back(timeData[i], proportionData[i]);
        myDataVector[i].printData();
    }

    // const vector<vector<double>> covarianceMat = readMatrixFromFile<double>("./Output/INFERENCE_SYNDATA/PROPOSAL_COVARIANCE_MATRIX/TL100_myPrior.txt");

    const vector<vector<double>> covarianceMat ({
            {priorParameters::logGenRate_proposalSD, 0.0, 0.0},
            {0.0, priorParameters::logMutationProb_proposalSD, 0.0},
            {0.0, 0.0, priorParameters::logDefThreshold_proposalSD} });

    const std::string sep = "\t";
    
    // --- PRIOR
    ModelParameters theta;
    ofstream priorFile("Output/INFERENCE_SYNDATA_SoS/PRIOR/myPrior.txt");
    priorFile << "log_genRate" << sep << "log_mutationProb" << sep << "log_defThreshold" << sep << "modelPrecision" << "\n";
    for (int i=0; i<10000; ++i) {
        theta.priorDraw();
        theta.printToFile(priorFile);
    }
    priorFile.close();
    
    // --- PARAMETER INFERENCE 
    // metropolisHastings myInference = metropolisHastings(nSim, timeData, proportionData);
    
    MetropolisGibbs myInference(nSim, myDataVector);

    myInference.setProposalCovarianceMatrix(covarianceMat);

    myInference.inferParameters(nOut);
    std::cout << myInference.getAcceptanceRatio() << std::endl;
    
    std::ofstream posteriorFile("Output/INFERENCE_SYNDATA_SoS/POST/TL100_myPrior_init.txt");
    myInference.printPosteriorToFile(posteriorFile);
    posteriorFile.close();

    
    /*
    // --- MODEL PREDICTION
    ofstream predictionFile("Output/Inference/modelPrediction.txt");
    myInference.modelPrediction(predictionFile);
    predictionFile.close();
   
    // --- PROPORTION PREDICTION
    const string thetaFile = "Output/Inference/posterior_gillespie.txt";
    ofstream file("Output/Inference/proportionPrediction_gillespie.txt");
    const double deltaT = constants::SECONDS_IN_YEAR / 2;
    const double tMax = constants::SECONDS_IN_YEAR * 65;
    proportionPrediction(file, thetaFile, 1000, 100, "\t", deltaT, tMax);
    file.close();
    */
}


/*
  // METROPOLIS-HASTINGS INFERENCE
    const vector<double> timeYears({35, 38, 40, 42, 45, 50, 60, 80});
    vector<double> timeSeconds(timeYears.size(), 0.0);
    for (size_t i=0; i<timeYears.size(); ++i) {
        timeSeconds[i] = timeYears[i]*secondsInYear;
    }
    
    const vector<vector<double>> proportionData = readMatrixFromFile<double>("./Output/syntheticData/synData_mat.txt");
    // const vector<vector<double>> covarianceMat = readMatrixFromFile<double>("./Output/Inference/PROPOSAL_COVARIANCE_MATRIX/proposalVariance_data10_Pr50.txt");
    
    const vector<vector<double>> covarianceMat({
        {priorParameters::logGenRate_proposalSD, 0.0, 0.0, 0.0},
        {0.0, priorParameters::logMutationProb_proposalSD, 0.0, 0.0},
        {0.0, 0.0, priorParameters::logDefThreshold_proposalSD, 0.0},
        {0.0, 0.0, 0.0, priorParameters::logModelError_proposalSD}
    });

    const std::string sep = "\t";

    // --- PRIOR
    modelParameters theta;
    ofstream priorFile("Output/Inference/PRIOR/prior1.txt");
    priorFile << "log_genRate" << sep << "log_mutationProb" << sep << "log_defThreshold" << sep << "log_modelError" << "\n";
    for (int i=0; i<10000; ++i) {
        theta.priorDraw();
        theta.printToFile(priorFile);
    }
    priorFile.close();

    // --- PARAMETER INFERENCE 
    
    metropolisHastings myInference = metropolisHastings(nSim, timeDataSeconds, proportionData);
    myInference.setProposalCovarianceMatrix(covarianceMat);
    myInference.inferParameters(nOut);
    std::cout << myInference.getAcceptanceRatio() << std::endl;

    std::ofstream posteriorFile("Output/Inference/POST/data100_TL100_Pr10.txt");
    myInference.printPosteriorToFile(posteriorFile);
    posteriorFile.close();

    // --- MODEL PREDICTION
    ofstream predictionFile("Output/Inference/modelPrediction.txt");
    myInference.modelPrediction(predictionFile);
    predictionFile.close();
   
    // --- PROPORTION PREDICTION
    const string thetaFile = "Output/Inference/posterior_gillespie.txt";
    ofstream file("Output/Inference/proportionPrediction_gillespie.txt");
    const double deltaT = constants::SECONDS_IN_YEAR / 2;
    const double tMax = constants::SECONDS_IN_YEAR * 65;
    proportionPrediction(file, thetaFile, 1000, 100, "\t", deltaT, tMax);
    file.close();
 */