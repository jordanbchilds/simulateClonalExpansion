//
//  readFiles.h
//  randomDrift_gillespie
//
//  Created by Jordan Childs on 09/07/2024.
//

#ifndef readFiles_h
#define readFiles_h

#include <sstream>

template<typename T>
std::vector<T> readVectorFromFile(const std::string& filename) {
    std::vector<T> vector;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return vector; // Return empty vector
    }
    std::string line;
    while (std::getline(file, line)) {
        /*
                T element;
        std::istringstream iss(line);
        T value;
        while (iss >> value) {
            vector.push_back(element);
        }
        */
        T element;
        std::istringstream iss(line);
        T value;
        while (iss >> value) {
            vector.push_back(value);
        }
    }
    file.close();
    return vector;
}

template<typename T>
std::vector<std::vector<T>> readMatrixFromFile(const std::string& filename) {

    std::vector<std::vector<T>> matrix;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return matrix; // Return empty matrix
    }
    std::string line;
    while (std::getline(file, line)) {
        std::vector<T> row;
        std::istringstream iss(line);
        T value;
        while (iss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }
    file.close();
    return matrix;
}


#endif /* readFiles_h */
