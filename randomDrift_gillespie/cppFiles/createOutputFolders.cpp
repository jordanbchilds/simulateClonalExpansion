
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace std;

int main () {
    cout << "Current path: " << std::__fs::filesystem::current_path() << endl;

    std::filesystem::create_directory("../Output");

    cout << "Directory created: " << "../Output" << "\n";
    return 0;
}
