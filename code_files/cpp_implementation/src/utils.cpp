#include "utils.hpp"
#include <fstream>

void writeCSV(const std::string& filename, const std::vector<std::vector<double>>& grid) {
    std::ofstream file(filename);
    for (const auto& row : grid) {
        for (size_t i = 0; i < row.size(); ++i)
            file << row[i] << (i + 1 < row.size() ? "," : "");
        file << "\n";
    }
    file.close();
}
