#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

// Helper functions
double calculateStd(const std::vector<double>& vec) {
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    double sq_sum = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0,
                                       [](double acc, double x) { return acc + x; },
                                       [mean](double x, double y) { return (x - mean) * (y - mean); });
    return std::sqrt(sq_sum / vec.size());
}

struct Data {
    std::vector<double> time;
    std::vector<double> high;
    std::vector<double> low;
    std::vector<double> close;
    int N;
};

std::time_t toEpoch(const std::string& date) {
    struct std::tm tm;
    std::istringstream ss(date);
    ss >> std::get_time(&tm, "%m/%d/%Y");
    std::time_t epoch = std::mktime(&tm);
    return epoch;
}

Data loadData(const std::string& filename) {
    Data data;
    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        throw std::runtime_error("File not found");
    }
    // Skipping header
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double t, h, l, c;
        iss >> t >> h >> l >> c; // Assuming CSV format is known and consistent
        data.time.push_back(t);
        data.high.push_back(h);
        data.low.push_back(l);
        data.close.push_back(c);
    }
    data.N = data.time.size();
    file.close();
    return data;
}

int main() {
    std::string dataFile = "/home/yui/Downloads/HO-5min.csv";

    std::time_t inSampleStart = toEpoch("10/02/2007");
    std::time_t inSampleEnd = toEpoch("10/01/2017");
    std::time_t outSampleStart = toEpoch("10/02/2017");
    std::time_t outSampleEnd = toEpoch("04/21/2023");

    int barsBack = 17001;
    double slpg = 65;
    double PV = 64000;
    double E0 = 200000;

    // Assuming MATLAB-like vector ranges
    std::vector<int> Length = {10000}; // Only one value in the MATLAB example
    std::vector<double> StopPct = {0.015}; // Only one value in the MATLAB example

    std::vector<std::vector<double>> resultInSample;
    std::vector<std::vector<double>> resultOutSample;

    Data data = loadData(dataFile);

    std::vector<double> E(data.N, E0);
    std::vector<double> trades(data.N, 0);
    std::vector<double> HH(data.N, std::numeric_limits<double>::lowest());
    std::vector<double> LL(data.N, std::numeric_limits<double>::max());

    for (int i = barsBack; i < data.N; i++) {
        for (int j = i - barsBack; j < i; j++) {
            HH[i] = std::max(HH[i], data.high[j]);
            LL[i] = std::min(LL[i], data.low[j]);
        }
    }

    double position = 0;
    for (int i = 1; i < data.N; i++) {
        double delta = 0;
        if (position == 0) {
            if (data.high[i] >= HH[i]) {
                position = 1;
                delta = -slpg / 2;
            } else if (data.low[i] <= LL[i]) {
                position = -1;
                delta = -slpg / 2;
            }
        } else if (position == 1 && data.low[i] <= LL[i]) {
            position = 0;
            delta = PV * (data.close[i] - data.close[i - 1]);
        } else if (position == -1 && data.high[i] >= HH[i]) {
            position = 0;
            delta = PV * (data.close[i - 1] - data.close[i]);
        }
        E[i] = E[i - 1] + delta;
    }

    // Print equity for demonstration
    std::cout << "Final equity: " << E.back() << std::endl;

        std::vector<double> PnL(data.N, 0.0);
    std::fill(PnL.begin(), PnL.begin() + barsBack, 0);  // equivalent to MATLAB zeros
    for (size_t k = barsBack; k < data.N; ++k) {
        PnL[k] = E[k] - E[k - 1]; // Calculate profit and loss
    }

    std::vector<std::vector<double>> resultsInSample(Length.size(), std::vector<double>(4));
    std::vector<std::vector<double>> resultsOutSample(Length.size(), std::vector<double>(4));

    // Loop over Length and StopPct not shown, assuming indices i and j are defined
    resultsInSample[i][0] = E[indInSample2] - E[indInSample1];
    resultsInSample[i][1] = *std::min_element(DD.begin() + indInSample1, DD.begin() + indInSample2);
    resultsInSample[i][2] = calculateStd(std::vector<double>(PnL.begin() + indInSample1, PnL.begin() + indInSample2));
    resultsInSample[i][3] = std::accumulate(trades.begin() + indInSample1, trades.begin() + indInSample2, 0);

    resultsOutSample[i][0] = E[indOutSample2] - E[indOutSample1];
    resultsOutSample[i][1] = *std::min_element(DD.begin() + indOutSample1, DD.begin() + indOutSample2);
    resultsOutSample[i][2] = calculateStd(std::vector<double>(PnL.begin() + indOutSample1, PnL.begin() + indOutSample2));
    resultsOutSample[i][3] = std::accumulate(trades.begin() + indOutSample1, trades.begin() + indOutSample2, 0);

    std::cout << StopPct[j] << ": in/out: " << resultsInSample[i][0] << "/" << resultsOutSample[i][0] << ", " 
              << resultsInSample[i][1] << "/" << resultsOutSample[i][1] << ", " 
              << resultsInSample[i][2] << "/" << resultsOutSample[i][2] << ", " 
              << resultsInSample[i][3] << "/" << resultsOutSample[i][3] << std::endl;

    return 0;
}
