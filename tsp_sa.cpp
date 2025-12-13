// Build: g++ -O3 -std=c++17 tsp_sa.cpp -o tsp_sa
// Run:   ./tsp_sa cities23.dat 1

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct City {
    std::string name;
    double lat_deg; // latitude
    double lon_deg; // longitude
};

static constexpr double EARTH_RADIUS_KM = 6371.0;

static inline double deg2rad(double x) { return x * M_PI / 180.0; }

double haversine_km(double lat1_deg, double lon1_deg, double lat2_deg, double lon2_deg) {
    const double lat1 = deg2rad(lat1_deg);
    const double lon1 = deg2rad(lon1_deg);
    const double lat2 = deg2rad(lat2_deg);
    const double lon2 = deg2rad(lon2_deg);

    const double dlat = lat2 - lat1;
    const double dlon = lon2 - lon1;

    const double s1 = std::sin(dlat / 2.0);
    const double s2 = std::sin(dlon / 2.0);

    const double a = s1 * s1 + std::cos(lat1) * std::cos(lat2) * (s2 * s2);
    const double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return EARTH_RADIUS_KM * c;
}

static inline std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    const auto b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

static inline std::string strip_quotes(std::string s) {
    s = trim(s);
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\''))) {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

// Parse: longitude latitude "City Name"
std::vector<City> read_cities_lonlat_quoted(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open file: " + path);

    std::vector<City> cities;
    std::string line;

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);

        double lon = 0.0, lat = 0.0;
        if (!(iss >> lon >> lat)) {
            throw std::runtime_error("Could not read lon/lat from line: " + line);
        }

        std::string rest;
        std::getline(iss, rest);
        rest = strip_quotes(trim(rest));
        if (rest.empty()) rest = "City" + std::to_string((int)cities.size());

        City c;
        c.lon_deg = lon;
        c.lat_deg = lat;
        c.name = rest;
        cities.push_back(c);
    }

    if (cities.size() < 3) throw std::runtime_error("Need at least 3 cities for TSP.");
    return cities;
}

std::vector<std::vector<double>> build_distance_matrix(const std::vector<City>& cities) {
    const int n = (int)cities.size();
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            const double dij = haversine_km(cities[i].lat_deg, cities[i].lon_deg,
                                            cities[j].lat_deg, cities[j].lon_deg);
            D[i][j] = dij;
            D[j][i] = dij;
        }
    }
    return D;
}

double tour_length(const std::vector<std::vector<double>>& D, const std::vector<int>& tour) {
    const int n = (int)tour.size();
    double s = 0.0;
    for (int k = 0; k < n; k++) {
        const int a = tour[k];
        const int b = tour[(k + 1) % n];
        s += D[a][b];
    }
    return s;
}

// 2-opt delta for reversing [i..j]
double delta_2opt(const std::vector<std::vector<double>>& D,
                  const std::vector<int>& tour, int i, int j) {
    const int n = (int)tour.size();
    const int a = tour[(i - 1 + n) % n];
    const int b = tour[i];
    const int c = tour[j];
    const int d = tour[(j + 1) % n];
    const double before = D[a][b] + D[c][d];
    const double after  = D[a][c] + D[b][d];
    return after - before;
}

void apply_2opt(std::vector<int>& tour, int i, int j) {
    std::reverse(tour.begin() + i, tour.begin() + j + 1);
}

struct SAConfig {
    // Melt stage (high-T Metropolis)
    double T_melt = 8000.0;
    int melt_steps = 30000;

    // Anneal schedule
    double T0 = 4000.0;
    double Tmin = 1e-2;
    double alpha = 0.995;     // T <- alpha*T
    int steps_per_T = 3000;   // Metropolis steps at each T
};

bool metropolis_accept(double dE, double T, std::mt19937& rng) {
    if (dE <= 0.0) return true;
    if (T <= 0.0) return false;
    std::uniform_real_distribution<double> U(0.0, 1.0);
    return U(rng) < std::exp(-dE / T);
}

std::pair<int,int> propose_2opt(int n, std::mt19937& rng) {
    std::uniform_int_distribution<int> U(0, n - 1);
    int i = U(rng), j = U(rng);
    if (i > j) std::swap(i, j);

    // avoid tiny/no-op
    if (j - i < 2) {
        std::uniform_int_distribution<int> Ui(0, n - 3);
        i = Ui(rng);
        std::uniform_int_distribution<int> Uj(i + 2, n - 1);
        j = Uj(rng);
    }
    // avoid full reversal (same cycle)
    if (i == 0 && j == n - 1) { i = 1; j = n - 2; }
    return {i, j};
}

std::pair<double,double> melt(const std::vector<std::vector<double>>& D,
                              std::vector<int>& tour,
                              const SAConfig& cfg,
                              std::mt19937& rng) {
    double L = tour_length(D, tour);
    int accept = 0;
    const int n = (int)tour.size();

    for (int s = 0; s < cfg.melt_steps; s++) {
        auto [i, j] = propose_2opt(n, rng);
        const double dL = delta_2opt(D, tour, i, j);
        if (metropolis_accept(dL, cfg.T_melt, rng)) {
            apply_2opt(tour, i, j);
            L += dL;
            accept++;
        }
    }
    return {L, (double)accept / std::max(1, cfg.melt_steps)};
}

struct SchedulePoint {
    double T;
    double Lbest;
    double Lcur;
};

void write_schedule_csv(const std::string& out_csv,
                        const std::vector<SchedulePoint>& sched) {
    std::ofstream f(out_csv);
    if (!f) throw std::runtime_error("Failed to write: " + out_csv);
    f << "T,Lbest,Lcur\n";
    f << std::fixed << std::setprecision(12);
    for (const auto& p : sched) {
        f << p.T << "," << p.Lbest << "," << p.Lcur << "\n";
    }
}

std::tuple<std::vector<int>, double, std::vector<SchedulePoint>>
anneal_with_schedule(const std::vector<std::vector<double>>& D,
                     const std::vector<int>& start,
                     const SAConfig& cfg,
                     std::mt19937& rng) {
    const int n = (int)start.size();
    std::vector<int> cur = start;
    double curL = tour_length(D, cur);

    std::vector<int> best = cur;
    double bestL = curL;

    std::vector<SchedulePoint> sched;
    sched.reserve(2048);

    double T = cfg.T0;
    while (T > cfg.Tmin) {
        for (int s = 0; s < cfg.steps_per_T; s++) {
            auto [i, j] = propose_2opt(n, rng);
            const double dL = delta_2opt(D, cur, i, j);
            if (metropolis_accept(dL, T, rng)) {
                apply_2opt(cur, i, j);
                curL += dL;
                if (curL < bestL) { bestL = curL; best = cur; }
            }
        }
        // record point after finishing this temperature block
        sched.push_back(SchedulePoint{T, bestL, curL});
        T *= cfg.alpha;
    }

    return {best, bestL, sched};
}

void write_route_lonlat(const std::string& out,
                        const std::vector<City>& cities,
                        const std::vector<int>& tour) {
    std::ofstream f(out);
    if (!f) throw std::runtime_error("Failed to write: " + out);
    f << std::fixed << std::setprecision(6);

    for (int idx : tour) {
        const auto& c = cities[idx];
        f << c.lon_deg << "\t" << c.lat_deg << "\t" << "\"" << c.name << "\"" << "\n";
    }
    // close the loop
    const auto& c0 = cities[tour[0]];
    f << c0.lon_deg << "\t" << c0.lat_deg << "\t" << "\"" << c0.name << "\"" << "\n";
}

static inline std::string tag_from_filename(const std::string& file) {
    // Simple heuristic for your requested tags
    if (file.find("150") != std::string::npos) return "150";
    if (file.find("1k")  != std::string::npos || file.find("1K") != std::string::npos ||
        file.find("1000")!= std::string::npos) return "1k";
    if (file.find("2k")  != std::string::npos || file.find("2K") != std::string::npos ||
        file.find("2000")!= std::string::npos) return "2k";
    if (file.find("23")  != std::string::npos) return "23";
    return "run";
}

int main(int argc, char** argv) {
    auto t_start = std::chrono::high_resolution_clock::now();

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " cities23.dat [seed]\n";
        return 1;
    }
    const std::string file = argv[1];
    unsigned int seed = 1;
    if (argc >= 3) seed = (unsigned int)std::stoul(argv[2]);

    SAConfig cfg; // tweak defaults here if desired

    try {
        std::mt19937 rng(seed);

        auto cities = read_cities_lonlat_quoted(file);
        const int n = (int)cities.size();
        std::cout << "Loaded " << n << " cities from " << file << "\n";

        auto D = build_distance_matrix(cities);

        // initial random tour (cycle start doesn't matter)
        std::vector<int> tour(n);
        for (int i = 0; i < n; i++) tour[i] = i;
        std::shuffle(tour.begin(), tour.end(), rng);

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Initial tour length: " << tour_length(D, tour) << " km\n";

        auto [Lm, acc] = melt(D, tour, cfg, rng);
        std::cout << "After melt: L=" << Lm << " km, accept~" << acc << "\n";

        auto [best_tour, bestL, sched] = anneal_with_schedule(D, tour, cfg, rng);
        std::cout << "Best found: L=" << bestL << " km\n";

        const std::string tag = tag_from_filename(file);
        write_route_lonlat("best_route_" + tag + ".dat", cities, best_tour);
        write_schedule_csv("schedule_" + tag + ".csv", sched);

        std::cout << "Wrote best_route_" << tag << ".dat\n";
        std::cout << "Wrote schedule_" << tag << ".csv\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }

    auto t_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = t_end - t_start;

    std::cout << "----------------------------------------\n";
    std::cout << "Total execution time: "
              << elapsed.count() << " seconds\n";
    std::cout << "----------------------------------------\n";

    return 0;

    return 0;
}
