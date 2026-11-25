// Standalone CSV visualization. Matplot++ requires C++17; quant_core stays C++11.
#include <matplot/matplot.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
// Each exported figure has exactly one axes. Avoid Matplot++'s unnecessary
// multiplot wrapper: gnuplot 6.0.6 warns about its inline data in multiplot mode.
class SingleAxesGnuplot : public matplot::backend::gnuplot {
public:
    void run_command(const std::string& command) override
    {
        if (command.find("set multiplot") == 0 || command.find("unset multiplot") == 0) return;
        matplot::backend::gnuplot::run_command(command);
    }
};

struct Row {
    std::string implementation;
    double threads, ms, minimum, maximum, speedup, efficiency;
};
struct Data {
    Row baseline;
    std::vector<Row> parallel;
    std::map<std::string, std::string> metadata;
};

// The benchmark writer emits unquoted, comma-separated scalar fields.
std::vector<std::string> split(std::string line)
{
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.find('"') != std::string::npos) throw std::runtime_error("Expected unquoted benchmark CSV fields");
    std::vector<std::string> fields;
    std::stringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    return fields;
}
double number(const std::string& value)
{
    std::size_t used = 0;
    const double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result)) throw std::runtime_error("Invalid CSV number: " + value);
    return result;
}
void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}
bool near(double a, double b) { return std::abs(a - b) <= 1e-8 * std::max(1.0, std::abs(b)); }
Data read_csv(const std::string& filename)
{
    std::ifstream input(filename);
    require(bool(input), "Cannot open CSV: " + filename);
    std::string line;
    require(bool(std::getline(input, line)), "Empty CSV");
    const auto header = split(line);
    std::map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < header.size(); ++i)
        require(columns.emplace(header[i], i).second, "Duplicate CSV column");
    Data data{};
    bool have_baseline = false;
    while (std::getline(input, line)) {
        if (line.empty() || line == "\r") continue;
        const auto fields = split(line);
        require(fields.size() == header.size(), "CSV row has wrong column count");
        const auto get = [&](const std::string& name) { return fields.at(columns.at(name)); };
        Row r{get("implementation"), number(get("threads")), number(get("median_ms")),
            number(get("min_ms")), number(get("max_ms")), number(get("speedup")), number(get("efficiency"))};
        require(r.threads >= 1 && std::floor(r.threads) == r.threads && r.minimum > 0 &&
                r.minimum <= r.ms && r.ms <= r.maximum && r.speedup > 0 && r.efficiency > 0,
                "Invalid runtime, thread count or ratio");
        for (const std::string name : {"paths", "steps", "repeats", "warmups", "s0", "sigma", "rate",
                                       "maturity", "strike", "base_seed", "build"}) {
            const auto inserted = data.metadata.emplace(name, get(name));
            require(inserted.second || inserted.first->second == get(name), "Mixed workloads in CSV: " + name);
        }
        if (r.implementation == "original_single") {
            require(!have_baseline && r.threads == 1, "Expected one serial baseline");
            data.baseline = r;
            have_baseline = true;
        } else {
            require(r.implementation == "parallel", "Unknown implementation");
            data.parallel.push_back(r);
        }
    }
    require(have_baseline && data.parallel.size() == 4, "Expected baseline and four parallel configurations");
    require(data.metadata.at("build") == "Release", "Expected Release benchmark results");
    for (const std::string name : {"paths", "steps", "repeats", "warmups"}) {
        const double n = number(data.metadata.at(name));
        require(n >= 1 && std::floor(n) == n, "Invalid workload count: " + name);
    }
    std::sort(data.parallel.begin(), data.parallel.end(), [](const Row& a, const Row& b) { return a.threads < b.threads; });
    const std::array<double, 4> threads = {{1, 2, 4, 8}};
    require(near(data.baseline.speedup, 1) && near(data.baseline.efficiency, 1), "Invalid baseline ratios");
    for (std::size_t i = 0; i < threads.size(); ++i) {
        const auto& r = data.parallel[i];
        require(r.threads == threads[i], "Expected unique thread counts 1, 2, 4, 8");
        require(near(r.speedup, data.baseline.ms / r.ms) && near(r.efficiency, r.speedup / r.threads),
                "CSV speedup/efficiency disagrees with serial baseline");
    }
    return data;
}

const std::array<float, 4> blue = {{0.0f, 0.0f, .447f, .698f}};
const std::array<float, 4> orange = {{0.0f, .835f, .369f, .0f}};
const std::array<float, 4> gray = {{0.0f, .45f, .48f, .52f}};
std::string fixed(double value, int decimals)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(decimals) << value;
    return out.str();
}
void draw(const Data& data, int metric, const std::filesystem::path& directory)
{
    using namespace matplot;
    const char* titles[] = {"Execution Time vs Threads", "Speedup vs Threads", "Parallel Efficiency vs Threads"};
    const char* ylabels[] = {"Median execution time (s)", "Speedup relative to original serial", "Parallel efficiency (%)"};
    const char* filenames[] = {"execution_time_vs_threads.png", "speedup_vs_threads.png", "parallel_efficiency_vs_threads.png"};
    std::vector<double> x, y;
    for (const auto& r : data.parallel) {
        x.push_back(r.threads);
        y.push_back(metric == 0 ? r.ms / 1000 : metric == 1 ? r.speedup : r.efficiency * 100);
    }
    auto fig = figure<SingleAxesGnuplot>(true); // Export only; no interactive window.
    fig->size(1600, 1000);
    fig->color("white");
    fig->font("Helvetica");
    fig->font_size(20);
    auto ax = fig->current_axes();
    ax->position({.14f, .16f, .81f, .66f});
    ax->font_size(20);
    ax->title_font_size_multiplier(1.4f);
    ax->title_font_weight("normal");
    const std::string caption = data.metadata.at("paths") + " paths | " + data.metadata.at("steps") +
        " steps | " + data.metadata.at("repeats") + " timed runs | " + data.metadata.at("warmups") + " warmup per configuration";
    ax->title("{/:Bold " + std::string(titles[metric]) + "}\\n{/Helvetica=16 " + caption + "}");
    ax->xlabel("Worker threads");
    ax->ylabel(ylabels[metric]);
    ax->x_axis().label_font_size(22);
    ax->y_axis().label_font_size(22);
    ax->xlim({.5, 8.5});
    ax->xticks(x);
    ax->xticklabels({"1", "2", "4", "8"});
    ax->y_grid(true);
    ax->x_grid(false);
    ax->grid_color({0, .88f, .90f, .92f});
    ax->grid_alpha(1.0f);
    ax->box(false);
    ax->hold(true);
    const double upper = metric == 0 ? std::max(data.baseline.ms / 1000, *std::max_element(y.begin(), y.end())) * 1.35
                         : metric == 1 ? std::max(8.0, *std::max_element(y.begin(), y.end())) * 1.15
                         : std::max(100.0, *std::max_element(y.begin(), y.end())) * 1.15;
    ax->ylim({0, upper});
    if (metric == 1) ax->yticks({0, 1, 2, 3, 4, 5, 6, 7, 8});
    if (metric == 2) ax->yticks({0, 20, 40, 60, 80, 100});
    auto measured = ax->plot(x, y, "-o");
    measured->color(blue).line_width(3).marker_size(10).marker_face_color(blue);
    std::string reference;
    if (metric == 0) {
        ax->plot(std::vector<double>{.5, 8.5}, std::vector<double>(2, data.baseline.ms / 1000), "--")
            ->color(orange).line_width(2);
        reference = "Original serial baseline (" + fixed(data.baseline.ms / 1000, 2) + " s)";
    } else if (metric == 1) {
        ax->plot(x, x, "--")->color(gray).line_width(2);
        reference = "Ideal linear speedup";
    } else {
        ax->plot(std::vector<double>{.5, 8.5}, std::vector<double>{100, 100}, "--")->color(gray).line_width(2);
        reference = "Ideal efficiency (100%)";
    }
    // Empty entries keep the four numeric annotations out of the legend.
    ax->legend(std::vector<std::string>{"Parallel engine (median)", reference, "", "", "", ""});
    auto key = ax->legend();
    key->location(metric == 0 ? legend::general_alignment::topright :
                  metric == 1 ? legend::general_alignment::topleft : legend::general_alignment::bottomleft);
    key->box(false);
    key->font_size(19);
    key->font_weight("normal");
    for (std::size_t i = 0; i < x.size(); ++i) {
        // Keep labels below runtime/efficiency points, above speedup points.
        const double offset = (metric == 1 ? .04 : -.055) * upper;
        auto label = ax->text(x[i], y[i] + offset, fixed(y[i], 2) + (metric == 0 ? " s" : metric == 1 ? "x" : "%"));
        label->font_size(18).alignment(labels::alignment::center).color(blue);
    }
    const auto output = directory / filenames[metric];
    require(fig->save(output.string()), "PNG export failed: " + output.string());
    std::cout << "Saved " << output << '\n';
}
}

int main(int argc, char** argv)
{
    try {
        if (argc != 3) {
            std::cerr << "Usage: plot_benchmark SUMMARY.csv OUTPUT_DIRECTORY\n";
            return 1;
        }
        const auto data = read_csv(argv[1]);
        // Use a portable decimal locale and a noninteractive terminal for the
        // child gnuplot processes. macOS gnuplot can crash on unavailable locales
        // such as C.UTF-8. These settings affect only this process and its children.
#ifdef _WIN32
        require(_putenv_s("LC_ALL", "C") == 0 && _putenv_s("GNUTERM", "pngcairo") == 0,
                "Cannot configure gnuplot environment");
#else
        require(setenv("LC_ALL", "C", 1) == 0 && setenv("GNUTERM", "pngcairo", 1) == 0,
                "Cannot configure gnuplot environment");
#endif
        std::filesystem::create_directories(argv[2]);
        for (int metric = 0; metric < 3; ++metric) draw(data, metric, argv[2]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Plot failed: " << error.what() << '\n';
        return 1;
    }
}
