#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

inline std::string csv_escape(const std::string &s)
{
    bool need_quotes = false;
    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            need_quotes = true;
            break;
        }
    }
    if (!need_quotes) {
        return s;
    }
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') {
            out.push_back('"');
        }
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

inline std::string csv_format_seconds(double sec)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(9) << sec;
    return oss.str();
}

inline bool file_is_nonempty(const std::string &path)
{
    struct stat st {};
    if (stat(path.c_str(), &st) != 0) {
        return false;
    }
    return st.st_size > 0;
}

class StepTimer {
public:
    struct Record {
        std::string name;
        double seconds;
    };

    void start()
    {
        _origin = _last = Clock::now();
        _records.clear();
    }

    double lap(const std::string &step)
    {
        const auto now = Clock::now();
        const double sec = std::chrono::duration<double>(now - _last).count();
        _records.push_back({step, sec});
        print_step(step, sec);
        _last = now;
        return sec;
    }

    //> Add a pre-measured interval without advancing the lap clock.
    void add(const std::string &step, double seconds)
    {
        if (seconds < 0.0) {
            return;
        }
        _records.push_back({step, seconds});
    }

    double elapsed() const
    {
        return std::chrono::duration<double>(Clock::now() - _origin).count();
    }

    double seconds_named(const std::string &step) const
    {
        double sum = 0.0;
        for (const auto &r : _records) {
            if (r.name == step) {
                sum += r.seconds;
            }
        }
        return sum;
    }

    const std::vector<Record> &records() const { return _records; }

    void summary() const
    {
        const double total = std::chrono::duration<double>(Clock::now() - _origin).count();
        std::cout << "\n========== timing summary ==========\n";
        for (const auto &r : _records) {
            const double pct = (total > 0.0) ? (100.0 * r.seconds / total) : 0.0;
            std::cout << "  " << std::left << std::setw(40) << r.name
                      << std::right << std::setw(12) << format(r.seconds)
                      << "  (" << std::fixed << std::setprecision(1) << pct << "%)\n";
        }
        std::cout << "  " << std::left << std::setw(40) << "total"
                  << std::right << std::setw(12) << format(total) << "\n";
        std::cout << "====================================\n";
    }

    //> Append one row per lap to a long-form CSV (creates header if file is new/empty).
    bool append_detail_csv(const std::string &path, const std::string &run_id,
                           const std::string &stage = "cpu_curvelet") const
    {
        if (path.empty()) {
            return true;
        }
        const bool write_header = !file_is_nonempty(path);
        std::ofstream out(path, std::ios::app);
        if (!out) {
            std::cerr << "Error: could not open timing detail CSV: " << path << std::endl;
            return false;
        }
        if (write_header) {
            out << "run_id,stage,detail,seconds\n";
        }
        for (const auto &r : _records) {
            out << csv_escape(run_id) << ','
                << csv_escape(stage) << ','
                << csv_escape(r.name) << ','
                << csv_format_seconds(r.seconds) << '\n';
        }
        return static_cast<bool>(out);
    }

private:
    using Clock = std::chrono::steady_clock;

    static std::string format(double sec)
    {
        std::ostringstream oss;
        oss << std::fixed;
        if (sec < 1.0)
            oss << std::setprecision(2) << (sec * 1000.0) << " ms";
        else
            oss << std::setprecision(4) << sec << " s";
        return oss.str();
    }

    static void print_step(const std::string &step, double sec)
    {
        std::cout << "[timer] " << step << ": " << format(sec) << std::endl;
    }

    Clock::time_point _origin{};
    Clock::time_point _last{};
    std::vector<Record> _records;
};

#endif  // TIMER_HPP
