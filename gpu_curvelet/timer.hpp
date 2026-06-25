#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class StepTimer {
public:
    void start()
    {
        _origin = _last = Clock::now();
        _records.clear();
    }

    double lap(const std::string& step)
    {
        const auto now = Clock::now();
        const double sec = std::chrono::duration<double>(now - _last).count();
        _records.push_back({step, sec});
        print_step(step, sec);
        _last = now;
        return sec;
    }

    void summary() const
    {
        const double total = std::chrono::duration<double>(Clock::now() - _origin).count();
        std::cout << "\n========== timing summary ==========\n";
        for (const auto& r : _records) {
            const double pct = (total > 0.0) ? (100.0 * r.seconds / total) : 0.0;
            std::cout << "  " << std::left << std::setw(40) << r.name
                      << std::right << std::setw(12) << format(r.seconds)
                      << "  (" << std::fixed << std::setprecision(1) << pct << "%)\n";
        }
        std::cout << "  " << std::left << std::setw(40) << "total"
                  << std::right << std::setw(12) << format(total) << "\n";
        std::cout << "====================================\n";
    }

private:
    using Clock = std::chrono::steady_clock;

    struct Record {
        std::string name;
        double seconds;
    };

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

    static void print_step(const std::string& step, double sec)
    {
        std::cout << "[timer] " << step << ": " << format(sec) << std::endl;
    }

    Clock::time_point _origin{};
    Clock::time_point _last{};
    std::vector<Record> _records;
};

#endif  // TIMER_HPP
