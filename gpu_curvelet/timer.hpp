#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//> High-level timing buckets for GPU / host pipeline stages.
enum class TimerCategory {
    MemoryAlloc,
    DataTransfer,
    Kernel,
    Thrust,
    Other,
    Count
};

inline const char *timer_category_label(TimerCategory cat)
{
    switch (cat) {
    case TimerCategory::MemoryAlloc:   return "memory_allocation_and_setup";
    case TimerCategory::DataTransfer:  return "data_transfer";
    case TimerCategory::Kernel:        return "kernel_execution";
    case TimerCategory::Thrust:        return "thrust_process";
    case TimerCategory::Other:         return "other";
    default:                           return "unknown";
    }
}

//> Accumulates lap timings from multiple functions into shared category totals.
class CategoryProfiler {
public:
    void set_title(const char *title)
    {
        _title = (title != nullptr) ? title : "";
    }

    void start()
    {
        _origin = _last = Clock::now();
        _details.clear();
        for (double &sec : _category_totals) {
            sec = 0.0;
        }
    }

    //> Elapsed since the previous lap (or start) is added to cat.
    double lap(TimerCategory cat, const char *detail = nullptr)
    {
        const auto now = Clock::now();
        const double sec = std::chrono::duration<double>(now - _last).count();
        add(cat, detail, sec);
        _last = now;
        return sec;
    }

    //> Add a pre-measured interval (e.g. from a scoped timer or CUDA event).
    void add(TimerCategory cat, const char *detail, double seconds)
    {
        if (seconds < 0.0) {
            return;
        }
        _category_totals[static_cast<size_t>(cat)] += seconds;
        if (detail != nullptr && detail[0] != '\0') {
            _details.push_back({cat, detail, seconds});
            print_step(cat, detail, seconds);
        }
    }

    double elapsed() const
    {
        return std::chrono::duration<double>(Clock::now() - _origin).count();
    }

    double category_total(TimerCategory cat) const
    {
        return _category_totals[static_cast<size_t>(cat)];
    }

    void summary() const
    {
        const double total = elapsed();
        const std::string heading = _title.empty()
            ? "timing by category"
            : (_title + " timing by category");

        std::cout << "\n========== " << heading << " ==========\n";
        for (size_t i = 0; i < static_cast<size_t>(TimerCategory::Count); ++i) {
            const auto cat = static_cast<TimerCategory>(i);
            const double sec = _category_totals[i];
            if (sec <= 0.0) {
                continue;
            }
            const double pct = (total > 0.0) ? (100.0 * sec / total) : 0.0;
            std::cout << "  " << std::left << std::setw(40) << timer_category_label(cat)
                      << std::right << std::setw(12) << format(sec)
                      << "  (" << std::fixed << std::setprecision(1) << pct << "%)\n";
        }
        std::cout << "  " << std::left << std::setw(40) << "total"
                  << std::right << std::setw(12) << format(total) << "\n";
        std::cout << "========================================\n";

        if (!_details.empty()) {
            const std::string detail_heading = _title.empty()
                ? "timing detail"
                : (_title + " timing detail");
            std::cout << "\n========== " << detail_heading << " ==========\n";
            for (const auto &r : _details) {
                const double pct = (total > 0.0) ? (100.0 * r.seconds / total) : 0.0;
                std::cout << "  [" << timer_category_label(r.category) << "] "
                          << std::left << std::setw(34) << r.detail
                          << std::right << std::setw(12) << format(r.seconds)
                          << "  (" << std::fixed << std::setprecision(1) << pct << "%)\n";
            }
            std::cout << "===================================\n";
        }
    }

private:
    using Clock = std::chrono::steady_clock;

    struct DetailRecord {
        TimerCategory category;
        const char *detail;
        double seconds;
    };

    static std::string format(double sec)
    {
        std::ostringstream oss;
        oss << std::fixed;
        if (sec < 1.0) {
            oss << std::setprecision(2) << (sec * 1000.0) << " ms";
        } else {
            oss << std::setprecision(4) << sec << " s";
        }
        return oss.str();
    }

    static void print_step(TimerCategory cat, const char *detail, double sec)
    {
        std::cout << "[timer] [" << timer_category_label(cat) << "] "
                  << detail << ": " << format(sec) << std::endl;
    }

    std::string _title;
    Clock::time_point _origin{};
    Clock::time_point _last{};
    double _category_totals[static_cast<size_t>(TimerCategory::Count)]{};
    std::vector<DetailRecord> _details;
};

//> RAII: on destruction, adds scope duration to the given category.
class ScopedCategoryTimer {
public:
    ScopedCategoryTimer(CategoryProfiler *profiler, TimerCategory cat, const char *detail)
        : _profiler(profiler), _cat(cat), _detail(detail), _start(Clock::now())
    {
    }

    ~ScopedCategoryTimer()
    {
        if (_profiler == nullptr) {
            return;
        }
        const auto now = Clock::now();
        const double sec = std::chrono::duration<double>(now - _start).count();
        _profiler->add(_cat, _detail, sec);
    }

private:
    using Clock = std::chrono::steady_clock;

    CategoryProfiler *_profiler;
    TimerCategory _cat;
    const char *_detail;
    Clock::time_point _start;
};

//> Legacy sequential timer (single function, string labels only). Prefer CategoryProfiler.
class StepTimer {
public:
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
        if (sec < 1.0) {
            oss << std::setprecision(2) << (sec * 1000.0) << " ms";
        } else {
            oss << std::setprecision(4) << sec << " s";
        }
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
