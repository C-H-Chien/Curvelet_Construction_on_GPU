#ifndef DATA_IO_HPP
#define DATA_IO_HPP

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

inline std::string edge_input_path(const std::string &filename)
{
    return "../test_files/" + filename;
}

inline bool is_blank_line(const std::string &line)
{
    return line.find_first_not_of(" \t\r\n") == std::string::npos;
}

template<typename T>
bool read_TO_edges_from_file(const std::string &filename, int values_per_edge, std::vector<T> &rd_data)
{
    const std::string path = edge_input_path(filename);
    std::ifstream in_file(path);
    if (!in_file) {
        std::cerr << "Error: input file does not exist: " << path << std::endl;
        return false;
    }

    std::cout << "reading data from file " << filename << std::endl;
    rd_data.clear();

    std::string line;
    int line_num = 0;
    int edge_count = 0;

    while (std::getline(in_file, line)) {
        line_num++;
        if (is_blank_line(line)) {
            std::cerr << "Error: blank line at line " << line_num << " in " << path << " (each line must contain exactly "
                      << values_per_edge << " numeric values)" << std::endl;
            rd_data.clear();
            return false;
        }

        std::istringstream iss(line);
        for (int j = 0; j < values_per_edge; j++) {
            T value;
            if (!(iss >> value)) {
                std::cerr << "Error: invalid edge file format in " << path << " at line " << line_num << ": expected "
                          << values_per_edge << " numeric values per line, found " << j << "." << std::endl;
                rd_data.clear();
                return false;
            }
            rd_data.push_back(value);
        }

        std::string remainder;
        std::getline(iss, remainder);
        if (remainder.find_first_not_of(" \t\r\n") != std::string::npos) {
            std::cerr << "Error: invalid edge file format in " << path << " at line " << line_num << ": expected exactly "
                      << values_per_edge << " numeric values per line, found extra data." << std::endl;
            rd_data.clear();
            return false;
        }

        edge_count++;
    }

    if (edge_count == 0) {
        std::cerr << "Error: input file contains no valid edges: " << path << std::endl;
        rd_data.clear();
        return false;
    }

    std::cout << "read " << edge_count << " edges from " << filename << std::endl;
    return true;
}

inline void write_int_array_to_file(std::string filename, int *wr_data, int first_dim, int second_dim)
{
    std::cout<<"writing data to a file "<<filename<<" ..."<<std::endl;
    const std::filesystem::path out_dir = "../outputs";
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "Failed to create output directory " << out_dir << ": " << ec.message() << std::endl;
        return;
    }

    std::string out_file_name = out_dir.string() + "/";
    out_file_name.append(filename);
	std::ofstream out_file;
    out_file.open(out_file_name);
    if ( !out_file.is_open() )
      std::cout<<"write data file cannot be opened!"<<std::endl;

	for (int i = 0; i < first_dim; i++) {
		for (int j = 0; j < second_dim; j++) {
            out_file << wr_data[i * second_dim + j] <<"\t";
		}
		out_file << "\n";
	}

    out_file.close();
}

inline void write_double_array_to_file(std::string filename, double *wr_data, int first_dim, int second_dim)
{
    std::cout<<"writing data to a file "<<filename<<" ..."<<std::endl;
    const std::filesystem::path out_dir = "../outputs";
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "Failed to create output directory " << out_dir << ": " << ec.message() << std::endl;
        return;
    }
    std::string out_file_name = out_dir.string() + "/";
    out_file_name.append(filename);
    std::ofstream out_file;
    out_file.open(out_file_name);
    if ( !out_file.is_open() )
      std::cout<<"write data file cannot be opened!"<<std::endl;

    for (int i = 0; i < first_dim; i++) {
        for (int j = 0; j < second_dim; j++) {
            out_file << wr_data[i * second_dim + j] <<"\t";
        }
        out_file << "\n";
    }

    out_file.close();
}

inline std::string chain_to_info_filename(const std::string &chain_file)
{
    const std::size_t pos = chain_file.find("chain");
    if (pos != std::string::npos) {
        std::string info_file = chain_file;
        info_file.replace(pos, 5, "info");
        return info_file;
    }
    return "info_" + chain_file;
}

#endif
