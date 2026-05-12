#pragma once

#include <string>

inline std::string trim(const std::string& s) {
    // trims leading and trailing spaces
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";

    size_t end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

inline bool is_numeric(const std::string& s) {
    std::string t = trim(s);

    return !t.empty() && t.find_first_not_of("+-0123456789.eE") == std::string::npos;
}

inline std::string remove_returns(const std::string& s) {
    std::string str;

    for (char c : s) {
        if (c != '\r' && c != '\n') {
            str += c;
        }
    }

    return str;
}
