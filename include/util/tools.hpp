// Copyright 2025-2026 Joseph Le
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "util/typedefs.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <string>

inline std::string trim(const std::string& s) {
    // trims leading and trailing spaces
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";

    size_t end = s.find_last_not_of(' ');
    return s.substr(start, end - start + 1);
}

inline bool is_numeric(const std::string& s, bool strict = false) {
    std::string t = trim(s);
    string list = strict ? "0123456789" : "+-0123456789.eE";

    return !t.empty() && t.find_first_not_of(list) == std::string::npos;
}

inline i32 last_numeric(const string& s, bool strict = false) {
    string t = trim(s);
    string list = strict ? "0123456789" : "+-0123456789.eE";
    auto idx = t.find_last_of(list);
    if (idx == string::npos) {
        return -1;
    }

    return static_cast<i32>(idx);
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

inline std::string make_lower(const std::string& s) {
    std::string str;

    for (char c : s) {
        str += std::tolower(c);
    }

    return str;
}

inline std::string make_upper(const std::string& s) {
    std::string str;

    for (char c : s) {
        str += std::toupper(c);
    }

    return str;
}

inline string space_to_underscore(const string str) {
    string out = str;
    for (char& c : out) {
        if (c == ' ') {
            c = '_';
        }
    }
    return out;
}

template <class T>
inline bool in_list(const T& x, std::initializer_list<T> options) {
    return std::find(options.begin(), options.end(), x) != options.end();
}

template <class T>
inline i32 index_of(const T& x, std::initializer_list<T> options) {
    auto it = std::find(options.begin(), options.end(), x);

    if (it == options.end()) {
        return -1;
    }

    return static_cast<i32>(std::distance(options.begin(), it));
}

inline string repeat_char(char c, i32 n) { return string(n, c); }
inline string repeat_string(const string& s, i32 n) {
    if (n <= 0) return "";

    string result;
    result.reserve(s.size() * n);

    for (i32 i = 0; i < n; i++) {
        result += s;
    }

    return result;
}
inline string repeat_char(const string& s, i32 n) { return repeat_string(s, n); }

template <class T>
inline string decimal_digits(T value, i32 digits) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(digits) << value;

    std::string s = out.str();
    size_t dot = s.find('.');

    return (dot == std::string::npos) ? "" : s.substr(dot + 1);
}

template <class T>
inline string to_string_precision(T value, i32 precision = 4) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;

    return out.str();
}
