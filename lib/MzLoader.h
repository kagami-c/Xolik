#pragma once

// MzLoader - an extreme lightweight mzML/mzXML loading library.
// 
// MzLoader aims at providing the essential functionality of loading *MS2* spectrum from most commonly used mzML and mxXML file.
// Comparing with other library, MzLoader merely focuses on loading core data (m/z, intensity, precursor, etc.), 
// which may be quite enough for most of the applications and experiments. API is simple and easy to learn (only one function for loading).
// All information in this header is enough for building your own experiment using this library.
// By using rapidxml as the xml parsing engine, it obtains extreme fast parsing speed. 
// Loading 588M file requires about total 2200ms for loading and parsing, which means that for loading total 20G data will only consume 100s
// (about 200M/s per second). This is already faster than usually used hard disk speed, so that the data loading will not be the bottleneck of 
// the analysis. And it is easy to understand, since the total original codes are about 500 lines.

// TODO: remove some useless fields, like total ion current
// TODO: cover all possible big endian data setting
// TODO: add filter, 1. charge filter, 2. mass filter, 3. ms level filter

#include <vector>
#include <memory>

class MzLoader {
public:
    typedef double Mass;
    typedef double Intensity;

    struct Spectrum {
        unsigned scan_num;
        unsigned ms_level;
        unsigned precursor_charge;
        double precursor_mz;
        double base_peak_mz;
        double base_peak_intensity;
        double total_ion_current;
        std::vector<std::pair<Mass, Intensity>> peaks;
    };

    MzLoader(const char* filename);
    ~MzLoader();

    // return whether next valid spectrum exists. 
    // parameter buffer is for output.
    bool LoadNext(Spectrum& buffer);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};
