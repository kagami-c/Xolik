#pragma once

#include <unordered_map>

// simply a struct storing all we need for a run of a search
struct Params {
    // data input and output, required
    std::string database_path = "";
    std::string mzfile_path = "";
    std::string output_path = "";

    // database digest
    std::string enzyme = "trypsin";
    unsigned max_miss_cleavage = 2;
    double min_allowed_mass = 1000;  // Da
    double max_allowed_mass = 5000;  // Da
    bool append_decoy = true;

    // mass match
    double ms1_tolerance = 50;  // ppm
    double ms2_tolerance = 0.5;  // Da

    // cross link
    char xlsite = 'K';
    double xlmass = 138.0680796;

    // search options
    bool use_LimXL_match = true;
    double threshold = 0.0001;  // Xcorr score threshold
    bool enable_rank = false;
    int rank = 1000;
    
    // evalue
    bool use_E_value = true;
    int histogram_size = 15000;

    // parallel
    bool enable_parallel = false;
    int thread = 4;

    // modifications
    std::unordered_map<char, double> fix_mods = { { 'C', +57.021464 } };
    std::unordered_map<char, double> var_mods = {};
};

const std::unordered_map<char, double> MassTable = {
    { 'G', 57.02147 },{ 'A', 71.03712 },{ 'S', 87.03203 },{ 'P', 97.05277 },
    { 'V', 99.06842 },{ 'T', 101.04768 },//{ 'C', 103.00919 + 57.021464 /* Fixed Mod on C */ },
    { 'C', 103.00919 },
    /*{ 'I', 113.08407 },*/{ 'L', 113.08407 },
    { 'N', 114.04293 },{ 'D', 115.02695 },{ 'Q', 128.05858 },
    { 'K', 128.09497 },{ 'E', 129.04260 },{ 'M', 131.04049 },{ 'H', 137.05891 },
    { 'F', 147.06842 },{ 'R', 156.10112 },{ 'Y', 163.06333 },{ 'W', 186.07932 }
};

// internal mass table
const double PROTON_MASS = 1.00727;
const double HYDROGEN_MASS = 1.00782;
const double OXYGEN_MASS = 15.99491;
const double WATER_MASS = OXYGEN_MASS + HYDROGEN_MASS + HYDROGEN_MASS;

// data interface, as the output of Search() routine
struct Record {
    unsigned spec_idx;
    double score;  // score = alpha_score + beta_score

    size_t alpha_idx;
    size_t alpha_link_site;
    size_t beta_idx;
    size_t beta_link_site;
    double alpha_score;
    double beta_score;
};
