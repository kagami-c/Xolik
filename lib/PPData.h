// Copyright (C) 2016

#pragma once

#include <memory>

class PPData {
public:
    struct Protein {
        const char* const name;
        const char* const sequence;
        const size_t sequence_length;

        Protein(const char* name, const char* sequence, size_t sequence_length);
    };

    struct Peptide {
        const char* sequence;
        size_t sequence_length;

        char n_term;
        char c_term;
        double mass;  // nonconst for adjusting modified amino acids

        const Protein* protein;
        size_t offset;  // offset in protein sequence

        Peptide(const Protein& protein, const char* compact_protein_sequence,
                size_t start_idx, size_t end_idx, double mass);
    };

    enum class EnzymeType { Trypsin };

    // ctors
    PPData(const char* filename, bool append_decoy, EnzymeType enzyme_type,
           unsigned max_miss_cleavage, double min_mass, double max_mass);
    PPData(const char* filename)
           : PPData(filename, false, EnzymeType::Trypsin, 0, 600.0, 5000.0) {}
    ~PPData();

    // access methods
    size_t size() const;
    const Peptide& operator[](const size_t index) const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};
