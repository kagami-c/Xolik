#include <fstream>
#include <string>
#include <chrono>
#include <tclap/CmdLine.h>
#include <MzLoader.h>
#include <PPData.h>
#include "Protocol.h"
#include "Fdr.h"

size_t FindProteinNameEnd(const char* name) {
    size_t end = 0;
    while (name[end] != '\0' && !isspace(name[end])) {
        ++end;
    }
    return end;
}

void WriteResults(const char* output_path, const std::vector<Record>& records,
                  const PPData& ppdata, const Params& params,
                  const std::vector<double>& q_values) {
    std::ofstream file(output_path);
    file << "ScanNum,Score,Peptide#1,LinkSite#1,Protein#1,Score#1,"
        "Peptide#2,LinkSite#2,Protein#2,Score#2,qValue\n";
    for (auto i = 0; i < records.size(); ++i) {
        const auto& record = records[i];
        const auto& alpha_peptide = ppdata[record.alpha_idx];
        const auto& beta_peptide = ppdata[record.beta_idx];
        const auto& alpha_protein = *alpha_peptide.protein;
        const auto& beta_protein = *beta_peptide.protein;
        file << record.spec_idx << ',' << record.score << ',' << alpha_peptide.n_term << '.'
             << std::string(alpha_peptide.sequence, alpha_peptide.sequence_length) << '.'
             << alpha_peptide.c_term << ',' << record.alpha_link_site << ','
             << std::string(alpha_protein.name, FindProteinNameEnd(alpha_protein.name)) << ','
             << record.alpha_score << ',' << beta_peptide.n_term << '.'
             << std::string(beta_peptide.sequence, beta_peptide.sequence_length) << '.'
             << beta_peptide.c_term << ',' << record.beta_link_site << ','
             << std::string(beta_protein.name, FindProteinNameEnd(beta_protein.name)) << ','
             << record.beta_score << ',' << q_values[i] << '\n';
    }
}

Params ParseArguments(int argc, const char** argv) {
    using StringArg = TCLAP::ValueArg<std::string>;
    using DoubleArg = TCLAP::ValueArg<double>;
    using IntArg = TCLAP::ValueArg<int>;
    TCLAP::CmdLine cmd("Xolik - A linear-time algorithm for searching cross-linked peptides", ' ', "beta");
    // TODO: add enzyme setting
//    std::vector<std::string> allowed_enzymes = { "trypsin" };
//    TCLAP::ValuesConstraint<std::string> enzyme_constraint(allowed_enzymes);
//    StringArg enzyme_arg("", "enzyme", "Enzyme for in silico digestion.", false, "trypsin", &enzyme_constraint, cmd);
    StringArg varmod_arg("", "varmod", "Variable modifications (Example: \"M+15.9949:S+79.96633\", use + or - to "
                         "connect AA with MASS, use : to separate multiple mods)", false, "", "PATTERN", cmd);
    StringArg fixmod_arg("", "fixmod", "Fixed modifications (Example: \"C+57.021464\", use + or - to connect AA with "
                         "MASS, use : to separate multiple mods)", false, "C+57.021464", "PATTERN", cmd);
    IntArg thread_arg("", "thread", "Number of threads for parallel computing, used together with --parallel", 
                                    false, 4, "INT", cmd);
    TCLAP::SwitchArg parallel_arg("", "parallel", "Enable parallel computing by multi-threading", cmd);
    IntArg histogram_size_arg("", "histogram_size",
                                            "Minimum data points required to build histogram for E-value estimation", 
                                            false, 15000, "SIZE", cmd);
    TCLAP::SwitchArg noevalue_arg("", "noevalue", 
                                  "Disable E-value estimation, E-value will be reported as -log10(evalue)", cmd);
    IntArg rank_arg("", "rank", "Rank threshold, used together with --enable_rank",
                                  false, 1000, "Rank", cmd);
    TCLAP::SwitchArg enable_rank_arg("", "enable_rank", "Enable rank-based filter on single peptide", cmd);
    DoubleArg thresh_arg("", "threshold", "Score threshold", false, 0.0001, "Xcorr", cmd);
    DoubleArg ms2_tol_arg("", "ms2tol", "MS2 tolerance", false, 0.5, "Da", cmd);
    DoubleArg ms1_tol_arg("", "ms1tol", "MS1 tolerance", false, 50, "PPM", cmd);
    DoubleArg xlmass_arg("", "xlmass", "Cross linker mass", false, 138.0680796, "Da", cmd);
    TCLAP::ValueArg<char> xlsite_arg("", "xlsite", "Cross-linked site", false, 'K', "AA", cmd);
    DoubleArg max_mass_arg("", "max", "Maximum peptide mass", false, 5000, "Da", cmd);
    DoubleArg min_mass_arg("", "min", "Minimum peptide mass", false, 1000, "Da", cmd);
    IntArg miss_cleavage_arg("", "miss", "Maximum miss cleavage", false, 2, "INT", cmd);
    TCLAP::SwitchArg disable_decoy_arg("", "disable_decoy", "Disable automatically decoy protein generation", cmd);
    TCLAP::SwitchArg disable_linear_arg("", "disable_linear", 
                                        "Disable linear-time algorithm, fallback to quadratic-time", cmd);
    StringArg output_path_arg("o", "output", "Output path", true, "", "FILE", cmd);
    StringArg mzfile_path_arg("s", "spectrum", "Mass spectrum file", true, "", "mzXML/mzML", cmd);
    StringArg database_path_arg("d", "database", "Database file", true, "", "FASTA", cmd);
    cmd.parse(argc, argv);

    Params params;
    params.database_path = database_path_arg.getValue();
    params.mzfile_path = mzfile_path_arg.getValue();
    params.output_path = output_path_arg.getValue();
    params.use_LimXL_match = disable_linear_arg.getValue() ? false : true;
    params.append_decoy = disable_decoy_arg.getValue() ? false : true;
    params.max_miss_cleavage = miss_cleavage_arg.getValue();
    params.min_allowed_mass = min_mass_arg.getValue();
    params.max_allowed_mass = max_mass_arg.getValue();
    params.xlsite = xlsite_arg.getValue();
    params.xlmass = xlmass_arg.getValue();
    params.ms1_tolerance = ms1_tol_arg.getValue();
    params.ms2_tolerance = ms2_tol_arg.getValue();
    params.threshold = thresh_arg.getValue();
    params.enable_rank = enable_rank_arg.getValue();
    params.rank = rank_arg.getValue();
    params.use_E_value = noevalue_arg.getValue() ? false : true;
    params.histogram_size = histogram_size_arg.getValue();
    params.enable_parallel = parallel_arg.getValue();
    params.thread = thread_arg.getValue();

    auto mod_pattern_parse = [](const std::string& arg, std::unordered_map<char, double>& map_out) {
        std::stringstream argstream(arg);
        char c;
        while (argstream >> c) {
            if (c == ':') {
                continue;
            } else if ('A' <= c && c <= 'Z' && c != 'X' && c != 'B' && c != 'Z' && c != 'J') {
                double mass_shift;
                if (argstream >> mass_shift) {
                    map_out[c] = mass_shift;
                } else {
                    std::cerr << "Parse mod pattern error: no mass shift value specified.\n";
                    exit(EXIT_FAILURE);
                }
            } else {
                std::cerr << "Parse mod pattern error: unknown amino acid " << c << ".\n";
                exit(EXIT_FAILURE);
            }
        }
    };
    mod_pattern_parse(fixmod_arg.getValue(), params.fix_mods);
    mod_pattern_parse(varmod_arg.getValue(), params.var_mods);
    
    return params;
}

void PrintSettings(const Params& params) {
    printf("Xolik ver.beta (%s, %s)\n", __DATE__, __TIME__);
    printf("Database:              %s\n", params.database_path.c_str());
    printf("Spectra:               %s\n", params.mzfile_path.c_str());
    printf("Output path:           %s\n", params.output_path.c_str());
    printf("Use linear-time mode:  %s\n", params.use_LimXL_match ? "true" : "false");
    printf("Append decoy:          %s\n", params.append_decoy ? "true" : "false");
    printf("Maximum miss cleavage: %d\n", params.max_miss_cleavage);
    printf("Minimum peptide mass:  %.2f Da\n", params.min_allowed_mass);
    printf("Maximum peptide mass:  %.2f Da\n", params.max_allowed_mass);
    printf("Cross-linked site:     %c\n", params.xlsite);
    printf("Cross linker mass:     %f Da\n", params.xlmass);
    printf("MS1 tolerance:         %.2f ppm\n", params.ms1_tolerance);
    printf("MS2 tolerance:         %.4f Da\n", params.ms2_tolerance);
    printf("Score threshold:       %.4f\n", params.threshold);
    printf("Enable rank filter:    %s\n", params.enable_rank ? "true" : "false");
    printf("Rank threshold:        %d\n", params.rank);
    printf("Use E-value:           %s\n", params.use_E_value ? "true" : "false");
    printf("Histogram size:        %d\n", params.histogram_size);
    printf("Parallel computing:    %s\n", params.enable_parallel ? "true" : "false");
    printf("Number of threads:     %d\n", params.thread);
    printf("Enzyme:                %s\n", params.enzyme.c_str());

    auto print_map = [](const std::unordered_map<char, double>& map) {
        if (map.empty()) {
            printf("None\n");
        } else {
            bool first = true;
            for (const auto& kv : map) {
                if (first) {
                    first = false;
                } else {
                    printf("                       ");
                }
                printf("%c %+.6f Da\n", kv.first, kv.second);
            }
        }
    };
    printf("Fixed mods:            ");
    print_map(params.fix_mods);
    printf("Variable mods:         ");
    print_map(params.var_mods);
}

int main(int argc, const char** argv) {
    Params params = ParseArguments(argc, argv);
    PrintSettings(params);

#pragma warning(disable: 4996)  // _CRT_SECURE_NO_WARNINGS
    auto start = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(start);
    printf(">> Begin at %s", std::ctime(&start_time));

    printf(">> Preparing data...\n");
    MzLoader loader(params.mzfile_path.c_str());
    PPData ppdata(params.database_path.c_str(), params.append_decoy,
                  PPData::EnzymeType::Trypsin, params.max_miss_cleavage,
                  params.min_allowed_mass, params.max_allowed_mass);

    printf(">> Searching...\n");
    auto results = Search(loader, ppdata, params);
    auto q_values = CalculateFDR(results, ppdata);
    WriteResults(params.output_path.c_str(), results, ppdata, params, q_values);

    auto end = std::chrono::system_clock::now();
    auto end_time = std::chrono::system_clock::to_time_t(end);
    printf(">> End at %s", std::ctime(&end_time));
    std::chrono::duration<double> elapsed = end - start;
    std::cout << ">> Elapsed: " << elapsed.count() << "s\n";

    return 0;
}
