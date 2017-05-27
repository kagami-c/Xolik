#include <fstream>
#include <string>
#include <tclap/CmdLine.h>
#include <MzLoader.h>
#include <PPData.h>
#include <Protocol.h>
#include <Fdr.h>
#include <chrono>

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
    TCLAP::CmdLine cmd("Xolik - A linear-time algorithm for searching cross-linked peptides", ' ', "beta");
    TCLAP::ValueArg<int> histogram_size_arg("", "histogram_size", "Minimum data points required to build histogram for E-value estimation", 
                                            false, 15000, "SIZE", cmd);
    TCLAP::SwitchArg enable_evalue_arg("", "use_evalue", "Estimate E-value for XCorr, reported as -log10(evalue)", cmd);
    TCLAP::ValueArg<int> rank_arg("", "rank", "Rank threshold, used together with --enable_rank", false, 1000, "Rank", cmd);
    TCLAP::SwitchArg enable_rank_arg("", "enable_rank", "Enable rank-based filter on single peptide", cmd);
    DoubleArg thresh_arg("", "threshold", "Score threshold", false, 0.0001, "Xcorr", cmd);
    DoubleArg ms2_tol_arg("", "ms2tol", "MS2 tolerance", false, 0.5, "Da", cmd);
    DoubleArg ms1_tol_arg("", "ms1tol", "MS1 tolerance", false, 50, "PPM", cmd);
    DoubleArg xlmass_arg("", "xlmass", "Cross linker mass", false, 138.0680796, "Da", cmd);
    TCLAP::ValueArg<char> xlsite_arg("", "xlsite", "Cross-linked site", false, 'K', "AA", cmd);
    DoubleArg max_mass_arg("", "max", "Maximum peptide mass", false, 5000, "Da", cmd);
    DoubleArg min_mass_arg("", "min", "Minimum peptide mass", false, 1000, "Da", cmd);
    TCLAP::ValueArg<int> miss_cleavage_arg("", "miss", "Maximum miss cleavage", false, 2, "INT", cmd);
    TCLAP::SwitchArg disable_decoy_arg("", "disable_decoy", "Disable automatically decoy protein generation", cmd);
    TCLAP::SwitchArg disable_linear_arg("", "disable_linear", "Disable linear-time algorithm, fallback to quadratic-time", cmd);
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
    params.use_E_value = enable_evalue_arg.getValue();
    params.histogram_size = histogram_size_arg.getValue();

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
