#include <fstream>
#include <gtest/gtest.h>

void Check(std::string reference_file, std::string test_file) {
	std::ifstream ref_file(reference_file);
	std::ifstream output_file(test_file);
	while (!output_file.eof() && !ref_file.eof()) {
		std::string output_line;
		output_file >> output_line;
		std::string ref_line;
		ref_file >> ref_line;
		EXPECT_EQ(ref_line, output_line);
	}
	EXPECT_TRUE(output_file.eof());
	EXPECT_TRUE(ref_file.eof());
}

auto database = "test\\data\\synth_plus_100_random.fasta";
auto mzfile = "test\\data\\XL_midK.mzXML";
auto reference = "test\\data\\reference.csv";

std::string prefix = "C:\\Users\\Jiaan\\Desktop\\Workspace\\LimXL\\";
auto cli_path = "x64\\Release\\cli.exe";
auto temp_output = "C:\\Users\\Jiaan\\Desktop\\temp.csv";

TEST(RegressionTest, LimXL) {
	std::string command = prefix + cli_path + " -o " + temp_output
		+ " -d " + prefix + database + " -s " + prefix + mzfile;
	system(command.c_str());
	Check(prefix + reference, temp_output);
}

TEST(RegressionTest, Naive) {
	std::string command = prefix + cli_path + " -o " + temp_output
		+ " -d " + prefix + database + " -s " + prefix + mzfile + " --disable_LimXL";
	system(command.c_str());
	Check(prefix + reference, temp_output);
}

TEST(RegressionTest, RankFilterFallback) {
	std::string command = prefix + cli_path + " -o " + temp_output
		+ " -d " + prefix + database + " -s " + prefix + mzfile + " --enable_rank --rank 100000";
	system(command.c_str());
	Check(prefix + reference, temp_output);
}