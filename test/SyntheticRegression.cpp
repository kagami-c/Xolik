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

TEST(SyntheticRegression, DISABLED_Xolik_lib1) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib1_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5";
	system(command.c_str());
	Check("..\\test\\data\\lib1.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Xolik_lib1_noevalue) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib1_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue";
    system(command.c_str());
    Check("..\\test\\data\\lib1_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, DISABLED_Xolik_lib2) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib2_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5";
    system(command.c_str());
    Check("..\\test\\data\\lib2.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Xolik_lib2_noevalue) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib2_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue";
    system(command.c_str());
    Check("..\\test\\data\\lib2_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, DISABLED_Xolik_lib3) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib3_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5";
    system(command.c_str());
    Check("..\\test\\data\\lib3.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Xolik_lib3_noevalue) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib3_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue";
    system(command.c_str());
    Check("..\\test\\data\\lib3_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, DISABLED_Xolik_lib1_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib1_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\lib1.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Xolik_lib1_noevalue_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib1_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\lib1_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, DISABLED_Xolik_lib2_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib2_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\lib2.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Xolik_lib2_noevalue_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib2_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\lib2_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, DISABLED_Xolik_lib3_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib3_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\lib3.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Xolik_lib3_noevalue_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib3_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\lib3_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Naive_lib1_noevalue) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib1_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --disable_linear";
    system(command.c_str());
    Check("..\\test\\data\\lib1_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Naive_lib2_noevalue) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib2_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --disable_linear";
    system(command.c_str());
    Check("..\\test\\data\\lib2_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Naive_lib3_noevalue) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib3_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --disable_linear";
    system(command.c_str());
    Check("..\\test\\data\\lib3_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Naive_lib1_noevalue_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib1_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --parallel --thread 4 --disable_linear";
    system(command.c_str());
    Check("..\\test\\data\\lib1_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Naive_lib2_noevalue_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib2_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --parallel --thread 4 --disable_linear";
    system(command.c_str());
    Check("..\\test\\data\\lib2_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(SyntheticRegression, Naive_lib3_noevalue_parallel) {
    const std::string command = "Release\\Xolik.exe "
        "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
        "-s C:\\Users\\Jiaan\\Desktop\\Experiments\\Exp2\\mxdb\\20111221_ananiav_DPDS_lib3_90min_CID35.mzXML "
        "-d C:\\Users\\Jiaan\\Desktop\\Experiments\\LimXL\\test\\data\\lib_disulfide_peptides.fasta "
        "--min 500 --xlsite C --xlmass -116.0430 --ms1tol 50 --ms2tol 0.5 --noevalue --parallel --thread 4 --disable_linear";
    system(command.c_str());
    Check("..\\test\\data\\lib3_noevalue.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

