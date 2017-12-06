#include <fstream>
#include <gtest/gtest.h>

void Check(std::string reference_file, std::string test_file);

TEST(RibosomeRegression, DESTRibo_092010exp2_LE_1and2rejected) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LE_1and2rejected.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LE_1and2rejected.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF10) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF10.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF10.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF1_101012172356) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF1_101012172356.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF1_101012172356.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DISABLED_DESTRibo_092010exp2_LEF2) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF2.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF2.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF3) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF3.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF3.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF4) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF4.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF4.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF5) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF5.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF5.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF6) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF6.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF6.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF7) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF7.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF7.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_LEF8) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF8.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF8.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DISABLED_DESTRibo_092010exp2_LEF9) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_LEF9.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_LEF9.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092010exp2_WD_1and2rejected) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092010exp2_WD_1and2rejected.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092010exp2_WD_1and2rejected.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LE_1and2rejected) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LE_1and2rejected.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LE_1and2rejected.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF10) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF10.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF10.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF1) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF1.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF1.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF2) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF2.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF2.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF3) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF3.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF3.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF4) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF4.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF4.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF5) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF5.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF5.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF6) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF6.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF6.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF7) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF7.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF7.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF8) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF8.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF8.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_LEF9) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_LEF9.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_LEF9.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_092210exp2_WD_1and2rejected) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_092210exp2_WD_1and2rejected.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_092210exp2_WD_1and2rejected.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF10wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF10wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF10wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF10wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF10wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF10wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF1wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF1wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF1wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF1wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF1wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF1wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF2wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF2wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF2wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF2wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF2wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF2wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF3wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF3wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF3wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF3wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF3wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF3wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF4wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF4wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF4wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF4wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF4wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF4wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF5wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF5wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF5wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF5wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF5wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF5wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF6wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF6wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF6wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF6wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF6wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF6wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF7wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF7wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF7wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF7wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF7wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF7wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF8wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF8wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF8wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF8wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF8wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF8wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF9wo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF9wo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF9wo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEF9wo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEF9wo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEF9wo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEwo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEwo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEwo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_LEwo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_LEwo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_LEwo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_noLEFwo2__071510) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_noLEFwo2+_071510.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_noLEFwo2+_071510.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

TEST(RibosomeRegression, DESTRibo_noLEFwo2__072710) {
    const std::string command = "Release\\Xolik.exe "
            "-d ..\\test\\data\\ribosome_ref\\ecoli30s50s.fasta "
            "-s D:\\Spectra\\PXD003381\\DESTRibo_noLEFwo2+_072710.mzXML "
            "-o C:\\Users\\Jiaan\\Desktop\\temp.csv "
            "--ms1tol 5 --ms2tol 0.02 --xlmass 136.10005 --min 500 --noevalue --parallel --thread 4";
    system(command.c_str());
    Check("..\\test\\data\\ribosome_ref\\DESTRibo_noLEFwo2+_072710.mzXML.output.csv", "C:\\Users\\Jiaan\\Desktop\\temp.csv");
}

