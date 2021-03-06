/*************************************************************
This is an implentation of greedy based tagSNP selection with fitted imputation modeling
Author: Dat T Nguyen <n.dat@outlook.com>
log:
Create: 12 June 2021
*************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <math.h>
#include <unordered_set>
#include <stdlib.h>

#include "utility/MyFunction.h"
#include "utility/SimpleArgsParser.h"
#include "utility/TimeLog.h"
#include "snp_effect_heap.cpp"
#include "find_snp_model.cpp"


using namespace std;

int main(int argc, char const *args[])
{
    // start time log
    TimeLog time_log;
    string syntax = "\nCommands:\n\tfind\t: find best tag SNP - part of model building step.\n\ttag\t: Run LmTag algorithm to find tag SNPs.\n";
    // command flags
    string find_mode = "find";
    string tag_mode = "tag";

    if(argc < 2){
        cout << syntax << endl;
        exit(1);
    }
    string command = args[1];
    // Argument processing
    SimpleArgsParser Agrs = SimpleArgsParser(argc, args);
    // find best tag SNP mode
    if(command.compare(find_mode) == 0)
    {
        if(argc < 3){
            cout << "\nUsage:\t./LmTag find\n\n\t--ld [ld file generated by Plink v1.9]\n\t--tag [tag SNP list generated by './build_naive_array.R']\n\t-o [output file name].\n" << endl;
            exit(1);
        }
        cout << "Running command:\tfind" << endl;
        // Agrs.find_flag(flag, defaut_value_string, error_mess);
        // PROGRAM WILL STOP IF defaut_value_string = "" and error_mess.length() > 0;
        string ld = Agrs.find_flag("--ld", "", "LD file generated by [Plink v1.9].");
        string tag = Agrs.find_flag("--tag", "", "generated by [./build_naive_array.R].");
        string out_file = Agrs.find_flag("-o", "find_out.txt");
        find_snp_model tem = find_snp_model(tag, ld, out_file);
    }else if (command.compare(tag_mode) == 0){
        if(argc < 3){
            cout << "\nUsage:\t./LmTag tag\n\n\t--ld_model [ld file generated by: 'fit_imputation_model.R']\n\t--eff [snp_effect data file]\n\t--vip [VIP snp data file]\n\t--exclude [excluded SNP data file]\n\t-k [k is number of top SNP considered in each round greedy]\n\t-o [output file name].\n" << endl;
            exit(1);
        } 
        cout << "Running command:\ttag" << endl;
        string ld_model = Agrs.find_flag("--ld_model", "", "LD file with imputation scores [fitted model].");
        string out_file = Agrs.find_flag("-o", "tagSNP_out.txt");
        string snp_effect = Agrs.find_flag("--eff");
        string snp_vip = Agrs.find_flag("--vip");
        string snp_exclude = Agrs.find_flag("--exclude");
        int K =  atoi(Agrs.find_flag("-k", "0").c_str());

        SNPs c_SNPs;
        c_SNPs.ld_scaning(ld_model);
        c_SNPs.read_ld_model(ld_model);
        c_SNPs.read_snp_effect(snp_effect);
        c_SNPs.read_snp_vip(snp_vip);
        // read exclude snps after => will over write VIP SNPs - 
        // noted: it is good to exclude bad SNP to be tag, if VIP snp is critical => can be add manually later
        c_SNPs.read_snp_exclude(snp_exclude);
        c_SNPs.tagging_with_heap(out_file, K);
        
    }else{
        cout << syntax << endl;
        exit(1);
    }
    time_log.print_clock();
    return 0;
};
// ./LmTag find --ld /media/datn/data/DatProjects/vn_array/4_paper/p_VN504/LmTag/chr10/10_ld_0.2.ld --tag /media/datn/data/DatProjects/vn_array/4_paper/p_VN504/LmTag/chr10/chr10_naive.txt -o test_find_tag.txt 
// ./LmTag tag --ld_model /media/datn/data/DatProjects/vn_array/4_paper/p_EUR/LmTag/chr10/ld_fitted_model.txt -o chr10_EUR_LmTag_v1.01.txt -k 1
