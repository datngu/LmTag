/************************************************************
This is an implentation of find best tagSNP for each tagged SNP to build model
Author: Dat T Nguyen <n.dat@outlook.com>
log:
Create: 12 March 2021
*************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <time.h>
#include <unordered_set>
#include <chrono>
#include <stdlib.h>
#include "utility/MyFunction.h"


using namespace std;

class find_snp_model
{
    private:
        struct snp {
            int snp = 0;
            int best_tag = 0;
            float r2 = 0.0;
        };      

        vector<snp> all_snp_info;
        vector<int> all_snp_list;
        vector<int> tagSNPs;
        vector<int> checking_tagSNPs;

        void read_tagSNP( string file_name, vector<int>& tagSNPs, vector<int>& checking_tagSNPs, int & N){
            tagSNPs.clear();
            checking_tagSNPs.resize(N+1,0);
            ifstream fi(file_name);
            string line;
            vector<string> tokens;
            int tagSNP;
            while(fi.good()){
                getline(fi, line);
                char delim = get_delim(line);
                tokenize(line, tokens, delim);
                //cout << tokens[0] << endl;
                tagSNP = atoi(tokens[1].c_str());
                N = (tagSNP > N) ? tagSNP : N;
                tagSNPs.push_back(tagSNP);
            }
            fi.close();
            for( auto x : tagSNPs){
                checking_tagSNPs[x] = 1;
            }
        }

        void read_ld( string file_name, vector<snp>& all_snp_info, vector<int>& all_snp_list, int & N){
            ifstream fi(file_name);
            all_snp_list.clear();
            all_snp_info.resize(N+1);
            vector<int> snp_check;
            snp_check.resize(N+1, 0);
            string line;
            vector<string> tokens;
            getline(fi, line); // skip 1st line
            char delim = get_delim(line);
            int snp_a = 0;
            int snp_b = 0;
            float r2 = 0.0;
            int max_pos = 0;
            while (fi.good())
            {
                getline(fi, line);
                tokenize(line, tokens, delim);
                if(tokens.size() < 7) continue;
                //cout << tokens.size() << endl;
                //cout << tokens[0] << endl;
                snp_a = atoi(tokens[1].c_str());
                snp_b = atoi(tokens[4].c_str());
                snp_check[snp_a] = 1;
                snp_check[snp_b] = 1;
                r2 = atof(tokens[6].c_str()) * atof(tokens[6].c_str());
                all_snp_info[snp_a].snp = snp_a;
                all_snp_info[snp_b].snp = snp_b;
                
                // checking the SNP pair belong to tagSNP list or not, only consider if SNP pair belong to tagSNPs.
                // snp_a checking
                if(checking_tagSNPs[snp_b] == 1){
                    if(all_snp_info[snp_a].r2 < r2){
                        all_snp_info[snp_a].r2 = r2;
                        all_snp_info[snp_a].best_tag = snp_b;
                    }
                }
                // snp_b checking
                if(checking_tagSNPs[snp_a] == 1){
                    if(all_snp_info[snp_b].r2 < r2){
                        all_snp_info[snp_b].r2 = r2;
                        all_snp_info[snp_b].best_tag = snp_a;
                    }
                }
            }
            fi.close();
            for(int i = 1; i < N+1; i++){
                if(snp_check[i] == 1) all_snp_list.push_back(i);
            }
        }  

        int find_max_snp( string file_name){
            ifstream fi(file_name);
            string line;
            vector<string> tokens;
            getline(fi, line); // skip 1st line
            int snp_a = 0;
            int snp_b = 0;
            int max_pos = 0;
            int tem;
            while (fi.good())
            {
                getline(fi, line);
                tokenize(line, tokens, ' ');
                if(tokens.size() < 7) continue;
                snp_a = atoi(tokens[1].c_str());
                snp_b = atoi(tokens[4].c_str());
                tem = snp_a >= snp_b ? snp_a : snp_b;
                max_pos = max_pos >= tem ? max_pos : tem; 
            }
            fi.close();
            return max_pos;
        }

        void write_results(string out_fn, vector<snp>& all_snp_info, vector<int>& all_snp_list){
            ofstream fo;
            fo.open (out_fn);
            fo << "SNP_position" << "\t" << "best_tagSNP_position" << "\t" << "LD_r2" << endl;
            for(auto x : all_snp_list){
                if(all_snp_info[x].best_tag > 0){
                fo << all_snp_info[x].snp << "\t" << all_snp_info[x].best_tag << "\t" << all_snp_info[x].r2 << endl;
                }
            }       

            fo.close();
        }

    public:
    
        find_snp_model(string tagSNP_file, string ld_file, string out_file ){
            int N = find_max_snp(ld_file);
            // cout << N << endl;
            read_tagSNP(tagSNP_file, tagSNPs, checking_tagSNPs, N);
            read_ld( ld_file, all_snp_info, all_snp_list, N);
            // print_result(all_snp_info, all_snp_list);
            write_results(out_file, all_snp_info, all_snp_list );
        }
};
