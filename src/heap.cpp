#ifndef MY_HEAP
#define MY_HEAP
/************************************************************
Author: Dat T Nguyen <n.dat@outlook.com>
log:
Create: 25 May 2021
*************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <math.h>
#include <unordered_set>
#include <stdlib.h>
#include <algorithm>

#include "utility/MyFunction.h"

using namespace std;

class SNPs
{
    private:
        int max_pos;
        vector<int> v_all_pos;
        vector<int> v_check_all_pos;

        struct snp_node
        {
            string chr;
            int pos;
            int degree;
            int effect_score;
            float sum_score;
            vector<int> adj_pos;
            vector<float> adj_scores;
        };
        


        snp_node node;
        vector<snp_node> heap;
        vector<snp_node> heap_2;
        vector<snp_node> top_K;
        map<int,int> pos2heap;
        int size;
        int size2;
        
        void snp_node_clear(snp_node& node){
            node.chr = "";
            node.pos = 0;
            node.degree = 0;
            node.effect_score = 0;
            node.sum_score = 0;
            node.adj_pos.clear();
            node.adj_scores.clear();
        }

        void max_heap(int i){
            int largest;
            int left = 2*i;
            int right = 2*i +1;

            if( left <= size and heap[left].sum_score > heap[i].sum_score){
                largest = left;
            }else{
                largest = i;
            }
        
            if(right <= size and heap[right].sum_score > heap[largest].sum_score) largest = right;

            if(largest != i){
                //swap(heap[i], heap[largest]);
                swap(heap[i], heap[largest]);
                pos2heap[ heap[i].pos ] = i;
                pos2heap[ heap[largest].pos ] = largest;
                max_heap(largest);
            }
        }

        void run_max_heap(){
            for (int i = size / 2; i >= 1; i--) {
                max_heap(i);
            }
        }

        void remove_heap(int heap_pos){
                if(size == 0) cout << "error: heap size is 0" << endl;
                snp_node tem = heap[heap_pos];
                heap[heap_pos] = heap[size];
                heap[size] = tem;
                pos2heap[ heap[heap_pos].pos ] = heap_pos;
                pos2heap[ heap[size].pos ] = size;
                size--;
        }

        void update_heap(int snp_pos, float imputed_score){
            heap[ pos2heap[snp_pos] ].sum_score = heap[ pos2heap[snp_pos] ].sum_score - imputed_score;
            heap[ pos2heap[snp_pos] ].degree --;
        }

    public:


        void ld_scaning(string file_name, float cutoff = 0.8){
            ifstream fi(file_name);
            cout << "\nScaning LD file: " << file_name << endl;
            string line;
            vector<string> tokens;
            getline(fi, line);
            char delim = get_delim(line);
            max_pos = 0;
            int N = 0;
            int pos1 = 0;
            int pos2 = 0;
            float score1 = 0.0;
            float score2 = 0.0;
            float ld_r2 = 0.0;
            v_check_all_pos.clear();
            v_all_pos.clear();
            while(fi.good()){
                getline(fi, line);
                tokenize(line, tokens, delim);
                //cout << line << endl;
                if( tokens.size() < 9 ) continue; // expected each line have 9 colunms
                pos1 = atoi(tokens[1].c_str());
                pos2 = atoi(tokens[4].c_str());
                ld_r2 = atof(tokens[6].c_str()) * atof(tokens[6].c_str());
                if(ld_r2 < cutoff) continue; // skip line if ld_2 < cutoff.
                max_pos = pos1 > pos2 ? pos1 : pos2;
                if( max_pos > N){
                    N = max_pos;
                    v_check_all_pos.resize(N+1, 0);
                }
                v_check_all_pos[pos1] = 1;
                v_check_all_pos[pos2] = 1;
            }
            fi.close();
            v_all_pos.push_back(-1); // dummy element
            for(int i = 0; i < v_check_all_pos.size(); i++){
                if(v_check_all_pos[i] == 1) v_all_pos.push_back(i);
            }
            cout << "\tScaning ... Done!" << endl;
            cout << "\tMax pos is: " << max_pos << endl;
            cout << "\tNumber of unique SNP position is: " << v_all_pos.size() - 1 << endl; 
        }

        void heap_init(){
            snp_node_clear(node);
            heap.clear();
            heap.resize(v_all_pos.size(), node);
            for(int i = 1; i < v_all_pos.size(); i++){
                pos2heap[ v_all_pos[i] ] = i;
                heap[i].pos = v_all_pos[i];
                //cout << v_all_pos[i] << "\tdeg :  " << heap[pos2heap[v_all_pos[i]]].degree << "\tsum_score : " << heap[pos2heap[v_all_pos[i]]].sum_score << "\tsize adj vec : " << heap[pos2heap[v_all_pos[i]]].adj_pos.size() << endl;
            }
            size = heap.size()-1;
            cout << "\tCreated heap! Heap size: " <<  size << endl;
        }

        void read_ld_model(string file_name, float cutoff = 0.8){
            ifstream fi(file_name);
            cout << "\nReading imputed scores to the program: " << file_name << endl;
            string line;
            vector<string> tokens;
            getline(fi, line);
            char delim = get_delim(line);
            int pos1 = 0;
            int pos2 = 0;
            float score1 = 0.0;
            float score2 = 0.0;
            float ld_r2 = 0.0;
            string chr;
            while(fi.good()){
                getline(fi, line);
                tokenize(line, tokens, delim);
                //cout << line << endl;
                if( tokens.size() < 9 ) continue; // expected each line have 9 colunms
                pos1 = atoi(tokens[1].c_str());
                pos2 = atoi(tokens[4].c_str());
                ld_r2 = atof(tokens[6].c_str()) * atof(tokens[6].c_str());
                score1 = atof(tokens[7].c_str());
                score2 = atof(tokens[8].c_str());
                chr = tokens[0];
                if(ld_r2 < cutoff) continue; // skip line if ld_2 < cutoff.
                // pos1
                heap[pos2heap[pos1]].chr = chr;
                heap[pos2heap[pos1]].degree++;
                heap[pos2heap[pos1]].adj_pos.push_back(pos2);
                heap[pos2heap[pos1]].adj_scores.push_back(score1);
                heap[pos2heap[pos1]].sum_score = heap[pos2heap[pos1]].sum_score + score1; 
                // pos2
                heap[pos2heap[pos2]].chr = chr;
                heap[pos2heap[pos2]].degree ++;
                heap[pos2heap[pos2]].adj_pos.push_back(pos1);
                heap[pos2heap[pos2]].adj_scores.push_back(score2);
                heap[pos2heap[pos2]].sum_score = heap[pos2heap[pos2]].sum_score + score2; 
            }
            fi.close();
            cout << "\tReading ... Done!" << endl;
        }


        void tagging_with_heap(string out_file){
            cout << "\nTagging... " << endl;
            run_max_heap(); // 1st run
            int tag_heap = 1;
            int k = 0;
            ofstream fo;
            fo.open(out_file);
            fo << "chr" << "\t" << "pos" << "\t" << "sum_score" << "\t" << "degree"   << endl;
            while(size > 0 and heap[tag_heap].sum_score > 0.7){
                //cout << size << "\t" << heap[tag_heap].degree << endl;
                fo << heap[tag_heap].chr << "\t" << heap[tag_heap].pos << "\t" << heap[tag_heap].sum_score << "\t" << heap[tag_heap].degree  << endl;
                snp_node n = heap[tag_heap];
                remove_heap(tag_heap);
                for(int i = 0; i < n.adj_pos.size(); i++){
                    int id = n.adj_pos[i];
                    if(v_check_all_pos[id] == 1){
                        v_check_all_pos[id] = 0;
                        snp_node m = heap[pos2heap[id]];
                        remove_heap(pos2heap[id]);
                        for(int j = 0; j < m.adj_pos.size(); j++){
                            int id2 = m.adj_pos[j];
                            update_heap(id2, m.adj_scores[j]);
                        }
                    }
                }
                run_max_heap();
                k++;
            }
            fo.close();
            cout << "\tTagging ... Done!" << endl;
            cout << "\tSelected: " << k << " tag SNPs" << endl;
            cout << "\tPlease check results at: " << out_file << endl;
            cout << endl;
        }
};

#endif