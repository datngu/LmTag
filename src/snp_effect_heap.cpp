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
            string id;
            int degree;
            float effect_score;
            float effect_score_orignal;
            float sum_score;
            vector<int> adj_pos;
            vector<float> adj_scores;
        };

        struct snp_node_2
        {
            int pos;
            float effect_score;
            //float effect_score_orignal;
            float sum_score;
        };

        snp_node node;
        snp_node_2 node_2;
        vector<snp_node> heap;
        //vector<snp_node_2> sub_heap;
        vector<snp_node_2> heap_2;
        vector<snp_node> top_K;
        map<int,int> pos2heap;
        int size;
        int size2;
        float min_effect_score = 0.0;
        float max_effect_score = 0.0;
        
        void snp_node_clear(snp_node& node){
            node.chr = "";
            node.pos = 0;
            node.degree = 0;
            node.effect_score = 0.0;
            node.effect_score_orignal = 0.0;
            node.sum_score = 0.0;
            node.adj_pos.clear();
            node.adj_scores.clear();
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
                swap(heap[i], heap[largest]);
                // maintain 2 heaps in the same time
                //swap(sub_heap[i], sub_heap[largest]);
                pos2heap[ heap[i].pos ] = i;
                pos2heap[ heap[largest].pos ] = largest;
                max_heap(largest);
            }
        }

        void max_heap_2(int i){
            int largest;
            int left = 2*i;
            int right = 2*i +1;

            if( left <= size2 and heap_2[left].sum_score > heap_2[i].sum_score){
                largest = left;
            }else{
                largest = i;
            }
        
            if(right <= size2 and heap_2[right].sum_score > heap_2[largest].sum_score) largest = right;

            if(largest != i){
                //swap(heap_2[i], heap_2[largest]);
                swap(heap_2[i], heap_2[largest]);
                max_heap_2(largest);
            }
        }

        void run_max_heap(){
            for (int i = size / 2; i >= 1; i--) {
                max_heap(i);
            }
        }

        void run_max_heap_2(){
            for (int i = size2 / 2; i >= 1; i--) {
                max_heap_2(i);
            }
        }

        snp_node_2 pop_max_heap_2(){
            if(size2 == 0) cout << "error: heap size is 0" << endl;
            snp_node_2 tem = heap_2[1];
            heap_2[1] = heap_2[size2];
            heap_2[size2] = tem;
            size2 --;
            run_max_heap_2();
            return tem;
        }

        void copy_heap_2(int n){
            heap_2.clear();
            heap_2.resize(n+1);
            for(int i = 0; i < n+1; i++){
                heap_2[i].pos = heap[i].pos;
                heap_2[i].effect_score = heap[i].effect_score;
                //heap_2[i].effect_score_orignal = heap[i].effect_score_orignal;
                heap_2[i].sum_score = heap[i].sum_score;
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

        int find_tag_snp(int K){ // return tagSNP genome position
            snp_node_2 tem;
            int res = 0;
            float max_eff = 0.0;
            if(K <= 1){
                res = heap[1].pos;
            }else{
                size2 = size;
                copy_heap_2(size2);
                for(int i = 0; i < K; i++){
                    tem = pop_max_heap_2();
                    if(i == 0){
                        res = tem.pos;
                        max_eff = tem.effect_score;
                    }else{
                        if(tem.effect_score > max_eff){
                            res = tem.pos;
                            max_eff = tem.effect_score;     
                        }
                    }
                }
            }
            return res;
        }


    public:

        void read_snp_effect( string file_name){
            ifstream fi(file_name);
            cout << "\nReading SNP effect file: " << file_name << endl;
            string line;
            float effect_score = 0.0;
            int tem_pos = 0;
            vector<string> tokens;
            getline(fi, line);
            char delim = get_delim(line);
            while(fi.good()){
                getline(fi, line);
                tokenize(line, tokens, delim);
                if( tokens.size() < 3 ) continue; // expected each line have 3 colunms: chr_id, pos, and effect score.
                tem_pos = atoi(tokens[1].c_str());
                //if(node.pos > max_pos) max_pos = node.pos;
                effect_score = atof(tokens[2].c_str());
                // find min and max effect score
                if(effect_score > max_effect_score) max_effect_score = effect_score;
                if(effect_score < min_effect_score) min_effect_score = effect_score;
                heap[pos2heap[tem_pos]].effect_score = effect_score;
                heap[pos2heap[tem_pos]].effect_score_orignal = effect_score;
            }
            fi.close();
            cout << "\tReading ... Done!" << endl;
            cout << "\tMax effect score: " << max_effect_score << endl;
            cout << "\tMin effect score: " << min_effect_score << endl;
        }

        void read_snp_vip( string file_name){
            ifstream fi(file_name);
            cout << "\nReading list of VIP SNPs: " << file_name << endl;
            string line;
            int tem_pos = 0;
            vector<string> tokens;
            getline(fi, line);
            char delim = get_delim(line);
            while(fi.good()){
                getline(fi, line);
                tokenize(line, tokens, delim);
                if( tokens.size() < 2 ) continue; // expected each line have 2 colunms: chr_id, pos
                tem_pos = atoi(tokens[1].c_str());
                // effect_score of VIP SNP is set with max_effect_score + 1 to prioritize it in tag SNP selection
                if(heap[1].chr.compare(tokens[0]) == 0) heap[pos2heap[tem_pos]].effect_score = max_effect_score + 1;
            }
            fi.close();
            cout << "\tReading ... Done!" << endl;
        }
  
        void read_snp_exclude( string file_name){
            ifstream fi(file_name);
            cout << "\nReading list of excluded SNPs: " << file_name << endl;
            string line;
            int tem_pos = 0;
            vector<string> tokens;
            getline(fi, line);
            char delim = get_delim(line);
            while(fi.good()){
                getline(fi, line);
                tokenize(line, tokens, delim);
                if( tokens.size() < 2 ) continue; // expected each line have 2 colunms: chr_id, pos
                tem_pos = atoi(tokens[1].c_str());
                // effect_score of excluded SNP is set with min_effect_score - 1 to remove it from tag SNP selection
                if(heap[1].chr.compare(tokens[0]) == 0 ) heap[pos2heap[tem_pos]].effect_score = min_effect_score - 1;
            }
            fi.close();
            cout << "\tReading ... Done!" << endl;
        }
    

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
            heap_init();  
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
            string id1 = "";
            string id2 = "";
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
                id1 = tokens[2];
                id2 = tokens[5];
                ld_r2 = atof(tokens[6].c_str()) * atof(tokens[6].c_str());
                score1 = atof(tokens[7].c_str());
                score2 = atof(tokens[8].c_str());
                chr = tokens[0];
                if(ld_r2 < cutoff) continue; // skip line if ld_2 < cutoff.
                // pos1
                heap[pos2heap[pos1]].chr = chr;
                heap[pos2heap[pos1]].id = id1;
                heap[pos2heap[pos1]].degree++;
                heap[pos2heap[pos1]].adj_pos.push_back(pos2);
                heap[pos2heap[pos1]].adj_scores.push_back(score1);
                heap[pos2heap[pos1]].sum_score = heap[pos2heap[pos1]].sum_score + score1; 
                // pos2
                heap[pos2heap[pos2]].chr = chr;
                heap[pos2heap[pos2]].id = id2;
                heap[pos2heap[pos2]].degree ++;
                heap[pos2heap[pos2]].adj_pos.push_back(pos1);
                heap[pos2heap[pos2]].adj_scores.push_back(score2);
                heap[pos2heap[pos2]].sum_score = heap[pos2heap[pos2]].sum_score + score2; 
            }
            fi.close();
            /*
            sub_heap.resize(heap.size());
            // copy sub_heap from heap
            for(int i = 0; i < heap.size(); i++){
                sub_heap[i].pos = heap[i].pos;
                sub_heap[i].sum_score = heap[i].sum_score;
                sub_heap[i].effect_score = heap[i].effect_score;
            }
            */
            cout << "\tReading ... Done!" << endl;
        }

        void tagging_with_heap(string out_file, int K = 1){
            cout << "\nTagging... " << endl;
            cout << "K setting: " << K << endl;
            run_max_heap(); // 1st run
            int tag_pos = find_tag_snp(K);
            int tag_heap = pos2heap[tag_pos];
            int k = 0;
            string flag  = "normal";
            if(heap[tag_heap].effect_score == min_effect_score - 1){
                flag = "excluded";
            }else if(heap[tag_heap].effect_score == max_effect_score + 1){
                flag = "vip";
            }else{
                flag  = "normal";
            };
            ofstream fo;
            fo.open(out_file);
            fo << "chr" << "\t" << "pos" << "\t" << "id" << "\t" << "sum_score" << "\t" << "degree" << "\t" << "effect_score" << "\t" << "flag" << "\t" << "tagged_pos"  << endl;
            while(size > 0 and heap[tag_heap].sum_score > 0.7){
                //cout << size << "\t" << heap[tag_heap].degree << endl;
                if(heap[tag_heap].effect_score == min_effect_score - 1){
                    flag = "excluded";
                }else if(heap[tag_heap].effect_score == max_effect_score + 1){
                    flag = "vip";
                }else{
                    flag  = "normal";
                };
                snp_node n = heap[tag_heap];
                string tagged_id = "";
                for(int i = 0; i < n.adj_pos.size(); i++){
                    int id = n.adj_pos[i];
                    if(tagged_id.length() <= 1){
                        tagged_id = tagged_id + to_string(id);
                    }else{
                        tagged_id = tagged_id + ";" + to_string(id);
                    }
                }
                fo << heap[tag_heap].chr << "\t" << heap[tag_heap].pos << "\t" << heap[tag_heap].id << "\t" << heap[tag_heap].sum_score << "\t" << heap[tag_heap].degree  <<  "\t" << heap[tag_heap].effect_score_orignal <<  "\t" << flag <<  "\t" << tagged_id << endl;
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
                tag_pos = find_tag_snp(K);
                tag_heap = pos2heap[tag_pos];
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