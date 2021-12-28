#ifndef AGRS_PARSER
#define AGRS_PARSER
#include <string>
#include <vector>

using namespace std;

class SimpleArgsParser
{
private:
    vector<string> all_agrs;
public:
    SimpleArgsParser(int argc, char const *argv[]){
        for(int i = 1; i < argc; i++){
            all_agrs.push_back(argv[i]);
        };
    };
    string find_flag(string flag, string defaut = "", string error = ""){
        string res = defaut;
        int is_found = 0;
        for(int i = 0; i < all_agrs.size(); i++){
            if( flag.compare(all_agrs[i]) == 0){
                res = all_agrs[i+1];
                is_found = 1;
                break;
            }
        }
        if(is_found == 0){
            cout << "\tCannot find argument: '" << flag << "'. " << endl;
            if(defaut.length() > 0){
                cout << "\tDefault values is used: " << defaut << endl;
            }else if(error.length() > 0){
                cout << "\tPROGRAM STOP!! =====>> MISSING : " << flag << ". " << error << endl;
                exit(1);
            }else{
                cout << "\tMissing argument: " << flag << " ??? It's okay =====>> CONTINUE..." << endl;
            }
        };
        if(is_found == 1) cout << "\tArgument: '" << flag << "' is specified as: '"<< res << "'" << endl;
        return res;
    };
};
#endif