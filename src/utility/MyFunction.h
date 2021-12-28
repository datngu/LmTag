#ifndef MY_FUNCTION
#define MY_FUNCTION
#include <string>
#include <vector>

using namespace std;

void tokenize(string &s, vector<string> &tokens, char delim = '\t')
{
    tokens.clear();
    int prev = -1;
    for (int i = 0; i <= s.length(); ++i)
    {
        if ((i == s.length()) || (s[i] == delim))
        {
            int tl = i - prev - 1;
            if (tl > 0) tokens.push_back(s.substr(prev+1, tl));
            prev = i;
        }
    }
}

char get_delim(string &s){
    int num_space = 0;
    int num_tab = 0;
    for(int i = 0; i < s.length(); i++){
        if(s[i] == ' ') num_space++;
        if(s[i] == '\t') num_tab++;
    }
    if(num_space > 0 and num_space > num_tab)
        return ' ';
    else
        return '\t'; // tab is defaut value
}

#endif