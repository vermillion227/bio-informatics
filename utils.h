#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>


using namespace std;

namespace utils {

bool cmp(pair<string, int> p1, pair<string,int> p2);

void ReversePattern(const string& ptn, string* reverse);

int SymbolToNumber(char c);

char NumberToSymbol(int n);

int PatternToNumber(const string& pattern);

void NumberToPattern(const int index, const int k, string* ptn);

void ComputeFrequencyArray(const string& s,
                           const int k,
                           vector<int> *freq);

void ComputeFrequencies(const string& s,
                        const int k,
                        map<string, int> *f);

void FindPos(const string& ptn, const string& genome,
             vector<int> *pos);

void Skew(const string& genome, const int i,
          vector<int> *s);

void FindMinims(const vector<int>& s,
                vector<int>* mins);

int HammingDistance(const string& a, const string& b);

void ApproxPatternMatch(const string& pattern,
                        const string& adn,
                        const int d,
                        vector<int>* pos);
int PatternCount(const string& adn, const string &pattern);

int ApproxPatternCount(const string& adn, const string& pattern, int d);

vector<string> Neighbors(const string& s, const int d);

void ComputeFrequenciesWithMismatch(const string& s,
                                    const int k,
                                    const int d,
                                    map<string, int> *f);

}  // namespace utils
