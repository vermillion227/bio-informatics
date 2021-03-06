#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

#include "utils.h"

using namespace std;

namespace utils {

bool cmp(pair<string, int> p1, pair<string,int> p2) {
  return p1.second < p2.second;
}

void ReversePattern(const string& ptn, string* reverse) {
  map<char, char> encode;
  encode['A'] = 'T';
  encode['G'] = 'C';
  encode['C'] = 'G';
  encode['T'] = 'A';

  for (int i = ptn.length() - 1; i >= 0; i--) {
    reverse->push_back(encode[ptn[i]]);
  }
}

int SymbolToNumber(char c) {
  switch(c) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      return -1;
  }
}

char NumberToSymbol(int n) {
  switch(n) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    default:
      return '-';
  }
}

int PatternToNumber(const string& pattern) {
  if (!pattern.size()) {
    return 0;
  }

  return 4 * PatternToNumber(pattern.substr(0, pattern.size() - 1)) +
      SymbolToNumber(pattern[pattern.size() - 1]);
}

void NumberToPattern(const int index, const int k, string* ptn) {
  if (k == 1) {
    ptn->push_back(NumberToSymbol(index));
    reverse(ptn->begin(), ptn->end());
    return;
  }

  int q = index / 4;
  int r = index % 4;
  ptn->push_back(NumberToSymbol(r));
  NumberToPattern(q, k - 1, ptn);
}

void ComputeFrequencyArray(const string& s,
                           const int k,
                           vector<int> *freq) {
  for (int i = 0; i <= s.size() - k; i++) {
    int index = PatternToNumber(s.substr(i, k));
    (*freq)[index]++;
  }
}

void ComputeFrequencies(const string& s,
                        const int k,
                        map<string, int> *f) {
  for (int i = 0; i <= s.length() - k; i++) {
    string pattern = s.substr(i, k);
    (*f)[pattern]++;
  }
}

void FindPos(const string& ptn, const string& genome,
             vector<int> *pos) {
  for (int i = 0; i <= genome.size() - ptn.size(); i++) {
    if (genome.substr(i, ptn.size()) == ptn) {
      pos->push_back(i);
    }
  }
}

void Skew(const string& genome, const int i,
          vector<int> *s) {
  s->push_back(0);
  for (int i = 1; i <= genome.length(); i++) {
    if (genome[i - 1] == 'C') {
        s->push_back((*s)[i - 1] - 1);
    } else if (genome[i - 1] == 'G') {
        s->push_back((*s)[i - 1] + 1);
    } else {
        s->push_back((*s)[i - 1]);
    }
  }
}

void FindMinims(const vector<int>& s,
                vector<int>* mins) {
  int m = s[0];
  for (int i = 1; i < s.size(); i++) {
    if (s[i] < m) {
      m = s[i];
    }
  }

  for (int i = 0; i < s.size(); i++) {
    if (s[i] == m) {
      mins->push_back(i);
    }
  }
}

int HammingDistance(const string& a, const string& b) {
  int d = 0;

  for (int i = 0; i < a.length(); i++) {
    if (a[i] != b[i]) {
      d++;
    }
  }

  return d;
}

int PatternCount(const string& adn, const string& pattern) {
  int count = 0;

  for (int i = 0; i <= adn.size() - pattern.length(); i++) {
    if (adn.substr(i, pattern.length()) == pattern) {
      count++;
    }
  }

  return count;
}

void ApproxPatternMatch(const string& pattern,
                        const string& adn,
                        const int d,
                        vector<int>* pos) {
  for (int i = 0; i <= adn.size() - pattern.size(); i++) {
    int dist = HammingDistance(adn.substr(i, pattern.length()), pattern);
    if (dist <= d) {
      pos->push_back(i);
    }
  }
}

int ApproxPatternCount(const string& adn, const string& pattern, int d) {
  int count = 0;

  for (int i = 0; i <= adn.size() - pattern.length(); i++) {
    int dist = HammingDistance(adn.substr(i, pattern.length()), pattern);
    if (dist <= d) {
      count++;
    }
  }

  return count;
}

vector<string> Neighbors(const string& s, const int d) {
  vector<string> neighs;
  if (d == 0) {
    neighs.push_back(s);
    return neighs;
  }

  if (s.length() == 1) {
    neighs.push_back("A");
    neighs.push_back("C");
    neighs.push_back("G");
    neighs.push_back("T");
    return neighs;
  }

  string suffix = s.substr(1, s.length() - 1);
    vector<string> suffix_neighs = Neighbors(suffix, d);
      for (string& n : suffix_neighs) {
        if (HammingDistance(n, suffix) < d) {
          string new_n = "";
          for (int i = 0; i < 4; i++) {
            new_n = new_n + NumberToSymbol(i);
            new_n = new_n + n;
            neighs.push_back(new_n);
            new_n.clear();
          }
        } else {
          string new_n = "";
          new_n = new_n + s[0];
          new_n = new_n + n;
          neighs.push_back(new_n);
        }
      }

  return neighs;
}

void ComputeFrequenciesWithMismatch(const string& s,
                                    const int k,
                                    const int d,
                                    map<string, int> *f) {
  for (int i = 0; i <= s.length() - k; i++) {
    string pattern = s.substr(i, k);
    vector<string> neighbors = Neighbors(pattern, d);
    for (string& neighbor : neighbors) {
      (*f)[neighbor]++;
    }
  }
}

}  // namespace utils
