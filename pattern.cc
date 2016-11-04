#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

#include "utils.h"

using namespace std;

void FrequentWords(const string& adn, const int k,
                   vector<string> *mfp) {
  map<string, int> counts;

  utils::ComputeFrequencies(adn, k, &counts);

  int max = (max_element(counts.begin(), counts.end(), utils::cmp))->second;

  for (auto& p : counts) {
    if (p.second == max) {
      mfp->push_back(p.first);
    }
  }
}

void FindClumpPatterns(const string& genome,
                       const int k, const int L, const int t,
                       vector<string> *fp) {
  map<string, int> frequencies;
  utils::ComputeFrequencies(genome.substr(0, L), k, &frequencies);

  for (auto p : frequencies) {
    if (p.second >= t) {
      fp->push_back(p.first);
    }
  }

  for (int i = 1; i <= genome.size() - L; i++) {
    string first_pattern = genome.substr(i - 1, k);
    string last_pattern = genome.substr(i + L - k, k);

    frequencies[first_pattern]--;
    frequencies[last_pattern]++;

    if (frequencies[last_pattern] >= t &&
        (find(fp->begin(), fp->end(), last_pattern) == fp->end())) {
      fp->push_back(last_pattern);
    }
  }
}

void FrequentWordsWithMismatch(const string& adn, const int k,
                               const int d,
                               vector<string> *mfp) {
  map<string, int> counts;

  utils::ComputeFrequenciesWithMismatch(adn, k, d, &counts);

  int max = (max_element(counts.begin(), counts.end(), utils::cmp))->second;

  for (auto& p : counts) {
    if (p.second == max) {
      mfp->push_back(p.first);
    }
  }
}

void FrequentWordsWithMismatchAndReverse(const string& adn, const int k,
                                         const int d,
                                         vector<string> *mfp) {
  map<string, int> counts;
  string rev;
  utils::ReversePattern(adn, &rev);

  utils::ComputeFrequenciesWithMismatch(adn, k, d, &counts);
  utils::ComputeFrequenciesWithMismatch(rev, k, d, &counts);

  int max = (max_element(counts.begin(), counts.end(), utils::cmp))->second;

  for (auto& p : counts) {
    if (p.second == max) {
      mfp->push_back(p.first);
    }
  }
}
