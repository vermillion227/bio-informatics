#include <iostream>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

#include "utils.h"

using namespace std;

string ProfileMostProbable(const string &text,
                           const int k,
                           const vector<vector<double> > &profile) {
  string pattern;
  double max_prob = 0;

  for (int i = 0; i < text.length() - k + 1; i++) {
    string aux = text.substr(i, k);
    double probability = 1;
    for (int j = 0; j < k; j++) {
      int index = utils::SymbolToNumber(aux[j]);
      probability *= profile[index][j];
    }

    if (probability > max_prob) {
      pattern = aux;
      max_prob = probability;
    }
  }

  return pattern;
}

vector<vector<double> > CreateProfile(const vector<string> &motifs) {
  int k = motifs[0].length();
  vector<vector<double> > profile;
  vector<double> line(k, 1);

  for (int i = 0; i < 4; i++) {
    profile.push_back(line);
  }

  for (const string &motif : motifs) {
    for (int i = 0; i < k; i++) {
      profile[utils::SymbolToNumber(motif[i])][i]++;
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < k; j++) {
      profile[i][j] /= (2 * motifs.size());
    }
  }

  return profile;
}

bool CommonPattern(const vector<string> &dna,
                   const string &pattern,
                   const int d) {
  vector<string> neighs = utils::Neighbors(pattern, d);
  for (const string &d : dna) {
    bool found = false;
    for (const string &p : neighs) {
      if (d.find(p) != string::npos) {
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }

  return true;
}

vector<string> GreedyMotifSearch(const vector<string> &dna,
                                 const int k,
                                 const int d) {
  vector<string> patterns;

  for (const string& current : dna) {
    for (int j = 0; j < current.length() - k + 1; j++) {
      string pattern = current.substr(j, k);
      vector<string> neighbors = utils::Neighbors(pattern, d);
      for (const string &pattern_ : neighbors) {
        if (CommonPattern(dna, pattern_, d)) {
          patterns.push_back(pattern_);
        }
      }
    }
  }

  sort(patterns.begin(), patterns.end());
  auto last = unique(patterns.begin(), patterns.end());
  patterns.erase(last, patterns.end());

  return patterns;
}

int main() {
  string t;
  int k, d;
  vector<string> dnas;

  cin  >> k >> d;

  while (cin >> t) {
    dnas.push_back(t);
  }

  vector<string> pp = GreedyMotifSearch(dnas, k, d);
  for (const string& p : pp) {
    cout << p << " ";
  }
  cout << endl;
}
