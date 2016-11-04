#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <string>
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

vector<string> MotifEnumeration(const vector<string> &dna,
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

int Distance(const vector<string> &dna, const string &ptn) {
  int d = 0;
  for (const string &seq : dna) {
    int min = INT_MAX;
    for (int i = 0; i < seq.length() - ptn.length() + 1; i++) {
      string current = seq.substr(i, ptn.length());
      int hamming = utils::HammingDistance(ptn, current);
      if (hamming < min) {
        min = hamming;
      }
    }
    d += min;
  }

  return d;
}

string MedianString(const vector<string> &dna, const int k) {
  int distance = INT_MAX;
  string median;
  for (int i = 0; i < (1 << (2 * k)); i++) {
    string aux;
    utils::NumberToPattern(i, k, &aux);
    int d = Distance(dna, aux);
    if (d < distance) {
      distance = d;
      median = aux;
    }
  }

  return median;
}

string GetConsensus(const vector<string> motifs) {
  vector<vector<double> > profile = CreateProfile(motifs);
  return ProfileMostProbable(motifs[0], motifs[0].size(), profile);
}

int Score(const vector<string> motifs) {
  string consensus = GetConsensus(motifs);
  return Distance(motifs, consensus);
}

vector<string> GreedyMotifSearch(const vector<string> dna,
                                 const int k,
                                 const int t) {
  vector<string> best_motifs;
  for (const string &seq : dna) {
    best_motifs.push_back(seq.substr(0, k));
  }

  string dna_0 = dna[0];
  for (int i = 0; i < dna_0.length() - k + 1; i++) {
    string motif_0 = dna_0.substr(i, k);
    vector<string> new_motifs;

    new_motifs.push_back(motif_0);
    for (int j = 1; j < t; j++) {
      vector<vector<double> > profile = CreateProfile(new_motifs);
      string new_motif = ProfileMostProbable(dna[j], k, profile);
      new_motifs.push_back(new_motif);
    }

    if (Score(new_motifs) < Score(best_motifs)) {
      best_motifs = new_motifs;
    }
  }

  return best_motifs;
}

int main() {
  int t;
  int k, d;
  vector<string> dnas;

  cin >> k >> t;

  for (int i = 0; i < t; i++) {
    string aux;
    cin >> aux;
    dnas.push_back(aux);
  }

  vector<string> motifs = GreedyMotifSearch(dnas, k, t);
  for (const string& m : motifs) {
    cout << m << endl;
  }

  return 0;
}
