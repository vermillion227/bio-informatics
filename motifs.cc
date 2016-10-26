#include <iostream>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

using namespace std;

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

string ProfileMostProbable(const string &text,
                           const int k,
                           const vector<vector<double> > &profile) {
  string pattern;
  double max_prob = 0;

  for (int i = 0; i < text.length() - k + 1; i++) {
    string aux = text.substr(i, k);
    double probability = 1;
    for (int j = 0; j < k; j++) {
      int index = SymbolToNumber(aux[j]);
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

  for (string &motif : motifs) {
    for (int i = 0; i < k; i++) {
      profile[SymbolToNumber(motif[i])][i]++;
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < k; j++) {
      profile /= (2 * motifs.length());
    }
  }

  return profile;
}

vector<string> GreedyMotifSearch(const vector<string> &dna,
                                 const int k,
                                 const int t) {
  // TODO:
}

int main() {
  string t;
  int k;
  vector<vector<double> > prof;

  cin >> t >> k;
  for (int i = 0; i < 4; i++) {
    vector <double> aux;
    for (int j = 0; j < k; j++) {
      double p;
      cin >> p;
      aux.push_back(p);
    }
    prof.push_back(aux);
    aux.clear();
  }

  cout << ProfileMostProbable(t, k, prof);
}
