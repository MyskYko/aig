#pragma once

#include <string>
#include <vector>
#include <random>

class aigman {
public:
  int nPis;
  int nPos;
  int nGates;
  int nObjs;
  bool fSorted;
  std::vector<int> vPos;
  std::vector<int> vObjs;
  std::vector<int> vValues;

  std::vector<bool> vDeads;
  std::vector<std::vector<int> > vvFanouts;

  std::vector<unsigned long long> vSims;

  std::vector<aigman> backup;

  aigman() {};
  aigman(int nPis, int nPos) : nPis(nPis), nPos(nPos) {
    nGates = 0;
    nObjs = nPis + 1;
    vObjs.resize(nObjs * 2);
    for(int i = 0; i < nPos; i++) {
      vPos.push_back(0);
    }
    fSorted = true;
  };
  aigman(const aigman & x) {
    nPis = x.nPis;
    nPos = x.nPos;
    nGates = x.nGates;
    nObjs = x.nObjs;
    fSorted = x.fSorted;
    vPos = x.vPos;
    vObjs = x.vObjs;
    vValues = x.vValues;
    vDeads = x.vDeads;
    vvFanouts = x.vvFanouts;
    vSims = x.vSims;
  }
  aigman &operator=(const aigman & x) {
    nPis = x.nPis;
    nPos = x.nPos;
    nGates = x.nGates;
    nObjs = x.nObjs;
    fSorted = x.fSorted;
    vPos = x.vPos;
    vObjs = x.vObjs;
    vValues = x.vValues;
    vDeads = x.vDeads;
    vvFanouts = x.vvFanouts;
    vSims = x.vSims;
    return *this;
  }
  void clear() {
    nPis = 0;
    nPos = 0;
    nGates = 0;
    nObjs = 0;
    fSorted = true;
    vPos.clear();
    vObjs.clear();
    vValues.clear();
    vDeads.clear();
    vvFanouts.clear();
    vSims.clear();
  }

  void read(std::string filename);
  void write(std::string filename);

  void getgates_rec(std::vector<int> & gates, int i);
  void getgates(std::vector<int> & gates);

  int getsim(int i);
  void simulate(std::vector<unsigned long long> const & inputs);

  void supportfanouts_rec(int i);
  void supportfanouts();
  void removenode(int i);
  void replacenode(int i, int j, bool prop = 1);

  int renumber_rec(int i, std::vector<int> & vObjsNew, int & nObjsNew);
  void renumber();

  void save(int i = 0);
  void load(int i = 0);

  int extract_rec(aigman * aig_new, int i);
  aigman * extract(std::vector<int> const & inputs, std::vector<int> const & outputs);

  void import(aigman * p, std::vector<int> const & inputs, std::vector<int> const & outputs);

  void getfocone_rec(std::vector<int> & gates, int i);
  void getfocone(std::vector<int> const & nodes, std::vector<int> & gates);

  void resimulate(std::vector<int> const & nodes, std::vector<unsigned long long> const & values);
};
