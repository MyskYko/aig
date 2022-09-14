#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>

#include "aig.hpp"

using namespace std;

unsigned char getnoneofch (ifstream & file)
{
  int ch = file.get();
  assert(ch != EOF);
  return ch;
}
unsigned decode (ifstream & file)
{
  unsigned x = 0, i = 0;
  unsigned char ch;
  while ((ch = getnoneofch (file)) & 0x80)
    x |= (ch & 0x7f) << (7 * i++);
  return x | (ch << (7 * i));
}
void encode (ofstream & file, unsigned x)
{
  unsigned char ch;
  while (x & ~0x7f)
    {
      ch = (x & 0x7f) | 0x80;
      file.put(ch);
      x >>= 7;
    }
  ch = x;
  file.put(ch);
}

void aigman::read(string filename) {
  ifstream f(filename);
  {
    string str;
    getline(f, str);
    stringstream ss(str);
    getline(ss, str, ' ');
    assert(str == "aig");
    getline(ss, str, ' ');
    nObjs = stoi(str);
    getline(ss, str, ' ');
    nPis = stoi(str);
    getline(ss, str, ' ');
    assert(str == "0");
    getline(ss, str, ' ');
    nPos = stoi(str);
    getline(ss, str, ' ');
    nGates = stoi(str);
    assert(nObjs == nGates + nPis);
  }

  nObjs++; // constant
  vPos.resize(nPos);
  vObjs.resize(nObjs * 2);

  {
    string str;
    for(int i = 0; i < nPos; i++) {
      getline(f, str);
      vPos[i] = stoi(str);
    }
  }

  for(int i = nPis + 1; i < nObjs; i++) {
    vObjs[i + i] = i + i - decode(f);
    vObjs[i + i + 1] = vObjs[i + i] - decode(f);
  }

  f.close();

  fSorted = true;
}

int aigman::renumber_rec(int i, std::vector<int> & vObjsNew, int & nObjsNew) {
  if(i <= nPis) return i;
  if(vValues[i]) return vValues[i];
  int l = renumber_rec(vObjs[i + i] >> 1, vObjsNew, nObjsNew);
  int r = renumber_rec(vObjs[i + i + 1] >> 1, vObjsNew, nObjsNew);
  vObjsNew[nObjsNew + nObjsNew] = (l << 1) + (vObjs[i + i] & 1);
  vObjsNew[nObjsNew + nObjsNew + 1] = (r << 1) + (vObjs[i + i + 1] & 1);
  vValues[i] = nObjsNew++;
  return vValues[i];
}

void aigman::renumber() {
  vDeads.clear();
  vValues.clear();
  vValues.resize(nObjs);
  int nObjsNew = nPis + 1;
  std::vector<int> vObjsNew((nObjsNew + nGates) * 2);
  for(int i = 0; i < nPos; i++) {
    int j = vPos[i] >> 1;
    vPos[i] = (renumber_rec(j, vObjsNew, nObjsNew) << 1) + (vPos[i] & 1);
  }
  vObjs = vObjsNew;
  nObjs = nObjsNew;
  assert(nObjs == 1 + nPis + nGates);
}

void aigman::write(string filename) {
  if(!fSorted) {
    renumber();
    fSorted = true;
  }
  ofstream f(filename);
  f << "aig " << nObjs - 1 << " " << nPis << " 0 " << nPos << " " << nObjs - nPis - 1 << endl;
  for(int i = 0; i < nPos; i++) {
    f << vPos[i] << endl;
  }
  for(int i = nPis + 1; i < nObjs; i++) {
    encode(f, i + i - vObjs[i + i]);
    encode(f, vObjs[i + i] - vObjs[i + i + 1]);
  }
  f.close();
}

int aigman::getsim(int i) {
  assert((i >> 1) < nObjs);
  if(i & 1) {
    return vSims[i >> 1] ^ 0xffffffffffffffffull;
  }
  return vSims[i >> 1];
}

void aigman::simulate(vector<unsigned long long> const & inputs) {
  assert(inputs.size() == (size_t)nPis);
  vSims.resize(nObjs);
  vSims[0] = 0; // constant
  for(int i = 0; i < nPis; i++) {
    vSims[i + 1] = inputs[i];
  }
  for(int i = nPis + 1; i < nObjs; i++) {
    vSims[i] = getsim(vObjs[i + i]) & getsim(vObjs[i + i + 1]);
  }
}

void aigman::supportfanouts_rec(int i) {
  if(i <= nPis || vValues[i]) return;
  {
    int j = vObjs[i + i] >> 1;
    vvFanouts[j].push_back(i);
    supportfanouts_rec(j);
  }
  {
    int j = vObjs[i + i + 1] >> 1;
    vvFanouts[j].push_back(i);
    supportfanouts_rec(j);
  }
  vValues[i] = 1;
}

void aigman::supportfanouts() {
  vValues.clear();
  vValues.resize(nObjs);
  vvFanouts.clear();
  vvFanouts.resize(nObjs);
  for(int i = 0; i < nPos; i++) {
    int j = vPos[i] >> 1;
    vvFanouts[j].push_back(- i - 1);
    supportfanouts_rec(j);
  }
}

void aigman::supportlevels() {
  if(vvFanouts.empty()) {
    supportfanouts();
  }
  vValues.clear();
  vValues.resize(nObjs);
  vLevels.resize(nObjs);
  std::vector<int> * cur = new std::vector<int>;
  std::vector<int> * next = new std::vector<int>;
  for(int i = 0; i < nPos; i++) {
    cur->push_back(vPos[i] >> 1);
  }
  int level = 0;
  while(!cur->empty()) {
    for(int i : *cur) {
      vValues[i] = 1;
      vLevels[i] = level;
      if(i <= nPis) continue;
      {
        int j = vObjs[i + i] >> 1;
        bool f = 1;
        for(int k : vvFanouts[j]) {
          if(k < 0) {
            continue;
          }
          if(!vValues[k]) {
            f = 0;
            break;
          }
        }
        if(f) {
          next->push_back(j);
        }
      }
      {
        int j = vObjs[i + i + 1] >> 1;
        bool f = 1;
        for(int k : vvFanouts[j]) {
          if(k < 0) {
            continue;
          }
          if(!vValues[k]) {
            f = 0;
            break;
          }
        }
        if(f) {
          next->push_back(j);
        }
      }
    }
    level++;
    auto tmp = cur;
    cur = next;
    next = tmp;
    next->clear();
  }
  delete cur;
  delete next;
}

void aigman::removenode(int i) {
  if(i <= nPis) return;
  if(vDeads.empty()) {
    vDeads.resize(nObjs);
  } else if(vDeads[i]) return;
  vDeads[i] = 1;
  nGates--;
  fSorted = false;
  if(vvFanouts.empty()) return;
  assert(vvFanouts[i].empty());
  {
    int j = vObjs[i + i] >> 1;
    auto it = find(vvFanouts[j].begin(), vvFanouts[j].end(), i);
    assert(it != vvFanouts[j].end());
    vvFanouts[j].erase(it);
    if(vvFanouts[j].empty()) {
      removenode(j);
    }
    else if(!vLevels.empty()) {
      vLevels[j] = 0;
      for(int k : vvFanouts[j]) {
        if(k < 0) continue;
        if(vLevels[j] < vLevels[k] + 1) {
          vLevels[j] = vLevels[k] + 1;
        }
      }
    }
  }
  {
    int j = vObjs[i + i + 1] >> 1;
    auto it = find(vvFanouts[j].begin(), vvFanouts[j].end(), i);
    assert(it != vvFanouts[j].end());
    vvFanouts[j].erase(it);
    if(vvFanouts[j].empty()) {
      removenode(j);
    }
    else if(!vLevels.empty()) {
      vLevels[j] = 0;
      for(int k : vvFanouts[j]) {
        if(k < 0) continue;
        if(vLevels[j] < vLevels[k] + 1) {
          vLevels[j] = vLevels[k] + 1;
        }
      }
    }
  }
}

void aigman::replacenode(int i, int j, bool prop) {
  assert(i >= 0);
  assert(j >= 0);
  if(vDeads.empty()) {
    vDeads.resize(nObjs);
  }
  assert(!vDeads[i]);
  assert(!vDeads[j >> 1]);
  if(vvFanouts.empty()) {
    supportfanouts();
  }
  std::vector<int> targets = vvFanouts[i];
  if(i == (j >> 1)) {
    if(j & 1) {
      for(int k : vvFanouts[i]) {
        if(k < 0) {
          vPos[- k - 1] ^= 1;
          continue;
        }
        if((vObjs[k + k] >> 1) == i) {
          vObjs[k + k] ^= 1;
        }
        if((vObjs[k + k + 1] >> 1) == i) {
          vObjs[k + k + 1] ^= 1;
        }
      }
    }
  }
  else {
    for(int k : vvFanouts[i]) {
      if(k < 0) {
        vPos[- k - 1] = j ^ (vPos[- k - 1] & 1);
      }
      else {
        if((vObjs[k + k] >> 1) == i) {
          vObjs[k + k] = j ^ (vObjs[k + k] & 1);
        }
        if((vObjs[k + k + 1] >> 1) == i) {
          vObjs[k + k + 1] = j ^ (vObjs[k + k + 1] & 1);
        }
        if(!vLevels.empty() && vLevels[j >> 1] < vLevels[k] + 1) {
          vLevels[j >> 1] = vLevels[k] + 1;
        }
      }
      vvFanouts[j >> 1].push_back(k);
    }
    vvFanouts[i].clear();
    removenode(i);
  }

  fSorted = false;
  if(!prop) return;

  for(int k : targets) {
    if(k < 0) {
      continue;
    }
    if(vDeads[k]) continue;
    if((vObjs[k + k] >> 1) == (vObjs[k + k + 1] >> 1)) {
      if((vObjs[k + k] & 1) == (vObjs[k + k + 1] & 1)) {
        replacenode(k, 1);
      }
      else {
        replacenode(k, 0);
      }
      continue;
    }
    if(!(vObjs[k + k] >> 1)) {
      if(vObjs[k + k] & 1) {
        replacenode(k, vObjs[k + k + 1]);
      }
      else {
        replacenode(k, 0);
      }
      continue;
    }
    if(!(vObjs[k + k + 1] >> 1)) {
      if(vObjs[k + k + 1] & 1) {
        replacenode(k, vObjs[k + k]);
      }
      else {
        replacenode(k, 0);
      }
    }
  }
}

void aigman::save(int i) {
  if(backup.size() <= (size_t)i) {
    backup.resize(i + 1);
  }
  backup[i] = *this;
}

void aigman::load(int i) {
  assert(backup.size() > (size_t)i);
  *this = backup[i];
}

int aigman::extract_rec(aigman * aig_new, int i) {
  if(vValues[i]) {
    return vValues[i];
  }
  int i0 = vObjs[i + i] >> 1;
  bool c0 = vObjs[i + i] & 1;
  i0 = extract_rec(aig_new, i0);
  int i1 = vObjs[i + i + 1] >> 1;
  bool c1 = vObjs[i + i + 1] & 1;
  i1 = extract_rec(aig_new, i1);
  aig_new->vObjs.push_back((i0 << 1) ^ c0);
  aig_new->vObjs.push_back((i1 << 1) ^ c1);
  vValues[i] = aig_new->nObjs++;
  aig_new->nGates++;
  return vValues[i];
}

aigman * aigman::extract(vector<int> const & inputs, vector<int> const & outputs) {
  vValues.clear();
  vValues.resize(nObjs);
  aigman * aig_new = new aigman(inputs.size(), 0);
  for(int i = 0; (size_t)i < inputs.size(); i++) {
    vValues[inputs[i]] = i + 1;
  }
  for(int i: outputs) {
    aig_new->vPos.push_back((extract_rec(aig_new, i >> 1) << 1) ^ (i & 1));
    aig_new->nPos++;
  }
  return aig_new;
}


void aigman::import(aigman * p, vector<int> const & inputs, vector<int> const & outputs) {
  if(vvFanouts.empty()) {
    supportfanouts();
  }
  p->vValues.clear();
  p->vValues.resize(p->nObjs);
  assert(inputs.size() == (size_t)p->nPis);
  for(int i = 0; i < p->nPis; i++) {
    p->vValues[i + 1] = inputs[i];
  }
  for(int i = p->nPis + 1; i < p->nObjs; i++) {
    int i0 = p->vObjs[i + i] >> 1;
    bool c0 = p->vObjs[i + i] & 1;
    i0 = p->vValues[i0];
    int i1 = p->vObjs[i + i + 1] >> 1;
    bool c1 = p->vObjs[i + i + 1] & 1;
    i1 = p->vValues[i1];
    vObjs.push_back((i0 << 1) ^ c0);
    vObjs.push_back((i1 << 1) ^ c1);
    vvFanouts[i0].push_back(nObjs);
    vvFanouts[i1].push_back(nObjs);
    p->vValues[i] = nObjs++;
    nGates++;
    vvFanouts.resize(nObjs);
  }
  vDeads.resize(nObjs);
  assert(outputs.size() <= (size_t)p->nPos);
  for(int i = 0; (size_t)i < outputs.size(); i++) {
    int j = (p->vValues[p->vPos[i] >> 1] << 1) ^ (p->vPos[i] & 1) ^ (outputs[i] & 1);
    replacenode(outputs[i] >> 1, j);
  }
}
