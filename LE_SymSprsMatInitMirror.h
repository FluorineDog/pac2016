#include "LE_SymSprsMatAux.h"
#include <atomic>
#include <cfloat>
#include <condition_variable>
#include <iostream>
#include <memory>
#include <mutex>
#include <pthread.h>
#include <thread>
// #include <hbwmalloc.h>
using std::cout;
using std::cerr;
using std::endl;

#include <cassert>
#include <map>
#include <set>
#include <vector>

#define dog_calloc(size, type)                                                 \
  (type *)aligned_alloc(64, ((size) * sizeof(type) + 63) & ~63);
#define dog_free(ptr) free(ptr);

constexpr int BLOCK = 528;
constexpr int THREAD_NUM = 0; // magic !! don't modify !!
typedef double __attribute((aligned(64))) aligned_double;
using su_t = SprsUMatRealStru;
void AdditionLU_SymbolicSymG(su_t *pFU){
  int *rs_u = pFU->uMax.rs_u;
  int *j_u = pFU->uMax.j_u;
	int iDim = pFU->uMax.iDim;
  std::vector<std::set<int>> refTable(iDim + 1);
  for (int i = 1; i < iDim + 1; ++i) {
    int kbeg = rs_u[i];
		int kend = rs_u[i + 1];
		for(int k = kbeg; k < kend; ++k){
			int j = j_u[k];
			refTable[j].insert(i);
		}
	}
  std::vector<std::set<int>> seqPart(iDim + 1);
	std::vector<std::vector<int>> paraPart(iDim + 1);
	for(int basei = 1; basei < iDim + 1; ++basei){
		int iend = std::min(basei + BLOCK, iDim);
    for (int i = basei; i > iend; i--) {
      const auto &ori = refTable[i];
      for (auto x : ori) {
        if (x > basei) {
          paraPart[i].push_back(x);
        } else {
          seqPart[i].insert(x);
          assert(x != i);
          seqPart[i].insert(seqPart[x].begin(), seqPart[x].end());
        }
      }
    }

	}
}