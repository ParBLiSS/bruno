/**
 * @file		quicktest.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


#include <vector>
#include <limits>
#include <cstdlib>
#include <unistd.h>

#include "sys/types.h"
#include "sys/sysinfo.h"

#include "common/Kmer.hpp"
#include "common/alphabets.hpp"
#include "common/AlphabetTraits.hpp"

#include "retired/kmer_index_element.hpp"
#include "io/sequence_iterator.hpp"
#include "common/sequence.hpp"
#include "io/fastq_loader.hpp"


using namespace std;
using namespace bliss::index;
using namespace bliss::io;


typedef bliss::Kmer<21, DNA, uint64_t> KmerType;


typedef KmerIndexElement<KmerType > KmerIndexType1;
typedef KmerIndexElementWithId<KmerType, bliss::io::FASTQSequenceId > KmerIndexType2;
typedef KmerIndexElementWithIdAndQuality<KmerType, bliss::io::FASTQSequenceId, float > KmerIndexType3;
typedef KmerIndexElementWithIdAndQuality<KmerType, bliss::io::FASTQSequenceId, double > KmerIndexType4;


void checkMemUsed(long long &phyMemUsed, long long &swapUsed, bool print) {
  //from http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
  struct sysinfo memInfo;

  sysinfo (&memInfo);
  phyMemUsed = (memInfo.totalram - memInfo.freeram) * memInfo.mem_unit;
  swapUsed = (memInfo.totalswap - memInfo.freeswap) * memInfo.mem_unit;
  if (print)
    printf("physical mem used %lld, swap used %lld\n", phyMemUsed, swapUsed);
}

void memUsedvsBaseline(const long long &phyMemUsed, const long long &swapUsed, long long &phyMemUsed2, long long &swapUsed2) {
  checkMemUsed(phyMemUsed2, swapUsed2, false);

  printf("physical mem used new %lld, swap used new %lld\n", phyMemUsed2 - phyMemUsed, swapUsed2 - swapUsed);
}


int main(int argc, char** argv) {
  srand(1);

  long long phyMemUsed, swapUsed;
  long long phyMemUsed2, swapUsed2;



  int size = 1000000;

  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElement copy.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType1 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      test.push_back(kmer);
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);

  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElement with move.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType1 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("Shared KmerIndexElement with move.  element size %lu\n", sizeof(KmerIndexType1));
    std::vector<KmerIndexType1 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    KmerIndexType1 kmer;
    for (int i = 0; i < size; ++i) {
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);





  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElementWithId.  element size %lu\n", sizeof(KmerIndexType2));
    std::vector<KmerIndexType2 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType2 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      kmer.id.file_pos = rand() % std::numeric_limits<uint64_t>::max();

      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElementWithIdAndQuality.  element size %lu\n", sizeof(KmerIndexType3));
    std::vector<KmerIndexType3 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType3 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      kmer.id.file_pos = rand() % std::numeric_limits<uint64_t>::max();
      kmer.qual = float(rand()) / float(RAND_MAX);

      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  printf("\n");
  checkMemUsed(phyMemUsed, swapUsed, true);
  {
    printf("KmerIndexElementWithIdAndQuality.  element size %lu\n", sizeof(KmerIndexType4));
    std::vector<KmerIndexType4 > test;
    printf("\tcreated vector\n");

    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

    for (int i = 0; i < size; ++i) {
      KmerIndexType4 kmer;
//      kmer.kmer = rand() % std::numeric_limits<uint64_t>::max();
      kmer.kmer.nextFromChar(rand() % (std::numeric_limits<uint8_t>::max() + 1));
      kmer.id.file_pos = rand() % std::numeric_limits<uint64_t>::max();
      kmer.qual = float(rand()) / float(RAND_MAX);

      test.push_back(std::move(kmer));
    }
    printf("\tpopulated vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);


    test.clear();
    printf("\tcleared vector %lu\n", test.size());
    memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);

  }
  printf("\tdeallocated vector\n");
  memUsedvsBaseline(phyMemUsed, swapUsed, phyMemUsed2, swapUsed2);
  sleep(2);


  return 0;
}
