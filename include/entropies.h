#ifndef __ENTROPIES_H_INCLUDED__ 
#define __ENTROPIES_H_INCLUDED__

#include <seqan/sequence.h>
#include <seqan/index.h>

#include "kernel_gen.h"

class EntropicProfiler
{
public:
  seqan::DnaString sequence;
  int lValue;
  int lowest_k;
  std::vector< seqan::String<unsigned> > kmerCounts;
  std::vector< seqan::Shape<seqan::Dna, seqan::SimpleShape> > shapes;
  KernelGen kernel_obj;

  EntropicProfiler();

  EntropicProfiler(seqan::DnaString const & seq,
                   KernelGen & kernel_obj,
  				         int const l, int const min_k = 1);

  int 
  checkAllPrefixes(seqan::Iterator<seqan::DnaString const, seqan::Rooted>::Type & start_it, 
                   int const min_length, int const max_length);

  int 
  countKmers();

  const seqan::String<unsigned> & 
  getCounts();

  void 
  printCounts();

  void 
  blurWithSuffixes(seqan::String<double> & entropies);

  void 
  blurWithPrefixes(seqan::String<double> & entropies);

  void 
  limbsBlur(seqan::String<double> & entropies);
  
  int
  blurWithReverseComplement(seqan::String<double> const & entropies_without_rc, 
                            seqan::String<double> & entropies_with_rc);
  void
  normalizeEntropies(seqan::String<double> & entropies);
  
  void 
  blurWithMidstrings(seqan::String<double> & entropies);
  
  //void
  //fullBlur(seqan::String<double> & entropies, double const phi);

  void 
  printEntropies(seqan::String<double> const & entropies);

  /******************** Arithmetic utilities ********************/

  double
  getMaxEntropy(seqan::String<double> const & entropies);

  double
  getAvgEntropy(seqan::String<double> const & entropies);

  double
  getStandDevEntropy(seqan::String<double> const & entropies, double const avgEntropy);

private:
  
};

#endif
