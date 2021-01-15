#ifndef __STATMODELS_H_INCLUDED__ 
#define __STATMODELS_H_INCLUDED__

#include <vector>
#include <seqan/sequence.h>
#include <seqan/alignment_free.h>
#include "kernel_gen.h"

class StatModels
{
public:
  seqan::StringSet<seqan::DnaString> training_set;
  int lValue;
  double kmerNumber;
  int seqLen;
  KernelGen kernel;

  StatModels();

  StatModels(seqan::StringSet<seqan::DnaString> & sequence_set, int l_value, int seq_len, KernelGen & kernel_obj);

  void
  getPairwiseBackgrounds(std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds,
                         int background_type);

  void 
  getExpectedCounts(std::vector< std::vector< seqan::String<double> > > & expectedCounts,
                    std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds);

  void
  getSimpleStdDevEntropies(std::vector< std::vector< seqan::String<double> > > & standDevEntropies,
                           std::vector< std::vector< seqan::String<double> > > & expectedCounts);

  void
  expected_entropies(seqan::String<double> & expected_entropies, seqan::MarkovModel<seqan::Dna> & pairwiseBackgrounds);

  void
  getExpectedEntropies(std::vector< std::vector< seqan::String<double> > > & expectedEntropies,
  					   std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds);

  void
  getStdDevEntropies(std::vector< std::vector< seqan::String<double> > > & standDevEntropies,
  					 std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds);
  void
  calculateGCovariance(double & covariance, seqan::DnaString const & w1, seqan::DnaString const & w2, seqan::MarkovModel<seqan::Dna> & bgModel);

  void
  getStdDevEntropies(std::vector< seqan::String<double> > & standDevEntropies,
                     std::vector< seqan::MarkovModel<seqan::Dna> > & bgModels);

   void
  entropyStdDev(seqan::String<double> & std_devs, seqan::MarkovModel<seqan::Dna> & bgModel);

  void
  getComplexStdDevEntropies(std::vector< std::vector< seqan::String<double> > > & standDevEntropies,
                     std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds);
private:
  double
  expected_count(int const hash, double const p_w);

  void
  expected_counts(seqan::String<double> & expected_counts, seqan::MarkovModel<seqan::Dna> bgModel);

  double
  suffixesProbabilities(int const hash, seqan::MarkovModel<seqan::Dna> & pairwiseBackgrounds);

  double
  expected_entropy(int const hash, seqan::MarkovModel<seqan::Dna> & pairwiseBackgrounds);

  double
  term1Cov(seqan::DnaString const & w1, seqan::DnaString const & w2, 
           seqan::MarkovModel<seqan::Dna> & bgModel,
           int const seqLen, double const p_w1, double const p_w2);

  double
  term23Cov(seqan::DnaString const & w1, seqan::DnaString const & w2, 
            seqan::MarkovModel<seqan::Dna> & bgModel,
            int const seqLen, int const d, double const p_w1, double const p_w2);

  double
  term4Cov(seqan::DnaString const & w1, seqan::DnaString const & w2,
           seqan::MarkovModel<seqan::Dna> & bgModel,
           int const seqLen, int const d, double const p_w1, double const p_w2);

  double
  term1234Cov(seqan::DnaString const & w1, seqan::DnaString const & w2, // l1 > l2
              seqan::MarkovModel<seqan::Dna> & bgModel,
              int const seqLen);

};
#endif
