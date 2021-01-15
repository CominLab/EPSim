#include <seqan/sequence.h>
#include <seqan/alignment_free.h>

#include "ep_statistics.h"

EPStatistics::EPStatistics() {}

EPStatistics::EPStatistics(int myBackgroundModel,
                           std::vector<double> & myKernel, int myL, int mySeqLen)
{
  backgroundModel = myBackgroundModel;
  kernel = myKernel;
  seqLen = mySeqLen;
  kmerNumber = pow(4, myL); // kernel.size();
  kmerProbabilities.resize(myL);
  for (int i = 0; i < kmerProbabilities.size(); i++)
    resize(kmerProbabilities[i], kmerNumber, 0);
}

/******************** Printing methods ********************/

void 
EPStatistics::print_bernoulli_model(seqan::String<double> const & nucleotideFrequencies)
{
  for(int i = 0; i < 4; ++i)
    std::cout << nucleotideFrequencies[i] << " ";
}

void 
EPStatistics::print_pairwise_bernoulli_model(std::vector< std::vector< seqan::String<double> > > const & pairwiseBackgroundFrequencies)
{
  for (int i = 0; i < pairwiseBackgroundFrequencies.size(); i++)
  {
    for (int j = 0; j < pairwiseBackgroundFrequencies[i].size(); j++)
    {
      print_bernoulli_model(pairwiseBackgroundFrequencies[i][j]);  // p(A) p(C) p(G) p(T)
      std::cout << "\t";
    }
    std::cout << std::endl;
  }
}

void
EPStatistics::printExpectedEntropies(seqan::String<double> const & expected_entropies)
{
  for (int i = 0; i < length(expected_entropies); i++)
  {
    std::cout << expected_entropies[i] << ":";
  }
  std::cout << std::endl;
}

/******************** Pairwise bernoulli model ********************/

void 
EPStatistics::estimate_bernoulli_model(seqan::DnaString const & seq,
	                                     seqan::String<double> & nucleotideFrequencies)
{
  seqan::String<unsigned> kmerCounts;
  countKmers(kmerCounts, nucleotideFrequencies, seq, 1);
}

void 
EPStatistics::estimate_pairwise_bernoulli_model(std::vector< std::vector< seqan::String<double> > > & pairwiseBackgroundFrequencies, 
                                                std::vector< seqan::String<double> > const & backgroundFrequencies)
{
  for (int i = 0; i < pairwiseBackgroundFrequencies.size(); i++)
  {
    for (int j = i; j < pairwiseBackgroundFrequencies[i].size(); j++)
    {
      seqan::String<double> model;
      resize(model, 4);
      for (int k = 0; k < 4; k++)
      {
        model[k] = (backgroundFrequencies[i][k] + backgroundFrequencies[j][k]) / 2.0;  // p(A) p(C) p(G) p(T)
      }
      pairwiseBackgroundFrequencies[i][j] = model;
      pairwiseBackgroundFrequencies[j][i] = model;
    }	
  }
}

/******************** Statistical properties of counts ********************/

double
EPStatistics::expected_count(int const hash)
{
  int lValue = kmerProbabilities.size();
  double weight = seqLen - lValue + 1; // = numero totale di occorrenze, metodo di seqan 
  double expected_count = weight * kmerProbabilities[lValue - 1][hash];
  return expected_count;
}

void
EPStatistics::expected_counts(seqan::String<double> & expected_counts)
{
  resize(expected_counts, kmerNumber);
  for (int i = 0; i < kmerNumber; i++)
  {
    expected_counts[i] = expected_count(i);
  }
}

/******************** Statistical properties of entropic profiles ********************/

void
EPStatistics::emitteted_probabilities()
{
  for (int i = 0; i < kmerProbabilities.size(); i++)
  {
    for (int j = 0; j < length(kmerProbabilities[i]); j++)
    {
      seqan::String<seqan::Dna> w;
      unhash(w, j, i + 1); // esult, hash, q
      kmerProbabilities[i][j] = backgroundModel->emittedProbability(w);
    }
  }
}

double
EPStatistics::suffixesProbabilities(int const hash) 
{
  int lValue = kmerProbabilities.size();
  double suffixesEntropy = 0.0;
  for (int j = 0; j < lValue; j++)
  {
    int subterm_index = hash % length(kmerProbabilities[j]);
    double weight = kernel[j] * (seqLen - j); // seqLen - j - 1 + 1
    suffixesEntropy += weight * kmerProbabilities[j][subterm_index];
  }
  return suffixesEntropy;
}

double
EPStatistics::expected_entropy(int const hash)
{
  double exp_entropy = suffixesProbabilities(hash);
  //exp_entropy += prefixesProbabilities(hash);
  return exp_entropy;
}

/* 
  Non è semplice come sembra, devo considerare:
  la probabilità della stringa
  la probabilità delle sottostringhe -> quali sottostringhe considerare dipende dal blur che ho applicato
  la probabilità del reverse complement
  la probabilità delle sottostringhe del reverse complement
*/
void
EPStatistics::expected_entropies(seqan::String<double> & expected_entropies)
{
  resize(expected_entropies, kmerNumber);
  for (int i = 0; i < kmerNumber; i++)
  {
    expected_entropies[i] = expected_entropy(i);
    
    //expected_entropies[i] += expected_entropy(kmerNumber - 1 - i); // rc
    //expected_entropies[i] = expected_entropies[i] / 2.0; // il diviso 2.0 non mi piace, non funziona per la varianza
  }
}

void
EPStatistics::expected_variances(seqan::String<double> & variances)
{
  for (int i = 0; i < kmerNumber; i++)
  {
    // prima di implementarla è meglio verificare la teoria
  }
}

/******************** Experimental methods ********************/

double
EPStatistics::prefixesProbabilities(int const hash) // experimental
{
  int lValue = kmerProbabilities.size();
  double prefixesEntropy = 0.0;
  for (int j = 0; j < lValue; j++)
  {
    int subterm_index = hash / (length(kmerProbabilities[lValue - 1]) / length(kmerProbabilities[j]));
    double weight = kernel[j] * (seqLen - j); // seqLen - j - 1 + 1
    prefixesEntropy += weight * kmerProbabilities[j][subterm_index];
  }
  return prefixesEntropy;
}