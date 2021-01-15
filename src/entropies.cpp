#include "entropies.h"

EntropicProfiler::EntropicProfiler() {}

EntropicProfiler::EntropicProfiler(seqan::DnaString const & seq,
                                   KernelGen & kernel_object, 
  				  				               int const l, int const min_k)
{
  // TODO: throw exception if min_k <= 0

  sequence = seq;
  lValue = l;
  lowest_k = min_k;
	
  kmerCounts.resize(lValue);
  shapes.resize(lValue);

  for (int i = 0; i < lValue; i++)
  {
  	int kmerNumber = pow(4, i + 1);
  	resize(kmerCounts[i], kmerNumber, 0);
  	resize(shapes[i], i + 1);
  }

  kernel_obj = kernel_object;
}

/******************** Printing methods ********************/

void 
EntropicProfiler::printCounts()
{
  for (int i = 0; i < kmerCounts.size(); i++)
  {
    for (int j = 0; j < length(kmerCounts[i]); j++)
    {
      std::cout << kmerCounts[i][j] << ":";
    }
    std::cout << std::endl;
  }
}

void 
EntropicProfiler::printEntropies(seqan::String<double> const & entropies)
{
  for (int i = 0; i < length(entropies); i++)
  {
    std::cout << entropies[i] << ":";
  }
  std::cout << std::endl;
}

/******************** Methods for the counts ********************/

int 
EntropicProfiler::checkAllPrefixes(seqan::Iterator<seqan::DnaString const, seqan::Rooted>::Type & start_it, 
								                   int const min_length, int const max_length)
{
  if (min_length <= 0)
    return -1;

  for (int i = min_length - 1; i < max_length; i++)
  {
    unsigned hashValue = hash(shapes[i], start_it);
    ++kmerCounts[i][hashValue];
  }
}

int 
EntropicProfiler::countKmers()
{
  seqan::Iterator<seqan::DnaString const, seqan::Rooted>::Type itSeq = begin(sequence);
  for (; itSeq <= end(sequence) - lValue; ++itSeq)
  {
    if (checkAllPrefixes(itSeq, lowest_k, lValue) == -1)
      return -1;
  }

  for (int i = 1; i < lValue; i++)
  {
  	if (checkAllPrefixes(itSeq, lowest_k, lValue - i) == -1)
      return -1;
  	++itSeq;
  }
}

const seqan::String<unsigned> &
EntropicProfiler::getCounts()
{
  return kmerCounts[lValue - 1];
}

/******************** Methods for the entropic profiles ********************/

void 
EntropicProfiler::blurWithSuffixes(seqan::String<double> & entropies) // m sperimentale
{
  clear(entropies);
  resize(entropies, pow(4, lValue), 0.0);

  for (int i = 0; i < length(entropies); i++)
  {
    for (int j = lowest_k - 1; j < lValue; j++)
    {
      int subterm_index = i % length(kmerCounts[j]);
      double weight = kernel_obj.kernel[j];
      entropies[i] += weight * kmerCounts[j][subterm_index];
      //std::cout << "j " << j << " subterm_index " << subterm_index << std::endl;
    }
  } 
}

/*
  non mi piace il / 2, probabilmente intriga nel calcolo della varianza
  le entropie finali sono simmetriche -> permette qualche ottimizzazione
*/
int 
EntropicProfiler::blurWithReverseComplement(seqan::String<double> const & entropies_without_rc, 
                                            seqan::String<double> & entropies_with_rc)
{
  clear(entropies_with_rc);
  resize(entropies_with_rc, pow(4, lValue), 0.0);

  if (length(entropies_without_rc) != length(entropies_with_rc))
  {
    std::cerr << "ERROR: cannot be length(entropies_without_rc) != length(entropies_with_rc)" << std::endl;
    return -1;
  }

  for (int i = 0; i < length(entropies_with_rc); i++)
    entropies_with_rc[i] = (entropies_without_rc[i] + entropies_without_rc[length(entropies_without_rc) - 1 - i]) / 2;
}

/*
  come cambia quando considero il reverse complement?
*/
void
EntropicProfiler::normalizeEntropies(seqan::String<double> & entropies)
{
  for (int i = 0; i < length(entropies); i++)
  {
    entropies[i] = entropies[i] / kernel_obj.norm_term;
  }
}

/******************** Arithmetic utilities ********************/

double
EntropicProfiler::getMaxEntropy(seqan::String<double> const & entropies) // or kmerCounts[lValue - 1]
{
  double maxEntropy = 0.0;
  for (int i = 0; i < length(entropies); i++)
  {
    if (entropies[i] > maxEntropy)
      maxEntropy = entropies[i];
  }
  std::cout << "maxEntropy " << maxEntropy << std::endl;
  return maxEntropy;
}

double
EntropicProfiler::getAvgEntropy(seqan::String<double> const & entropies) // or kmerCounts[lValue - 1]
{
  double avgEntropy = 0.0;
  for (int i = 0; i < length(entropies); i++)
    avgEntropy += entropies[i];
  avgEntropy = avgEntropy / length(entropies);
  std::cout << "avgEntropy " << avgEntropy << std::endl;
  return avgEntropy;
}

double
EntropicProfiler::getStandDevEntropy(seqan::String<double> const & entropies, double const avgEntropy) // or kmerCounts[lValue - 1]
{
  double standDevEntropy = 0.0;
  for (int i = 0; i < length(entropies); i++)
    standDevEntropy += pow(entropies[i] - avgEntropy, 2.0);
  standDevEntropy = sqrt(standDevEntropy / (length(entropies) - 1));
  std::cout << "standDevEntropy " << standDevEntropy << std::endl;
  return standDevEntropy;
}

/******************** Experimental methods ********************/

void
EntropicProfiler::blurWithPrefixes(seqan::String<double> & entropies)
{
  clear(entropies);
  resize(entropies, pow(4, lValue), 0.0);

  for (int i = 0; i < length(entropies); i++)
  {
    for (int j = lowest_k - 1; j < lValue; j++)
    {
      int subterm_index = i / (length(kmerCounts[lValue - 1]) / length(kmerCounts[j]));
      double weight = kernel_obj.kernel[j];
      entropies[i] += weight * kmerCounts[j][subterm_index];
    }
  }
}

void
EntropicProfiler::limbsBlur(seqan::String<double> & entropies) // m sperimentale
{
  clear(entropies);
  resize(entropies, pow(4, lValue), 0.0);

  seqan::String<double> prefix_entropies;
  blurWithPrefixes(prefix_entropies);

  seqan::String<double> suffix_entropies;
  blurWithSuffixes(suffix_entropies);

  for (int i = 0; i < length(entropies); i++)
  {
  	entropies[i] = (prefix_entropies[i] + suffix_entropies[i]) / 2.0;
  }
}

void 
EntropicProfiler::blurWithMidstrings(seqan::String<double> & entropies) // m sperimentale
{
  clear(entropies);
  resize(entropies, pow(4, lValue), 0.0);

  for (int i = 0; i < length(entropies); i++)
  {
  	for (int j = lValue - 1; j >= 0; j = j - 2)
  	{
  	  int subterm_index = i / (length(kmerCounts[lValue - 1]) / length(kmerCounts[j]));
  	  double weight = kernel_obj.kernel[j];
  	  entropies[i] += weight * kmerCounts[j][subterm_index];
  	}
  }	
}

//void
//fullBlur(seqan::String<double> & entropies, double const phi) {}

/*

double
EntropicProfiler::getArithmeticAverageEntropy()
{
	return ;
}

double
EntropicProfiler::getStandardDeviationEntropy()
{
	return ;
}

double
EntropicProfiler::getSuffixCountsSumList()
{
	return ;
}

*/