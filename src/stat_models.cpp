#include "stat_models.h"

#include <seqan/sequence.h>
#include <seqan/alignment_free.h>

StatModels::StatModels() {}

StatModels::StatModels(seqan::StringSet<seqan::DnaString> & sequence_set, int l_value, int seq_len, KernelGen & kernel_obj)
{
  training_set = sequence_set;
  lValue = l_value;
  kmerNumber = pow(4, lValue);
  seqLen = seq_len;
  kernel = kernel_obj;
}

void
StatModels::getPairwiseBackgrounds(std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds,
                                   int background_type)
{
  pairwiseBackgrounds.resize(length(training_set));
  seqan::Iterator< seqan::StringSet<seqan::DnaString> >::Type it1;
  for (it1 = begin(training_set); it1 != end(training_set); ++it1) 
  {
    seqan::Iterator< seqan::StringSet<seqan::DnaString> >::Type it2;
    for (it2 = begin(training_set); it2 != end(training_set); ++it2) 
    {
      seqan::DnaString seq1seq2;
      append(seq1seq2, value(it1));
      append(seq1seq2, value(it2));
      seqan::StringSet<seqan::DnaString> bgSequences;
      stringToStringSet(bgSequences, seq1seq2);

      seqan::MarkovModel<seqan::Dna> backgroundModel(background_type);
      buildMarkovModel(backgroundModel, bgSequences);
      
      pairwiseBackgrounds[it1 - begin(training_set)].push_back(backgroundModel);
    }
  }
}

double
StatModels::expected_count(int const hash, double const p_w)
{
  // in alternativa seqLen - lValue + 1 = numero totale di occorrenze, metodo di seqan 
  double expected_count = (seqLen - lValue + 1) * p_w;
  return expected_count;
}

void
StatModels::expected_counts(seqan::String<double> & expected_counts, seqan::MarkovModel<seqan::Dna> bgModel)
{
  resize(expected_counts, kmerNumber);
  for (int i = 0; i < kmerNumber; i++)
  {
    seqan::String<seqan::Dna> w;
    unhash(w, i, lValue); // result, hash, q COLLO DI BOTTIGLIA
    double p_w = bgModel.emittedProbability(w);
    expected_counts[i] = expected_count(i, p_w);
  }
}

void
StatModels::getExpectedCounts(std::vector< std::vector< seqan::String<double> > > & expectedCounts,
                              std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds)
{
  for (int i = 0; i < pairwiseBackgrounds.size(); i++)
  {
    expectedCounts.resize(length(training_set));
    for (int j = i; j < pairwiseBackgrounds[i].size(); j++) // 0 al posto di i
    {
      expectedCounts[i].resize(length(training_set));
      expected_counts(expectedCounts[i][j], pairwiseBackgrounds[i][j]);
      
      expectedCounts[j].resize(length(training_set));
      expectedCounts[j][i] = expectedCounts[i][j];
    }
  }
}

/*
  Approssimazione valida approssimando la binomiale con poisson (non vero per gli EP perchè covarianza non trascurabile),
  il che implica media = varianza,
  dev_stand = sqrt(varianza) = sqrt(media)
*/
void
StatModels::getSimpleStdDevEntropies(std::vector< std::vector< seqan::String<double> > > & standDevEntropies,
                                     std::vector< std::vector< seqan::String<double> > > & expectedCounts)
{
  for (int i = 0; i < expectedCounts.size(); i++)
  {
    standDevEntropies.resize(length(training_set));
    for (int j = 0; j < expectedCounts[i].size(); j++)
    {
      standDevEntropies[i].resize(length(training_set));
      for (int k = 0; k < length(expectedCounts[i][j]); k++)
      {
        resize(standDevEntropies[i][j], kmerNumber);
        standDevEntropies[i][j][k] = sqrt(expectedCounts[i][j][k]);
      }
    }
  }
}

double
StatModels::suffixesProbabilities(int const hash,
                                  seqan::MarkovModel<seqan::Dna> & pairwiseBackgrounds)
{
  double suffixesEntropy = 0.0;
  for (int j = 0; j < lValue; j++) // G TG TGG
  {
  	//std::cout << "loop suffixesProbabilities" << std::endl;

    double weight = kernel.kernel[j] * (seqLen - j); // seqLen - j - 1 + 1
                                                     // perché j = 0 significa k = 1
                                                     // peso * numero medio di occorrenze

    if (weight < 0.1) // per non fare calcoli inutili, occhio che anche le probabilità sono basse!
      continue;

    seqan::String<seqan::Dna> w;
    unhash(w, hash, j + 1); // result, hash, q
    double p_w = pairwiseBackgrounds.emittedProbability(w);

    //std::cout << "w p_w " << w << " " << p_w << std::endl;

    //std::cout << "kernel[j] (seqLen - j)" << kernel.kernel[j] << " " << seqLen - j << std::endl;
    //std::cout << "weight " << weight << std::endl;
    //std::cout << "che cosa scartiamo? " << weight * p_w << std::endl;

    suffixesEntropy += weight * p_w;
  }

  return suffixesEntropy;
}

double
StatModels::expected_entropy(int const hash,
                             seqan::MarkovModel<seqan::Dna> & pairwiseBackgrounds)
{
  double exp_entropy = suffixesProbabilities(hash, pairwiseBackgrounds) / kernel.norm_term;
  //std::cout << "suffixesEntropy " << exp_entropy << std::endl;
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
StatModels::expected_entropies(seqan::String<double> & expected_entropies,
                               seqan::MarkovModel<seqan::Dna> & pairwiseBackgrounds)
{
  resize(expected_entropies, kmerNumber);
  for (int i = 0; i < kmerNumber; i++)
  {
  	//std::cout << "pre_expected_entropies[i]" << std::endl;
    expected_entropies[i] = expected_entropy(i, pairwiseBackgrounds);
    
    //expected_entropies[i] += expected_entropy(kmerNumber - 1 - i); // rc
    //expected_entropies[i] = expected_entropies[i] / 2.0; // il diviso 2.0 non mi piace, non funziona per la varianza
  }
}

// da spezzare in due
void
StatModels::getExpectedEntropies(std::vector< std::vector< seqan::String<double> > > & expectedEntropies,
                                 std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds)
{
  for (int i = 0; i < pairwiseBackgrounds.size(); i++)
  {
    expectedEntropies.resize(length(training_set));
    for (int j = i; j < pairwiseBackgrounds[i].size(); j++)
    {
      expectedEntropies[i].resize(length(training_set));
      expected_entropies(expectedEntropies[i][j], pairwiseBackgrounds[i][j]);

      expectedEntropies[j].resize(length(training_set));
      expectedEntropies[j][i] = expectedEntropies[i][j];
    }
  }
}

double
StatModels::term1Cov(seqan::DnaString const & w1, seqan::DnaString const & w2, 
                     seqan::MarkovModel<seqan::Dna> & bgModel,
                     int const seqLen, double const p_w1, double const p_w2)
{
  double covariance = 0.0;

  int l1 = length(w1);
  int l2 = length(w2);

  int minl12 = std::min(l1,l2);
  double mult_term = seqLen - minl12 + 1;

  //double p_w1 = emittedProbability(bgModel, w1);
  //double p_w2 = emittedProbability(bgModel, w2);

  double cov_part2 = p_w1 * p_w2;
  covariance += mult_term * (p_w1 - cov_part2); // d = 0 primo termine

  return covariance;
}

double
StatModels::term23Cov(seqan::DnaString const & w1, seqan::DnaString const & w2, 
                      seqan::MarkovModel<seqan::Dna> & bgModel,
                      int const seqLen, int const d, double const p_w1, double const p_w2)
{
  double covariance = 0.0;

  int l1 = length(w1);
  int l2 = length(w2);

  int minl12 = std::min(l1,l2);
  double mult_term = seqLen - minl12 + 1;

  //double p_w1 = emittedProbability(bgModel, w1);
  //double p_w2 = emittedProbability(bgModel, w2);

  if (d <= l1 - l2) // in pratica non serve, sarebbe gestito dal clump
  {
    //std::cout << "no clump d " << d << std::endl;
    covariance += (mult_term - d) * p_w1; // secondo termine
  }
  else
  {
    // w1 ACGA
    // w2 --GAC
    // prefisso della più lunga + la più corta
    seqan::DnaString clump = prefix(w1, d);
    append(clump, w2);
    //std::cout << "clump " << clump << " d " << d << std::endl;
    double p_clump = emittedProbability(bgModel, clump);
    covariance += (mult_term - d) * p_clump; // terzo termine
  }

  double cov_part2 = p_w1 * p_w2;
  covariance = covariance - cov_part2;

  return covariance;
}

double
StatModels::term4Cov(seqan::DnaString const & w1, seqan::DnaString const & w2,
                     seqan::MarkovModel<seqan::Dna> & bgModel,
                     int const seqLen, int const d, double const p_w1, double const p_w2)
{
  double covariance = 0.0;

  int l1 = length(w1);
  int l2 = length(w2);

  int minl12 = std::min(l1,l2);
  double mult_term = seqLen - minl12 + 1;

  //double p_w1 = emittedProbability(bgModel, w1);
  //double p_w2 = emittedProbability(bgModel, w2);

  // w1 -ACGA
  // w2 GAC
  // prefisso della più corta (G) + la più lunga (ACGA)
  seqan::DnaString clump = prefix(w2, d);
  append(clump, w1);
  //std::cout << "clump " << clump << " d " << d << std::endl;
  double p_clump = emittedProbability(bgModel, clump);
  covariance += (mult_term - d) * p_clump; // quarto termine (value(d) + 1)

  double cov_part2 = p_w1 * p_w2;
  covariance = covariance - cov_part2;

  return covariance;
}

double
StatModels::term1234Cov(seqan::DnaString const & w1, seqan::DnaString const & w2, // l1 > l2
                        seqan::MarkovModel<seqan::Dna> & bgModel,
                        int const seqLen)
{
  double covariance = 0.0;

  // da fare una volta sola perché sono operazioni costose
  double p_w1 = emittedProbability(bgModel, w1);
  double p_w2 = emittedProbability(bgModel, w2);

  covariance += term1Cov(w1, w2, bgModel, seqLen, p_w1, p_w2);

  seqan::String<int> periodicity1;
  calculatePeriodicity(periodicity1, w1, w2);
  seqan::Iterator<seqan::String<int>, seqan::Rooted>::Type d;
  for (d = begin(periodicity1); d < end(periodicity1); ++d)
  {
    covariance += term23Cov(w1, w2, bgModel, seqLen, value(d), p_w1, p_w2);
  }

  seqan::String<int> periodicity2;
  calculatePeriodicity(periodicity2, w2, w1);
  for (d = begin(periodicity2); d < end(periodicity2); ++d)
  {
    covariance += term4Cov(w1, w2, bgModel, seqLen, value(d), p_w1, p_w2);
  }

  return covariance;
}

void
StatModels::calculateGCovariance(double & covariance, seqan::DnaString const & w1, seqan::DnaString const & w2, seqan::MarkovModel<seqan::Dna> & bgModel)
{
  covariance = 0.0;

  if (w1 == w2)
  {
    calculateVariance(covariance, w1, bgModel, seqLen);
    return;
  }

  int l1 = length(w1);
  int l2 = length(w2);

  if (l1 == l2)
  {
  	calculateCovariance(covariance, w1, w2, bgModel, seqLen);
  	return;
  }

  //std::cout << "new formula" << std::endl;

  if (l1 > l2)
    covariance = term1234Cov(w1, w2, bgModel, seqLen);
  else
    covariance = term1234Cov(w2, w1, bgModel, seqLen);

  //std::cout << "covariance " << covariance << std::endl;
}

void
StatModels::entropyStdDev(seqan::String<double> & std_devs, seqan::MarkovModel<seqan::Dna> & bgModel)
{
  resize(std_devs, kmerNumber);
  for (int i = 0; i < kmerNumber; i++) // itero su ogni parola
  {
  	std_devs[i] = 0.0;
  	for (int j = 0; j < lValue; j++)
  	{
      if (kernel.kernel[j] < 0.1) // per evitare calcoli inutili
        continue;

  	  seqan::String<seqan::Dna> w1;
  	  unhash(w1, i, j + 1); // result, hash, q
  	  //std::cout << "w1 " << w1 << std::endl;

  	  for (int k = 0; k < lValue; k++)
  	  {
        if (kernel.kernel[k] < 0.1) // per evitare calcoli inutili
          continue;

  	    seqan::String<seqan::Dna> w2;
  	    unhash(w2, i, k + 1); // result, hash, q
  	    //std::cout << "w2 " << w2 << std::endl;

        double covariance = 0.0;
        calculateGCovariance(covariance, w1, w2, bgModel);
        //std::cout << "covariance " << covariance << std::endl;

        double weight = kernel.kernel[j] * kernel.kernel[k];
        //std::cout << "weight " << kernel.kernel[j] << " " << kernel.kernel[k] << std::endl;
  	    
        std_devs[i] += weight * covariance;
  	    //std::cout << "variance not normalized" << std_devs[i] << std::endl;
      }  
  	}
	  
    // normalize and calculate stand dev = sqrt(variance)
    std_devs[i] = sqrt(std_devs[i] / pow(kernel.norm_term, 2.0));

   //return;
  }
}

void
StatModels::getStdDevEntropies(std::vector< seqan::String<double> > & standDevEntropies,
                               std::vector< seqan::MarkovModel<seqan::Dna> > & bgModels)
{
  standDevEntropies.resize(length(training_set));
  for (int i = 0; i < bgModels.size(); i++)
  {
    entropyStdDev(standDevEntropies[i], bgModels[i]);
  }
}

// lentissima
void
StatModels::getComplexStdDevEntropies(std::vector< std::vector< seqan::String<double> > > & standDevEntropies,
                               std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > & pairwiseBackgrounds)
{
  for (int i = 0; i < pairwiseBackgrounds.size(); i++)
  {
    standDevEntropies.resize(length(training_set));
    for (int j = i; j < pairwiseBackgrounds[i].size(); j++)
    {
      standDevEntropies[i].resize(length(training_set));
      entropyStdDev(standDevEntropies[i][j], pairwiseBackgrounds[i][j]);

      standDevEntropies[j].resize(length(training_set));
      standDevEntropies[j][i] = standDevEntropies[i][j];
    }
  }
}