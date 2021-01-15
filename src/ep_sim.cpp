/*
 * ep_sim.cpp
 *
 * Testando il programma con i conteggi
 * D2 e D2Star con m0 sono identici ad ALF
 * D2Star con m1 varia di poco, pare un errore di approssimazione
 *
 * Created by Morris Antonello on 13/03/14
 */

#include <string>
#include <iostream>
#include <map>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/alignment_free.h>
#include <seqan/basic.h>

#include "entropies.h"
#include "stat_models.h"
#include "score_filler.h"
#include "kernel_gen.h"

/* Encapsulate user settings */
struct ep_sim_options
{
  seqan::CharString sequence_filename;
  seqan::CharString output_filename;
  int l_value;
  double phi_value;
  int background_type;
  bool revComplFlag;
  std::string measure_type;
  bool debugWithCounts;
  bool approximateStandDev;
  bool fullPairwiseStandDev;
  ep_sim_options()
  {}
};

/*
 * Parse the command line
 *
 * Params:
 *  ep_sim_options & options: struct which stores user settings
 *  int argc
 *  char const ** argv
 */
seqan::ArgumentParser::ParseResult
parseCommandLine(ep_sim_options & options, int argc, char const ** argv)
{
  // Setup ArgumentParser
  seqan::ArgumentParser parser("ep_sim");
  
  // Set short description, version, and date
  setShortDescription(parser, "ep_sim");
  setVersion(parser, "1.0");
  setDate(parser, "April 2014");

  // Define usage line and long description
  addUsageLine(parser,
             "\\fIsequence_filename\\fP \\fIoutput_filename\\fP [\\fIOPTIONS\\fP]");
  addDescription(parser,
               "./ep_sim ../dataset/mod4-RANDOM3-len100-p1-m1.fasta modep1_random_sep.txt -L 5 -P 0.85 -B 1 -D D2 ");

  // Two arguments required
  addArgument(parser, seqan::ArgParseArgument(
    seqan::ArgParseArgument::INPUTFILE, "sequence_filename"));
  addArgument(parser, seqan::ArgParseArgument(
    seqan::ArgParseArgument::OUTPUTFILE, "output_filename"));

  // Define options
  addOption(parser, seqan::ArgParseOption(
    "L", "l_value", "Set l_value.", 
    seqan::ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "l_value", "5");
  setMinValue(parser, "l_value", "1");

  addOption(parser, seqan::ArgParseOption(
    "P", "phi_value", "Set phi_value.", 
    seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "phi_value", "0.5");

  addOption(parser, seqan::ArgParseOption(
    "B", "background_type", "Set background_type.", 
    seqan::ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "background_type", "1");

  addOption(parser, seqan::ArgParseOption(
    "D", "measure_type", "Select measure_type.", 
    seqan::ArgParseArgument::STRING, "STRING"));
  setDefaultValue(parser, "measure_type", "D2");

  addOption(parser, seqan::ArgParseOption(
    "R", "reverse_complement", "Take into account reverse_complement."));
  
  addOption(parser, seqan::ArgParseOption(
    "C", "debug_with_counts", "Debug with counts."));

  addOption(parser, seqan::ArgParseOption(
    "A", "approximate_stand_dev", "Approximate standard deviation."));

  addOption(parser, seqan::ArgParseOption(
    "F", "full_pair_stand_dev", "Full pairwise standard deviation."));

  // Parse command line
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // Only extract  options if the program will continue after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

  // Extract arguments values
  getArgumentValue(options.sequence_filename, parser, 0);
  getArgumentValue(options.output_filename, parser, 1);

  // Extract option values
  getOptionValue(options.l_value, parser, "l_value");
  getOptionValue(options.phi_value, parser, "phi_value");
  getOptionValue(options.background_type, parser, "background_type");
  getOptionValue(options.measure_type, parser, "measure_type");

  options.revComplFlag = isSet(parser, "reverse_complement");
  options.debugWithCounts = isSet(parser, "debug_with_counts");
  options.approximateStandDev = isSet(parser, "approximate_stand_dev");
  options.fullPairwiseStandDev = isSet(parser, "full_pair_stand_dev");
    
  return seqan::ArgumentParser::PARSE_OK;
}

/*
 * Print out user settings
 *
 * Params:
 *  ep_sim_options & options: struct which stores user settings
 */
void
print_options(ep_sim_options const & options)
{
  // Print out options
  std::cout << "L " << options.l_value << std::endl;
  std::cout << "phi " << options.phi_value << std::endl;
  std::cout << "background_type " << options.background_type << std::endl;
  std::cout << "revComplFlag " << options.revComplFlag << std::endl;
  std::cout << "debugWithCounts " << options.revComplFlag << std::endl;
  std::cout << "approximateStandDev " << options.approximateStandDev << std::endl;
  std::cout << "fullPairwiseStandDev " << options.fullPairwiseStandDev << std::endl;
}

void 
coutVector(std::vector<double> v) 
{
  std::cout << "[ ";
  for(int i=0; i<v.size(); i++) 
  {
	std::cout << v[i];
	if(i<v.size()-1)
	{
		std::cout << ",";
	}
	std::cout << " ";
  }
std::cout << "]" << std::endl;
}

void
coutString(seqan::String<double> const & entropies)
{
  for (int i = 0; i < length(entropies); i++)
  {
    std::cout << entropies[i] << ":";
  }
  std::cout << std::endl;
}

int 
main(int argc, char const ** argv)
{
  /* Parse the command line */

  ep_sim_options options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
  
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  print_options(options);

  bool printTiming = true;
  double startTime = seqan::sysTime();

  /* Read sequences file (FASTA/FASTQ) */

  seqan::StringSet<seqan::DnaString> training_set;
  int seq_len;

  seqan::SequenceStream seqStream(toCString(options.sequence_filename));
  if (!isGood(seqStream))
  {
    std::cerr << "ERROR: could not open file " << options.sequence_filename << std::endl;
    return -1;
  }

  while (!atEnd(seqStream))
  {

    seqan::CharString id; 
    seqan::DnaString seq;
    if (readRecord(id, seq, seqStream) != 0)
    {
      std::cerr << "ERROR: could not read from " << options.sequence_filename << std::endl;
      return -1;
    }
    std::cout << id << '\t' << seq << std::endl;
    seq_len = length(seq);

    appendValue(training_set, seq);
    std::cout << "number of sequences is " << length(training_set) << std::endl;
    std::cout << "total sequences length is " << lengthSum(training_set) << std::endl;
  }

  if (printTiming)
    std::cerr << "Time to read sequences file " << seqan::sysTime() - startTime << " s." << std::endl;
  startTime = seqan::sysTime();

  /* Compute entropies record by record */

  std::vector< seqan::String<double> > all_entropies(length(training_set));
  KernelGen kernel(options.phi_value / 10.0, options.l_value); // / 10.0 utile per gli script
  std::vector< seqan::MarkovModel<seqan::Dna> > singleBackgrounds; // for N2

  seqan::Iterator< seqan::StringSet<seqan::DnaString> >::Type it;
  for (it = begin(training_set); it != end(training_set); ++it) 
  {
    // TODO check sequence is not empty

  	EntropicProfiler ep(value(it), kernel, options.l_value);
    ep.countKmers();
    //ep.printCounts();

    if (options.debugWithCounts)
    {
      all_entropies[it - begin(training_set)] = ep.getCounts();
    }
    else
    {
      seqan::String<double> entropies;

      std::cout << "blurWithSuffixes" << std::endl;
      ep.blurWithSuffixes(entropies);
      ep.normalizeEntropies(entropies); // NORMALIZZAZIONE
      //ep.printEntropies(entropies);

      // Valori aritmetici calcolati sui valori normalizzati
   	  double maxEntropy = ep.getMaxEntropy(entropies);
   	  double avgEntropy = ep.getAvgEntropy(entropies);
   	  double standDevEntropy = ep.getStandDevEntropy(entropies, avgEntropy);

   	  if (options.measure_type.compare("AR") == 0) // arithmetic
  	  {
  	  	for (int i = 0; i < length(entropies); i++)
  	  	  entropies[i] =  (entropies[i] - avgEntropy) / standDevEntropy;

  	  }
  	  else if (options.measure_type.compare("ARM") == 0) // arithmetic max
  	  {
  	  	// Giustificazione per il massimo
  	  	// http://statistics.about.com/od/Descriptive-Statistics/a/Range-Rule-For-Standard-Deviation.htm
  	  	if (maxEntropy > pow(10.0, -10) && maxEntropy < pow(10.0, 10))
  	  	{
  	  	  for (int i = 0; i < length(entropies); i++)
  	  	    entropies[i] =  entropies[i] / maxEntropy;
  	  	}
  	  	else
  	  	{
  	  	  std::cerr << "Max cannot be 0 or infinity " << maxEntropy <<  std::endl;
  	  	  return -1;
  	  	}
      }
      else if (options.measure_type.compare("D2*") == 0
               && !options.debugWithCounts
               && !options.approximateStandDev
               && !options.fullPairwiseStandDev) // neighbourhood of suffixes
      {
        // La probabilità di A in Bernoulli e M1 per la stessa sequenza è diversa
        // perché in Bernoulli / seq_len mentre in M1 / (seq_len - 1) 

        std::cout << "Computing single background model" << std::endl;
        seqan::MarkovModel<seqan::Dna> model(options.background_type);
        
        // modo1 con stringset
        seqan::StringSet<seqan::DnaString> bgSeq;
        stringToStringSet(bgSeq, value(it));
        buildMarkovModel(model, bgSeq);

        // modo2 senza stringset, meno pulito per via dell'ultimo parametro
        //seqan::String<unsigned> kmerCounts;
        //countKmers(kmerCounts, model, value(it), 2);
        
        singleBackgrounds.push_back(model);

        /*std::cout << "stationaryDistribution " << model.stationaryDistribution << std::endl;
        std::cout << "model.transition " << model.transition << std::endl;

        double p_w = emittedProbability(model, "A");
        std::cout << "pr A " << p_w << std::endl;
        p_w = model.emittedProbability("C");
        std::cout << "pr C " << p_w << std::endl;
        p_w = model.emittedProbability("G");
        std::cout << "pr G " << p_w << std::endl;
        p_w = model.emittedProbability("T");
        std::cout << "pr T " << p_w << std::endl;*/
      }

	    all_entropies[it - begin(training_set)] = entropies; //entropies_with_rc;
    }

    //std::cout << "blurWithPrefixes" << std::endl; 
    //ep.blurWithPrefixes(entropies);
    //ep.normalizeEntropies(entropies);
    //ep.printEntropies(entropies);

    //std::cout << "limbsBlur" << std::endl;
    //ep.limbsBlur(entropies);
    //ep.normalizeEntropies(entropies);
    //ep.printEntropies(entropies);
    
    //std::cout << "blurWithReverseComplement" << std::endl;
    //seqan::String<double> entropies_with_rc;
    //ep.blurWithReverseComplement(entropies, entropies_with_rc);
    ////ep.normalizeEntropies(entropies_with_rc);
    //ep.printEntropies(entropies_with_rc);
  }

  if (printTiming)
    std::cerr << "Time to compute entropies and backgrounds record by record " << seqan::sysTime() - startTime << " s." << std::endl;
  startTime = seqan::sysTime();

  std::vector< std::vector< seqan::String<double> > > expectedEntropies;
  std::vector< std::vector< seqan::String<double> > > standDevEntropies;
  if (options.measure_type.compare("D2*") == 0)
  {
  	std::vector< std::vector< seqan::MarkovModel<seqan::Dna> > > pairwiseBackgrounds;
  	StatModels sm(training_set, options.l_value, seq_len, kernel);

    if (options.debugWithCounts)
    {
      sm.getPairwiseBackgrounds(pairwiseBackgrounds, options.background_type);
      sm.getExpectedCounts(expectedEntropies, pairwiseBackgrounds);
      sm.getSimpleStdDevEntropies(standDevEntropies, expectedEntropies);
    }
    else if (options.approximateStandDev)
    {
      sm.getPairwiseBackgrounds(pairwiseBackgrounds, options.background_type);
  	  sm.getExpectedEntropies(expectedEntropies, pairwiseBackgrounds);
  	  sm.getSimpleStdDevEntropies(standDevEntropies, expectedEntropies);
    }
    else if (options.fullPairwiseStandDev)
    {
      sm.getPairwiseBackgrounds(pairwiseBackgrounds, options.background_type);
      sm.getExpectedEntropies(expectedEntropies, pairwiseBackgrounds);
      sm.getComplexStdDevEntropies(standDevEntropies, pairwiseBackgrounds);
    }
    else
    {
      for (int i = 0; i < all_entropies.size(); i++) // itero sulle sequenze
      {
        // fissata una sequenza,
        // media e deviazione standard dell'entropia
        // per ogni parola
        // TODO: incapsulo

        /*double covariance = 0;
        sm.calculateGCovariance(covariance, "AAAA", "CAA", singleBackgrounds[i]);
        std::cout << "covariance" << covariance << std::endl;
        covariance = 0;
        sm.calculateGCovariance(covariance, "AA", "AAA", singleBackgrounds[i]);
        std::cout << "covariance" << covariance << std::endl;
        covariance = 0;
        sm.calculateGCovariance(covariance, "AA", "AA", singleBackgrounds[i]);
        std::cout << "covariance" << covariance << std::endl;
        covariance = 0;
        sm.calculateGCovariance(covariance, "AA", "A", singleBackgrounds[i]);
        std::cout << "covariance" << covariance << std::endl;
        covariance = 0;
        sm.calculateGCovariance(covariance, "A", "A", singleBackgrounds[i]);
        std::cout << "covariance" << covariance << std::endl;

        return 0;*/

        seqan::String<double> expEntropies;
        seqan::String<double> stdDevEntropies;
        sm.expected_entropies(expEntropies, singleBackgrounds[i]);
        sm.entropyStdDev(stdDevEntropies, singleBackgrounds[i]);

        //std::cout << "expEntropies for N2-style standardization" << std::endl;
        //coutString(expEntropies);
        //std::cout << "stdDevEntropies for N2-style standardization" << std::endl;
        //coutString(stdDevEntropies);

        for (int j = 0; j < length(all_entropies[i]); j++) // no divisione per zero
        {
          if ((stdDevEntropies[j] > pow(10.0, -10)) && (stdDevEntropies[j] < pow(10.0, 10)))
          {
            all_entropies[i][j] = (all_entropies[i][j] - expEntropies[j]) / stdDevEntropies[j];
          }
          else
          {
            all_entropies[i][j] = 0.0;
          }
        }

        //std::cout << "After N2-style standardization" << std::endl;
        //coutString(all_entropies[i]);
      }
    }
  }

  if (printTiming)
    std::cerr << "Time to compute expectations and variances " << seqan::sysTime() - startTime << " s." << std::endl;
  startTime = seqan::sysTime();


  for (int i = 0; i < all_entropies.size(); i++)
  {
    std::cout << "entropies " << all_entropies[i] << std::endl;
    coutString(all_entropies[i]);
  }

  if (options.revComplFlag) // non gestisce eventuali postprocessing, ok se d2 o e2
  {
    for (int i = 0; i < all_entropies.size(); i++)
    {
      for (int j = 0; j < length(all_entropies[i]) / 2; j++)
      {
        int last = length(all_entropies[i]) - j -  1;
        if (all_entropies[i][j] > all_entropies[i][last])
        {
          all_entropies[i][last] = all_entropies[i][j];
        }
        else
        {
          all_entropies[i][j] = all_entropies[i][last];
        }
      }
    }
  }

  if (printTiming)
    std::cerr << "Time to consider reverse complement " << seqan::sysTime() - startTime << " s." << std::endl;
  startTime = seqan::sysTime();

  std::cout << "first" << std::endl;
  coutString(all_entropies[0]);
  std::cout << "last" << std::endl;
  coutString(all_entropies[all_entropies.size() - 1]);

  /* Fill score matrix */
  ScoreFiller sf(all_entropies);
  sf.fill_score_matrix(options.measure_type, expectedEntropies, standDevEntropies);
  sf.write_score_matrix(toCString(options.output_filename));

  if (printTiming)
    std::cerr << "Time to fill score matrix " << seqan::sysTime() - startTime << " s." << std::endl;

  return 0;
}