/*
 * real_seq.cpp
 *
 * Created by Morris Antonello on 03/06/14
 */

#include <string>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <stdlib.h> // srand rand

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int 
main(int argc, char ** argv)
{

  if (argc < 3)
  {
  	std::cerr << "Not enough args" << std::endl;
    return -1;
  }

  char * starting_background = argv[1];
  char * output_seq = argv[2];
  double probability = atof(argv[3]);
  std::vector < std::string > implanted_words; 

  for (int i = 4; i < argc; i++)
  {
  	implanted_words.push_back(argv[i]);
  }

  std::cout << "vuoi impiantare tot parole " << implanted_words.size() << std::endl;
  std::cout << "con probabilità " << probability << std::endl;


  // open file to read

  seqan::SequenceStream inSeqStream(starting_background);
  if (!isGood(inSeqStream))
  {
    std::cerr << "ERROR: could not open file " << starting_background << std::endl;
    return -1;
  }

  // open file to write

  seqan::SequenceStream outSeqStream(output_seq, seqan::SequenceStream::WRITE);
  if (!isGood(outSeqStream))
  {
      std::cerr << "ERROR: Could not open the file.\n";
      return 1;
  }

  // read record by record
  // modify it
  // write it

  while (!atEnd(inSeqStream))
  {

    // read
    seqan::CharString id; 
    seqan::CharString seq;
    if (readRecord(id, seq, inSeqStream) != 0)
    {
      std::cerr << "ERROR: could not read from " << starting_background << std::endl;
      return -1;
    }
    //std::cout << id << '\t' << seq << std::endl;


    // write
    id += "with implanted_words";
    std::string modified_seq = toCString(seq);

    srand ( time(NULL) );

    int current_word = 0;
    for (int i = 0; i < modified_seq.length() - implanted_words[current_word].length() + 1; i++)
    {
      double random_num = fRand(0, 1);
      //std::cout << "random_num is " << random_num << std::endl;
      if (random_num < probability)
      {
        modified_seq.replace(i, implanted_words[current_word].length(), implanted_words[current_word]);

        //std::cout << "Ho inserito la parola " << implanted_words[current_word] << std::endl;

        i = i + implanted_words[current_word].length() - 1; // -1 per l'i++ // no overlapping

        current_word++;
        if (current_word >= implanted_words.size())
          current_word = 0;
      }
    }

    if (writeRecord(outSeqStream, id, modified_seq) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }

  }

  return 0;
}