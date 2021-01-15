/*
 * seq_splitter.cpp
 *
 * Created by Morris Antonello on 26/06/14
 */


#include <string>
#include <iostream>
#include <map>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/alignment_free.h>
#include <seqan/basic.h>

#include <sstream>

#include <iostream>
#include <fstream>

int 
main(int argc, char ** argv)
{
  char * seq_filename;
  char * output_filename;
  int num_seq_to_extract;
  int sub_length;
  if (argc == 5)
  {
  	seq_filename = argv[1];
  	output_filename = argv[2];
  	num_seq_to_extract = atoi(argv[3]);
  	sub_length = atoi(argv[4]);
  }
  else
  {
  	std::cerr << "Check arguments" << std::endl;
  	return -1;
  }

  seqan::SequenceStream seqStream(seq_filename);
  if (!isGood(seqStream))
  {
    std::cerr << "ERROR: could not open file " << seq_filename << std::endl;
    return -1;
  }

  std::ofstream myfile;
  myfile.open (output_filename	);

  int extracted_strings = 0;
  while (!atEnd(seqStream) && extracted_strings < num_seq_to_extract)
  {
    seqan::CharString id; 
    seqan::CharString seq;
    if (readRecord(id, seq, seqStream) != 0)
    {
      std::cerr << "ERROR: could not read from " << seq_filename << std::endl;
      return -1;
    }
    //std::cout << id << '\t' << seq << std::endl;
    int seq_len = length(seq);
    std::cout << "seq_len " << seq_len << std::endl;

    char * char_seq = toCString(seq);
    std::string str_seq(char_seq);
    int start_index = 0;
    while (seq_len > sub_length && start_index < seq_len && extracted_strings < num_seq_to_extract)
    {
      std::stringstream sub_id;
      sub_id << ">number " << extracted_strings << " start_index " << start_index << " seq_len " << seq_len;
      std::string sub_str = str_seq.substr(start_index, sub_length);
      std::transform(sub_str.begin(), sub_str.end(), sub_str.begin(), toupper);
      myfile << sub_id.str();
      myfile << "\n";
      myfile << sub_str;
      myfile << "\n";
      //appendValue(output_ids, sub_id.str());
      //appendValue(output_set, sub_str);
      start_index += sub_length;
      extracted_strings++;
    }
  }

  myfile.close();

  return 0;
}