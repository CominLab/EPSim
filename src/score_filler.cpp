#include "score_filler.h"

#include <iostream>
#include <fstream>

ScoreFiller::ScoreFiller() {}

ScoreFiller::ScoreFiller(std::vector< seqan::String<double> > const & all_ep) 
{
  all_entropies = all_ep;

  int number_of_records = all_entropies.size();
  score_matrix.resize(number_of_records);
  for (int i = 0; i < number_of_records; i++)
  	score_matrix[i].resize(number_of_records);
}

double 
ScoreFiller::measure_distance(double const & ep1, double const & ep2, std::string const & measure_type)
{
  if (measure_type.compare("D2") == 0 || measure_type.compare("D2*") == 0 || measure_type.compare("AR") == 0 || measure_type.compare("ARM") == 0) 
  {
    return ep1 * ep2;
  }
  else
  {
    std::cerr << "Unknown measure" << std::endl;
  }
}

void 
ScoreFiller::fill_score_matrix(std::string const & measure_type,
                               std::vector< std::vector< seqan::String<double> > > const & expectedEntropies,
                               std::vector< std::vector< seqan::String<double> > > const & standDevEntropies) 
{
  int number_of_records = all_entropies.size();
  int number_of_motifs = length(all_entropies[0]);

  for (int i = 0; i < number_of_motifs; i++)
  {
    for (int j = 0; j < number_of_records; j++)
    {
   	  for (int k = 0; k < number_of_records; k++)
      {
        double ep1 = all_entropies[j][i];
        double ep2 = all_entropies[k][i];
        if (measure_type.compare("D2*") == 0)
        {
          if (expectedEntropies.size() > 0 && standDevEntropies.size() > 0) // cioè se standardizzazione approssimata
          {
            // Test if variance is larger than 0 and smaller than inf before dividing
            if ((standDevEntropies[j][k][i] > pow(10.0, -10)) && (standDevEntropies[j][k][i] < pow(10.0, 10)))
            {
              ep1 = (ep1 - expectedEntropies[j][k][i]) / standDevEntropies[j][k][i]; 
              ep2 = (ep2 - expectedEntropies[j][k][i]) / standDevEntropies[j][k][i];
            }
            else
            {
              ep1 = ep2 = 0.0;
            }
          }
        }
  	    score_matrix[j][k] += measure_distance(ep1, ep2, measure_type);
      }
    }
  }
}

void 
ScoreFiller::write_score_matrix(std::string const & output_filename) 
{
  std::ofstream output_file(output_filename.c_str());
  for (int i = 0; i < score_matrix.size(); i++)
  {
    for (int j = 0; j < score_matrix[i].size(); j++)
      output_file << score_matrix[i][j] << "\t";
    output_file << std::endl;
  }
  output_file.close();
}