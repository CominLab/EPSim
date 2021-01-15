#ifndef __SCOREFILLER_H_INCLUDED__ 
#define __SCOREFILLER_H_INCLUDED__

#include <vector>
#include <string>
#include <seqan/sequence.h>

class ScoreFiller
{
public:
  std::vector< seqan::String<double> > all_entropies;
  std::vector< std::vector<double> > score_matrix;

  ScoreFiller();

  ScoreFiller(std::vector< seqan::String<double> > const & all_ep);

  double 
  measure_distance(const double & ep1, const double & ep2, std::string const & measure_type);

  void 
  fill_score_matrix(std::string const & measure_type,
                    std::vector< std::vector< seqan::String<double> > > const & expectedEntropies,
                    std::vector< std::vector< seqan::String<double> > > const & standDevEntropies);

  void 
  write_score_matrix(std::string const & output_filename);

};

#endif