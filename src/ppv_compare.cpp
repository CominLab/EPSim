#include <dirent.h> // getdir opendir readdir
#include <errno.h>
#include <vector>
#include <string> // std::string, std::to_string
#include <sstream> // stringstream
#include <iostream> // cout cerr
#include <fstream> // ifstream
#include <numeric> // accumulate
#include <stdlib.h> // atoi
#include <map>

// Get files in dir
int getdir (std::string dir, std::vector<std::string> &files)
{

  unsigned char isFile =0x8;
  DIR *dp;
  struct dirent *dirp;

  if((dp  = opendir(dir.c_str())) == NULL) // open dir
  {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }
  
  while(dirp=readdir(dp)) // read dir
  {
    if ( dirp->d_type == isFile) // if isFile
    {
      files.push_back(std::string(dirp->d_name));
    }
  }
  
  closedir(dp); // close
  
  return 0;

}

void extract_params(std::string current_file, std::string & delimiter, std::string & k, std::string & length, std::string & l)
{
  int j = 0;
  size_t pos_delimiter = 0;
  
  while ((pos_delimiter = current_file.find(delimiter)) != std::string::npos)
  {
    std::string token = current_file.substr(0, pos_delimiter); // substring between two delimiters (or till or from)

    if (j == 0)
      k = token;
    else if (j == 1)
      length = token;

    current_file.erase(0, pos_delimiter + delimiter.length());
    j++;
  }
  l = current_file;
}

int get_avg(std::string & current_pathname, double & avg)
{
  std::ifstream current_filestream(current_pathname.c_str());
  if (!current_filestream.is_open())
  {
    std::cerr << "Unable to open file: " << current_pathname;
    return 1;
  }

  std::string tmp; // "Performance:"
  double perf_value;
  std::vector<double> performances;
  while(current_filestream >> tmp >> perf_value)
  {
    performances.push_back(perf_value);
    std::cout << tmp << " " << perf_value << std::endl; 
  }
  current_filestream.close();
    
  double sum = std::accumulate(performances.begin(), performances.end(), 0.0);
  avg = sum / performances.size();
}

int main( int argc, char** argv )
{

  /*
   * Get files in directories
   */

  std::vector<std::string> directories;
  std::cout << argc << std::endl;
  for (int i = 1; i < argc; i++)
  {
    //std::cout << argv[i] << std::endl;
    directories.push_back(argv[i]);//("./performances_d2_AGCCA/"); //("./performances_d2star_AGCCA/");
  }

  std::map < std::string, std::map <std::string, std::map <std::string, std::map <int, double> > > > k_length_method_l_avg;

  for (int j = 0; j < directories.size(); j++)
  {
  
    std::vector<std::string> files;

    getdir(directories[j], files);

    // print out filenames in dir
    for (int i = 0; i < files.size(); i++)
    {
      std::cout << files[i] << std::endl;
    }

    /*
     * Get params from filename
     */

    std::string delimiter = "-";

    // For each file
    for (int i = 0; i < files.size(); i++)
    {
      std::string k; // 1 or 2
      std::string length; // 100 150 250 500 1000 if k = 1 or 250 500 1000 2000 3000 4000 5000 10000 if k = 2
      std::string l;
      extract_params(files[i], delimiter, k, length, l);

      std::cout << "method is " << directories[j]
                << "k is " << k
                << " - length is " << length
                << " - l is " << l << std::endl;

      /*
       * Open file
       * get avg
       */

      std::string current_pathname = directories[j] + files[i];
      double avg = 0.0;
      get_avg(current_pathname, avg);

      std::cout << "Avg is " << avg << std::endl;
      
      // l is int instead of string because:
      // 1) the elements in a map are always sorted by its key following a specific strict weak ordering criterion 
      // indicated by its internal comparison object (of type Compare).
      // 2) l entries in .dat files must be ordered following int criterion instead of string, otherwise plot is invalid
      k_length_method_l_avg[k][length][directories[j]][atoi(l.c_str())] = avg; 
    }

  }

  /*
   * Write .dat
   */

  std::map < std::string, std::map <std::string, std::map <std::string, std::map <int, double> > > >::iterator it_k;
  for (it_k = k_length_method_l_avg.begin(); it_k != k_length_method_l_avg.end(); ++it_k)
  {

    std::map <std::string, std::map <std::string, std::map <int, double> > >::iterator it_length;
    for (it_length = it_k->second.begin(); it_length != it_k->second.end(); ++it_length)
    {

       std::string output_filename = "../dat/k" + it_k->first + "-length" + it_length->first + ".dat";
       std::ofstream datafile (output_filename.c_str());
      
       if (datafile.is_open())
       {
         
         datafile << "# k " << it_k->first << " length " << it_length->first << std::endl;

         int method_index = 0;
         std::map <std::string, std::map <int, double> >::iterator it_method;
         for (it_method = it_length->second.begin(); it_method != it_length->second.end(); ++it_method)
         {

           datafile << "# method " << it_method->first << std::endl;
           datafile << "# l " << "ppv " << "color" << std::endl;

           std::map <int, double>::iterator it_l;
           for (it_l = it_method->second.begin(); it_l != it_method->second.end(); ++it_l)
           {
             datafile << it_l->first  << " " << it_l->second << " " << method_index/*it_method->first*/ << std::endl;
           }
           datafile << std::endl;
           method_index++;
         }

         datafile.close();
       
       }
       else 
       {
         std::cerr << "Unable to open output file " << output_filename << std::endl;
       }
    
    }
  }

  return 0;
}
