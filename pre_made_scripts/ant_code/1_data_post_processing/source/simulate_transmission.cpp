//to compile with option -std=c++11 

using namespace std;

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>

int get_tag_index( vector <int> taglist, int tag){
  int taglist_size (taglist.size());
  int tested_tag(-1);int tag_index (-1);int idx (-1);
  while (tested_tag!=tag){
    idx = idx+1;
    tested_tag = taglist[idx];
  }
  return(idx);
};

struct event{
	double absolute_contamination_time; 
	double relative_contamination_time; 
	int tag; 
	int contaminated_by; 
	int nb_susceptible; 
	int nb_contaminated; 
	double	initial_load;
	double final_load;
	bool infectious;
};

// [[Rcpp::export]]
DataFrame simulate_transmission(DataFrame i_table, DataFrame ant_list, double t0){
  // Fixed input parameters 
  double attenuation = 0.00066;	
  double load_threshold = 0.00024;
  double p_transm = 0.39;
  
  //define variables from input
  vector <int> taglist       = ant_list["tag"];
  vector <int> status        = ant_list["status"];
  vector <double> load       = ant_list["load"];
  int tag_count (taglist.size());
    
  vector <int> Tag1          = i_table["Tag1"];
  vector <int> Tag2          = i_table["Tag2"];
  vector <int> Startframe    = i_table["Startframe"];
  vector <int> Stopframe     = i_table["Stopframe"];
  vector <double> Starttime  = i_table["Starttime"];
  vector <double> Stoptime   = i_table["Stoptime"];
  int table_size (Tag1.size());

	//declare variables that will be used in the loop
	double time;int index_source; int index_recipient;double proba_transfer; int transfer_duration;double random;double x0; double xtot; double amount_transferred;int frame_number;double threshold;
	int index1; int index2; double load1; double load2;
	
	// build random number generator
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();//set a seed linked with the time at which the program is runs; ensures random number sequence will always be differemt
  	std::default_random_engine generator (seed); //set the random number generator using the desited seed
  	std::uniform_real_distribution<double> distribution(0.0,1.0);	//define that you want a uniform diestribution between 0 and 1

	// // print parameter values
	// cout  << "attenuation = " << attenuation << "; load_threshold = " << load_threshold << "; p_transm = " << p_transm << endl;

	/////////////////////////////////////////////////////////
	//perform checks
	/////////////////////////////////////////////////////////		
//  	//print first 2 lines of i_table
// 	cout << "first line of interaction table: " << endl;
// 	cout << Tag1[0] << "," << Tag2[0] << "," << Startframe[0] << "," << Stopframe[0] << "," << Starttime[0] << "," << Stoptime[0] << endl;
// 	cout << "last line of interaction table: " << endl;
// 	cout << Tag1[table_size-1] << "," << Tag2[table_size-1] << "," << Startframe[table_size-1] << "," << Stopframe[table_size-1] << "," << Starttime[table_size-1] << "," << Stoptime[table_size-1] << endl;
// 	
// 	//checking if the get tag function works
// 	cout << "checking get_tag_index function" << endl;
// 	index1 = get_tag_index(taglist,Tag1[0]); index2 = get_tag_index(taglist,Tag2[0]);
// 	cout << "first line: tag1 = " << taglist[index1] << ", infection status = "<< status[index1] <<", load = "<< load[index1]
// 	     << "; tag2 = " << taglist[index2] << ", infection status = "<< status[index2] <<", load = "<< load[index2] << endl;
// 	index1 = get_tag_index(taglist,Tag1[table_size-1]); index2 = get_tag_index(taglist,Tag2[table_size-1]);
// 	cout << "last line: tag1 = " << taglist[index1] << ", infection status = "<< status[index1] <<", load = "<< load[index1]
//       << "; tag2 = " << taglist[index2] << ", infection status = "<< status[index2] <<", load = "<< load[index2] << endl;

	/////////////////////////////////////////////////////////	
	//initialise ant counts 
	/////////////////////////////////////////////////////////
	int contaminated(0);
	int susceptible(0);

	for (int t1 (0); t1<tag_count; t1++){
	  if (status[t1]==0){//untreated workers
	    susceptible=susceptible+1;
	  }
	  if (status[t1]==1){//treated workers
	    contaminated                             =contaminated+1;
	  }
	}
	// cout << "At the beginning: " << contaminated << " infected ants; " << susceptible << " susceptible ants." << endl;
	
	/////////////////////////////////////////////////////////	
	// initialise contamination events table
	/////////////////////////////////////////////////////////
	vector<event> events;
	for (int t1 (0); t1<tag_count; t1++){
	  if (status[t1]==1){//treated workers
	    event infection;
	    infection.absolute_contamination_time = t0; 
	    infection.relative_contamination_time = 0;
	    infection.tag                         = taglist[t1]; 
	    infection.contaminated_by             = -1; 
	    infection.nb_susceptible              = susceptible; 
	    infection.nb_contaminated                 = contaminated; 
	    infection.initial_load                = 1;
	    infection.final_load                  = -1;
	    infection.infectious                  = 1;
	    events.push_back(infection);
	  }
	}
	
	/////////////////////////////////////////////////////////	
	// Loop over interactions
	/////////////////////////////////////////////////////////
	for (int i(0); i < table_size;i++){ //read all interaction lines
		amount_transferred = 0;
	  
		//get frame number
		frame_number = Startframe[i];

		// get the tag index of each ant involved in the interaction
		index1 = get_tag_index(taglist,Tag1[i]); index2 = get_tag_index(taglist,Tag2[i]);
		
		// get the current load of each ant involved in the interaction
		load1 = load[index1];		load2 = load[index2];
	
		// comparison: only performs interaction trial if at least one ant is infectious; i.e. above the load_threshold; and the 2 ants have different load
		if (((load1>=load_threshold)|(load2>=load_threshold))&&(load1!=load2)){
			time = Starttime[i]; //get time of beginning of interactions
			
			//determine which is the source ant and which is the recipient ant
			if (load1>load2){
				index_source=index1;
				index_recipient=index2;
			}else{
				index_source=index2;
				index_recipient=index1;
			}

			//stochastic transmission event:depending on the concentration of the source ant (the larger the amount of infectious propagules, the higher the probability of transfer)
			// Define transmission threshold
			threshold=p_transm*load[index_source];

			// Draw a random number between 0 and 1
			random = distribution(generator);
			
			// Check if a stochastic transmission event should occur (i.e., if random < threshold); if so, perform the rest
			if (random <= threshold){
				// Determine the duration of the interaction in frames
				transfer_duration = Stopframe[i] - Startframe[i] +1;
	    	// Determine the amount transferred depending on the transfer duration
     		 		x0 = load[index_recipient];
     	 			xtot = x0 +  load[index_source];
			      double amount_transferred = ((pow((1-2*attenuation),transfer_duration))*(x0-(xtot/2))+(xtot/2)) - x0;
			      
			// If recipient ant was NOT contaminated yet, add a new event to table events
				if (status[index_recipient]!=1){
					// update status so that ant is now listed as contaminated
					status[index_recipient] = 1;
				  //update number of ants in each category
				  contaminated = contaminated + 1;
				  susceptible = susceptible -1 ;
				  //create infection event
				  event infection;
					infection.absolute_contamination_time   = time; 
					infection.relative_contamination_time   = time - t0;
					infection.tag                           = taglist[index_recipient]; 
					infection.contaminated_by               = taglist [index_source];
					infection.nb_susceptible                = susceptible; 
					infection.nb_contaminated               = contaminated; 
					infection.initial_load                  = amount_transferred;
					infection.final_load                    = -1;
					if (infection.initial_load>=load_threshold){
						infection.infectious                  =1;
					}else{
						infection.infectious                  =0;
					}
					events.push_back(infection);
				}else{
				// if ant was already contaminated, and if her current load is above load threshold, then retrieve event and declare her infectious
					if ((load[index_recipient] + amount_transferred)>=load_threshold){//check if new load higher than threshold
						// find index of the right line
						for (int j(0); j< events.size(); j++){//for each line of the table
						  int tag(events[j].tag);//get tag of interest
						  if (tag==taglist[index_recipient]){
						    events[j].infectious =1;
						  }
						}
					}
				}
				//4/ in any case, update load vector with new load values
				load[index_source] = load[index_source] -amount_transferred ;
				load[index_recipient] = load[index_recipient] +amount_transferred ;
			}
		}
	}//end of transmission loop

	//Prepare dataframe columns for output
	vector <double> absolute_contamination_time (events.size());
	vector <double> relative_contamination_time (events.size());
	vector <int>    tags                        (events.size());
	vector <int>    contaminated_by             (events.size());
	vector <int>    nb_susceptible              (events.size());
	vector <int>    nb_contaminated             (events.size());
	vector <double> initial_load                (events.size());
	vector <double> final_load                  (events.size());
	vector <bool>   infectious                  (events.size());
	
	
	//fill in final load column in events table as well as output columns
	for (int j(0); j< events.size(); j++){//for each line of the table
	  //update final load in events
	  int tag(events[j].tag);//get tag of interest
	  int index = get_tag_index(taglist,tag);
	  events[j].final_load = load[index];
	  // fill in column
	  absolute_contamination_time [j]=events[j].absolute_contamination_time;
	  relative_contamination_time [j]=events[j].relative_contamination_time;
	  tags                        [j]=events[j].tag;
	  contaminated_by             [j]=events[j].contaminated_by;
	  nb_susceptible              [j]=events[j].nb_susceptible;
	  nb_contaminated             [j]=events[j].nb_contaminated;
	  initial_load                [j]=events[j].initial_load;
	  final_load                  [j]=events[j].final_load;
	  infectious                  [j]=events[j].infectious;
	}
return DataFrame::create(_["absolute_contamination_time"]= absolute_contamination_time,
                         _["relative_contamination_time"]= relative_contamination_time,
                         _["tag"]                        = tags,
                         _["contaminated_by"]            = contaminated_by,
                         _["nb_susceptible"]             = nb_susceptible,
                         _["nb_contaminated"]            = nb_contaminated,
                         _["initial_load"]               = initial_load,
                         _["final_load"]                 = final_load,
                         _["infectious"]                 = infectious
                         );
}