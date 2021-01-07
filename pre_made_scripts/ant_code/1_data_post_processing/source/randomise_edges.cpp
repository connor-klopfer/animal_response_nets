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


// [[Rcpp::export]]
DataFrame  randomise_edges(DataFrame i_table){
  /////////////////////////////////////////////////////////		
  //define columns
  vector <int> Tag1       = i_table["Tag1"];
  vector <int> Tag2       = i_table["Tag2"];
  vector <int> Startframe = i_table["Startframe"];
  vector <int> Stopframe  = i_table["Stopframe"];

  /////////////////////////////////////////////////////////
  //get size of interaction table
  int table_size;
  table_size = Tag1.size();
  
  //initialise resampled data
  vector<int> resampled;
  //initialise tag id variables
  int i; int ii; int j; int jj;
  int new_i; int new_ii; int new_j; int new_jj;
  // initialise new index variable
  int index;
  
  double min_time; double max_time; double min_ttime; double max_ttime;
  // build random number generator
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();//set a seed linked with the time at which the program is runs; ensures random number sequence will always be differemt
  std::default_random_engine generator (seed); //set the random number generator using the desired seed
  
  /////////////////////////////////////////////////////////
  // Open input files
  
  // initialise random number generators	
  std::uniform_int_distribution<int> uniform_integer_distribution(0,(table_size-1));
  std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
  
  /////////////////////////////////////////////////////////	
  // randomisation loop	
  /////////////////////////////////////////////////////////
  for (int edge(0); edge < table_size;edge++){ //read all interaction lines
    
    //check if edge has already been resampled
    bool sampled (0);
    if (resampled.size()>0){
      for (int j(0); j < resampled.size(); j++){
        if (resampled[j]==edge){
          sampled = 1;
        }
      }
    }
    //only perform resampling if edge has not already been resampled
    if (!sampled){
      // get tags 
      i = 	Tag1[edge];
      j = 	Tag2[edge]; 
      // get min and max time of the interaction
      min_time = Startframe[edge] ;
      max_time = Stopframe[edge] ;
      //initialise condition variable
      bool condition (0);
      while(!condition){// while condition not fulfilled
        index = uniform_integer_distribution(generator);//index = index of the other edge with which we are going to draw
        // get tag indices
        ii = Tag1[index];
        jj = Tag2[index];
        
        // get min and max time of the interaction
        min_ttime = Startframe[index];
        max_ttime = Stopframe[index];
        
        // now draw a random number between 0 and 1 to determine which switch you want to do
        double rand (uniform_distribution(generator));
        if (rand<0.5){
          //cout << "case1" << endl;//with a probability of 1/2, replace (i,j) and (ii,jj) by (i,jj) and (ii,j)
          new_i = i;new_ii = ii;
          new_j = jj; new_jj = j;
        }else{
          //cout << "case2" << endl; //with a probability of 1/2, replace (i,j) and (ii,jj) by (i,ii) and (j,jj)
          new_i = i;new_ii = j;
          new_j = ii;new_jj = jj;
        }
        
        // check that one is not creating any multiple edge
        // for this, check that there are no overlapping interaction with edge or index that already involves the same two parners
        bool multiple_edge(0);
        for (int k(0); k < table_size;k++){
          //interactions overlapping with edge
          if (((Startframe[k]>=min_time)&&(Startframe[k]<=max_time))||((Stopframe[k]>=min_time)&&(Stopframe[k]<=max_time))){
            if (((Tag1[k]==new_i)&&(Tag2[k]==new_j))||((Tag1[k]==new_j)&&(Tag2[k]==new_i))){
              multiple_edge=1;								
            }
          }
          //interactions overlapping with index
          if (((Startframe[k]>=min_ttime)&&(Startframe[k]<=max_ttime))||((Stopframe[k]>=min_ttime)&&(Stopframe[k]<=max_ttime))){
            if (((Tag1[k]==new_ii)&&(Tag2[k]==new_jj))||((Tag1[k]==new_jj)&&(Tag2[k]==new_ii))){
              multiple_edge=1;								
            }
          }
        }
        // now test whether rearrangement is valid
        if 	(
            (edge!=index)//condition 1: have drawn something different
          &&
            (new_i!=new_j)// condition 2: no self edge on edge
          &&
            (new_ii!=new_jj)// condition 3: no self edge on index
          &&
            (!multiple_edge)// condition 4: have not created a multiple edge
        ){
          condition = 1;
        }
      }//while condition
      //once condition is met, make the change into the table
      
      ////first updates new tag indices
      Tag1[edge] = new_i;
      Tag2[edge] = new_j;	
      Tag1[index] = new_ii;
      Tag2[index] = new_jj;	
      
      ////second update resampled
      resampled.push_back(edge);
      resampled.push_back(index);				
    }//if(!sampled)
  }// for edge
  return DataFrame::create(_["Tag1"]= Tag1, _["Tag2"]= Tag2);
}
