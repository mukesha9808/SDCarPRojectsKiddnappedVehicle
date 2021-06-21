/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

//Random generator engine
std::default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  /* Choosing number of particles 
  Error in position 0.3m, 0.3m, 0.01 rad
  desired accuracy  0.01, 0.01, 0.01 rad
  taking twice the particle in each dimension to calculate desired numer of particles*/
  num_particles = (std[0]/0.05)*(std[1]/0.05)*(std[2]/0.005);  // TODO: Set the number of particles
 
  //Initialise particles
  for(int i=0; i< num_particles; ++i){
    Particle prtcl;
    // randon gaussian distrubation for gps input as first data
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);
  
    //Define particles
    prtcl.id=i;
    prtcl.x=dist_x(gen);
    prtcl.y=dist_y(gen);
    prtcl.theta=dist_theta(gen);
    prtcl.weight=1;
  
    particles.push_back(prtcl);
  }
  //Initialisation complete
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  //Predict next position of each particle
  for(unsigned int i=0; i< particles.size(); ++i){
    
    double temp_x;
    double temp_y;
    double temp_theta=particles[i].theta +(yaw_rate*delta_t);
    
    //Calculate position states
    if (fabs(yaw_rate) > 1.0e-5){
      temp_x=particles[i].x + ((sin(temp_theta)-sin(particles[i].theta))*velocity/yaw_rate);
      temp_y=particles[i].y + ((cos(particles[i].theta)-cos(temp_theta))*velocity/yaw_rate);
//std::cout << "I was here 21" << std::endl;
    } else{
      temp_x=particles[i].x + (velocity*delta_t*cos(particles[i].theta));
      temp_y=particles[i].y + (velocity*delta_t*sin(particles[i].theta));
 		//std::cout << "I was here 20" << std::endl;
    }
      
    //Adding gaussian noise to position vector
    std::normal_distribution<double> dist_x(temp_x, std_pos[0]);
    std::normal_distribution<double> dist_y(temp_y, std_pos[1]);
    std::normal_distribution<double> dist_theta(temp_theta, std_pos[2]);

    particles[i].x=dist_x(gen);
    particles[i].y=dist_y(gen);
    particles[i].theta=dist_theta(gen) ;     
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  /* To find associated data first we need to find all observation near particular landmark
  and then we have to find nearest neighbor among those observation.
  
  My solution  first run through all observation and find closest landmark to this observation
  I have created two array of vector distance and index size of array is same as landmark.
  Once closest landmark is found, I push index of observation and distance to landmark at array 
  position of landmark. This way all landmarks will have observation corresponding to this landmark.
  
  Now I traverse through each landmark index in array and find nearest observation to this landmark 
  and asssociate observation to landmark
  */
  
  //Local valid associated observations
  vector<LandmarkObs> validobs;
  
  // Variable to store shortest distance from observation and index of observation
  vector<int> obsIndex[predicted.size()];
  vector<double> obsDist[predicted.size()];
  
  //Loop through each observation to find nearest landmark
  for(unsigned int i=0; i< observations.size();++i) {
    //Reset value for distance and index for finding sortest distance
    double distance=100;
    unsigned int closest_landmark=predicted.size();
    
    // Observation coordinates
    double obs_x=observations[i].x;
    double obs_y=observations[i].y;

    //Find closest landmark
    for(unsigned int j=0; j< predicted.size();++j) {
      //Landmark coordinates
      double landmark_x=predicted[j].x;
      double landmark_y=predicted[j].y;
      
      //ditance from landmark
      double temp_dist=dist(landmark_x, landmark_y, obs_x, obs_y);
      
      //Check if landmark closest
      if(temp_dist < distance) {
        distance=temp_dist;
        closest_landmark=j;
      }
    }   
    //push observations to closest landmark vector array
    if(closest_landmark < predicted.size()) {    
      obsIndex[closest_landmark].push_back(i);
      obsDist[closest_landmark].push_back(distance);
    }
  }
  
  //Traverse through array to find nearest neighbor to each landmark
  for (unsigned int k=0; k< predicted.size(); ++k) {
    
    //Proceed only if observations collected for given landmark
    if (!obsIndex[k].empty()) {
      vector<int> associtionIndex;
      vector<double> associtionDist;
      
      //Reset value for distance and index for finding sortest distance
      double distance=100;
      unsigned int index=observations.size();
      
      //Get indeces off of observations around this landmark
      associtionIndex=obsIndex[k];
      associtionDist=obsDist[k];
      
      //Loop through indeices to find nearest neighbor
      for(unsigned int z=0; z< associtionIndex.size(); ++z) {
        if(associtionDist[z] < distance) {
          distance=associtionDist[z];
          index=associtionIndex[z];
        }        
      }
      
      //Assign association
      if(index < observations.size()) {
        LandmarkObs  foundAssociation;

        //Index of associated landmark
        foundAssociation.id=k;
        //Coordinate of associated observations
        foundAssociation.x=observations[index].x;
        foundAssociation.y=observations[index].y;

        //Push the observations to vecotor
        validobs.push_back(foundAssociation);
      }
    }
  }
  //Update association
  observations=validobs;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
 
  //Store landmarks Predicted observations 
  vector<LandmarkObs> predicted_obs;
  for (unsigned int i=0; i <map_landmarks.landmark_list.size(); ++i) {
    LandmarkObs obs;
    obs.id=map_landmarks.landmark_list[i].id_i;
    obs.x=map_landmarks.landmark_list[i].x_f;
    obs.y=map_landmarks.landmark_list[i].y_f;
    predicted_obs.push_back(obs);
  }
     
  vector<double> weights_temp;
  
  //Process each particle for calculating weigth
  for (unsigned int p=0; p< particles.size(); ++p) {
    
    //Vector to store processed particles to save calculation effore
    vector<LandmarkObs> processPartcle;
    
    //Flag to indicate weight calculation should bre performed
    bool notProcessed=true;
    
    if (!processPartcle.empty()) {      
      for(unsigned int j=0; j < processPartcle.size(); ++j) {
        if ((processPartcle[j].x == particles[p].x) && (processPartcle[j].y == particles[p].y)) {
          notProcessed=false;
          break;
        }
      }
    }
    
    if (notProcessed) {
      
      //Variable for setting association
      vector<int> prtclassociations;
      vector<double> prtclsense_x;
      vector<double> prtclsense_y;
      
      //Initialise weight for this particle
      particles[p].weight=1;

      //Add noise to sensed measurement data
      vector<LandmarkObs> prtclObs_mapcoord;
      
      //loop through each observation  
      for(unsigned int k=0; k < observations.size(); ++k) {
        //Placeholder to store observations in map coordinate
        LandmarkObs Obs_mapcoord;
        
        //Process the observation if within range only
        if(observations[k].x < sensor_range && observations[k].y <sensor_range) {
        
          //Obesrvation in vehicle coordinates
          double x_c=observations[k].x;
          double y_c=observations[k].y;
          
          //Transform coordinates to map frame
          Obs_mapcoord.x= (x_c*cos(particles[p].theta)) - (y_c*sin(particles[p].theta)) + particles[p].x;
          Obs_mapcoord.y= (x_c*sin(particles[p].theta)) + (y_c*cos(particles[p].theta)) + particles[p].y;
          
          //Push global coordinate of obsevation into vector
          prtclObs_mapcoord.push_back(Obs_mapcoord);
        }
      }
      
      //Associate observations to landmarks where each predicted observations has index of transformed observation
      dataAssociation(predicted_obs,prtclObs_mapcoord);
      
      //Calculate weight of particle through each valid observation
      for(unsigned int z=0; z < prtclObs_mapcoord.size(); ++z) {
        //Observation's coordinates in map frame
        double obs_x=prtclObs_mapcoord[z].x;
        double obs_y=prtclObs_mapcoord[z].y;
        
        //Associated landmarks coordinates
        int index=prtclObs_mapcoord[z].id;
        double mu_x=predicted_obs[index].x;
        double mu_y=predicted_obs[index].y;
        
        //Push data for setting association 
        prtclassociations.push_back(predicted_obs[index].id);
        prtclsense_x.push_back(obs_x);
        prtclsense_y.push_back(obs_y);
        
        // calculate normalization term
        double gauss_norm;
        gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
        
        // calculate exponent
        double exponent;
        exponent = (pow(obs_x - mu_x, 2) / (2 * pow(std_landmark[0], 2)))
                   + (pow(obs_y - mu_y, 2) / (2 * pow(std_landmark[1], 2)));
        
        // calculate weight using normalization terms and exponent
        double weight;
        weight = gauss_norm * exp(-exponent);
        
        //Final particle weight
        particles[p].weight=particles[p].weight *weight;
      }
      
      //Push each particle's final weight into vector
      weights_temp.push_back(particles[p].weight);
      
      //Set association data for each particles
      SetAssociations(particles[p],prtclassociations,prtclsense_x,prtclsense_y);
    }
    weights=weights_temp;
  }  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  //Descrete distribution based on weights
  std::discrete_distribution<int> dist_wght(weights.begin(), weights.end());
  
  //Resample particles
  // Set of current particles
  std::vector<Particle> prtcl_resample;
  for(unsigned int n=0; n< particles.size(); ++n) {
    int new_partIndex=dist_wght(gen);
    prtcl_resample.push_back(particles[new_partIndex]);
  }
  
  //Update re-sampled particles
  particles=prtcl_resample;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  
  return s;
}