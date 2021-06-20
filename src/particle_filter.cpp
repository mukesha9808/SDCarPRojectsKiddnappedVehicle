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
  num_particles = (std[0]/0.01)*(std[1]/0.01);  // TODO: Set the number of particles

  // randon gaussian distrubation for gps input as first data
  std::normal_distribution<double> dist_x(0, std[0]);
  std::normal_distribution<double> dist_y(0, std[1]);
  std::normal_distribution<double> dist_theta(0, std[2]);
  
  //Initialise particles
  for(int i=0; i< num_particles; ++i){
    Particle prtcl;
   	prtcl.id=i;
    prtcl.x=x + dist_x(gen);
    prtcl.y=y + dist_y(gen);
    prtcl.theta=theta + dist_theta(gen);
    prtcl.weight=1;
  	
    particles.push_back(prtcl);
  }

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
  // randon gaussian distrubation for noise
  
  //predict particles
  for(unsigned int i=0; i< particles.size(); ++i){

    std::normal_distribution<double> dist_theta(yaw_rate, std_pos[2]);
    double temp_x;
    double temp_y;
    double noisy_yawrate=  dist_theta(gen);
    double theta1=particles[i].theta +(noisy_yawrate*delta_t); //+ theta_noise;

    
    if (fabs(noisy_yawrate) > 1.0e-5){
      temp_x=particles[i].x + ((sin(theta1)-sin(particles[i].theta))*velocity/noisy_yawrate);
      temp_y=particles[i].y + ((cos(particles[i].theta)-cos(theta1))*velocity/noisy_yawrate);
//std::cout << "I was here 21" << std::endl;
    } else{
      temp_x=particles[i].x + (velocity*delta_t*cos(particles[i].theta));
      temp_y=particles[i].y + (velocity*delta_t*sin(particles[i].theta));
 		//std::cout << "I was here 20" << std::endl;
    }
     
    particles[i].theta=theta1 ; 
    //std::cout << "I was here 19" << std::endl;
    std::normal_distribution<double> dist_x(temp_x, std_pos[0]);
    std::normal_distribution<double> dist_y(temp_y, std_pos[1]);
    
    particles[i].x=dist_x(gen);
    particles[i].y=dist_y(gen);
    
  }
  
  //std::cout << "theta" << particles[0].theta  <<  std::endl;
  //std::cout << "noise   ("<< ((sin(particles[0].theta +(yaw_rate*delta_t) + theta_noise)-sin(particles[0].theta))*velocity/yaw_rate)<< " ," << ((cos(particles[0].theta)-cos(particles[0].theta +(yaw_rate*delta_t) + theta_noise))*velocity/yaw_rate)<<"," << theta_noise <<")" << std::endl;
  //std::cout << "vel, thetha   ("<< velocity << " ," << yaw_rate<<"," << delta_t<<")" << std::endl;
  //std::cout << "xy vale   ("<< particles[0].x << " ," << particles[0].y << "," << particles[0].theta <<")" << std::endl;

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
  // Variable to store shortest distance from observation and index of observation
  //Local valid associated observations
  vector<LandmarkObs> validobs;
  
   // Variable to store shortest distance from observation and index of observation
  vector<int> obsIndex[predicted.size()];
  vector<double> obsDist[predicted.size()];
  
  //Loop through each observation to find nearest landmark
  for(unsigned int i=0; i< observations.size();++i) {
    
    // Observation coordinates
    double obs_x=observations[i].x;
    double obs_y=observations[i].y;
    
    double distance=100;
    int closest_landmark=predicted.size();
    //Find closest landmark
    for(unsigned int j=0; j< predicted.size();++j) {
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
    
    //Collect all the observations to particular landmark
    if(closest_landmark < predicted.size()) {
      
      obsIndex[closest_landmark].push_back(i);
      obsDist[closest_landmark].push_back(distance);
    }
  }
  
  for (unsigned int k=0; k< predicted.size(); ++k) {
    if (obsIndex[k].size() > 0) {
      vector<int> associtionIndex;
      vector<double> associtionDist;
      
      associtionIndex=obsIndex[k];
      associtionDist=obsDist[k];
      
      double distance=100;
      int index=observations.size();
      
      for(unsigned int z=0; z< associtionIndex.size(); ++z) {
        if(associtionDist[z] < distance) {
          distance=associtionDist[z];
          index=associtionIndex[z];
        }        
      }
      
      if(index < observations.size()) {
        LandmarkObs  foundAssociation;
        foundAssociation.id=k;
        foundAssociation.x=observations[index].x;
        foundAssociation.y=observations[index].y;
        //std::cout << "Observation " << index << " coord  (" << observations[index].x << "," << observations[index].y << ")" << std::endl;
        //std::cout << "Landmark " << k  << " coord  (" << predicted[k].x << "," << predicted[k].y << ")" << std::endl;
        validobs.push_back(foundAssociation);
      }
    }
  }
  observations=validobs;
 //std::cout << "I was here 111";
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
  //radom noise in measurement
  
  std::normal_distribution<double> dist_x(0, std_landmark[0]);
  std::normal_distribution<double> dist_y(0, std_landmark[1]);
   
   
  //Store land marks Predicted observations 
  vector<LandmarkObs> predicted_obs;
  for (unsigned int i=0; i <map_landmarks.landmark_list.size(); ++i) {
    LandmarkObs obs;
    obs.x=map_landmarks.landmark_list[i].x_f;
    obs.y=map_landmarks.landmark_list[i].y_f;
    //std::cout << "I was here 11" << std::endl;
    predicted_obs.push_back(obs);
  }
  
   
  vector<double> weights_temp;
  
  //Process each particle for calculating weigth
  for (unsigned int p=0; p< particles.size(); ++p) {
 	
    
    particles[p].weight=1;
    
    vector<int> prtclassociations;
  vector<double> prtclsense_x;
  vector<double> prtclsense_y;
     
    
    
    //Add noise to sensed measurement data
    vector<LandmarkObs> prtclObs_mapcoord;
    for(unsigned int k=0; k < observations.size(); ++k) {
      LandmarkObs Obs_mapcoord;
     //std::cout << "Obs "<< k <<"   cord is   ("<< observations[k].x << "," <<observations[k].y << ")" << std::endl;
      if(observations[k].x < sensor_range && observations[k].y <sensor_range) {
        //Add gaussian noise to measuremnet
        double x_c=observations[k].x;// +dist_x(gen);
        double y_c=observations[k].y;// +dist_y(gen);
        
    	
        //Transform coordinates to map frame
        Obs_mapcoord.x= (x_c*cos(particles[p].theta)) - (y_c*sin(particles[p].theta)) + particles[p].x;
        Obs_mapcoord.y= (x_c*sin(particles[p].theta)) + (y_c*cos(particles[p].theta)) + particles[p].y;
        
		//std::cout << "I was here 9" << std::endl;
        prtclObs_mapcoord.push_back(Obs_mapcoord);
      }
    }
    
    //std::cout << "I was here 8" << std::endl;
    //Associate observations to landmarks where each predicted observations has index of transformed observation
    dataAssociation(predicted_obs,prtclObs_mapcoord);
     //std::cout << "I was here 7" << std::endl; 
    
    //Calculate weight
    for(unsigned int z=0; z < prtclObs_mapcoord.size(); ++z) {
      double obs_x=prtclObs_mapcoord[z].x;
      double obs_y=prtclObs_mapcoord[z].y;

      //Associated observation as x,y
      int index=prtclObs_mapcoord[z].id;
      
      prtclassociations.push_back(index);
      prtclsense_x.push_back(obs_x);
      prtclsense_y.push_back(obs_y);
      
      
      //take landmark observation coordinates as mu for multi variate disrbituion
      double mu_x=predicted_obs[index].x;
      double mu_y=predicted_obs[index].y;
		//std::cout << "myxy  (" << mu_x<<",  "<< mu_y  <<")"<< std::endl;
      //std::cout << "obsxy  (" << obs_x<<",  "<< obs_y <<")"<< std::endl;
      
      // calculate normalization term
      double gauss_norm;
      gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      //std::cout << "I was here 6"<< gauss_norm << std::endl;
      // calculate exponent
      double exponent;
      exponent = (pow(obs_x - mu_x, 2) / (2 * pow(std_landmark[0], 2)))
               + (pow(obs_y - mu_y, 2) / (2 * pow(std_landmark[1], 2)));
      //std::cout << "I was here 7"<< exponent << std::endl;
      // calculate weight using normalization terms and exponent
      double weight;
      weight = gauss_norm * exp(-exponent);
      //std::cout << "exp" << exponent<< std::endl;
		//std::cout << "weight " << p <<" of obs  "<< z  <<"th"<< "  is  "<<weight << ")"  << std::endl;
      	//std::cout << "weight " << p <<" after obs  "<< z << "  is  "<<particles[p].weight << ")"  << std::endl;
      //Final particle weight
      particles[p].weight=particles[p].weight *weight;
      }
    
    weights_temp.push_back(particles[p].weight);

      //std::cout << "weight of   " << p <<"  is  " <<particles[p].weight << ")"  << std::endl;
    SetAssociations(particles[p],prtclassociations,prtclsense_x,prtclsense_y);
    //std::cout << "I was here 5  "<< std::endl;
  }
  weights=weights_temp;
 //std::cout << "I was here 4" << std::endl;
  
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
  //std::cout << "I was here 31   "<< weights.size() << std::endl;

  //std::cout << "all weights (" ;
  
  //for(int check=0; check < weights.size(); ++check) {
  //  std::cout << weights[check] <<", ";
  //}
 // std::cout << std::endl;
  
  //Resample particles
  // Set of current particles
  std::vector<Particle> prtcl_resample;
  for(unsigned int n=0; n< particles.size(); ++n) {
    int new_partIndex=dist_wght(gen);
    //std::cout << "gnerated    " << new_partIndex << std::endl;
    //std::cout << "I was here 32" << new_partIndex << std::endl;
    prtcl_resample.push_back(particles[new_partIndex]);
  
  }
  particles=prtcl_resample;
  //std::cout << "I was here 3" << std::endl;
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
  //std::cout << "I was here 2" << std::endl;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  //std::cout << "I was here 1" << std::endl;
  
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
  
  //std::cout << "I was here 0" << std::endl;
  return s;
}