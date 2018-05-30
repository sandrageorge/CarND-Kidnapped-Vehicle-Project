/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Create normal distributions for x, y and theta.
	normal_distribution<double> dist_x(0.0, std[0]);
	normal_distribution<double> dist_y(0.0, std[1]);
	normal_distribution<double> dist_theta(0.0, std[2]);

	num_particles = 500;
	particles.resize(num_particles); // Resize the `particles` vector to fit desired number of particles
	weights.resize(num_particles);
	double init_weight = 1.0;///num_particles;

	for (unsigned int i = 0; i < particles.size(); i++){

		particles[i].id = i;
		particles[i].x = x + dist_x(gen);
		particles[i].y = y + dist_y(gen);
		particles[i].theta = theta + dist_theta(gen);
		particles[i].weight = init_weight;
	}	

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);

	for (unsigned int i = 0; i < particles.size(); i++){

        if (fabs(yaw_rate) < 0.00001){

            particles[i].x += velocity * cos(particles[i].theta) * delta_t + dist_x(gen);
            particles[i].y += velocity * sin(particles[i].theta) * delta_t + dist_y(gen);
			particles[i].theta += dist_theta(gen);
        }
        else{

            const double theta = particles[i].theta + yaw_rate * delta_t;
            particles[i].x += velocity / yaw_rate * (sin(theta) - sin(particles[i].theta)) + dist_x(gen);
            particles[i].y += velocity / yaw_rate * (-cos(theta) + cos(particles[i].theta)) + dist_y(gen);
            particles[i].theta = theta + dist_theta(gen);
        }
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	// loop over all observations
	for (unsigned int i = 0; i < observations.size(); i++) {
		LandmarkObs obs = observations[i];
		// set initial minimum distance to the largest possible value
		double min_dist = numeric_limits<double>::max(); //max = 1.7976931348623157e+308
		int landmark_id = -1;
		
		// loop over all predicted groudntruth landmarks
		for (unsigned int j = 0; j < predicted.size(); j++) {
			LandmarkObs pred = predicted[j];
			
			// calculate distance between observation and predicted groundtruth landmark
			double curr_dist = dist(obs.x, obs.y, pred.x, pred.y);

			// search for predicted landmark with minimum distance to current observation
			if (curr_dist < min_dist) {
				min_dist = curr_dist;
				landmark_id = pred.id;
			}
		}
		// set the observation's id to the nearest predicted landmark's id
		observations[i].id = landmark_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// loop over all particles
	for (unsigned int i = 0; i < particles.size(); i++) {

		// double p_x = particles[i].x;
    	// double p_y = particles[i].y;
    	// double p_theta = particles[i].theta;
		Particle p = particles[i];

		// filter out the map landmark locations predicted to be outside sensor range of the particle
	    vector<LandmarkObs> predictions;
		// loop over all landmarks
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			float lm_x = map_landmarks.landmark_list[j].x_f;
      		float lm_y = map_landmarks.landmark_list[j].y_f;
      		int lm_id = map_landmarks.landmark_list[j].id_i;

			if (dist(lm_x, lm_y, p.x, p.y) <= sensor_range) {
        		predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
      		}
		}

		// transform observations from vehicle reference frame to map reference frame
		vector<LandmarkObs> t_observations;
    	for (unsigned int j = 0; j < observations.size(); j++) {
			double t_x = cos(p.theta) * observations[j].x - sin(p.theta) * observations[j].y + p.x;
			double t_y = sin(p.theta) * observations[j].x + cos(p.theta) * observations[j].y + p.y;
			t_observations.push_back(LandmarkObs{ observations[j].id, t_x, t_y });
		}

	    dataAssociation(predictions, t_observations);

		particles[i].weight = 1.0;

		for (unsigned int j = 0; j < t_observations.size(); j++) {
			LandmarkObs t_obs = t_observations[j];
			LandmarkObs pred;

			int pred_id = t_obs.id;

			for (unsigned int k = 0; k < predictions.size(); k++) {
				if (predictions[k].id == pred_id) {
				pred.x = predictions[k].x;
				pred.y = predictions[k].y;
				}
			}

			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			particles[i].weight *= ( 1 / (2 * M_PI * std_x * std_y ) ) * exp( -( pow( pred.x - t_obs.x, 2 ) / ( 2 * pow( std_x, 2 ) ) + 
												( pow( pred.y - t_obs.y, 2 ) / ( 2  *pow( std_y, 2 ) ) ) ) );

		}	
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> resampled_particles;

	for (int i = 0; i < num_particles; i++) {
		weights[i] = particles[i].weight;
	}

	double max_weight = *max_element(weights.begin(), weights.end());

	uniform_real_distribution<double> unirealdist(0.0, 2.0 * max_weight);

	uniform_int_distribution<int> uniintdist(0, num_particles-1);

	auto index = uniintdist(gen);

	double beta = 0.0;

	for (int i = 0; i < num_particles; i++) {
		beta += unirealdist(gen);
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampled_particles.push_back(particles[index]);
	}

	particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
