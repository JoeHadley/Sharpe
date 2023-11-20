#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
using namespace std;

// Class for weapons. Internal variables are firingSpeed, firingHeight and shotSpread
class weapon{
   
    public:
    
    double firingSpeed;
    double firingHeight;
    double shotSpread;

    // Constructor sets the weapon's stats
    weapon(double v0, double h1, double sigma){
        firingSpeed = v0;
        firingHeight = h1;
        // Shot Spread is given in degrees, but internally we want to use radians:
        shotSpread = sigma*M_PI/180;
    }



    // Function to find the optimal angle for the weapon at a given distance with a given target 
    double findOptimalAngle(double d, double yt, double g){

        // Calculate the optimal angle. See Joe's notes for working out that this is the correct formula
        double tanTheta =  firingSpeed*firingSpeed/(g*d) - sqrt( pow(firingSpeed,4)/pow(g*d,2)  -  (1 + 2*pow(firingSpeed,2)*(yt-firingHeight)/(g*pow(d,2)))  );
        double thetaC = atan(tanTheta);
        return thetaC;
    }
};

void rnum(double& input, opt_angle){ //gaussian distribution for the angle with standard deviation  3 degrees (does this need to be in radians?)
      const double mean = opt_angle;
      const double standdev = 3;
      random_device rand_dev;
      mt19937 gen(rand_dev());
      normal_distribution<double> distr(mean, standdev);
      input = distr(gen);
}


// Function that runs a hypothetical fireShot() function 15 times and gathers the outputs in a vector 
vector<double> doTrial(weapon gun, double distance, double targetLowerBound, double targetUpperBound){

    int shotNumber = 15;
    vector<double> shotLocations = {};


    // For each shot in the trial, get a y-value at the distance
    for (int i = 0; i<shotNumber; i++){
        // Placeholder for fireShot() function:
        /*
        double y_value = fireShot(weapon gun, double distance);
        shotLocations.push_back(y_value);
        */
    }
    
}


// A function that takes the array of shot locations from a trial and counts how many hit the target.
int countHits(vector<double> shotArray,double targetLowerBound, double targetUpperBound){
    
    int hitCount = 0;
    for (int i = 0; i<shotArray.size(); i++){

        if (shotArray[i] <= targetUpperBound and shotArray[i] >= targetLowerBound){
            hitCount++;
        }
    }   
    return hitCount;
}






int main(){

    // Constants that won't change for this project:
    const double targetUpperBound = 2.3;
    const double targetLowerBound = 0.3;
    const double targetCentre = (targetUpperBound + targetLowerBound)/2;
    const double g = 9.81;


    // Constants that may change:
    double distance = 100;


    // Declare weapon objects with their firing speed, firing height, and shot spread in degrees. The weapon class constructor will convert to radians.
    weapon musket(450,0.7, 3.0);
    weapon rifle(600, 0.3, 1.0);


    // Find optimal angle
    double thetaC = musket.findOptimalAngle(distance, targetCentre, g);
    cout << thetaC;


    return 0;
}
