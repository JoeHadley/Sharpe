#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
using namespace std;


template <typename T>
class WritePyData {
public:
    std::string filename = "m_temp.py"; // temporary name of the data file

    WritePyData(const std::string& x, bool initialise = true)
        : filename(x)
    {
        std::ofstream myfile;
        if (initialise) {
            myfile.open(filename);
            myfile << "import numpy as np" << "\n" << "\n";
        }
        myfile.close();
    }

    void write_out_vector(const std::string& array_name, const std::vector<T>& array)
    {
        std::ofstream myfile;
        myfile.open(filename, std::ios_base::app);

        myfile << array_name << " = np.array((" << "\n";

        for (auto i = 0; i < array.size(); ++i) {
            if (i != (array.size() - 1)) {
                myfile << array.at(i) << ", " << "\n";
            } else {
                myfile << array.at(i) << "\n";
            }
        }

        myfile << "))" << "\n";
        myfile.close();
    }
};

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





// Function that runs a hypothetical fireShot() function 15 times and gathers the outputs in a vector 
/*vector<double> doTrial(weapon gun, double distance, double targetLowerBound, double targetUpperBound){

    int shotNumber = 15;
    vector<double> shotLocations = {};


    // For each shot in the trial, get a y-value at the distance
    for (int i = 0; i<shotNumber; i++){
        // Placeholder for fireShot() function:
        /*
        double y_value = fireShot(weapon gun, double distance);
        shotLocations.push_back(y_value);
        
    }
    
}*/

//code for random number distribution

double random_number(){


random_device rand_dev;
mt19937 generator(rand_dev());
normal_distribution<double> distr(0, 0.5*M_PI/(double)180);
return distr(generator);
}


// finding the height of the hit

double actualHeight( double inputAngle)
{   
    
    double height;
    double initial_height = 0.7; //dont know syntax to call this from the class
    double g =9.81;
    double v_0 = 450;
    double d =100;

    height = initial_height - g*d*d/(2*v_0*v_0) + d*tan(inputAngle) -g*d*d*tan(inputAngle)*tan(inputAngle)/(2*v_0*v_0);  
    return height;
}
    



// A function that takes the array of shot locations from a trial and counts how many hit the target.
int countHits(vector<double> shotArray,double targetLowerBound, double targetUpperBound){
    
    vector <int> hit;
    for (int i = 0; i<shotArray.size(); i++){

        if (shotArray[i] <= targetUpperBound and shotArray[i] >= targetLowerBound){
            hit.push_back(i);
        }
    }   
    return hit.size();
}
 double doTrial(double thetaC){   

    vector <double> trial; 
    

    for(int i=0; i <=15; i++)
    {    
         double shootingAngle = thetaC + random_number();         
         trial.push_back(actualHeight(shootingAngle));                

    }
    return countHits(trial,0.3,2.3);  
    
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

    //finding the shooting angle

   

    

    //setting up the shot array test

   

    
    vector <int>  trialHits;

    for(int i=0; i<= 1000; i++)
    {
        trialHits.push_back(doTrial(thetaC));


    }

  
   

   /* vector <int> size_quartile;
    

    size_quartile.push_back( countHits(ShotArray,0.3,0.7));
    size_quartile.push_back( countHits(ShotArray,0.7,1.1));
    size_quartile.push_back( countHits(ShotArray,1.1,1.5));
    size_quartile.push_back( countHits(ShotArray,1.5,1.9));
    size_quartile.push_back( countHits(ShotArray,1.9,2.3));
    */
    

    WritePyData<int> writer("hitCounts.py");

    string name = "hitCounts";

    writer.write_out_vector("my_array", trialHits );    
    



    return 0;
}
