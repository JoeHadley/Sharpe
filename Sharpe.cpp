#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include <ctime>
using namespace std;

const double g = 9.81;


mt19937 initializeRandomGenerator() {
    unsigned seed = static_cast<unsigned>(time(nullptr));
    mt19937 rng(seed);
    return rng;
}

void linspace(vector<double> &interpVector, double startValue, double endValue, int numPoints) {
    interpVector.clear();

    if (numPoints == 1){
        interpVector.push_back(startValue);
    }
    else {
        double step = (endValue - startValue) / (numPoints - 1);
        for (int n = 0; n < numPoints; n++) {
            interpVector.push_back(startValue + step * n);
        }
    }
}

// Target class to hold distances
class target{

    public:
    const double upper;
    const double lower;
    const double centre;
    double distance;

    target(double lowerBound, double upperBound, double distance) : lower(lowerBound), upper(upperBound), centre((lowerBound + upperBound)/2), distance(distance) {}

};

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
            myfile << "import numpy as np" << "\n" << "import matplotlib.pyplot as plt" << "\n";
        }
        myfile.close();
    }




    void write_out_hist_vector(const std::string& array_name, const std::vector<T>& array)
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
        myfile << "bins = np.unique(" << array_name<<")"<<"\n";
        myfile << "plt.hist("<<array_name<<",bins)" << "\n" << "plt.show()"  << "\n";
        myfile.close();
    }


    void write_out_comparison_vectors(const std::string& indepVar_name, const std::vector<T>& indepVar,const vector<std::string>& depVar_names, const vector<std::vector<T>>& depVars,const vector<std::string>& depVarError_names, const vector<std::vector<T>>& depVarErrors)
    {
        std::ofstream myfile;
        myfile.open(filename, std::ios_base::app);

        // Write the independent variable
        myfile << indepVar_name << " = np.array((" << "\n";
        for (auto i = 0; i < indepVar.size(); ++i) {
            if (i != (indepVar.size() - 1)) {
                myfile << indepVar.at(i) << ", " << "\n";
            } else {
                myfile << indepVar.at(i) << "\n";
            }
        }
        myfile << "))" << "\n";

        // Write each dependent variable
        for (int j = 0 ; j < depVar_names.size(); j++){



            myfile << depVar_names[j] << " = np.array((" << "\n";
            for (auto i = 0; i < depVars[j].size(); ++i) {
                if (i != (depVars[j].size() - 1)) {
                    myfile << depVars[j].at(i) << ", " << "\n";
                } else {
                    myfile << depVars[j].at(i) << "\n";
                }
            }
            myfile << "))" << "\n";


            myfile << depVarError_names[j] << " = np.array((" << "\n";
            for (auto i = 0; i < depVarErrors[j].size(); ++i) {
                if (i != (depVarErrors[j].size() - 1)) {
                    myfile << depVarErrors[j].at(i) << ", " << "\n";
                } else {
                    myfile << depVarErrors[j].at(i) << "\n";
                }
            }
            myfile << "))" << "\n";
        
        
        
        
            myfile << "plt.errorbar(" << indepVar_name << "," << depVar_names[j] << "," << depVarError_names[j] << ")"<< "\n";
        
        
        
        }

        myfile << "plt.show()"  << "\n";
        myfile.close();
    }
};

// Class for weapons. Internal variables are firingSpeed, firingHeight and shotSpread. Aim angle is set using takeAim();
class weapon{

    public:
    string name;
    double firingSpeed;
    double firingHeight;
    double shotSpread;
    double firingRate;
    double aimAngle;

    // Each weapon comes pre-loaded with a normal distribution
    normal_distribution<double> distribution;

    // Constructor sets the weapon's stats
    weapon(string name, double firingSpeed, double firingHeight, double shotSpread, double firingRate)
        : aimAngle(0), name(name), firingSpeed(firingSpeed), firingHeight(firingHeight), shotSpread(shotSpread), firingRate(firingRate), distribution(0, shotSpread) {}

    // Function to change internally stored angle
    void takeAim(target target){
        aimAngle = findOptimalAngle(target);
    }

    double sampleDistribution(mt19937 &rng){

        double sample = distribution(rng);

        return sample;
    }


    // Function to find the optimal angle for the weapon at a given distance with a given target
    double findOptimalAngle(target target){

        // Calculate the optimal angle. See Joe's notes for working out that this is the correct formula
        double tanTheta =  firingSpeed*firingSpeed/(g*target.distance) - sqrt( pow(firingSpeed,4)/pow(g*target.distance,2)  -  (1 + 2*pow(firingSpeed,2)*(target.centre-firingHeight)/(g*pow(target.distance,2)))  );
        double optimalAngle = atan(tanTheta);
        return optimalAngle;
    }

    // finding the height of the hit
    double actualHeight(weapon gun, target target, double inputAngle){   
        
        double height;
        double initial_height = gun.firingHeight;
        double v_0 = gun.firingSpeed;
        double d = target.distance;

        height = initial_height - g*d*d/(2*v_0*v_0) + d*tan(inputAngle) -g*d*d*tan(inputAngle)*tan(inputAngle)/(2*v_0*v_0);
    return height;
}

    // Returns the y value when the bullet is at the target, for a pre-aimed weapon
    double fireShot(target target, mt19937 &rng){

        
        double kick = distribution(rng);
        double angle = aimAngle + kick;
        
        //cout << "Aim angle: " << aimAngle << ", kick: " << kick << ", firing angle: " << angle << endl;

        // Calculate the y value when bullet is at target's x value using information from target and gun
        double yAtDistance = firingHeight + target.distance*tan(angle) - g*target.distance*target.distance*(1+tan(angle)*tan(angle))/(2*firingSpeed*firingSpeed);

        return yAtDistance;

    };


    // Function to find the optimal angle for the weapon at a given distance with a given target 
    double findOptimalAngle(double d, double yt){

        // Calculate the optimal angle. See Joe's notes for working out that this is the correct formula
        double tanTheta =  firingSpeed*firingSpeed/(g*d) - sqrt( pow(firingSpeed,4)/pow(g*d,2)  -  (1 + 2*pow(firingSpeed,2)*(yt-firingHeight)/(g*pow(d,2)))  );
        double thetaC = atan(tanTheta);
        return thetaC;
    }
};


// Changes the distance of the target object, and aims each gun at the new target. Doesn't work right now!
void changeDistance(vector<weapon> &guns, target target,double distance){

    target.distance = distance;
    for (int i = 0; i < guns.size(); i++){
        guns[i].takeAim(target);
    }

}



double average(vector<double> Q, bool showResult = false){ //code for finding the average of a vector
   int N = Q.size();
   double sum = 0;
   
   for (int i = 0; i < N; i++)
    {
        sum += Q[i];
    }
    double av = sum/N;

    if (showResult){
        cout << "average is " << av << endl;
    }
    
    return av;
}

double stand_dev(vector<double> Q, double mean, bool showResult = false){ //code for finding the standard deviation of a vector

   int N = Q.size();
   double var = 0;
   for (int i = 0; i < N; i++)
    {
       var = var += pow(Q[i]-mean,2)/N;
    }
   double sd = sqrt(var);
    if (showResult){
        cout << "standard deviation is " << sd << endl;
    }

   
   return sd;
   
}
double stand_dev(vector<double> Q, bool showResult = false){ //code for finding the standard deviation of a vector

    double mean = average(Q);

   int N = Q.size();
   double var = 0;
   for (int i = 0; i < N; i++)
    {
       var = var += pow(Q[i]-mean,2)/N;
    }
   double sd = sqrt(var);
    if (showResult){
        cout << "standard deviation is " << sd << endl;
    }

   
   return sd;
   
}





// A function that takes the array of shot locations from a trial and counts how many hit the target.
int countHits(vector<double> shotArray, target target){
    
    vector <int> hit;
    for (int i = 0; i<shotArray.size(); i++){

        if (shotArray[i] <= target.upper and shotArray[i] >= target.lower){
            hit.push_back(i);
        }
    }   
    return hit.size();
}


vector<double>  doTrial(weapon gun, target target, double timeLength, int &hitCount, mt19937&rng){


    // Use a number of shots based on each gun's firing rate and the time in a trial
    int shotNumber = floor(timeLength*gun.firingRate);
    //cout << gun.name << ": " << shotNumber << " shots" << endl;
    

    vector<double> shotLocations(shotNumber);

    // For each shot in the trial, get a y-value at the distance using fireShot()
    for (int i = 0; i<shotNumber; i++){

        shotLocations[i] = gun.fireShot(target,rng);



        if (shotLocations[i] <= target.upper and shotLocations[i] >= target.lower){
            hitCount++;
        }

    }
    return shotLocations;
}






int main(){


    mt19937 rng = initializeRandomGenerator();

    // Constants that won't change for this project:
    const double targetUpperBound = 2.3;
    const double targetLowerBound = 0.3;
    double distance = 100;


    target target(targetLowerBound,targetUpperBound,distance);

    // Declare weapon objects with their firing speed, firing height, and shot spread in degrees. The weapon class constructor will convert to radians.
    weapon musket("Musket",450,0.7,2*M_PI/180,3);
    weapon rifle("Rifle",600,0.3,1*M_PI/180,2);
    weapon sniper("Sniper rifle", 1000,0.9,0.1*M_PI/180,2.5);
    weapon baseball("Phalanx of baseball pitchers", 44,1,0*M_PI/180,400 );


    // Find optimal angle
    double thetaC = musket.findOptimalAngle(target);
    cout << thetaC << endl;

   
    int trialLength = 5; // Minutes
    int trialNumber = 1000000;



    
    vector <int>  musketTrialHits(trialNumber);
    vector <double>  musketTrialAccuracy(trialNumber);


    vector <int>  rifleTrialHits(trialNumber);
    vector <double>  rifleTrialAccuracy(trialNumber);


    vector<double> shots;
    int musketHitCount = 0;
    int rifleHitCount = 0;


    musket.takeAim(target);
    rifle.takeAim(target);
    sniper.takeAim(target);
    baseball.takeAim(target);

    vector<double> distances = {};
    linspace(distances,1,100,100);

    vector <double>  musketDistanceAccuracies(distances.size());
    vector <double>  musketDistanceErrors(distances.size());


    vector <double>  rifleDistanceAccuracies(distances.size());
    vector <double>  rifleDistanceErrors(distances.size());


    for (int j = 0; j < distances.size(); j++ ){

        cout << distances[j] <<"/" << "100" << endl;

        // Function to change where the target is, and where the weapons in the array are aiming (Doesn't work right now!)
        //changeDistance({&musket},target,distances[j]);

        target.distance = distances[j];
        musket.takeAim(target);
        rifle.takeAim(target);


        // Run a number of 5-minute trials at that distance 
        for (int i = 0; i < trialNumber; i++){





            vector<double> musketShots = doTrial(musket,target,trialLength,musketHitCount,rng);
            //cout << "Musket: " << musketHitCount << endl;
            musketTrialHits[i] = musketHitCount;
            //musketTrialAccuracy[i] = double(musketHitCount)/double(floor(trialLength*musket.firingRate));
            musketTrialAccuracy[i] = double(musketHitCount);
            
            musketHitCount = 0;


            // Something about the rifle is returning no results
            
            vector<double> rifleShots = doTrial(rifle,target,trialLength,rifleHitCount,rng);
            //cout << "Rifle: " << rifleHitCount << endl;
            rifleTrialHits[i] = rifleHitCount;
            //rifleTrialAccuracy[i] = double(rifleHitCount)/double(floor(trialLength*rifle.firingRate));
            rifleTrialAccuracy[i] = double(rifleHitCount);            
            rifleHitCount = 0;


        }

        musketDistanceAccuracies[j] = average(musketTrialAccuracy);
        musketDistanceErrors[j] = stand_dev(musketTrialAccuracy);
        

        rifleDistanceAccuracies[j] = average(rifleTrialAccuracy);
        rifleDistanceErrors[j] = stand_dev(rifleTrialAccuracy);


    }




    for (int i = 0; i < distances.size(); i++){
        cout << "Distance " << distances[i] << "m, " << musket.name << ": " << 100*musketDistanceAccuracies[i] << " +- " << 100*musketDistanceErrors[i] <<"%" << ", " << rifle.name << ": " << 100*rifleDistanceAccuracies[i] << " +- " << 100*rifleDistanceErrors[i] <<  endl;
    }
  


    WritePyData<int> writer("hitCounts.py");

    string name = "hitCounts";

    writer.write_out_hist_vector("my_array", musketTrialHits );




    WritePyData<double> writer2("weaponComparison.py");



    vector<string> varNames = {"musketAccuracies","rifleAccuracies"};
    vector<string> varErrorNames = {"musketErrors","rifleErrors"};
    vector<vector<double>> overallAccuracies = {musketDistanceAccuracies,rifleDistanceAccuracies};
    vector<vector<double>> overallErrors = {musketDistanceErrors,rifleDistanceErrors};
    


    writer2.write_out_comparison_vectors("distances", distances,varNames,overallAccuracies,varErrorNames,overallErrors);
    



    return 0;
}
