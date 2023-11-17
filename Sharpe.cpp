#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <time.h>
#include <ctime>

#include <bits/stdc++.h>
using namespace std;


class weapon{
    public:
    
    double firingSpeed;
    double firingHeight;

    weapon(double v0, double h1){
        firingSpeed = v0;
        firingHeight = h1;
    }


};

double findOptimalAngle(double v0, double g, double d, double yt, double h1){

    double thetaC = atan( pow(v0,2)/(g*d) - sqrt( pow(v0,4)/pow(g*d,2)  -  (1 + 2*pow(v0,2)/(g*pow(d,2)*(yt-h1)))  ));

    return thetaC;

}



int main(){

    double g = 9.81;
    double d = 100;
    double yt = 1.3;
    

    weapon musket(450,0.7);

    double v0 = musket.firingSpeed;
    double h1 = musket.firingHeight;




    double thetaC = findOptimalAngle(v0,g,d,yt,h1);


    cout << thetaC;




    return 0;
}