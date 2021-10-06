using namespace std;
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#define pi 3.141592653589793 

double e = pow(2, 0.5);

double mp = pow(0.2, 0.5);
double me = mp / 4;
double ma = mp;
//double ma = me / (pow(3, 2.0 / 3) * pow(pi, 4.0 / 3));
double G = 1;
double h = 1;

double eps = (8 * me * e * e) / (h * h * pow(3, 2.0 / 3) * pow(pi, 1.0 / 3));
double eta = 4 * pi * e * e * 2 * ma / (h * h);


double f(double t, double ue, double ua, double sa1, double sa2, double sa3) {
    if (t == 0)
    {
        return -sa2 * (-sa2 / ua) - eta * (4 - G * pow(ma, 2) / pow(e, 2)) * ua * ua * ua + eta * (2 + G * me * ma / pow(e, 2)) * pow(ue, 3.0 / 2) * ua;
    }
    else if (ua == 0) {
        return 0;
    }
    else
    {
        return -sa3 * (4.0 / t - 2 * sa1 / (ua)) - sa2 * (-8.0 * sa1 / (t * ua) - sa2 / (ua)+2 * sa1 * sa1 / (ua * ua)) - sa1 * (4 * sa1 * sa1 / (t * ua * ua)) - eta * (4 - G * pow(ma, 2) / pow(e, 2)) * ua * ua * ua + eta * (2 + G * me * ma / pow(e, 2)) * pow(ue, 3.0 / 2) * ua;
        //return -sa3 * (4.0 / t - 2 * sa1 / (ua)) - sa2 * (-8.0 * sa1 / (t * ua) - sa2 / (ua)+2 * sa1 * sa1 / (ua * ua)) + sa1 * (4 * sa1 * sa1 / (t * ua * ua)) + eta * (4 - G * pow(ma, 2) / pow(e, 2)) * ua * ua * ua - eta * (2 + G * me * ma / pow(e, 2)) * pow(ue, 3.0 / 2) * ua;
    }
}
double g(double t, double ue, double ua, double se) {
    if (ue == 0)
    {
        return 0;
    }
    else if (t == 0) {
        return -(2 + G * me * ma / pow(e, 2)) * eps * ua * ua + eps * (1 - G * pow(me, 2) / pow(e, 2)) * pow(ue, 3.0 / 2);
    }

    else
    {
        return -(2.0 / t) * se - (2 + G * me * ma / pow(e, 2)) * eps * ua * ua + eps * (1 - G * pow(me, 2) / pow(e, 2)) * pow(ue, 3.0 / 2) ;
    }
}


int main() {
    //cout << eta << " " << eps << "\n";
    ofstream logfile("logfileEuler.txt");
    //Inputs
    //'time'
    double t;
    // function values. we convert the second order system to a first order system we need four variables.
    //The v's track the densities, the w's their derivatives.

    double ue;
    double se;
    double nue;
    double nse;

    double ua;
    double nua;
    double sa1;
    double sa2;
    double sa3;
    double nsa3;
    double nsa2;
    double nsa1;

    //Begin iteration
    // 1,1 -0.193
    //1,0.9, -0.1170092
    // 1, 1.361, 0
    t = 0;
    ue = pow(0.19, 2.0/3);
    ua = 1;
    double step = 0.00000001;
    se = 0;
    sa1 =0;
    sa2 = 0;
    sa3 = 0;
    int flag = 0;
    int i = 0;
    double Ne;
    double Na;

    do {

        if (ue + step * se <= 0) { ue = 0; se = 0; }
        else if (ua + step * sa1 <= 0) { ua = 0; sa1 = 0; sa2 = 0; sa3 = 0; }
        if (isnan(ua)) { break; }

        //Na+=(step*ua*ua*0.5+step*0.5*pow(ua+step*sa1,2))*t*t;
        //Ne+=(step*pow(ue, 3.0/2)*0.5+step*0.5*pow(ue+step*se,3.0/2))*t*t;

        nsa3 = sa3 + step * f(t, ue, ua, sa1, sa2, sa3);
        nsa2 = sa2 + step * sa3;
        nsa1 = sa1 + step * sa2;
        nse = se + step * g(t, ue, ua, se);

        nua = ua + step * sa1;
        nue = ue + step * se;
        sa3 = nsa3;
        sa2 = nsa2;
        sa1 = nsa1;
        ua = nua;
        se = nse;
        ue = nue;

        if (ua <= 0) { ua = 0; sa1 = 0; sa2 = 0; sa3 = 0;}


        t = t + step;
        i += 1;

        if (logfile.is_open() && i % 100000 == 0 && t>0)
            logfile << t << " " << ua << " " << pow(ue, 3.0 / 2)  << " " << se << " " << sa1 << "\n";
    } while (ue < 100  && i < 1000000000);
    //cout << pow(ue, 3.0 / 2) << " "<< se << "\n";
}