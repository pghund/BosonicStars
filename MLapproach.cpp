using namespace std;
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include<string>
#include<iomanip>
#define pi 3.141592653589793


//double lambdaa = -2.903651591728209;
double lambdaa = -2;

double DescentParameter = 0.1;
double Tol = 0.05;
double e = pow(2, 0.5);
double mp = pow(0.2, 0.5);
double me = mp / 4;
double ma = mp;
double G = 1;
double h = 1;
double eps = (8 * me * e * e) / (h * h * pow(3, 2.0 / 3) * pow(pi, 1.0 / 3));
double eta = 4 * pi * e * e * 2 * ma / (h * h);
double FGrid[8001][8001];

string Estr = "ElectricGrid1C.txt";
string Gstr = "GravityGrid1C.txt";
string elecstr = "iteration1electronC.txt";
string alphastr = "iteration1alphaC.txt";
string gradstr = "gradient1.txt";


double f(double t, double ue, double ua, double sa1, double sa2, double sa3) {
    if (t == 0) {
        return -sa2 * (-sa2 / ua) - eta * (4 - G * ma * ma / (e * e)) * ua * ua * ua + eta * (2 + G * me * ma / (e * e)) * pow(ue, 3.0 / 2) * ua;
    }
    else if (ua == 0) {
        return 0;
    }

    else {
        return -sa3 * (4.0 / t - 2 * sa1 / (ua)) - sa2 * (-8.0 * sa1 / (t * ua) - sa2 / (ua)+2 * sa1 * sa1 / (ua * ua)) - sa1 * (4 * sa1 * sa1 / (t * ua * ua)) - eta * (4 - G * ma * ma / (e * e)) * ua * ua * ua + eta * (2 + G * me * ma / (e * e)) * pow(ue, 3.0 / 2) * ua;
    }
}

double g(double t, double ue, double ua, double se) {
    if (ue == 0) {
        return 0;
    }
    else if (t == 0) {
        return -(2 + G * me * ma / (e * e)) * eps * ua * ua + eps * (1 - G * me * me / (e * e)) * pow(ue, 3.0 / 2);
    }
    else {
        return -(2.0 / t) * se - (2 + G * me * ma / (e * e)) * eps * ua * ua + eps * (1 - G * me * me / (e * e)) * pow(ue, 3.0 / 2);
    }
}

double  F(double r1, double r2, int n) {
    if (r1 < r2) {
        return pow(r1, n + 2) / pow(r2, n + 1);
    }

    else if (r2 < r1) {
        return pow(r2, n) / pow(r1, n - 1);
    }
    else {
        return r1;
    }

}

int main() {

    ofstream logfile;
    double t = 0;
    double ue = 1;
    double ua = 1.0 / pow(2, 0.5);
    double step = 0.002;
    double r1;
    double sa1;
    double nsa1;
    double nua;
    int flag=0;

    
    int i = 0;
    double alphaGrid[8001];
    double electronGrid[8001];
    alphaGrid[0] = 1.0 / pow(2, 0.5);
    electronGrid[0] = 1;
    for (int j = 1;j < 8001; j++) {
        t = j * step;
        if (4 * sin(t / 4) / t < 0) {
            flag = 1;
        }
        if (flag == 0){
            electronGrid[j] = pow(4 * sin(t / 4) / t, 2.0 / 3);
            alphaGrid[j] = pow(2 * sin(t / 4) / t, 1.0 / 2);
        }
        else {
            electronGrid[j] = 0;
            alphaGrid[j] = 0;
        }

    }

    for (int j = 1;j < 8001; j++) {
        r1 = j * step;
        for (int i = 1;i < 8001; i++) {
            FGrid[j][i] = F(r1, i*step, 0);
        }

    }
    double sum;

    logfile.open(alphastr);
    for (int j = 0;j < 8001; j++) {
        logfile << j << ": ";
        logfile << j * step << " " << alphaGrid[j] * alphaGrid[j] << "\n";

    }
    logfile.close();

    logfile.open(elecstr);
    for (int j = 0;j < 8001; j++) {
        logfile << j << ": ";
        logfile << j * step << " " << pow(electronGrid[j], 3.0 / 2) << "\n";

    }
    logfile.close();

    //at this point should have initial guess, so we can start the iteration
    //we first need to make the electric and gravitational grids
    double D_G;
    double D_E;
    double r2;
    double GravityGrid[8001];
    double ElectricGrid[8001];
    double mGrid[8001];
    for (int i = 0; i < 8001; i++) {
        if (i == 0 || i == 8000) {
            mGrid[i] = 7;
        }
        else if (remainder(i, 4) == 0) {
            mGrid[i] = 14;
        }

        else if (remainder(i, 4) == 2) {
            mGrid[i] = 12;
        }

        else if (remainder(i, 2) == 1 || remainder(i, 2) == -1) {
            mGrid[i] = 32;
        }
    }


    for (int j = 0; j < 8001; j++) {
        D_G = 0;
        D_E = 0;
        r1 = j * step;
        for (int i = 0; i < 8001; i++) {
            r2 = i * step;

            if (i == 0 || i == 8000) {
                D_G += 7 * (me * pow(electronGrid[i], (3.0 / 2)) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                D_E += 7 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
                // an example of the ordering D_E += 7 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * F(r2, r1, 0);
            }

            else if (remainder(i, 4) == 0) {
                D_G += 14 * (me * pow(electronGrid[i], 3.0 / 2) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                D_E += 14 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
            }

            else if (remainder(i, 4) == 2) {
                D_G += 12 * (me * pow(electronGrid[i], 3.0 / 2) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                D_E += 12 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
            }

            else if (remainder(i, 2) == 1 || remainder(i, 2) == -1) {
                D_G += 32 * (me * pow(electronGrid[i], 3.0 / 2) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                D_E += 32 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
            }

        }
        GravityGrid[j] = 2 * pi * 2*step / 45 * D_G;
        ElectricGrid[j] = 2 * pi * 2*step / 45 * D_E;
    }

    logfile.open(Gstr);
    for (int j = 0;j < 8001; j++) {
        logfile << j << ": ";
        logfile << j * step << " " << GravityGrid[j] << "\n";

    }
    logfile.close();

    logfile.open(Estr);
    for (int j = 0;j < 8001; j++) {
        logfile << j << ": ";
        logfile << j * step << " " << ElectricGrid[j] << "\n";

    }
    logfile.close();

    //now we should have the grids calculated.So we calculate the new densities
    double newelectronGrid[8001];
    double newalphaGrid[8001];
    double lambdae;

    lambdae = me * GravityGrid[0] + e * e * ElectricGrid[0] - 1.0 / (2 * me / (h * h * pow(3, 2.0 / 3) * pow(pi, 4.0 / 3)));
    //cout << lambdae << "\n";
    flag = 0;
    for (int i = 0; i < 8001; i++) {
        if (flag == 0) {
            newelectronGrid[i] =  max((-lambdae + me * GravityGrid[i] + e * e * ElectricGrid[i]) * (2 * me / (h * h * pow(3, 2.0 / 3) * pow(pi, 4.0 / 3))), 0.0);
            if (newelectronGrid[i] == 0) {
                flag = 1;
            }
        }
        else {
            newelectronGrid[i] = 0;
        }
    }

    ua = 1.0 / pow(2, 0.5);
    newalphaGrid[0] = 1.0 / pow(2, 0.5);
    i = 1;
    sa1 = step * ua * (2 * ma / (h * h) * (lambdaa - ma * GravityGrid[0] + 2 * e * e * ElectricGrid[0]));
    t = step;
    //now for the alpha particles, we need to use an ODE solver
    // if we were to convert this to RK4, we would need to use h=0.002
    while (i < 8001) {
        if (ua + step * sa1 <= 0) {
            ua = 0;
            sa1 = 0;
        }
        nsa1 = sa1 + step * (-(2.0 / t) * sa1 + ua * 2 * ma / (h * h) * (lambdaa - ma * GravityGrid[i] + 2 * e * e * ElectricGrid[i]));
        nua = ua + step * sa1;
        ua = nua;
        sa1 = nsa1;
        if (ua <= 0) {
            ua = 0;
            sa1 = 0;
        }
        newalphaGrid[i] = ua;
        t = t + step;
        i += 1;
    }
    
    double du;
    double newdw;
    double dw;
    double newdwp;
    double dwp;
    int adjustedi;
    double prefactor1 = 2 * pi * 2 * step / 45;
    double prefactor2 = 2 * me * pow(3, 2.0 / 3) * pow(pi, 4.0 / 3) / (h * h);
    double prefactor3 = 2 * ma / (h * h);
    // at this point, we need to compute the gradient of the loss function so we have 8001 entries for both densities, and 1 entry for lambdae. is then gradient (electron, alpha, lambda)
    // in computing the gradient, I have removed mGrid since this artificially makes some terms more important than others.
    double Gradient[16003];
    for (int i = 0; i < 16003; i++) {
        sum = 0;
        if (i < 8001) {
            // this is for the first 8001 electron entries
            for (int j = 0; j < 8001; j++) {
                //first let us compute duj/dui=du and dwj/dui=dw
                if (newalphaGrid[j] == 0) {
                    du = 0;
                }
                else {
                    du = prefactor1*prefactor2* pow(electronGrid[i], 1.0 / 2)*(3.0/2)*mGrid[i] *(-me * me *FGrid[i][0] + e * e* FGrid[i][0] +me * me  * FGrid[i][j] - e * e * FGrid[i][j]);
                    
                }
                if (j == 0) {
                    newdw = 0;
                    newdwp = 0;
                }
                else if (j == 1) {
                    newdw = dw + step * dwp;
                    newdwp = step*alphaGrid[j - 1] * 2 * ma / (h * h) * (-ma * 2 * pi * 2 * step / 45 * (3.0 / 2)  * me * pow(electronGrid[i], 1.0 / 2) * FGrid[i][j-1] + 2 * e * e * 2 * 2 * pi * step / 45 * (-3.0 / 2)  * pow(electronGrid[i], 1.0 / 2) * FGrid[i][j-1]);

                }
                else {
                    newdw = dw + step * dwp;
                    newdwp = dwp * (1 + 2.0 / (j - 1)) + step*prefactor3*( dw *(lambdaa - ma * GravityGrid[j - 1] + 2 * e * e * ElectricGrid[j - 1]) + alphaGrid[j - 1] * pow(electronGrid[i], 1.0 / 2)*mGrid[i] * FGrid[i][j-1] *(3.0/2)*prefactor1*(-ma * me * -2*e*e));
                }

                if (i == j) {
                    sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (1 - du) - 2 * (alphaGrid[j] - newalphaGrid[j]) * newdw;
                }
                else {
                    sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (- du) - 2 * (alphaGrid[j] - newalphaGrid[j]) * newdw;
                }
                dw = newdw;
                dwp = newdwp;
            }
            Gradient[i] = sum*step;
        }
        else if (i < 16002) {
            adjustedi = i - 8001;
            //this is for the next 8001 alpha entries
            for (int j = 0; j < 8001; j++) {
                if (newalphaGrid[j] == 0) {
                    du = 0;
                }
                else {
                    du = prefactor2*prefactor1*alphaGrid[adjustedi]*mGrid[i]*(-me  * 2 * ma  * FGrid[adjustedi][0] - e * e * 4 * FGrid[adjustedi][0] + me * 2 * ma * FGrid[adjustedi][j] + e * e * 4 * FGrid[adjustedi][j]);
                }
                if (j == 0) {
                    newdw = 0;
                    newdwp = 0;
                }
                else if (j == 1) {
                    newdw = dw + step * dwp;
                    newdwp = step* alphaGrid[j - 1] * 2 * ma / (h * h) * (-ma * 2 * pi * 2 * step / 45 * 2  * ma * alphaGrid[adjustedi] * FGrid[adjustedi][j-1] + 2 * e * e * 2 * 2 * pi * step / 45 * (4)  * alphaGrid[adjustedi] * FGrid[adjustedi][j-1]);

                }
                else {
                    newdw = dw + step * dwp;
                    
                    newdwp = dwp * (1 + 2.0 / (j - 1)) + step*prefactor3*(dw  * (lambdaa - ma * GravityGrid[j - 1] + 2 * e * e * ElectricGrid[j - 1]) + alphaGrid[j - 1]  *prefactor1*alphaGrid[adjustedi]* mGrid[i]*FGrid[adjustedi][j-1] *(-ma  * 2  * ma  + 8 * e * e ));
                }

                if (adjustedi == j) {
                    sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (- du) + 2 * (alphaGrid[j] - newalphaGrid[j]) *(1 -newdw);
                }
                else {
                    sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (-du) + 2 * (alphaGrid[j] - newalphaGrid[j]) *( -newdw);
                }
                dw = newdw;
                dwp = newdwp;      
            }
            Gradient[i] = sum*step;
        }
        else {
            //this is for the lambdaa entry
            for (int j = 0; j < 8001; j++) {
                if (j == 0) {
                    newdw = 0;
                    newdwp = 0;
                }
                else if (j == 1) {
                    newdw = dw + step * dwp;
                    newdwp = step*2 * ma / (h * h) * (dw * lambdaa + alphaGrid[j - 1]);

                }
                else {
                    newdw = dw + step * dwp;
                    newdwp = dwp * (1 + 2.0 / (j - 1)) + step*prefactor3 * (dw *( lambdaa-ma*GravityGrid[j-1]+2*e*e*ElectricGrid[j-1]) + alphaGrid[j - 1]);
                }

                sum += 2 * (alphaGrid[j] - newalphaGrid[j]) * (-newdw);

                dw = newdw;
                dwp = newdwp;
            }
        
            Gradient[i] = sum*step;
        }
    }
    
    // seems like I need to normalize the gradient
    double most=0;
    int loc=0;
    for (int j = 0; j < 16002; j++) {
        if (abs(Gradient[j]) > most) {
            most = abs(Gradient[j]);
            loc = j;
        }

    }
    logfile.open(gradstr);
    for (int j = 0; j < 16003; j++) {
        logfile << Gradient[j] << "\n";

    }
    logfile.close();

    cout << most << " " << loc << "\n";

    for (int j = 0; j < 16002; j++) {
        Gradient[j] = Gradient[j]/most;
    }
    double error = 0;
    for (int j = 0; j < 8001; j++) {
        error += (alphaGrid[j] - newalphaGrid[j]) * (alphaGrid[j] - newalphaGrid[j]) + (electronGrid[j] - newelectronGrid[j]) * (electronGrid[j] - newelectronGrid[j]);
    }
    cout << error*0.002<<"\n";
    
    double intermediateelectronGrid[8001];
    int iter = 1;
    char Charcount;
    double olderror;
    while (iter < 10 && error>0.1) {
        olderror = error;
        iter += 1;
        flag = 0;
        Charcount = '0' + iter;
        alphastr[9] = Charcount;
        elecstr[9] = Charcount;
        Estr[12] = Charcount;
        Gstr[11] = Charcount;
        gradstr[8] = Charcount;
        cout << iter << "\n";
        for (int j = 0; j < 8001; j++) {
            intermediateelectronGrid[j] =electronGrid[j]-DescentParameter*Gradient[j];
            //electronGrid[j] -= DescentParameter * Gradient[j];
            alphaGrid[j] -= DescentParameter*Gradient[j+8001] ;
            if (intermediateelectronGrid[j] < 0.000001) {
                intermediateelectronGrid[j] = 0;
            }
            if (alphaGrid[j] < 0.000001) {
                alphaGrid[j] = 0;
            }
        }
        
        for (int j = 0; j < 8001; j++) {
            if (j < 20) {
                electronGrid[j] = intermediateelectronGrid[j];
            }
            else if(j<7500) {
                sum = 0;
                for (int i = -15; i < 16; i++) {
                    sum += intermediateelectronGrid[j + i];
                }
                electronGrid[j] = sum / 31;
            }
            else {
                electronGrid[j] = intermediateelectronGrid[j];
            }
        }
        
        //lambdaa -= DescentParameter*0.1 * Gradient[16002];
        //cout << lambdaa << "\n";
        
        logfile.open(alphastr);
        for (int j = 0;j < 8001; j++) {
            logfile << j << ": ";
            logfile << j * step << " " << alphaGrid[j] * alphaGrid[j] << "\n";

        }
        logfile.close();

        logfile.open(elecstr);
        for (int j = 0;j < 8001; j++) {
            logfile << j << ": ";
            logfile << j * step << " " << pow(electronGrid[j], 3.0 / 2) << "\n";

        }
        logfile.close();
        

        for (int j = 0; j < 8001; j++) {
            D_G = 0;
            D_E = 0;
            r1 = j * step;
            for (int i = 0; i < 8001; i++) {
                r2 = i * step;

                if (i == 0 || i == 8000) {
                    D_G += 7 * (me * pow(electronGrid[i], (3.0 / 2)) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                    D_E += 7 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
                    // an example of the ordering D_E += 7 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * F(r2, r1, 0);
                }

                else if (remainder(i, 4) == 0) {
                    D_G += 14 * (me * pow(electronGrid[i], 3.0 / 2) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                    D_E += 14 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
                }

                else if (remainder(i, 4) == 2) {
                    D_G += 12 * (me * pow(electronGrid[i], 3.0 / 2) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                    D_E += 12 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
                }

                else if (remainder(i, 2) == 1 || remainder(i, 2) == -1) {
                    D_G += 32 * (me * pow(electronGrid[i], 3.0 / 2) + ma * pow(alphaGrid[i], 2)) * FGrid[i][j];
                    D_E += 32 * (2 * pow(alphaGrid[i], 2) - pow(electronGrid[i], 3.0 / 2)) * FGrid[i][j];
                }

            }
            GravityGrid[j] =  2 * pi * 2*step / 45 * D_G ;
            ElectricGrid[j] =  2 * pi * 2*step / 45 * D_E;
        }
        /*
        logfile.open(Gstr);
        for (int j = 0;j < 8001; j++) {
            logfile << j << ": ";
            logfile << j * step << " " << GravityGrid[j] << "\n";

        }
        logfile.close();

        logfile.open(Estr);
        for (int j = 0;j < 8001; j++) {
            logfile << j << ": ";
            logfile << j * step << " " << ElectricGrid[j] << "\n";

        }
        logfile.close();
        */


        lambdae = me * GravityGrid[0] + e * e * ElectricGrid[0] - 1.0 / (2 * me / (h * h * pow(3, 2.0 / 3) * pow(pi, 4.0 / 3)));
        for (int i = 0; i < 8001; i++) {
            if (flag == 0) {
                newelectronGrid[i] =  max((-lambdae + me * GravityGrid[i] + e * e * ElectricGrid[i]) * (2 * me / (h * h * pow(3, 2.0 / 3) * pow(pi, 4.0 / 3))), 0.0);
                if (newelectronGrid[i] == 0) {
                    flag = 1;
                }
            }
            else {
                newelectronGrid[i] = 0;
            }
        }


        ua = 1.0 / pow(2, 0.5);
        newalphaGrid[0] = 1.0 / pow(2, 0.5);
        i = 1;
        sa1 = step * ua * (2 * ma / (h * h) * (lambdaa - ma * GravityGrid[0] + 2 * e * e * ElectricGrid[0]));
        t = step;
        //now for the alpha particles, we need to use an ODE solver

        while (i < 8001) {
            if (ua + step * sa1 <= 0) {
                ua = 0;
                sa1 = 0;
            }
            nsa1 = sa1 + step * (-(2.0 / t) * sa1 + ua * 2 * ma / (h * h) * (lambdaa - ma * GravityGrid[i] + 2 * e * e * ElectricGrid[i]));
            nua = ua + step * sa1;
            ua = nua;
            sa1 = nsa1;
            if (ua <= 0) {
                ua = 0;
                sa1 = 0;
            }
            newalphaGrid[i] =  ua;
            t = t + step;
            i += 1;
        }
        for (int i = 0; i < 16003; i++) {
            sum = 0;
            if (i < 8001) {
                // this is for the first 8001 electron entries
                for (int j = 0; j < 8001; j++) {
                    //first let us compute duj/dui=du and dwj/dui=dw
                    if (newalphaGrid[j] == 0) {
                        du = 0;
                    }
                    else {
                        du = prefactor1 * prefactor2  * pow(electronGrid[i], 1.0 / 2) * (3.0 / 2) *mGrid[i]* (-me * me * FGrid[i][0] + e * e * FGrid[i][0] + me * me * FGrid[i][j] - e * e * FGrid[i][j]);

                    }
                    if (j == 0) {
                        newdw = 0;
                        newdwp = 0;
                    }
                    else if (j == 1) {
                        newdw = dw + step * dwp;
                        newdwp = step* alphaGrid[j - 1] * 2 * ma / (h * h) * (-ma * 2 * pi * 2 * step / 45 * (3.0 / 2)  * me * pow(electronGrid[i], 1.0 / 2) * FGrid[i][j - 1] + 2 * e * e * 2 * 2 * pi * step / 45 * (-3.0 / 2) * pow(electronGrid[i], 1.0 / 2) * FGrid[i][j - 1]);

                    }
                    else {
                        newdw = dw + step * dwp;
                        newdwp = dwp * (1 + 2.0 / (j - 1)) + step* prefactor3 * (dw * (lambdaa - ma * GravityGrid[j - 1] + 2 * e * e * ElectricGrid[j - 1]) + alphaGrid[j - 1]  * pow(electronGrid[i], 1.0 / 2) * mGrid[i]*FGrid[i][j - 1] * (3.0 / 2) * prefactor1 * (-ma * me * -2 * e * e));
                    }

                    if (i == j) {
                        sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (1 - du) - 2 * (alphaGrid[j] - newalphaGrid[j]) * newdw;
                    }
                    else {
                        sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (-du) - 2 * (alphaGrid[j] - newalphaGrid[j]) * newdw;
                    }
                    dw = newdw;
                    dwp = newdwp;
                }
                Gradient[i] = sum*step;
            }
            else if (i < 16002) {
                adjustedi = i - 8001;
                //this is for the next 8001 alpha entries
                for (int j = 0; j < 8001; j++) {
                    if (newalphaGrid[j] == 0) {
                        du = 0;
                    }
                    else {
                        du = prefactor2 * prefactor1 * alphaGrid[adjustedi] *mGrid[i] * (-me * 2 * ma * FGrid[adjustedi][0] - e * e * 4 * FGrid[adjustedi][0] + me * 2 * ma * FGrid[adjustedi][j] + e * e * 4 * FGrid[adjustedi][j]);
                    }
                    if (j == 0) {
                        newdw = 0;
                        newdwp = 0;
                    }
                    else if (j == 1) {
                        newdw = dw + step * dwp;
                        newdwp = step*alphaGrid[j - 1] * 2 * ma / (h * h) * (-ma * 2 * pi * 2 * step / 45 * 2 * ma * alphaGrid[adjustedi] * FGrid[adjustedi][j - 1] + 2 * e * e * 2 * 2 * pi * step / 45 * (4)  * alphaGrid[adjustedi] * FGrid[adjustedi][j - 1]);

                    }
                    else {
                        newdw = dw + step * dwp;
                        newdwp = dwp * (1 + 2.0 / (j - 1)) + step*prefactor3 * (dw * (lambdaa - ma * GravityGrid[j - 1] + 2 * e * e * ElectricGrid[j - 1]) + alphaGrid[j - 1] * prefactor1  * alphaGrid[adjustedi]*mGrid[i] * FGrid[adjustedi][j - 1] * (-ma * 2 * ma + 8 * e * e));
                    }

                    if (adjustedi == j) {
                        sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (-du) + 2 * (alphaGrid[j] - newalphaGrid[j]) * (1 - newdw);
                    }
                    else {
                        sum += 2 * (electronGrid[j] - newelectronGrid[j]) * (-du) + 2 * (alphaGrid[j] - newalphaGrid[j]) * (-newdw);
                    }
                    dw = newdw;
                    dwp = newdwp;
                }
                Gradient[i] = sum*step;
            }

            else {
                //this is for the lambdaa entry
                for (int j = 0; j < 8001; j++) {
                    if (j == 0) {
                        newdw = 0;
                        newdwp = 0;
                    }
                    else if (j == 1) {
                        newdw = dw + step * dwp;
                        newdwp = step*2 * ma / (h * h) * (dw * lambdaa + alphaGrid[j - 1]);

                    }
                    else {
                        newdw = dw + step * dwp;
                        newdwp = dwp * (1 + 2.0 / (j - 1)) + step * prefactor3 * (dw * (lambdaa - ma * GravityGrid[j - 1] + 2 * e * e * ElectricGrid[j - 1]) + alphaGrid[j - 1]);
                    }

                    sum += 2 * (alphaGrid[j] - newalphaGrid[j]) * (-newdw);

                    dw = newdw;
                    dwp = newdwp;
                }

                Gradient[i] = sum*step;
            }
        }
        most = 0;
        loc = 0;
        
        logfile.open(gradstr);
        for (int j = 0; j < 16002; j++) {
            logfile << Gradient[j] << "\n";
            if (abs(Gradient[j]) > most) {
                most = abs(Gradient[j]);
                loc = j;
            }

        }
        //logfile.close();
        

        cout << most << " " << loc << "\n";
        error = 0;
        for (int j = 0; j < 16002; j++) {
            Gradient[j] = Gradient[j] / most;
        }


        for (int j = 0; j < 8001; j++) {
            error += (alphaGrid[j] - newalphaGrid[j]) * (alphaGrid[j] - newalphaGrid[j]) + (electronGrid[j] - newelectronGrid[j]) * (electronGrid[j] - newelectronGrid[j]);
        }
        cout << error*0.002 << "\n";
        if (error > olderror) {
            DescentParameter = DescentParameter / 2;
            cout << "half-step" << "\n";
        }

    }
}


