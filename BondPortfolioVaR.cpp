#include <iostream>
#include <cmath>
#include <numeric>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

using namespace std;

/* INPUT. As in examples 1 & 2 in A Primer on Bond Portfolio Value at Risk - Smith, i.e.
          Bond data and correlation matrix obtain from historical data.

   OUTPUT. The output.csv file containing modified durations and single VaRs for each bond, as 
           well as the Undiversified and Bond Portfolio VaR for the full bond Portfolio, using 
           Yield or Price historical data. 

*/

// Helper functions

// Get bonds from a csv file named bonds.csv 
vector<string> loadCSV(string filename) {

	// File pointer
	fstream fin;

	// Open the file
	fin.open(filename, ios::in);

	vector<string> data;
	string line, word;

    // Extract data from the file
	while (getline(fin, line)) {
		stringstream s(line);
		while (getline(s, word, ',')) {
			data.push_back(word);
		}
	}
    return data;
}

// Write ouput to a csv file named output.csv
void writeOutputToCSV(double individualVaRsFromYield[], double individualVaRsFromPrice[], double MD[], double UVaRYield, double PVaRYield, double UVaRPrice, double PVaRPrice, int numberOfBounds) {
	
    // file pointer
	fstream fout;

	// opens an existing csv file or creates a new file.
	fout.open("output.csv", ios::out);

    // Insert the data to file
    fout << "Single Bond VaR from Yield, Single Bond VaR from Price, Modified Duration\n";
    for(int i = 0; i< numberOfBounds; i++){
        fout << individualVaRsFromYield[i] << ", " << individualVaRsFromPrice[i] << ", " << MD[i] << "\n";
    }
    fout << "\n" 
    << "Undiversified VaR for the full Bond Portfolio from Yield, Bond Portfolio VaR from Yield" << "\n"
    << UVaRYield << ", " << PVaRYield
    << "\n\n"
    << "Undiversified VaR for the full Bond Portfolio from Price, Bond Portfolio VaR from Price" << "\n"
    << UVaRPrice << ", " << PVaRPrice;
}

// Compute the modified duration of a bond
double getModifiedDuration(double yield, double couponRate, int yearsToMaturity) {
    double macaulayDuration = 0;
    double modifiedDuration = 0;
    if(couponRate == 0){ // improve accuracy by simplifying formula if zero-coupon bond
        macaulayDuration = yearsToMaturity;
        modifiedDuration = (double) macaulayDuration/(1+yield);
    } else { // coupon bond general formula
        macaulayDuration = (double) (1+yield)/yield - (double)(1+yield + yearsToMaturity*(couponRate - yield))/(couponRate*(pow(1+yield, yearsToMaturity)-1)+yield);
        modifiedDuration = (double) macaulayDuration/(1+yield);
    }
    return modifiedDuration;
}

// Compute VaR of a bond from yield
double getSingleBondVaRfromYield(double MV, double MD, double yield, double yieldVolatility, double zScore) {
    double bondVaR = MV*MD*yield*yieldVolatility*zScore;
    return bondVaR;
}

// Compute VaR of a bond from price
double getSingleBondVaRfromPrice(double SD, double zScore) {
    double bondVaR = SD*zScore;
    return bondVaR;
}

// Compute the Undiversified VaR of the bond portfolio
double getUndiversifiedVaR(double individualBondVaRs[], int numberOfBonds){
    double UVaR = 0;
    for(int i=0; i<numberOfBonds; i++) {
        UVaR+= individualBondVaRs[i];
    }
    return UVaR;
}

// Compute the Bond Portfolio VaR of our bond portfolio
double getBondPortfolioVaR(double individualBondVaRs[], double *corMat[], int numberOfBonds) {
    double bondPortfolioVaR=0;
    for(int i=0;i < numberOfBonds; i++){
        for(int j=0;j < numberOfBonds; j++){
            bondPortfolioVaR+= individualBondVaRs[i]*individualBondVaRs[j]*corMat[i][j];
        }
    }
    return sqrt(bondPortfolioVaR);
}

/* The main function calls the helper functions above to read the INPUT from the csv files "bonds.csv" and "corrMatrix.csv",
 to produce the output, and to store it in "output.csv".
*/
int main(void) 
{

    // Get input from csv files "bonds.csv" and "corrMatrix.csv"
    vector<string> data = loadCSV("bonds.csv");
    vector<string> corrMatrix = loadCSV("corrMatrix.csv");

    double zScore = 2.326;
    int numberOfBonds = data.size()/11-1; // positive integer

    int maturity[numberOfBonds];; // years
    double parValue[numberOfBonds]; // dollars
    double parPrice[numberOfBonds]; // percentages
    double marketValue[numberOfBonds]; // dollars 
    double yield[numberOfBonds]; // =2*(yield per period); unit of years 
    double yieldVolatility[numberOfBonds]; // percentages
    double sharePrice[numberOfBonds]; // dollars
    double numberOfShares[numberOfBonds]; // positive integer
    double shareOfMV[numberOfBonds]; // percentage
    double dailySD[numberOfBonds]; // dollars
    double couponRate[numberOfBonds]; // percentages
    double modifiedDuration[numberOfBonds] ; // years 
    double *correlationMatrix[numberOfBonds]; // matrix of correlations of changes in yield
    for(int i = 0; i < numberOfBonds; i++) {
        correlationMatrix[i] = new double[numberOfBonds];
    }
    
    // Populate data variables
    for(int i = 0; i < numberOfBonds; i++) {
        maturity[i]= stod(data[11+i*11]);
        parValue[i]= stod(data[11+i*11+1]);
        parPrice[i]= stod(data[11+i*11+2]);
        marketValue[i]= stod(data[11+i*11+3]);
        yield[i]= stod(data[11+i*11+4]);
        yieldVolatility[i]= stod(data[11+i*11+5]);
        sharePrice[i]= stod(data[11+i*11+6]);
        numberOfShares[i]= stod(data[11+i*11+7]);
        shareOfMV[i]= stod(data[11+i*11+8]);
        dailySD[i]= stod(data[11+i*11+9]);
        couponRate[i]= stod(data[11+i*11+10]);
        for(int j = 0; j < numberOfBonds; j++) {
            correlationMatrix[i][j] = stod(corrMatrix[i*numberOfBonds+j]);
        }
    }

    // Preliminary computations 

    double individualBondVaRsFromYield[numberOfBonds];
    double individualBondVaRsFromPrice[numberOfBonds];

    // Compute modified durations of bonds
    for(int i=0; i< numberOfBonds; i++){
        modifiedDuration[i] = getModifiedDuration(yield[i], couponRate[i], maturity[i]);
    }

    // Compute individual bond VaR from yield and from price
    for(int i=0; i< numberOfBonds; i++){
        individualBondVaRsFromYield[i]=getSingleBondVaRfromYield(marketValue[i], modifiedDuration[i], yield[i], yieldVolatility[i], zScore);
        individualBondVaRsFromPrice[i]=getSingleBondVaRfromPrice(dailySD[i], zScore);
    }

    // Store OUPUT in csv file "output.csv"

    double UVaRYield = getUndiversifiedVaR(individualBondVaRsFromYield, numberOfBonds);
    double PVaRYield = getBondPortfolioVaR(individualBondVaRsFromYield, correlationMatrix, numberOfBonds);
    double UVaRPrice = getUndiversifiedVaR(individualBondVaRsFromPrice, numberOfBonds);
    double PVaRPrice = getBondPortfolioVaR(individualBondVaRsFromPrice, correlationMatrix, numberOfBonds);

    writeOutputToCSV(individualBondVaRsFromYield, individualBondVaRsFromPrice, modifiedDuration, UVaRYield, PVaRYield, UVaRPrice, PVaRPrice, numberOfBonds);
}