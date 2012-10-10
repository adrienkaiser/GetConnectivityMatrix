
/*
	./GetMatrix-build$ GetMatrix --connecMapsFile ../ConnecMaps.csv --labelMap /home/akaiser/Networking/Freesurfer/TestPenelope/subjects/Penelope/mri/aparc+aseg_resampled.nrrd --matrixFile ../Matrix.txt
*/
/* std classes */
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <sstream> // to convert int to std::string
#include <math.h> // for the sqrt and pow
#include <fstream> // to open a file

/* GenerateCLP */
#include "GetMatrixCLP.h" //generated when ccmake

/*itk classes*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include <itksys/SystemTools.hxx> // for GetParentDirectory() and GetFilenameLastExtension()

  /////////////////////////////////////////
 //           MAIN FUNCTION             //
/////////////////////////////////////////

int main (int argc, char *argv[])
{
	PARSE_ARGS; //thanks to this line, we can use the variables entered in command line as variables of the program
	// std::string : csv file ConnecMapsFile, map image LabelMap, output file MatrixFile


/* Test if all files are here */
	if( ConnecMapsFile.empty() )
	{
		std::cout<<"| The CSV file containing the connectivity maps for each label is missing."<<std::endl;
		return 0;
	}
	if( LabelMap.empty() )
	{
		std::cout<<"| The Label Map image is missing."<<std::endl;
		return 0;
	}
	if( MatrixFile.empty() )
	{
		std::cout<<"| The output Matrix file is missing. It will be created in the current work directory as \'./ConnectivityMatrix.txt\'."<<std::endl;
		MatrixFile = "./ConnectivityMatrix.txt";
	}

/* Test if files are readable */
	if( access(ConnecMapsFile.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
	{
		std::cout<<"| The CSV file containing the connectivity maps for each label is not readable : \'"<< ConnecMapsFile <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	if( itksys::SystemTools::GetFilenameLastExtension( ConnecMapsFile.c_str() ) != ".csv" ) // Test if the file is a csv
	{
		std::cout<<"| The given file containing the connectivity maps for each label is not a CSV file : \'"<< ConnecMapsFile <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	if( access(LabelMap.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
	{
		std::cout<<"| The Label Map image is not readable : \'"<< LabelMap <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	std::string MatrixFileFolder = itksys::SystemTools::GetParentDirectory( MatrixFile.c_str() );
	if( access(MatrixFileFolder.c_str(), W_OK) != 0 ) // Test if the folder is writable => unistd::access() returns 0 if W(write)_OK
	{
		std::cout<<"| The containing folder given for the output Matrix file is not writable : \'"<< MatrixFileFolder <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}

/* Test if the connectivity maps in the csv files are readable and count them */
	std::ifstream file (ConnecMapsFile.c_str() , std::ios::in);  // opening in reading
	if(! file) // error while opening
	{
		std::cout<<"| Error opening the file \'"<< ConnecMapsFile <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	else std::cout<<"| Successfully loaded : \'"<< ConnecMapsFile <<"\'."<<std::endl;

	int NbConnectMaps=0;
	std::string ConnecMap;
	while( std::getline ( file, ConnecMap, ',' ) ) 
	{
		if( access(ConnecMap.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
		{
			std::cout<<"| The following Connectivity map image is not readable : \'"<< ConnecMap <<"\'."<<std::endl;
			std::cout<<"=> ABORT <="<<std::endl;
			return -1;
		}
		NbConnectMaps ++ ;
	}
	file.close();

	std::cout<<"| "<< NbConnectMaps <<" Connectivity maps found."<<std::endl;

/* Open the files */

/*
	std::ofstream file (ConnecMapsFile.c_str() , std::ios::out | std::ios::trunc);  // opening in writing with erasing the open file
	if(! file) // error while opening
	{
		std::cout<<"| File can not be opened: Nothing will be saved"<<std::endl;
		return 0; // nothing will be saved
	}

	file << "Number of nodes in the network: "  << std::endl;

	file.close();
*/

/* End of Main function */

	return 0;
}

