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
#include "itkImageRegionIterator.h"
#include <itksys/SystemTools.hxx> // for GetParentDirectory() and GetFilenameLastExtension()

  /////////////////////////////////////////
 //          LABEL STRUCTURE            //
/////////////////////////////////////////

struct Label
{
	int LabelID;
	std::vector< double > Values; // contains all the values in each label
};

  /////////////////////////////////////////
 //           ITK IMAGE TYPES           //
/////////////////////////////////////////
	typedef itk::Image < double , 3 > ImageType; //itk type for image
	typedef itk::ImageFileReader <ImageType> ReaderType; //itk reader class to open an image
	typedef itk::ImageRegionIterator< ImageType > IteratorImageType;

  /////////////////////////////////////////
 //            TEST FILES               //
/////////////////////////////////////////

int testFiles (std::string ConnecMapsFile, std::string LabelMap, std::string MatrixFile) // csv file ConnecMapsFile, map image LabelMap, output file MatrixFile : returns -1 if problem, otherwise returns the nb of connec maps
{
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

/* Test if the connectivity maps in the csv file are readable and count them */
	std::ifstream CSVfileStream (ConnecMapsFile.c_str() , std::ios::in);  // opening in reading
	if(! CSVfileStream) // error while opening
	{
		std::cout<<"| Error opening the Connectivity maps file \'"<< ConnecMapsFile <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	else std::cout<<"| Successfully loaded : \'"<< ConnecMapsFile <<"\'."<<std::endl;

	int NbConnectMaps=0;
	std::string ConnecMap;
	while( std::getline ( CSVfileStream, ConnecMap, ',' ) ) 
	{
		if( access(ConnecMap.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
		{
			std::cout<<"| The following Connectivity map image is not readable : \'"<< ConnecMap <<"\'."<<std::endl;
			std::cout<<"=> ABORT <="<<std::endl;
			//return -1;
		}
		NbConnectMaps ++ ;
	}
	CSVfileStream.close();

	std::cout<<"| "<< NbConnectMaps <<" Connectivity maps found."<<std::endl;
	return NbConnectMaps;

}

  /////////////////////////////////////////
 //           GET NB LABELS             //
/////////////////////////////////////////

int GetNbLabels(std::string LabelMap) // returns -1 if failed to load image
{
/* Open the Label Map */
	ImageType::Pointer LabelMapImage;
	try
	{
		ReaderType::Pointer reader=ReaderType::New();
		reader->SetFileName( LabelMap );
		reader->Update();
		LabelMapImage = reader->GetOutput();
	}
	catch( itk::ExceptionObject exception )
	{
		std::cerr << exception << std::endl ;
		return -1 ;
	}
	std::cout<<"| Successfully loaded image : \'"<< LabelMap <<"\'."<<std::endl;

/* Main Loop : test all the voxels and count labels */
	std::vector< int > Labels; // contains all the labels already found : NO duplicate
	IteratorImageType LabelIter ( LabelMapImage , LabelMapImage->GetLargestPossibleRegion() );

	for( LabelIter.GoToBegin() ; !LabelIter.IsAtEnd() ; ++LabelIter )
	{
		int value = LabelIter.Get();
		if( value!=0 && std::find(Labels.begin(), Labels.end(), value) == Labels.end() ) Labels.push_back(value);
	}

/* Return */
	return Labels.size();
}

  /////////////////////////////////////////
 //            COMPARE SIZE             //
/////////////////////////////////////////

int CompareSize(std::string Image1, std::string Image2) // returns 1 if same size, 0 if different
{
/* Variables definitions */
	ReaderType::Pointer reader1=ReaderType::New();
	ImageType::RegionType region1;
	reader1->SetFileName( Image1 );
	reader1->UpdateOutputInformation(); // get the informations in the header
	region1 = reader1->GetOutput()->GetLargestPossibleRegion();

	ReaderType::Pointer reader2=ReaderType::New();
	ImageType::RegionType region2;
	reader2->SetFileName( Image2 );
	reader2->UpdateOutputInformation(); // get the informations in the header
	region2 = reader2->GetOutput()->GetLargestPossibleRegion();

/* Test sizes */
	if( (int)region1.GetSize()[0] != (int)region2.GetSize()[0] ) return 0; // x coordinate
	if( (int)region1.GetSize()[1] != (int)region2.GetSize()[1] ) return 0; // y coordinate
	if( (int)region1.GetSize()[2] != (int)region2.GetSize()[2] ) return 0; // z coordinate

	return 1;
}

  /////////////////////////////////////////
 //      FOR SORT : COMPARE LABELS      //
/////////////////////////////////////////

bool CompareLabels (const Label &a, const Label &b) // if = -> return 0 | if a<b return <0 | if a>b return >0
{
	return  a.LabelID < b.LabelID ;
}
	
  /////////////////////////////////////////
 //         METRICS COMPUTATIONS        //
/////////////////////////////////////////

std::vector< double > GetMeans (std::vector< Label > Labels)
{
	std::vector< double > Means;

	for(unsigned int i=0;i<Labels.size();i++) // for all labels in Labels
	{
		double Mean;
		double Sum=0;

		for(unsigned int j=0;j<Labels[i].Values.size();j++) Sum = Sum + Labels[i].Values[j]; // sum all the values in the label
		Mean = Sum / Labels[i].Values.size(); // Mean = sum of the values for all the voxels in the label divided by nb of voxels in the label

		Means.push_back(Mean);
	}

	return Means;
}

  /////////////////////////////////////////
 //           MAIN FUNCTION             //
/////////////////////////////////////////

int main (int argc, char *argv[])
{
	PARSE_ARGS; //thanks to this line, we can use the variables entered in command line as variables of the program
	// std::string : csv file ConnecMapsFile, map image LabelMap, output file MatrixFile

/* Test the files given by user */
	int NbConnectMaps = testFiles(ConnecMapsFile, LabelMap, MatrixFile);
	if(NbConnectMaps==-1) return -1; // info display in the function

/* Test if right number of labels */
	if( GetNbLabels(LabelMap) != NbConnectMaps )
	{
		std::cout<<"| The number of labels in Label Map \'"<< LabelMap <<"\' is not the same as the number of connectivity maps in the csv file \'"<< ConnecMapsFile <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1; // nothing will be saved
	}

/* Open the matrix file */
	std::ofstream MatrixfileStream (MatrixFile.c_str() , std::ios::out | std::ios::trunc);  // opening in writing with erasing the open file
	if(! MatrixfileStream) // error while opening
	{
		std::cout<<"| Error creating the matrix file \'"<< MatrixfileStream <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	else std::cout<<"| Successfully created : \'"<< MatrixFile <<"\'."<<std::endl;

/* Open the connect maps CSV file */
	std::ifstream CSVfileStream (ConnecMapsFile.c_str() , std::ios::in);  // opening in reading
	if(! CSVfileStream) // error while opening
	{
		std::cout<<"| Error opening the Connectivity maps file \'"<< ConnecMapsFile <<"\'."<<std::endl;
		std::cout<<"=> ABORT <="<<std::endl;
		return -1;
	}
	else std::cout<<"| Successfully loaded : \'"<< ConnecMapsFile <<"\'."<<std::endl;

	std::string ConnecMap; // will contain the paths to the connec maps (updated in the main FOR loop)

/* Open the Label Map with ITK */
	ImageType::Pointer LabelMapImage;
	try
	{
		ReaderType::Pointer LabelMapReader=ReaderType::New();
		LabelMapReader->SetFileName( LabelMap );
		LabelMapReader->Update();
		LabelMapImage = LabelMapReader->GetOutput();
	}
	catch( itk::ExceptionObject exception )
	{
		std::cerr << exception << std::endl ;
		return -1 ;
	}
	std::cout<<"| Successfully loaded image : \'"<< LabelMap <<"\'."<<std::endl;

	IteratorImageType LabelIter ( LabelMapImage , LabelMapImage->GetLargestPossibleRegion() );

/* Vector Declaration */
	std::vector< Label > Labels; // int LabelID, std::vector< double > Values
	Labels.resize( NbConnectMaps ); // prior test has been done to test that there is the same nb of connect maps in the csv and labels in the label map

/* Connect Maps ITK variables */
	ImageType::Pointer ConnecMapImage;
	ReaderType::Pointer ConnecMapReader = ReaderType::New();

/* Main FOR loop */
	for(int i=0;i<NbConnectMaps;i++)
	{
		// Clear the vector
		Labels.clear();

		// Get the path to the connect map
		std::getline ( CSVfileStream, ConnecMap, ',' );

		// Test if connect map and Label map have the same size
		if( CompareSize(ConnecMap, LabelMap)== 0 ) // returns 1 if same size, 0 if different
		{
			std::cout<<"| The connectivity map \'"<< ConnecMap <<"\' has not the same size as the label map \'"<< LabelMap <<"\'."<<std::endl;
			std::cout<<"=> ABORT <="<<std::endl;
			return -1;
		}

		// Open the connectivity map with ITK
		try
		{
			ConnecMapReader->SetFileName( ConnecMap );
			ConnecMapReader->Update();
			ConnecMapImage = ConnecMapReader->GetOutput();
		}
		catch( itk::ExceptionObject exception )
		{
			std::cerr << exception << std::endl ;
			return -1 ;
		}
		std::cout<<"| Successfully loaded image : \'"<< ConnecMap <<"\'."<<std::endl;
		IteratorImageType ConnecIter ( ConnecMapImage , ConnecMapImage->GetLargestPossibleRegion() );

		// Build the vector by getting all the points in the image
		ConnecIter.GoToBegin();
		for( LabelIter.GoToBegin() ; !LabelIter.IsAtEnd() ; ++LabelIter )
		{
			// Get the values in the image
			int Labelvalue = LabelIter.Get();
			double connecValue = ConnecIter.Get();

			if(Labelvalue!=0)
			{
				// Find the Label and test if already in vector
				int VectorIndex=Labels.size();
				for(unsigned int j=0;j<Labels.size();j++) if(Labels[j].LabelID == Labelvalue) VectorIndex=j;

				if( (unsigned int)VectorIndex == Labels.size() ) // if label not yet in vector
				{
					VectorIndex = Labels.size(); // here the label size corresponds to the future index for the new label because we will do a push_back
					Label NewLabel;
					NewLabel.LabelID = Labelvalue;
					NewLabel.Values.clear();
					Labels.push_back(NewLabel);
				}

				// Update vector
				Labels[VectorIndex].Values.push_back( connecValue );
			}

			++ConnecIter;
		}

		// Sort the labels vector by label :
		std::sort ( Labels.begin(), Labels.end(), CompareLabels ); // CompareLabels() compare the LabelID in each struct Label and returns a bool

		// DEBUG: Display the Labels vector
//		for(unsigned int j=0;j<Labels.size();j++) std::cout<< Labels[j].LabelID << "\t" << Labels[j].Values.size() << std::endl; // display the label ID and the nb of voxels in this label
//		std::cout<<"=========================================================="<<std::endl;

		// Build the matrix row
		std::vector< double > Row;

		Row = GetMeans(Labels); // Compute the mean connectivity in each label

		/* Write the matrix file */
		for(unsigned int j=0;j<Row.size();j++) MatrixfileStream << Row[j]  << " " ;
		MatrixfileStream << std::endl;
	}

	std::cout<<"| Successfully computed Matrix : \'"<< MatrixFile <<"\'."<<std::endl;

/* Close the matrix file */
	MatrixfileStream.close();

/* End of Main function */
	return 0;
}

