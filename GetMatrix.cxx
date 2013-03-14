/*
GetMatrix --connecMapsFile ../ConnecMaps.csv --labelMap /home/akaiser/Networking/Freesurfer/TestPenelope/subjects/Penelope/mri/aparc+aseg_resampled_noWMnoCSF.nrrd --matrixFile ../FSL_Matrix.txt --matrixMetric Mean --connecMapsFileIndex 1 --useRegionBoundary --normalizeFSLConnec

GetMatrix --connecMapsFile ../CostMaps.csv --labelMap /home/akaiser/Networking/Freesurfer/TestPenelope/subjects/Penelope/mri/aparc+aseg_resampled_noWMnoCSF.nrrd --matrixFile ../Cost_Matrix.txt --matrixMetric Mean --connecMapsFileIndex 1 --isCostMap --useRegionBoundary

GetMatrix --tractsFiles /home/akaiser/Networking/ukf/tracts.vtk,/home/akaiser/Networking/ukf/tractsSecondTensor.vtk --labelMap /home/akaiser/Networking/Freesurfer/TestPenelope/subjects/Penelope/mri/aparc+aseg_resampled_noWMnoCSF.nrrd --matrixFile ../Tracts_Matrix.txt --FA ~/Networking/from_Utah/Data/b1000_fa.nrrd
*/

/* std classes */
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <fstream> // to open a file
#include <sstream>

/* GenerateCLP */
#include "GetMatrixCLP.h" //generated when ccmake

/*itk classes*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkVTKImageIO.h"
#include <itksys/SystemTools.hxx> // for GetParentDirectory() and GetFilenameLastExtension()

/*vtk classes*/
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkPoints.h>

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
typedef itk::Image < double , 3 >                     ImageType; //itk type for image
typedef itk::ImageFileReader <ImageType>              ReaderType; //itk reader class to open an image
typedef itk::ImageRegionIterator< ImageType >         IteratorImageType;
typedef itk::ConstNeighborhoodIterator< ImageType >   NeighborhoodIteratorType;
typedef itk::VTKImageIO                               ImageIOType;

  /////////////////////////////////////////
 //           OPEN ITK IMAGE            //
/////////////////////////////////////////
ImageType::Pointer OpenITKimage (std::string ImagePath)
{
  ImageType::Pointer itkImage;
  try
  {
    ReaderType::Pointer reader=ReaderType::New();
    reader->SetFileName( ImagePath );
    reader->Update();
    itkImage = reader->GetOutput();
  }
  catch( itk::ExceptionObject exception )
  {
    std::cerr << exception << std::endl ;
    return NULL;
  }

  return itkImage;
}

  /////////////////////////////////////////
 //            TEST FILES               //
/////////////////////////////////////////
int testFiles (  std::string ConnecMapsFile,
    int ConnecMapsFileIndex,
    std::string LabelMap,
    std::string MatrixFile,
    std::vector< std::string > & TractsFiles, // reference to the tracts files vector so we can pop_back the files that are not ok (not readable or not vtk)
    std::string & FA) // reference to the FA image so we can empty it if error -> empty = not used
// csv file ConnecMapsFile, map image LabelMap, output file MatrixFile, vtk tracts files vector TractsFiles, image FA : returns -1 if problem, 0 if no connectivity maps (use tracts), otherwise returns the nb of connec maps found
{
/* Test if all files are here */
  if( ConnecMapsFile.empty() && TractsFiles.empty() ) // TractsFiles is a vector of string
  {
    std::cout<<"| Please give either a csv file containing paths to connectivity maps or tract file(s) containing the full brain tractography"<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }
  if( LabelMap.empty() )
  {
    std::cout<<"| The Label Map image is missing."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }

/* Test if files are readable */
  if( !ConnecMapsFile.empty() && access(ConnecMapsFile.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
  {
    std::cout<<"| The CSV file containing the connectivity maps for each label is not readable : \'"<< ConnecMapsFile <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }
  if( !ConnecMapsFile.empty() && itksys::SystemTools::GetFilenameLastExtension( ConnecMapsFile.c_str() ) != ".csv" ) // Test if the file is a csv
  {
    std::cout<<"| The given file containing the connectivity maps for each label is not a CSV file : \'"<< ConnecMapsFile <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }
  if( ConnecMapsFile.empty() ) // if tracts files used, remove the ones that are bad and test the FA
  {
    std::vector< int > filesToRemoveFromVector; // indexes of files that are bad
    for(unsigned int i=0;i<TractsFiles.size();i++) // for all the tracts files
    {
      if( access(TractsFiles[i].c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
      {
        std::cout<<"| This tracts file containing the full brain tractography is not readable : \'"<< TractsFiles[i] <<"\', it will not be processed."<<std::endl;
        filesToRemoveFromVector.insert (filesToRemoveFromVector.begin(),i); // insert at the begining so the last indexes will be removed first
      }
      else if( itksys::SystemTools::GetFilenameLastExtension( TractsFiles[i].c_str() ) != ".vtk" ) // Test if the file is a vtk image file
      {
        std::cout<<"| This tracts file containing the full brain tractography is not a VTK image file : \'"<< TractsFiles[i] <<"\', it will not be processed."<<std::endl;
        filesToRemoveFromVector.insert (filesToRemoveFromVector.begin(),i); // insert at the begining so the last indexes will be removed first
      }
    }
    for(unsigned int i=0 ; i<filesToRemoveFromVector.size(); i++) TractsFiles.erase( TractsFiles.begin() + filesToRemoveFromVector[i] ); // (the values are sorted already because of the order they are inserted) we remove the last files first, so the indexes of the other files don't change and we can still use the right indexes in the vector

    if( TractsFiles.empty() )
    {
      std::cout<<"| None of the given tracts files can be processed."<<std::endl;
      std::cout<<"=> ABORT <="<<std::endl;
      return -1;
    }

    if( !FA.empty() && access(FA.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
    {
      std::cout<<"| The FA image is not readable : \'"<< FA <<"\', it will not be used."<<std::endl;
      FA=""; // empty the string so it is not used
    }
  }
  if( access(LabelMap.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
  {
    std::cout<<"| The Label Map image is not readable : \'"<< LabelMap <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }
  std::string MatrixFileFolder = itksys::SystemTools::GetParentDirectory( MatrixFile.c_str() );
  if ( MatrixFileFolder.empty() ) MatrixFileFolder=".";
  if( access(MatrixFileFolder.c_str(), W_OK) != 0 ) // Test if the folder is writable => unistd::access() returns 0 if W(write)_OK
  {
    std::cout<<"| The containing folder given for the output Matrix file is not writable : \'"<< MatrixFileFolder <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }

/* Test if the connectivity maps in the csv file are readable and count them */
  if( !ConnecMapsFile.empty() )
  {
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
    while( std::getline ( CSVfileStream, ConnecMap ) ) 
    {
      std::stringstream strstr(ConnecMap);
      for(int j=0;j<ConnecMapsFileIndex;j++) std::getline(strstr, ConnecMap, ',');

      if( access(ConnecMap.c_str(), R_OK) != 0 ) // Test if the file is readable => unistd::access() returns 0 if R(read)_OK
      {
        std::cout<<"| The following Connectivity map image is not readable : \'"<< ConnecMap <<"\'."<<std::endl;
        std::cout<<"=> ABORT <="<<std::endl;
        return -1;
      }
      NbConnectMaps ++ ;
    }
    CSVfileStream.close();

    std::cout<<"| "<< NbConnectMaps <<" Connectivity maps found."<<std::endl;
    return NbConnectMaps;
  }

  return 0; // use tracts
}

  /////////////////////////////////////////
 //             GET LABELS              //
/////////////////////////////////////////

std::vector< int > GetLabels(std::string LabelMap) // returns empty vector if failed to load image
{
/* Labels Vector */
  std::vector< int > LabelsVector; // contains all the labels already found : NO duplicate

/* Open the Label Map */
  ImageType::Pointer LabelMapImage = OpenITKimage(LabelMap);
  if ( LabelMapImage ) std::cout<<"| Successfully loaded image to find labels: \'"<< LabelMap <<"\'."<<std::endl; //if !=NULL
  else return LabelsVector; // empty at this point

/* Main Loop : test all the voxels and count labels */
  IteratorImageType LabelIter ( LabelMapImage , LabelMapImage->GetLargestPossibleRegion() );

  for( LabelIter.GoToBegin() ; !LabelIter.IsAtEnd() ; ++LabelIter )
  {
    int value = LabelIter.Get();
    if( value!=0 && std::find(LabelsVector.begin(), LabelsVector.end(), value) == LabelsVector.end() ) LabelsVector.push_back(value);
  }

/* Return */
  std::sort ( LabelsVector.begin(), LabelsVector.end() );
  return LabelsVector;
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
 //            TEST BOUNDARY            //
/////////////////////////////////////////

bool voxelBelongToBoundary ( NeighborhoodIteratorType LabelIter) // looks in the block neighborhood and returns true if a voxel from another label is found (not outside the brain), otherwise false
{
  int Center = LabelIter.GetCenterPixel();

  for(unsigned int i=0;i<LabelIter.Size();i++) // for all the neighbors
  {
    int Neighbor = LabelIter.GetPixel(i);

    if( Neighbor!=0 && Neighbor != Center ) return true; // if Neighbor==0 =  not considered as boundary
  }

  return false;
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

std::vector< double > GetMinima (std::vector< Label > Labels)
{
  std::vector< double > Minima;

  for(unsigned int i=0;i<Labels.size();i++) // for all labels in Labels
  {
    double Min=1000000;

    for(unsigned int j=0;j<Labels[i].Values.size();j++) if( Labels[i].Values[j] < Min ) Min = Labels[i].Values[j]; // search the min value for all the values in the label

    Minima.push_back(Min);
  }

  return Minima;
}

std::vector< double > GetMaxima (std::vector< Label > Labels)
{
  std::vector< double > Maxima;

  for(unsigned int i=0;i<Labels.size();i++) // for all labels in Labels
  {
    double Max=0;

    for(unsigned int j=0;j<Labels[i].Values.size();j++) if( Labels[i].Values[j] > Max ) Max = Labels[i].Values[j]; // search the max value for all the values in the label

    Maxima.push_back(Max);
  }

  return Maxima;
}


std::vector< double > GetQuantiles (std::vector< Label > Labels, double Quantile) // double Quantile = 0.1|0.5|0.9 // the Quantile 0.5 (50%) is the median
{
  std::vector< double > Quantiles;

  for(unsigned int i=0;i<Labels.size();i++) // for all labels in Labels
  {
    // Sort the values in the label
    std::sort ( Labels[i].Values.begin(), Labels[i].Values.end() );

    Quantiles.push_back( Labels[i].Values[ Labels[i].Values.size()*Quantile ] ); // The Quantile is a value from the set of values, the value at [Quantile]%
  }

  return Quantiles;
}

  /////////////////////////////////////////
 //       MATRIX WITH CONNEC MAPS       //
/////////////////////////////////////////

std::vector< std::vector< double > > ComputeMatrixWithConnec (  std::string LabelMap, 
                                                                int NbConnectMaps, 
                                                                std::string ConnecMapsFile, 
                                                                int ConnecMapsFileIndex, 
                                                                NeighborhoodIteratorType LabelIter,
                                                                bool useRegionBoundary,
                                                                bool useOnlyReachedVoxels,
                                                                std::string MatrixMetric,
                                                                double Quantile,
                                                                bool isCostMap,
                                                                bool normalizeFSLConnec) // returns empty vector if fail, matrix if ok
{ 
/* Exit Fail */
  std::vector< std::vector< double > > EmptyVector;

/* Test if right number of labels */
  std::vector< int > LabelsInLabelMap =  GetLabels(LabelMap);
  int NbLabels = LabelsInLabelMap.size();
  if( NbLabels != NbConnectMaps )
  {
    std::cout<<"| The number of labels in Label Map \'"<< LabelMap <<"\' ("<< NbLabels <<") is not the same as the number of connectivity maps in the csv file \'"<< ConnecMapsFile <<"\' ("<< NbConnectMaps <<")."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return EmptyVector;
  }

/* Open the connect maps CSV file */
  std::ifstream CSVfileStream (ConnecMapsFile.c_str() , std::ios::in);  // opening in reading
  if(! CSVfileStream) // error while opening
  {
    std::cout<<"| Error opening the Connectivity maps file \'"<< ConnecMapsFile <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return EmptyVector;
  }
  else std::cout<<"| Successfully loaded : \'"<< ConnecMapsFile <<"\'."<<std::endl;

  std::string ConnecMap; // will contain the paths to the connec maps (updated in the main FOR loop)

/* Vector Declaration */
  std::vector< Label > Labels; // int LabelID, std::vector< double > Values
  Labels.resize( NbConnectMaps ); // prior test has been done to test that there is the same nb of connect maps in the csv and labels in the label map

  std::vector< std::vector< double > > ConnecMatrix;

/* Main FOR loop */
  double MinConnec=10000000;
  double MaxConnec=0;

  for(int i=0;i<NbConnectMaps;i++)
  {
    // Clear the vector
    Labels.clear();

    // Get the path to the connect map
    std::getline ( CSVfileStream, ConnecMap ); // ConnecMap if the line
    std::stringstream strstr(ConnecMap); // ConnecMap if the line
    for(int j=0;j<ConnecMapsFileIndex;j++) std::getline(strstr, ConnecMap, ','); // splitting the line by ','

    // Test if connect map and Label map have the same size
    if( CompareSize(ConnecMap, LabelMap)== 0 ) // returns 1 if same size, 0 if different
    {
      std::cout<<"| The connectivity map \'"<< ConnecMap <<"\' has not the same size as the label map \'"<< LabelMap <<"\'."<<std::endl;
      std::cout<<"=> ABORT <="<<std::endl;
      return EmptyVector;
    }

    // Open the connectivity map with ITK
    ImageType::Pointer ConnecMapImage = OpenITKimage(ConnecMap);
    if ( ConnecMapImage ) std::cout<<"| ["<< i+1 <<"\t/ "<< NbConnectMaps <<"\t] Successfully loaded image : \'"<< ConnecMap <<"\'"; //if !=NULL
    else return EmptyVector;
    
    IteratorImageType ConnecIter ( ConnecMapImage , ConnecMapImage->GetLargestPossibleRegion() );

    // Get the total number of tracts used in FSL (the 'waytotal' file created by PROBTRACKX gives this nb) => normalization
    double FSLNbTracts;
    if(normalizeFSLConnec)
    {
      // Open the 'waytotal' file
      std::string FSLwaytotalFile = itksys::SystemTools::GetFilenamePath( ConnecMap ) + "/waytotal";
      std::ifstream waytotalFileStream (FSLwaytotalFile.c_str() , std::ios::in);  // opening in reading

      // Get the number
      waytotalFileStream>>FSLNbTracts;

      std::cout<<" \t("<<FSLNbTracts<<" \ttracts used to compute this connectivity map)."<<std::endl;
    }
    else std::cout<<"."<<std::endl; //display

    // Build the vector by getting all the points in the image
    for( LabelIter.GoToBegin(), ConnecIter.GoToBegin() ; !LabelIter.IsAtEnd() ; ++LabelIter, ++ConnecIter)
    {
      // Get the values in the image
      int Labelvalue = LabelIter.GetCenterPixel(); // !! LabelIter is a NeighborhoodIteratorType
      double connecValue = (double) ConnecIter.Get();

      if( Labelvalue!=0 && connecValue>=0 ) // take into account only the values that belong to a label and the positive values (negative values are not supposed to be in region of interests)
      {
        // Value processing
        if( isCostMap ) connecValue = 1 / ( 1 + connecValue ) ; // done here to avoid connecValue = -1 and /0
        else if( normalizeFSLConnec && FSLNbTracts!=0 ) connecValue = connecValue / FSLNbTracts ;

        // test Boundary
        bool boundaryOK = true ;
        if( useRegionBoundary ) if( ! voxelBelongToBoundary( LabelIter ) ) boundaryOK = false;

        // test reached voxel
        bool ReachedOK = true ;
        if( useOnlyReachedVoxels ) if( connecValue < 0.0001 ) ReachedOK = false; // voxels which connectivity is under 0.0001 are supposed not reached

        if( boundaryOK && ReachedOK )
        {
//DEBUG : INDEXES OF SELECTED VOXELS//  std::cout<<"Voxel : "<< LabelIter.GetIndex()[0] <<":" << LabelIter.GetIndex()[1] <<":" << LabelIter.GetIndex()[2] <<std::endl;

          // Find the Label and test if already in vector
          int VectorIndex=Labels.size();
          for(unsigned int j=0;j<Labels.size();j++) if(Labels[j].LabelID == Labelvalue) VectorIndex=j;

          // Create Label element if label not yet in vector
          if( (unsigned int)VectorIndex == Labels.size() )
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
      }
    } // connec image done

    // Check if all labels have been found and fill vector with missing ones
    for(unsigned int j=0;j<LabelsInLabelMap.size();j++) // theoratical labels in given label map
    {
      bool LabelMissing=true;
      for(unsigned int k=0;k<Labels.size();k++) // search all the labels vector to find the jth label in map
      {
        if(Labels[k].LabelID==LabelsInLabelMap[j]) LabelMissing = false; // label found in vector
      }
      if(LabelMissing) // if value not in the vector (so not reached in the connec map) => create one with the value -1 as connec => turned to minimum connec value when normalizing
      {
//std::cout<<"Connec map #"<<i<<" | Label missing : "<<LabelsInLabelMap[j]<<std::endl;
        Label LabelNotFound;
        LabelNotFound.LabelID = LabelsInLabelMap[j];
        LabelNotFound.Values.push_back( -1 );
        Labels.push_back(LabelNotFound);
      }
    }

    // Sort the labels vector by labelID :
    std::sort ( Labels.begin(), Labels.end(), CompareLabels ); // CompareLabels() compare the LabelID in each struct Label and returns a bool

///////////////////////////////////
/*
if(i==15) // to write out the region sizes in a file
{
  std::ofstream fstream ("regionSizes" , std::ios::out | std::ios::trunc);  // opening in writing with erasing the open file
  if( fstream) // error while opening
  {
    for(unsigned int j=0;j<Labels.size();j++) fstream <<Labels[j].Values.size() << std::endl;
    fstream.close();
  }
}
*/
//////////////////////////////

    // DEBUG: Display the Labels vector
//    for(unsigned int j=0;j<Labels.size();j++) std::cout<< Labels[j].LabelID << "\t" << Labels[j].Values.size() << std::endl; // display the label ID and the nb of voxels in this label
//    std::cout<<"=========================================================="<<std::endl;

    // Build the matrix row
    std::vector< double > Row;

    // Compute the metric
    if(MatrixMetric == "Mean") Row = GetMeans(Labels); // Returns the mean connectivity in each label
    else if(MatrixMetric == "Minimum") Row = GetMinima(Labels); // Returns the min connectivity in each label
    else if(MatrixMetric == "Maximum") Row = GetMaxima(Labels); // Returns the max connectivity in each label
    else if(MatrixMetric == "Quantile") Row = GetQuantiles(Labels,Quantile); // Returns the given quantile connectivity in each label

    // Get the minimum and maximum values in the matrix
    for(unsigned int j=0;j<Row.size();j++)
    {
//      if( isCostMap ) Row[j] = Row[j] * Labels[i].Values.size(); // if cost map, normalize with the nb of voxels in source label

      if( i!=(int)j && Row[j]!=-1 ) // values on the diag are not real connec values (usually 0)
      {
         if( Row[j] < MinConnec ) MinConnec = Row[j];
         if( Row[j] > MaxConnec ) MaxConnec = Row[j];
      }
    }

    ConnecMatrix.push_back( Row );
  }

  // Display min max
  std::cout<<"| The range of the connectivy values in the matrix is: "<<MinConnec<<" -> "<<MaxConnec<<std::endl;

  // Normalize values
  for(unsigned int i=0;i<ConnecMatrix.size();i++)
  {
    for(unsigned int j=0;j<ConnecMatrix.size();j++)
    {
      if( i == j ) ConnecMatrix[i][j] = MinConnec ; // Minimum value on the diagonal -> allow the diag value to be in the range of all values and not completely different
      else if( ConnecMatrix[i][j]==-1 ) ConnecMatrix[i][j] = MinConnec ; // Minimum value for the non reached regions
//      ConnecMatrix[i][j] = ( ConnecMatrix[i][j] - MinConnec ) / ( MaxConnec - MinConnec ) ; // normalize the values between 0 and 1
    }
  }

  return ConnecMatrix; // EXIT_OK

}

  /////////////////////////////////////////
 //         MATRIX WITH TRACTS          //
/////////////////////////////////////////

ImageType::Pointer ImageLPStoRAS (ImageType::Pointer Image)
{
/* LPS to RAS = flip x and y axis + change origin : x=x+Nx , y=y+Ny
Nx, Ny = nb of voxels on x and y axis
(x)    (-1  0  0 Nx )   (x)
(y) =  ( 0 -1  0 Ny ) x (y)
(z)    ( 0  0  1  0 )   (z)
(1)ras ( 0  0  0  1 )   (1)lps

    { Xras = Xlps + Nx
<=> { Yras = Ylps + Ny
    { Zras = Zlps

=> get the values to add to the origin : Nx and Ny :

N1, N2 = nb of voxels on 2 1st directions | D1,2,3 = directions | spx,y,z = spacings
(Ox)    (Ox)     (spx 0  0 )   (D1x D2x D3x)   (N1)
(Oy) =  (Oy)  +  ( 0 spy 0 ) x (D1y D2y D3y) x (N2)
(Oz)ras (Oz)lps  ( 0  0 spz)   (D1z D2z D3z)   (0 )

    { Oxras = Oxlps + spx * (D1x.N1 + D2x.N2)
<=> { Oyras = Oylps + spy * (D1y.N1 + D2y.N2)
    { Ozras = Ozlps + 0

    { Oxras = Oxlps + Nx, with Nx = spx * (D1x.N1 + D2x.N2)
<=> { Oyras = Oylps + Ny, with Ny = spy * (D1y.N1 + D2y.N2)
    { Ozras = Ozlps + 0
*/
  ImageType::SizeType size = Image->GetLargestPossibleRegion().GetSize(); // get N1 and N2
  ImageType::SpacingType spacing = Image->GetSpacing(); // get spacing
  ImageType::DirectionType Direction = Image->GetDirection(); // Direction[row][col]
  ImageType::PointType Origin = Image->GetOrigin();

  double Nx = spacing[0] * (Direction[0][0] * size[0] + Direction[0][1] * size[1]);
  double Ny = spacing[1] * (Direction[1][0] * size[0] + Direction[1][1] * size[1]);

  // change origin : x=x+Nx , y=y+Ny, with Nx and Ny in real values, not index
  Origin[0] = Origin[0] + Nx;
  Origin[1] = Origin[1] + Ny;
  Image->SetOrigin( Origin );

  // flip 2first axis : change 2 first columns of direction matrix
  Direction[0][0] = -1 * Direction[0][0];
  Direction[1][0] = -1 * Direction[1][0];
  Direction[2][0] = -1 * Direction[2][0];
  Direction[0][1] = -1 * Direction[0][1];
  Direction[1][1] = -1 * Direction[1][1];
  Direction[2][1] = -1 * Direction[2][1];
  Image->SetDirection( Direction );

  return Image;
}

std::vector< std::vector< double > > ComputeMatrixWithTracts(  std::string LabelMap,
                                                               ImageType::Pointer LabelMapImage,
                                                               NeighborhoodIteratorType LabelIter,
                                                               std::vector< std::string > TractsFiles,
                                                               std::string FA) // returns empty vector if fail, matrix if ok
{
/* Exit Fail Vector */
  std::vector< std::vector< double > > EmptyVector;

/* Vector Declaration */
  std::vector< std::vector< double > > ConnecMatrix;
  std::vector< int > Labels = GetLabels(LabelMap); // Contains all the labels, sorted
  int NbLabels = Labels.size();
  if( NbLabels==0 ) return EmptyVector;

  ConnecMatrix.resize(NbLabels);
  for(int i=0;i<NbLabels;i++) ConnecMatrix[i].resize(NbLabels);

/* Open FA image with ITK */
  ImageType::Pointer FAimage; // is here to avoid compilation error (variable missing)
  if( !FA.empty() )
  {
    FAimage = OpenITKimage(FA);
    if ( FAimage ) std::cout<<"| Successfully loaded FA image : \'"<< FA <<"\'."<<std::endl; //if !=NULL
    else return EmptyVector; // empty at this point

    FAimage = ImageLPStoRAS( FAimage );
  }

/* For all the tract files */
  int NbFibersInMatrix=0; // contains the nb of fibers in matrix, so the fibers linking 2 "non 0" regions

  if( TractsFiles.size() >= 2 ) std::cout<<"| "<< TractsFiles.size() << " tracts files will be processed"<<std::endl;
  else std::cout<<"| 1 tracts file will be processed"<<std::endl; // for the 's' in the end of 'file'

  for(unsigned int tf=0;tf<TractsFiles.size();tf++)
  {
    std::cout<<"| ["<< tf+1 <<"/"<< TractsFiles.size() << "] tracts file"<<std::endl;

/* Open Tracts File */ // open file: fvlight::FiberViewerLightGUI.cxx line 384
    vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
    reader->SetFileName( TractsFiles[tf].c_str() );
    std::cout<<"| Opening vtk tracts file \'"<< TractsFiles[tf] <<"\'..."<<std::endl;
    reader->Update(); // takes a "long" time
    vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

    if( PolyData==NULL )
    {
      std::cout<<"| Error loading vtk file : \'"<< TractsFiles[tf] <<"\' -> Ignoring this tract file."<<std::endl;
      if( TractsFiles.size()<2 ) // if the only tract file fails, then abort
      {
        std::cout<<"=> ABORT <="<<std::endl;
        return EmptyVector;
      }
    }
    else
    {
      std::cout<<"| Successfully loaded fibers vtk image : \'"<< TractsFiles[tf] <<"\'."<<std::endl;

      int NbTracts = PolyData->GetNumberOfCells();
      std::cout<<"| "<< NbTracts <<" fibers found in the file \'"<< TractsFiles[tf] <<"\'."<<std::endl;

/* Main FOR loop */
      for(int i=0;i<NbTracts;i++) // for all the tracts
      {
        // Display nb of tracts done
        std::cout<<"\r| ["<< i+1 <<"/"<< NbTracts << "] fibers read";

        // Get the RAS (Right, Anterior, Superior) coordinates of the begin and end point of the tract
        vtkSmartPointer<vtkCell> Cell = PolyData->GetCell( i );
        double * FiberBounds = Cell->GetBounds(); // (xmin,xmax,ymin,ymax,zmin,zmax)
// std::cout<<"("<< FiberBounds[0] <<" , "<< FiberBounds[2] <<" , "<< FiberBounds[4] <<") -> ("<< FiberBounds[1] <<" , "<< FiberBounds[3] <<" , "<< FiberBounds[5] <<")"<<std::endl; // info display : begin and end points 

        // Convert the RAS coordinates into index coordinates (ITK: physical to index) : bool TransformPhysicalPointToIndex (const Point< TCoordRep, VImageDimension > &point, IndexType &index) const 
        ImageType::IndexType BeginVoxelIndex;
        ImageType::IndexType EndVoxelIndex;
        ImageType::PointType PhysicalBeginPoint;
        ImageType::PointType PhysicalEndPoint;

        PhysicalBeginPoint[0] = FiberBounds[0];
        PhysicalBeginPoint[1] = FiberBounds[2];
        PhysicalBeginPoint[2] = FiberBounds[4];
        PhysicalEndPoint[0] = FiberBounds[1];
        PhysicalEndPoint[1] = FiberBounds[3];
        PhysicalEndPoint[2] = FiberBounds[5];

        bool BeginIsInside = LabelMapImage -> TransformPhysicalPointToIndex( PhysicalBeginPoint, BeginVoxelIndex ); // true if index inside the image, false if outside
        bool EndIsInside = LabelMapImage -> TransformPhysicalPointToIndex( PhysicalEndPoint, EndVoxelIndex ); // true if index inside the image, false if outside

        if( !BeginIsInside || !EndIsInside )
        {
          if(!BeginIsInside) std::cout<<" : This fiber begins outside the label map: physical coordinates = ("<< PhysicalBeginPoint[0] <<" , "<< PhysicalBeginPoint[1] <<" , "<< PhysicalBeginPoint[2] <<").";
          if(!EndIsInside) std::cout<<" : This fiber ends outside the label map: physical coordinates = ("<< PhysicalEndPoint[0] <<" , "<< PhysicalEndPoint[1] <<" , "<< PhysicalEndPoint[2] <<").";
          std::cout<<" It will not be taken into account."<<std::endl;
        }
        else
        {
          // Find the corresponding labels
          LabelIter.SetLocation( BeginVoxelIndex );
          int BeginLabel=LabelIter.GetCenterPixel();

          LabelIter.SetLocation( EndVoxelIndex );
          int EndLabel=LabelIter.GetCenterPixel();

//std::cout<<BeginLabel<<" -> "<< EndLabel<<std::endl;

          // Add a number in the matrix case
          if( BeginLabel!=0 && EndLabel!=0)
          {
            int BeginMatrixIndex = std::distance( Labels.begin(), std::find(Labels.begin(), Labels.end(), BeginLabel) ); // distance to convert std::iterator to int
            int EndMatrixIndex = std::distance( Labels.begin(), std::find(Labels.begin(), Labels.end(), EndLabel) );

            if( FA.empty() ) ConnecMatrix[ BeginMatrixIndex ][ EndMatrixIndex ] = ConnecMatrix[ BeginMatrixIndex ][ EndMatrixIndex ] + 1 ; // add 1 for each tract that links the 2 regions 
            else // compute the mean FA along the tract
            {
              double FAsum = 0;
              vtkSmartPointer<vtkPoints> CellPoints = Cell->GetPoints();
              int NbPoints = CellPoints->GetNumberOfPoints(); // nb of points in the fiber

              for(int j=0;j<NbPoints;j++)
              {
                ImageType::PointType PhysicalPoint =  CellPoints->GetPoint( j ); // GetPoint() returns a double[3] with the PHYSICAL points
                ImageType::IndexType PointIndex;

                bool PointIsInside = FAimage -> TransformPhysicalPointToIndex( PhysicalPoint, PointIndex ); // true if index inside the image, false if outside
                if(!PointIsInside)  std::cout<<"| This point is outside the FA image: physical coordinates = ("<< PhysicalPoint[0] <<" , "<< PhysicalPoint[1] <<" , "<< PhysicalPoint[2] <<")"<<std::endl;
                else
                {
                  IteratorImageType FAimageIter ( FAimage , FAimage->GetLargestPossibleRegion() );
                  FAimageIter.SetIndex( PointIndex );
                  double FAvalue = FAimageIter.Get();
                  FAsum = FAsum + FAvalue ;
                }
              }

              ConnecMatrix[ BeginMatrixIndex ][ EndMatrixIndex ] = ConnecMatrix[ BeginMatrixIndex ][ EndMatrixIndex ] + FAsum/NbPoints ; // add the mean FA for each tract that links the 2 regions
            }

            NbFibersInMatrix ++ ; // contains the nb of fibers in matrix, so the fibers linking 2 "non 0" regions
          }
        }// else of if( !BeginIsInside || !EndIsInside )

      } // for(int i=0;i<NbTracts;i++) // for all the tracts

      std::cout<< std::endl; // for the number of fibers read

    } // else of if( PolyData==NULL )

  } // for(unsigned int tf=0;tf<TractsFiles.size();tf++)

  std::cout<< "| "<< NbFibersInMatrix << " fibers used to build the matrix (the others begin or end outside of the regions)."<<std::endl;

/* Normalize and Get the minimum and maximum values in the matrix */
/*  double MinConnec=10000000;
  double MaxConnec=0;
  for(unsigned int i=0;i<ConnecMatrix.size();i++)
  {
    for(unsigned int j=0;j<ConnecMatrix.size();j++)
    {
      if( i!=j) // values on the diag are not real connec values
      {
        ConnecMatrix[i][j] = ConnecMatrix[i][j] / NbFibersInMatrix ;

        // Threshold
        double threshold=50;
        if( FA.empty() ) threshold=0.01344;
        if( ConnecMatrix[i][j] > threshold ) ConnecMatrix[i][j]=threshold;

        if( ConnecMatrix[i][j] < MinConnec ) MinConnec = ConnecMatrix[i][j];
        if( ConnecMatrix[i][j] > MaxConnec ) MaxConnec = ConnecMatrix[i][j];
      }
    }
  }
  std::cout<<"| The range of the connectivy values in the matrix is: "<<MinConnec<<" -> "<<MaxConnec<<std::endl;
*/
/* Normalize */
/*  for(unsigned int i=0;i<ConnecMatrix.size();i++)
  {
    for(unsigned int j=0;j<ConnecMatrix.size();j++)
    {
      if( i==j ) ConnecMatrix[i][j] = MinConnec; // min on the diagonal
      ConnecMatrix[i][j] = ( ConnecMatrix[i][j] - MinConnec ) / ( MaxConnec - MinConnec ) ; // normalize the values between 0 and 1
    }
  }
*/
  return ConnecMatrix; // EXIT_OK
}

  /////////////////////////////////////////
 //           MAIN FUNCTION             //
/////////////////////////////////////////

int main (int argc, char *argv[])
{
  PARSE_ARGS; //thanks to this line, we can use the variables entered in command line directly as variables of the program
/* Generate CLP arguments :
csv file  ConnecMapsFile
int    ConnecMapsFileIndex
map image  LabelMap
image    FA
output file  MatrixFile
std::string  MatrixMetric = "Mean|Minimum|Maximum|Median|Quantile"
bool    useRegionBoundary
bool    useOnlyReachedVoxels
bool    isCostMap
bool    normalizeFSLConnec
vtk file  TractsFiles
double    Quantile
*/

/* Test the files given by user */
  if( MatrixFile.empty() )
  {
    std::cout<<"| The output Matrix file is missing. It will be created in the current work directory as \'./ConnectivityMatrix.txt\'."<<std::endl;
    MatrixFile = "./ConnectivityMatrix.txt";
  }
  int NbConnectMaps = testFiles(ConnecMapsFile, ConnecMapsFileIndex, LabelMap, MatrixFile, TractsFiles, FA); // returns -1 if problem, 0 if no connectivity maps (use tracts), otherwise returns the nb of connec maps found
  if(NbConnectMaps==-1) return -1; // info display in the function

/* Display Infos */
  if( !ConnecMapsFile.empty() )
  {
    std::cout<<"| The matrix will be computed using the metric: "<< MatrixMetric; // info display
    if(MatrixMetric=="Quantile") std::cout<<" = "<< Quantile<< std::endl; // info display
    else std::cout<< std::endl;
    if(useRegionBoundary) std::cout<<"| Only the boundary voxels will be taken into account."<< std::endl ; // info display
    if( isCostMap ) std::cout<<"| As the given connectivity map images are cost images, the final values will be inverted with x -> 1/(1+x)."<< std::endl ; // info display
  }

/* Open the Label Map with ITK */
  ImageType::Pointer LabelMapImage = OpenITKimage(LabelMap);
  if ( LabelMapImage ) std::cout<<"| Successfully loaded image to compute matrix: \'"<< LabelMap <<"\'."<<std::endl; //if !=NULL
  else return -1; // empty at this point

  // If tracts used, convert label map image coordinate system from LPS (nrrd, itk) to RAS (vtk) to compare with the tracts points
  if( ConnecMapsFile.empty() ) LabelMapImage = ImageLPStoRAS( LabelMapImage ); // if using tracts

  // Neighborhood iterator for label map
  NeighborhoodIteratorType::RadiusType Neighborhoodradius;
  Neighborhoodradius.Fill(1); // block neighborhood 
  NeighborhoodIteratorType LabelIter ( Neighborhoodradius, LabelMapImage , LabelMapImage->GetLargestPossibleRegion() );

/* Get the matrix */ // functions return empty vector if fail, connec matrix if OK
  std::vector< std::vector< double > > ConnecMatrix;

  if( !ConnecMapsFile.empty() ) ConnecMatrix = ComputeMatrixWithConnec ( LabelMap, NbConnectMaps, ConnecMapsFile, ConnecMapsFileIndex, LabelIter, useRegionBoundary, useOnlyReachedVoxels, MatrixMetric, Quantile, isCostMap, normalizeFSLConnec );
  else ConnecMatrix = ComputeMatrixWithTracts ( LabelMap, LabelMapImage, LabelIter, TractsFiles, FA );

  if( ConnecMatrix.empty() ) return -1;

/* Write the matrix file */
  std::ofstream MatrixfileStream (MatrixFile.c_str() , std::ios::out | std::ios::trunc);  // opening in writing with erasing the open file
  if(! MatrixfileStream) // error while opening
  {
    std::cout<<"| Error creating the matrix file \'"<< MatrixFile <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }
  else std::cout<<"| Successfully created : \'"<< MatrixFile <<"\'."<<std::endl;

  for(unsigned int i=0;i<ConnecMatrix.size();i++)
  {
    for(unsigned int j=0;j<ConnecMatrix.size();j++) MatrixfileStream << ConnecMatrix[i][j]  << "," ;

    MatrixfileStream << std::endl;
  }

  std::cout<<"| Successfully computed Matrix : \'"<< MatrixFile <<"\'."<<std::endl;

  MatrixfileStream.close();

/* Create Labels File */
  std::string LabelsFile = itksys::SystemTools::GetFilenamePath( MatrixFile ) + "/Matrix_Labels.txt";
  std::ofstream LabelsfileStream (LabelsFile.c_str() , std::ios::out | std::ios::trunc);  // opening in writing with erasing the open file
  if(! LabelsfileStream) // error while opening
  {
    std::cout<<"| Error creating the labels file \'"<< LabelsFile <<"\'."<<std::endl;
    std::cout<<"=> ABORT <="<<std::endl;
    return -1;
  }
  else std::cout<<"| Successfully created : \'"<< LabelsFile <<"\'."<<std::endl;

  // Put the labels into the file
  std::vector< int > Labels = GetLabels(LabelMap);
  for(unsigned int i=0;i<Labels.size();i++) LabelsfileStream << Labels[i] << std::endl;

  // Close the file
  LabelsfileStream.close();

/* End of Main function */
  return 0;
}

