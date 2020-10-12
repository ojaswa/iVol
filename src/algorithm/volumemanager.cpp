#include "volumemanager.h"

#include <vtkImageData.h>
#include <vtkImageGradient.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkInformation.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageAccumulate.h>
#include <vtkImageInterpolator.h>
#include <vtkImageShiftScale.h>

//Clustering

// Template declarations
template void VolumeManager::computeLH<unsigned char>();
template void VolumeManager::computeLH<unsigned short>();


VolumeManager::VolumeManager()
{
    imageData = vtkSmartPointer<vtkImageData>::New();
}

void VolumeManager::setData(vtkImageData* idata)
{
    imageData->ShallowCopy(idata);
}

vtkImageData* VolumeManager::getData()
{
    return imageData;
}

vtkImageData* VolumeManager::getLH()
{
    return imageLH;
}

template <class T>
void VolumeManager::computeLH()
{
    //Compute the gradient of gaussian smoothed volume
    vtkImageGaussianSmooth *gaussianFilter  = vtkImageGaussianSmooth::New();
    gaussianFilter->SetInputData(imageData);
    gaussianFilter->SetDimensionality(3);
    gaussianFilter->SetStandardDeviations(1.0, 1.0, 1.0);
    gaussianFilter->Update();

    vtkImageGradient *gradient = vtkImageGradient::New();
    gradient->SetInputData(gaussianFilter->GetOutput());
    gradient->SetDimensionality(3);
    gradient->HandleBoundariesOn();
    gradient->Update();

    //Label voxels for LH computation
    vtkImageData *labelImage = computeGradientLabelImageLH<T>(gaussianFilter->GetOutput());

    //Create LH array
    int dim[3], type;
    double spa[3], ori[3];
    imageData->GetDimensions(dim);
    imageData->GetSpacing(spa);
    imageData->GetOrigin(ori);
    type = imageData->GetScalarType();

    imageLH = vtkSmartPointer<vtkImageData>::New();
    imageLH->SetDimensions(dim);
    imageLH->SetSpacing(spa);
    imageLH->SetOrigin(ori);
    imageLH->AllocateScalars(type, 2);

    //Compute LH histogram
    T *imageDataPtr = static_cast<T*> (gaussianFilter->GetOutput()->GetScalarPointer());
    T *imageLHPtr = static_cast<T*> (imageLH->GetScalarPointer());
    unsigned char *labelImagePtr = static_cast<unsigned char*> (labelImage->GetScalarPointer());

    vtkImageInterpolator *gradientInterpolator = vtkImageInterpolator::New(); //For linear interpolation of normals
    gradientInterpolator->Initialize(gradient->GetOutput());
    gradientInterpolator->SetOutValue(0.0);
    gradientInterpolator->SetInterpolationModeToLinear();
    gradientInterpolator->Update();

    vtkImageInterpolator *imageInterpolator = vtkImageInterpolator::New();
    imageInterpolator->Initialize(gaussianFilter->GetOutput());
    imageInterpolator->SetOutValue(0.0);
    imageInterpolator->SetInterpolationModeToLinear();
    imageInterpolator->Update();

    T value, lh[2];
    double delta_t = 1.0;
    unsigned char label;
    for (int z=0; z<dim[2]; z++)
        for(int y=0; y<dim[1]; y++)
            for(int x=0; x<dim[0]; x++)
            {
                //Seed position is (x, y, z) with following gradient
                value = *(imageDataPtr + x + (y + z*dim[1])*dim[0]);
                label = *(labelImagePtr + x + (y + z*dim[1])*dim[0]);
                if (label == 0) { //LH = Voxel value, for inner voxels
                    *(imageLHPtr + (x + (y + z*dim[1])*dim[0])*2) = value;
                    *(imageLHPtr + (x + (y + z*dim[1])*dim[0])*2 + 1) = value;
                } else { //Trace along the normal using Runge Kutta
                    //Trace along +ve direction
                    lh[0] = (T)gradientTrack(imageInterpolator, gradientInterpolator, delta_t, x, y, z);

                    //Trace along -ve direction
                    lh[1] = (T)gradientTrack(imageInterpolator, gradientInterpolator, -delta_t, x, y, z);

                    //Assign values
                    *(imageLHPtr + (x + (y + z*dim[1])*dim[0])*2) = (lh[0] < lh[1])?lh[0]:lh[1]; //L
                    *(imageLHPtr + (x + (y + z*dim[1])*dim[0])*2 + 1) = (lh[0] < lh[1])?lh[1]:lh[0];//H
                }
            }

    /* //Testing code
    std::ofstream stream("LH.bin", std::ios::binary );
    long elems = dim[0]*dim[1]*dim[2]*2;
    stream.write(reinterpret_cast<const char*>(imageLH->GetScalarPointer()), elems*sizeof(T));
    stream.close();
    fprintf(stderr, "Wrote binary LH file: %s\n", "LH.bin");
*/

    //Create a density image of bivariate histogram (for display purpose)
    imageLH->GetDimensions(dim);
    double valuesRange[2];
    imageData->GetScalarRange(valuesRange);
    double scale = 255.0/(valuesRange[1] - valuesRange[0]);
    vtkSmartPointer<vtkImageShiftScale> scaler = vtkSmartPointer<vtkImageShiftScale>::New();
    scaler->SetOutputScalarTypeToUnsignedChar();
    scaler->SetScale(scale);
    scaler->SetShift(-valuesRange[0]);
    scaler->SetInputData(imageLH);
    scaler->Update();
    unsigned char *scaledLHPtr = static_cast<unsigned char*> (scaler->GetOutput()->GetScalarPointer());

    imageLHDensityPlot = vtkSmartPointer<vtkImageData>::New();
    imageLHDensityPlot->SetDimensions(256, 256, 1); // Is a 2D image of counts
    imageLHDensityPlot->SetSpacing(1, 1, 1);
    imageLHDensityPlot->SetOrigin(0, 0, 0);
    imageLHDensityPlot->AllocateScalars(VTK_DOUBLE, 1);

    double *densityPtr = static_cast<double*>(imageLHDensityPlot->GetScalarPointer());
    long nelems = dim[0]*dim[1]*dim[2];
    int L, H;
    long ndens = 256*256;
    double freq;
    for(long i=0; i<ndens; i++)
        *(densityPtr + i) = 0.0; // Initialize with 1 so that taking log is not problematic.
    for(long i=0; i<nelems; i++) {
        L = *(scaledLHPtr + i*2);
        H = *(scaledLHPtr + i*2 + 1);
        freq = *(densityPtr + L + H*256);
        *(densityPtr + L + H*256) = freq + 1.0;
    }
    for(long i=0; i<ndens; i++) //Take log for better color contrast
        *(densityPtr + i) = logf(1.0 + *(densityPtr + i));

    //Free memory
    labelImage->Delete();
    imageInterpolator->ReleaseData();
    imageInterpolator->Delete();
    gradientInterpolator->ReleaseData();
    gradientInterpolator->Delete();
    gradient->Delete();
    gaussianFilter->Delete();
}

double VolumeManager::gradientTrack(vtkImageInterpolator *image, vtkImageInterpolator *gradient, double delta_t, int x, int y, int z)
{
    double pos[3];
    double value, value_old;
    int sign, sign_old;
    int dim[3];

    imageData->GetDimensions(dim);
    pos[0] = x; pos[1] = y; pos[2] = z;

#if LH_INTERPOLATE_IMAGE_WHILE_RK_TACKING
    image->SetInterpolationModeToLinear();
#else
    image->SetInterpolationModeToNearest();
#endif
    image->Update();
    image->Interpolate(pos, &value);

    value_old = value;    
    rungeKuttaEstimate(gradient, delta_t, pos);//New position/time is returned in pos, t
    image->Interpolate(pos, &value);
    sign_old = value_old - value;
    sign = sign_old;
    //fprintf(stderr, "<P:%c", (sign>0)?'+':(sign<0)?'-':'0');
    while(sign*sign_old > 0) { // If change of sign (product -ve) or no change in value (product 0), stop.
        value_old = value;
        rungeKuttaEstimate(gradient, delta_t, pos);//New position/time is returned in pos, t
        image->Interpolate(pos, &value);
        sign_old = sign;
        sign = value_old - value;
        //fprintf(stderr, "%c", (sign>0)?'+':(sign<0)?'-':'0');
    }
    //fprintf(stderr, ">\n");
    return value_old;
}

void VolumeManager::rungeKuttaEstimate(vtkImageInterpolator *gradient, double delta_t, double *pos)
{
    //RK2 estimate of position in a gradient field
    double p[3];
    double grad1[3], grad2[3];
    double gmag;
    double direction[3];

    //Get gradient at initial position
    p[0] = pos[0]; p[1] = pos[1]; p[2] = pos[2];
    gradient->Interpolate(p, grad1);
    gmag = sqrt(grad1[0]*grad1[0] + grad1[1]*grad1[1] + grad1[2]*grad1[2]);
    //Get direction vector from gradient
    direction[0] = grad1[0]/(gmag+TINY);
    direction[1] = grad1[1]/(gmag+TINY);
    direction[2] = grad1[2]/(gmag+TINY);
    //Alter position in this direction
    p[0] += delta_t*direction[0];
    p[1] += delta_t*direction[1];
    p[2] += delta_t*direction[2];

    //Compute grad2 at p
    gradient->Interpolate(p, grad2);
    //Compute average of the two gradients
    grad2[0] = (grad1[0] + grad2[0])/2.0;
    grad2[1] = (grad1[1] + grad2[1])/2.0;
    grad2[2] = (grad1[2] + grad2[2])/2.0;
    gmag = sqrt(grad2[0]*grad2[0] + grad2[1]*grad2[1] + grad2[2]*grad2[2]);
    //Get direction vector from gradient
    direction[0] = grad2[0]/(gmag+TINY);
    direction[1] = grad2[1]/(gmag+TINY);
    direction[2] = grad2[2]/(gmag+TINY);
    //Alter origin position in this direction (to be returned to caller)
    pos[0] += delta_t*direction[0];
    pos[1] += delta_t*direction[1];
    pos[2] += delta_t*direction[2];

#if !LH_INTERPOLATE_IMAGE_WHILE_RK_TACKING
    //Convert to integer and bounds check.
    pos[0] = int(pos[0]);
    pos[1] = int(pos[1]);
    pos[2] = int(pos[2]);
    int dim[3];
    imageData->GetDimensions(dim);
    if(pos[0] < 0) pos[0] = 0; if(pos[0] >= dim[0]) pos[0] = dim[0]-1;
    if(pos[1] < 0) pos[1] = 0; if(pos[1] >= dim[1]) pos[1] = dim[1]-1;
    if(pos[2] < 0) pos[2] = 0; if(pos[2] >= dim[2]) pos[2] = dim[2]-1;
#endif
}

template <class T>
vtkImageData* VolumeManager::computeGradientLabelImageLH(vtkImageData *smoothedImage)
{
    vtkImageGradientMagnitude *gmag = vtkImageGradientMagnitude::New(); // This is the same data type as the input image
    gmag->SetDimensionality(3);
    gmag->SetInputData(smoothedImage);
    gmag->Update();

    double epsilon = computeEpsilon(gmag);

    int dim[3];
    double spa[3], ori[3];
    smoothedImage->GetDimensions(dim);
    smoothedImage->GetSpacing(spa);
    smoothedImage->GetOrigin(ori);

    vtkImageData *labelImage = vtkImageData::New();
    labelImage->SetDimensions(dim);
    labelImage->SetSpacing(spa);
    labelImage->SetOrigin(ori);
    labelImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

    T *gmagPtr = static_cast<T*> (gmag->GetOutput()->GetScalarPointer());
    T *imagePtr = static_cast<T*> (smoothedImage->GetScalarPointer());
    unsigned char *labelImagePtr = static_cast<unsigned char*> (labelImage->GetScalarPointer());

    T gradmag, value, neib;
    int sum;
    for (int z=0; z<dim[2]; z++)
        for(int y=0; y<dim[1]; y++)
            for(int x=0; x<dim[0]; x++)
            {
#if LH_BASIC_THRESHOLDING
                gradmag = *(gmagPtr + x + (y + z*dim[1])*dim[0]);
                *(labelImagePtr + x + (y + z*dim[1])*dim[0]) = (gradmag > epsilon)?255:0;
#else
                if(x == 0 || x == (dim[0]-1) || y == 0 || y == (dim[1]-1) || z == 0 || z == (dim[2]-1)) {
                    *(labelImagePtr + x + (y + z*dim[1])*dim[0]) = 0;
                    continue;
                }
                value = *(imagePtr + x + (y + z*dim[1])*dim[0]);
                sum = 0;
                for (int k=-1; k<2; k++) for (int j=-1; j<2; j++) for (int i=-1; i<2; i++) {
                    if (k==0 && j == 0 && i==0) continue;
                    neib = *(imagePtr + (x+i) + ((y+j) + (z+k)*dim[1])*dim[0]);
                    sum += (abs(neib - value) <= epsilon)?1:0;
                }
                *(labelImagePtr + x + (y + z*dim[1])*dim[0]) = (sum >= 23)?0:255;
#endif
            }

    gmag->Delete();
    return labelImage;
}

double VolumeManager::computeEpsilon(vtkImageGradientMagnitude *gmag)
{
    vtkImageAccumulate *accum = vtkImageAccumulate::New();
    accum->SetInputData(gmag->GetOutput());
    double range[2];
    gmag->GetOutput()->GetScalarRange(range);
    fprintf(stderr, "@VolumeManager::computeEpsilon: Image magnitude range: [%f, %f]\n", range[0], range[1]);
    int irange[] = {int(range[0]), int(range[1])};
    accum->SetComponentExtent(irange[0], irange[1], 0, 0, 0, 0);
    accum->SetComponentOrigin(0, 0, 0);
    accum->SetComponentSpacing(1, 0, 0);
    accum->Update();
    int dims_accum[3], freq;
    accum->GetOutput()->GetDimensions(dims_accum);
    int *cumulativeHist = new int[dims_accum[0]];
    int old = 0;
    for(vtkIdType bin = 0; bin < dims_accum[0]; ++bin) //Compute cumulative histogram
    {
        freq = *(static_cast<int*>(accum->GetOutput()->GetScalarPointer(bin, 0, 0)));
        cumulativeHist[bin] = old + freq;
        old = cumulativeHist[bin];
        //fprintf(stderr, "%d -> %d ->%d\n", bin + irange[0], freq, cumulativeHist[bin]);
    }
    float targetFreq = LH_CUM_HISTOGRAM_FRACTION*cumulativeHist[dims_accum[0]-1];
    double epsilon = 0.0;
    for(vtkIdType bin = dims_accum[0]-1; bin >=0 ; bin--){
        if(cumulativeHist[bin] < targetFreq) {
            epsilon = bin + irange[0];
            break;
        }
    }

    fprintf(stderr, "@VolumeManager::computeEpsilon: epsilon: (%f, %f)\n", targetFreq, epsilon);

    accum->Delete();
    delete []cumulativeHist;

    return epsilon;
}

void VolumeManager::writeLH()
{
    //Write binary dump of data
    int dim[3];
    imageLH->GetDimensions(dim);
    long elems = dim[0]*dim[1]*dim[2]*2;

    std::ofstream stream("LH.bin", std::ios::binary );
    stream.write(reinterpret_cast<const char*>(imageLH->GetScalarPointer()), elems*imageLH->GetScalarSize());
    stream.close();

    //Write header information
    stream.open("LH.hdr");
    stream <<"Volume size: "<<dim[0] <<" "<<dim[1]<<" "<<dim[2]<<"\nData type: "<<imageLH->GetScalarTypeAsString();
    stream.close();
    fprintf(stderr, "Wrote binary LH file: LH.bin, and header LH.hdr\n");
}

void VolumeManager::clusteringKMeans()
{
    //TODO:
}
