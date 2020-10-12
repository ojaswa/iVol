#ifndef VOLUMEMANAGER_H
#define VOLUMEMANAGER_H

#include <vtkSmartPointer.h>

#define TINY 1e-12
#define LH_CUM_HISTOGRAM_FRACTION 0.9
#define LH_BASIC_THRESHOLDING 0
#define LH_INTERPOLATE_IMAGE_WHILE_RK_TACKING 1

class vtkImageGradientMagnitude;
class vtkImageData;
class vtkImageInterpolator;

class VolumeManager
{
public:
    VolumeManager();
    void setData(vtkImageData*);
    vtkImageData* getData();
    vtkImageData* getLH();
    vtkImageData* getLHDensity() { return imageLHDensityPlot; }
    template<class T> void computeLH();
    void writeLH();
    void clusteringKMeans();

private:
    //Data
    vtkSmartPointer<vtkImageData> imageData;
    vtkSmartPointer<vtkImageData> imageLH; // 32-bit LH for 16-bit volume, 16-bit LH for 8-bit volume. Stores LH pair per voxel: LS bytes store L value, MS byts store H value
    vtkSmartPointer<vtkImageData> imageLHDensityPlot;

    //Functions
    double computeEpsilon(vtkImageGradientMagnitude *gmag);
    template<class T> vtkImageData* computeGradientLabelImageLH(vtkImageData *smoothedImage);
    double gradientTrack(vtkImageInterpolator *image, vtkImageInterpolator *gradient, double delta_t, int x, int y, int z);
    void rungeKuttaEstimate(vtkImageInterpolator *gradient, double delta_t, double *pos);
};

// Extern template declarations
extern template void VolumeManager::computeLH<unsigned char>();
extern template void VolumeManager::computeLH<unsigned short>();

#endif // VOLUMEMANAGER_H
