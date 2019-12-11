#include "Core/Utility/Stat.h"
#include "Core/Utility/Type.h"
#include "Core/Utility/Value.h"
#include "Core/Utility/ValueArray.h"
#include <kvs/Message>
#include <kvs/StructuredVolumeObject>
#include <kvs/StructuredVolumeImporter>
#include <kvs/PolygonObject>
#include <kvs/OrthoSlice>
#include <kvs/HydrogenVolumeData>
#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/ScreenCaptureEvent>
#include <kvs/Stat>
#include <kvs/ImageObject>
#include <kvs/ImageRenderer>
#include <kvs/Endian>
#include <kvs/ColorImage>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <kvs/RayCastingRenderer>
#include <chrono>

static const int NX = 301;
static const int NY = 301;
static const int NZ = 50;
static const int SIZE = NX*NY*NZ;

kvs::ValueArray<float> createCoords(int nx, int ny, int nz);

//static const kvs::ValueArray<float> coords = createCoords(NX, NY, NZ);

enum Parameter {
  U,   // 東西風 : X-wind component (m s-1)
  V,   // 南北風 : Y-wind component (m s-1)
  W,   // 鉛直風 : Z-wind component (m s-1)
  T,   // 気温   : Temperature (K)
  P,   // 気圧   : Pressure (Pa)
  QV,  // 水蒸気混合比 : Water vapor mixing ratio (kg kg-1) 
  QC,  // 雲水混合比   : Cloud water mixing ratio (kg kg-1)
  QR,  // 雨混合比     : Rain mixing ratio (kg kg-1)
  QCI, // 雲氷混合比   : Cloud ice mixing ratio (kg kg-1)
  QS,  // 雪混合比     : Snow mixing ratio (kg kg-1)
  QG   // あられ混合比 : Graupel mixing ratio (kg kg-1)
};

kvs::ValueArray<float> createCoords(int nx, int ny, int nz) {
  // x方向は1.09545294622*10^-3°間隔
  std::vector<float> x(nx);
  float radius_earth = 6360000.0; // [m]
  // 基準地点=(東経133.590905 , 北緯33.51538902)
  float latitude = 33.5153;
  float longitude_interval = 0.00109545294622;
  float radius_point = radius_earth * cos(latitude/360.0);
  float circumference_point = radius_point * 2.0 * 3.141592;

  for (int i = 0; i < nx; ++i) {
    x[i] = (i-nx/2.) * circumference_point * longitude_interval / 360.0;
  }
  // y方向は8.99279260651*10^-4°間隔
  float latitude_interval = 0.000899279260;
  std::vector<float> y(ny);
  for (int i = 0; i < ny; ++i) {
    y[i] = (i-ny/2.) * radius_earth * 2.0 * 3.141592 * latitude_interval / 360.0;
  }

  // z方向 : (単位はメートル[m])
  //   1   -20.0000  2    20.0000  3    60.0000  4   118.0000  5   194.0000
  //   6   288.0000  7   400.0000  8   530.0000  9   678.0000 10   844.0000
  //  11  1028.0000 12  1230.0000 13  1450.0000 14  1688.0000 15  1944.0000
  //  16  2218.0000 17  2510.0000 18  2820.0000 19  3148.0000 20  3494.0000
  //  21  3858.0000 22  4240.0000 23  4640.0000 24  5058.0000 25  5494.0000
  //  26  5948.0000 27  6420.0000 28  6910.0000 29  7418.0000 30  7944.0000
  //  31  8488.0000 32  9050.0000 33  9630.0000 34 10228.0000 35 10844.0000
  //  36 11478.0000 37 12130.0000 38 12800.0000 39 13488.0000 40 14194.0000
  //  41 14918.0000 42 15660.0000 43 16420.0000 44 17198.0000 45 17994.0000
  //  46 18808.0000 47 19640.0000 48 20490.0000 49 21358.0000 50 22244.0000
  std::vector<float> z = {
        -20.0000  ,    20.0000  ,    60.0000  ,   118.0000  ,   194.0000
     ,   288.0000  ,   400.0000  ,   530.0000  ,   678.0000 ,   844.0000
    ,  1028.0000 , 1230.0000 ,  1450.0000 ,  1688.0000 ,  1944.0000
    ,  2218.0000 , 2510.0000 ,  2820.0000 ,  3148.0000 ,  3494.0000
    ,  3858.0000 , 4240.0000 ,  4640.0000 ,  5058.0000 ,  5494.0000
    ,  5948.0000 , 6420.0000 ,  6910.0000 ,  7418.0000 ,  7944.0000
    ,  8488.0000 , 9050.0000 ,  9630.0000 , 10228.0000 , 10844.0000
    , 11478.0000 ,12130.0000 , 12800.0000 , 13488.0000 , 14194.0000
    , 14918.0000 ,15660.0000 , 16420.0000 , 17198.0000 , 17994.0000
    , 18808.0000 ,19640.0000 , 20490.0000 , 21358.0000 , 22244.0000
  };

  // normalize
  // 倍率は適当
  auto norm = std::max(x[nx-1], std::max(y[ny-1], z[nz-1]));
  for (int i = 0; i < nx; ++i) {
    x[i] /= norm/301.0;
  }
  for (int i = 0; i < ny; ++i) {
    y[i] /= norm/301.0;
  }
  for (int i = 0; i < nz; ++i) {
    z[i] /= norm/301.0;
  }
  
  x.insert(x.end(),y.begin(), y.end());
  x.insert(x.end(),z.begin(), z.end());
  return kvs::ValueArray<float>(x);
}

kvs::ValueArray<float> loadValueArray(std::ifstream &ifs, int size) {
  kvs::ValueArray<float> array(size);
  ifs.read((char*)array.data(), size*4);
  
  // convert endian
  if (kvs::Endian::IsLittle()) {
    for (auto&& i: array) {
      kvs::Endian::Swap(&i);
    }
  }
  return array;
}

kvs::StructuredVolumeObject *load(std::ifstream &ifs, kvs::ValueArray<float> array) {
  kvs::ValueArray<float> kvs_value = loadValueArray(ifs, SIZE);
  
  kvs::StructuredVolumeObject *vol = new kvs::StructuredVolumeObject();
  vol->setGridTypeToUniform();
  vol->setVeclen(1);
  vol->setResolution(kvs::Vector3ui(NX, NY, NZ));
  vol->setValues(array);
  vol->updateMinMaxValues();

  //vol->setGridTypeToRectilinear();
  //vol->setCoords(coords.clone());
  //vol->updateMinMaxCoords();
  return vol;
}

kvs::StructuredVolumeObject *loadData(std::string filename, Parameter p) {
  std::ifstream ifs (filename, std::ios::in | std::ios::binary);
  if (ifs.fail()) {
    std::cerr << "file not found" << std::endl;
    return nullptr;
  }

  ifs.seekg(SIZE*4*p, std::ios_base::beg);
  kvs::StructuredVolumeObject* vol = load(ifs, loadValueArray(ifs, SIZE));
   
  ifs.close();
  return vol;
}

// https://qiita.com/Ushio/items/f5630d87f55c7afa984e
class Kahan {
public:
    Kahan() {}
    Kahan(double value) {}

    void add(double x) {
        double y = x - _c;
        double t = _sum + y;
        _c = (t - _sum) - y;
        _sum = t;
    }
    void operator=(double x) {
        _sum = x;
        _c = 0.0;
    }
    void operator+=(double x) {
        add(x);
    }
    operator double() const {
        return _sum;
    }
private:
    double _sum = 0.0;
    double _c = 0.0;
};

class OnlineVariance {
public:
    void addSample(double x) {
        _n++;
        double delta = x - _mean;
        _mean += delta / _n;
        double delta2 = x - _mean;
        _M2 += delta * delta2;
    }
    double variance() const {
        // return _M2 / (_n - 1);
        return _M2 / _n;
    }
    double avarage() const {
        return _mean;
    }
private:
    Kahan _M2 = 0.0;
    Kahan _mean = 0.0;
    int _n = 0;
};

class OnlineArrayVariance {
public:
    OnlineArrayVariance(size_t s) :_size(s) {
      _M2 = std::vector<Kahan>(_size);
      _mean = std::vector<Kahan>(_size);
    };
    void addArray(kvs::ValueArray<float> array) {
      int i = 0;
      _n++;
      for (auto&& a: array) {
        addSample(a, i++);
      }
    }
    void addSample(double x, int i) {
        double delta = x - _mean[i];
        _mean[i] += delta / _n;
        double delta2 = x - _mean[i];
        _M2[i] += delta * delta2;
    }
    kvs::ValueArray<float> variance() const {
        // return _M2 / (_n - 1);
      kvs::ValueArray<float> v(_size);
        for (size_t i = 0; i < _size; ++i) {
          v[i] = _M2[i] / _n;
        }
        return v;
    }
    kvs::ValueArray<float> avarage() const {
      kvs::ValueArray<float> v(_size);
      for (size_t i = 0; i < _size; ++i) {
        v[i] = _mean[i];
      }
      return v;
    }
private:
    std::vector<Kahan> _M2;
    std::vector<Kahan> _mean;
    size_t _size;
    int _n = 0;
};

//int main() {
//  OnlineArrayVariance a(5);
//  std::vector<double> v1 = {1, 2, 3, 4, 5};
//  std::vector<double> v2 = {2, 2, 3, 4, 5};
//  std::vector<double> v3 = {3, 2, 3, 4, 5};
//  std::vector<double> v4 = {4, 2, 3, 4, 5};
//  a.addArray(v1);
//  a.addArray(v2);
//  a.addArray(v3);
//  a.addArray(v4);
//  auto average = a.avarage();
//  auto variance = a.variance();
//  for (int i = 0; i < 5; ++i) {
//    std::cout << i << ": average: " << average[i] << " " << "variance: " << variance[i] << std::endl;
//  }
//
//}

int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );

    kvs::StructuredVolumeObject* volume = loadData("ensemble_data/gs0030.bin", Parameter::QV);

    kvs::StructuredVolumeObject* volume2 = loadData("ensemble_data/gs0030.bin", Parameter::V);
    std::cout << kvs::Stat::Corr(volume->values().asValueArray<float>(), volume2->values().asValueArray<float>()) << std::endl;

    OnlineArrayVariance oav(SIZE);
   // for (int i = 0; i < 20; ++i) {
   //   std::string renban_name = make_renban(i);
   //   kvs::StructuredVolumeObject* vol = loadData(renban_name, Parameter::V);
   //   oav.addArray(vol->values().asValueArray<float>());
   //   delete vol;
   // }

    oav.addArray(volume->values().asValueArray<float>());
    oav.addArray(volume2->values().asValueArray<float>());
    kvs::StructuredVolumeObject *average_vol = new kvs::StructuredVolumeObject();
    average_vol->setGridTypeToUniform();
    average_vol->setVeclen(1);
    average_vol->setResolution(kvs::Vector3ui(NX, NY, NZ));
    average_vol->setValues(oav.avarage());
    average_vol->updateMinMaxValues();

    kvs::StructuredVolumeObject *variance_vol = new kvs::StructuredVolumeObject();
    variance_vol->setGridTypeToUniform();
    variance_vol->setVeclen(1);
    variance_vol->setResolution(kvs::Vector3ui(NX, NY, NZ));
    variance_vol->setValues(oav.variance());
    variance_vol->updateMinMaxValues();

    if ( !volume )
    {
        kvsMessageError( "Cannot create a structured volume object." );
        return( false );
    }

    // 指定されたz軸断面の位置
    const float p = 15;

    const kvs::OrthoSlice::AlignedAxis a = kvs::OrthoSlice::ZAxis;
    const kvs::TransferFunction t( 256 );
    kvs::PolygonObject* object = new kvs::OrthoSlice( variance_vol, p, a, t );
    if ( !object )
    {
        kvsMessageError( "Cannot create a polygon object." );
        delete volume;
        return( false );
    }

    kvs::glsl::RayCastingRenderer* ren = new kvs::glsl::RayCastingRenderer();
    ren->enableShading();

    kvs::glut::Screen screen( &app );
    //screen.registerObject( object );
    screen.registerObject( variance_vol, ren );
    screen.setGeometry( 0, 0, 512, 512 );
    screen.setTitle( "kvs::OrthoSlice" );

    screen.create();
    screen.paintEvent();
    screen.paintEvent();
    kvs::ColorImage image = screen.scene()->camera()->snapshot();
    image.write("image.bmp");

    //kvs::ScreenCaptureEvent capture_event;
    ////capture_event.setFilename("file1");
    //screen.addEvent(&capture_event);
    //screen.show();

    return( app.run() );
    //app.quit();
    delete volume;
}
