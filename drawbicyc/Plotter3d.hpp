#pragma once
#include "Plotter2d.hpp"
#include "common.hpp"

// 3d version
template<>
class Plotter_t<3> : public PlotterCommon 
{
  public:
  using arr_t = blitz::Array<float, 3>;
  blitz::Array<int, 2> k_i;
  blitz::thirdIndex LastIndex;

  protected:
  using parent_t = PlotterCommon;
  hsize_t n[3];
  enum {x, y, z};
  arr_t tmp;
  blitz::Range yrange;

  public:

  auto h5load(
    const string &file, 
    const string &dataset
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    parent_t::h5load(file, dataset);
    this->h5s.getSimpleExtentDims(n, NULL);
  
    hsize_t 
      cnt[3] = { n[x],  n[y],  n[z] }, 
      off[3] = { 0,     0,     0    };
    this->h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[3] = {
      hsize_t(tmp.extent(0)), 
      hsize_t(tmp.extent(1)), 
      hsize_t(tmp.extent(2)) 
    };
    this->h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);

    return blitz::safeToReturn(tmp + 0);
  }

  auto h5load_timestep(
    const string &file, 
    const string &dataset,
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    return h5load(timestep_file, dataset);
  }

  auto horizontal_mean(
    const arr_t &data
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 1>() + 0))
  {
    using namespace blitz::tensor;
    auto tmp = blitz::mean(data(i,k,j), k);
    blitz::Array<float, 2> mean2d(tmp);
    auto tmp2 = blitz::mean(mean2d(j,i), j);
    blitz::Array<float, 1> mean(tmp2);
    return blitz::safeToReturn(mean + 0);
  }

  auto horizontal_sum(
    const arr_t &data
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 1>() + 0))
  {
    using namespace blitz::tensor;
    auto tmp = blitz::sum(data(i,k,j), k);
    blitz::Array<float, 2> mean2d(tmp);
    auto tmp2 = blitz::sum(mean2d(j,i), j);
    blitz::Array<float, 1> mean(tmp2);
    return blitz::safeToReturn(mean + 0);
  }

  template <class gp_t, class data_t>
  void plot(gp_t &gp, const data_t &data)
  {
  //  throw std::runtime_error("3d fields plotting doesn't work yet");
    blitz::Array<float, 3> tmp3d(data);
    using namespace blitz::tensor;
    // select a slize in second dimension to average over
    auto tmp3dslice = tmp3d(blitz::Range::all(), yrange, blitz::Range::all());
    auto tmp2d = blitz::mean(tmp3dslice(i,k,j), k); // mean over second dimension
    blitz::Array<float, 2> tmp(tmp2d);

    gp << "set xrange [0:" << tmp.extent(0)-1 << "]\n";
    gp << "set yrange [0:" << tmp.extent(1)-1 << "]\n";
    gp << "splot '-' binary" << gp.binfmt(tmp.transpose(blitz::secondDim, blitz::firstDim)) << " scan=yx origin=(0,0,0) with image failsafe notitle\n";
    gp.sendBinary(tmp);
  }

  //ctor
  Plotter_t(const string &file):
    parent_t(file)
  {

    po::options_description opts3d("profile plotting options");
    opts3d.add_options()
      ("y0", po::value<int>()->default_value(0) , "index of first cell in y dimension to vaerage over")
      ("y1", po::value<int>()->default_value(0) , "index of last cell in y dimension to vaerage over")
    ;
    po::variables_map vm; 
    handle_opts(opts3d, vm);
    int y0 = vm["y0"].as<int>();
    int y1 = vm["y1"].as<int>();
    if(y1!=0 && y0!=0)
      yrange = blitz::Range(y0, y1);
    else
      yrange = blitz::Range::all();

    // read number of timesteps
    this->h5f.openDataSet("T").getSpace().getSimpleExtentDims(n, NULL);
    this->map["t"] = n[0];

    // read number of cells
    this->h5f.openDataSet("X").getSpace().getSimpleExtentDims(n, NULL); // X gives cell-border coordinates (+1)
    this->map["x"] = n[0]-1;
    this->map["y"] = n[1]-1;
    this->map["z"] = n[2]-1;
    tmp.resize(n[0], n[1], n[2]);
    k_i.resize(n[0]-1, n[1]-1);

    // read dx,dy,dz
    h5load(file + "/const.h5", "X");
    this->map["dx"] = tmp(1,0,0) - tmp(0,0,0);
    h5load(file + "/const.h5", "Y");
    this->map["dy"] = tmp(0,1,0) - tmp(0,0,0);
    h5load(file + "/const.h5", "Z");
    this->map["dz"] = tmp(0,0,1) - tmp(0,0,0);
    this->CellVol = this->map["dx"] * this->map["dy"] * this->map["dz"];
    this->DomainSurf = this->map["dx"] * this->map["dy"] * this->map["x"] * this->map["y"];


    // other dataset are of the size x*z, resize tmp
    tmp.resize(n[0]-1, n[1]-1, n[2]-1);
  }
};

