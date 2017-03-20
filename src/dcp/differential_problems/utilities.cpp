#include <dcp/differential_problems/utilities.h> //"utilities.h"
#include <stdio.h>

#ifdef SprittlesShikhmurzaev
std::string savepath ("/u/laureandi/ifumagalli/dcp_test_output/SprittlesShikhmurzaev/");
#else
//std::string savepath ("/u/laureandi/ifumagalli/dcp_test_output/impostoVertCos/imposto/");
//std::string savepath ("/u/laureandi/ifumagalli/dcp_test_output/vertVel/");
//std::string savepath ("/u/laureandi/ifumagalli/dcp_test_output/vertVelFineDrho/");
//std::string savepath ("/u/laureandi/ifumagalli/dcp_test_output/poiseuilleImpostoBla/");
//std::string savepath ("/u/dati/numer/ifumagalli/dcp_test_output/okPoiseuilleLargoAngoloTer/");
//std::string savepath ("/u/laureandi/ifumagalli/dcp_test_output/impostoBla/");
std::string savepath ("/u/dati/numer/ifumagalli/dcp_test_output/tmp/");
#endif

/*#if defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev) && !defined(Yamamoto)
double        lx (13.6),
              ly (4*lx),
              ustar (0.02),
              yxratio (ly/lx);
std::size_t   nx (10),
              ny ((std::size_t) nx*yxratio);
double        TP_EPS (lx/nx+6e-16),
              re(0.25*lx * 0.81/1.95),
              st(0),
              ca(1.0e-2);
double        betaVal(36);

#elif defined(SprittlesShikhmurzaev) && !defined(GerbeauLelievre) && !defined(Yamamoto)
double        lx (1),
              ly (2),
              ustar (0.02),
              yxratio (ly/lx);
std::size_t   nx (20),
              ny ((std::size_t) nx*yxratio);
double        TP_EPS (lx/nx+6e-16),
              re(100),
              st(0),
              ca(1.0e-2);
double        betaVal(36);

#elif defined(Yamamoto) && !defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev)
  #ifdef vanMourik
  double  lx (1e-2),
  #else
double        lx (0.92e-3),
  #endif
              ly (lx),
  #ifdef vanMourik
                ustar (0.23),
  #else
              ustar (0.1),
  #endif
              yxratio (ly/lx);
std::size_t   nx (16),
  #ifdef vanMourik
              ny ((std::size_t) nx*yxratio);
  #else
              ny ((std::size_t) 2*nx*yxratio);
  #endif
double        TP_EPS (lx/nx+6e-16),
              rho (1115),
              re(0),
              st(0),
              ca(0);
double        betaVal(36),
              cosThetaSVal (cos(69.8*3.14159265/180.0));

#endif*/

bool compare (std::tuple<std::size_t,double,double> p1, std::tuple<std::size_t,double,double> p2)
// p = [idx, xcoord, ycoord]
{
    return ( ( std::get<2>(p1) < std::get<2>(p2) ) || ( std::get<2>(p1) == std::get<2>(p2) && std::get<1>(p1) < std::get<1>(p2) ) );
}
bool compareDofs (std::tuple<std::size_t,double,double,std::size_t> p1, std::tuple<std::size_t,double,double,std::size_t> p2)
// p = [idx, xcoord, ycoord, component]
{
    return ( ( std::get<2>(p1) < std::get<2>(p2) ) || ( std::get<2>(p1) == std::get<2>(p2) && std::get<1>(p1) < std::get<1>(p2) ) 
             || ( std::get<2>(p1) == std::get<2>(p2) && std::get<1>(p1) == std::get<1>(p2) && std::get<3>(p1) < std::get<3>(p2) ) );
}
bool compareFirst (std::tuple<std::size_t,double,double> p1, std::tuple<std::size_t,double,double> p2)
{
    return ( std::get<0>(p1) < std::get<0>(p2) );
}
std::size_t getFirst (const std::tuple<std::size_t,double,double> t)
{
    return std::get<0>(t);
}
std::size_t getFirstDofs (const std::tuple<std::size_t,double,double,std::size_t> t)
{
    return std::get<0>(t);
}
std::pair<double,double> getDofsCoords (const std::tuple<std::size_t,double,double,std::size_t> t)
{
    return std::make_pair(std::get<1>(t), std::get<2>(t));
}

void dcp::print2csv (const dolfin::Function & fun, const std::string & valname, const std::string & filename, const std::vector<dolfin::la_index>::const_iterator dofsBegin, const std::vector<dolfin::la_index>::const_iterator dofsEnd, const std::vector<double> & dofsCoords)
{
    const dolfin::GenericVector & vec (* fun.vector());
    std::ofstream file (filename.c_str(), std::ofstream::out);
    file.precision (15);
    file << "x,y,z,"+valname+"," << std::endl;
    file << std::scientific;
      // NB keep std::scientific, otherwise importing csv in ParaView may yield (mysterious) errors
    for (std::vector<dolfin::la_index>::const_iterator it (dofsBegin); it!=dofsEnd; ++it)
      file << dofsCoords[2*(*it)] << "," << dofsCoords[2*(*it)+1] << ",0," << vec[*it] << "," << std::endl;
    file.close();
}
void dcp::print2csv (const dolfin::Function & fun, std::string & filename, const std::vector<dolfin::la_index>::const_iterator dofsBegin, const std::vector<dolfin::la_index>::const_iterator dofsEnd, const std::vector<double> & dofsCoords)
{
    const dolfin::GenericVector & vec (* fun.vector());
    std::ofstream file (filename.c_str(), std::ofstream::out);
    file.precision (15);
    file << "x,y,z,scalar," << std::endl;
    file << std::scientific;
      // NB keep std::scientific, otherwise importing csv in ParaView may yield (mysterious) errors
    for (std::vector<dolfin::la_index>::const_iterator it (dofsBegin); it!=dofsEnd; ++it)
      file << dofsCoords[2*(*it)] << "," << dofsCoords[2*(*it)+1] << ",0," << vec[*it] << "," << std::endl;
    file.close();
}
void dcp::print2csv (const dolfin::Function & fun, std::string filename, const std::vector<std::size_t>::const_iterator dofsBegin, const std::vector<std::size_t>::const_iterator dofsEnd, const std::vector<double> & dofsCoords)
{
    const dolfin::GenericVector & vec (* fun.vector());
    std::ofstream file (filename.c_str(), std::ofstream::out);
    file.precision (15);
    file << "x,y,z,scalar" << std::endl;
    for (std::vector<std::size_t>::const_iterator it (dofsBegin); it!=dofsEnd; ++it)
      file << dofsCoords[2*(*it)] << "," << dofsCoords[2*(*it)+1] << ",0," << vec[*it] << std::endl;
    file.close();
}
void dcp::print2csv (const dolfin::Function & fun, std::string filename, const std::vector<std::size_t>::const_iterator dofsBegin1, const std::vector<std::size_t>::const_iterator dofsBegin2, const std::vector<std::size_t>::const_iterator dofsEnd, const std::vector<double> & dofsCoords)
{
    const dolfin::GenericVector & vec (* fun.vector());
    std::ofstream file (filename.c_str(), std::ofstream::out);
    file.precision (15);
    file << "x coord, y coord, z coord, uX, uY" << std::endl;
    std::vector<std::size_t>::const_iterator it2 (dofsBegin2);
    for (std::vector<std::size_t>::const_iterator it1 (dofsBegin1); it1!=dofsEnd; ++it1, ++it2)
      file << dofsCoords[2*(*it1)] << ", " << dofsCoords[2*(*it1)+1] << ", 0, " << vec[*it1] << ", " << vec[*it2] << std::endl;
    file.close();
}

/*class NormalAtVertices : public dolfin::Expression
{
public:

  NormalAtVertices (std::shared_ptr<const dolfin::Mesh> mesh) :
    dolfin::Expression (2),
    mesh_ (mesh)
  {}

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x,
            const ufc::cell& ufc_cell) const
  {
    dolfin_assert(ufc_cell.local_facet >= 0);

    dolfin::Cell cell(mesh, ufc_cell.index);
    dolfin::Point n = cell.normal(ufc_cell.local_facet);

    values[0] = n[0];
    values[1] = n[1];
  }

 private:

   std::shared_ptr<const dolfin::Mesh> mesh_;

 };*/
