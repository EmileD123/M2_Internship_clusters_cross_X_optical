#pragma once

/* ---------------------------------------------------------------   
   Classes  AstroMatch  
     Association d'objets sur la base de coordonnes spheriques 
     typiquement des sources celestes  avec coordonnes (ra,dec) 
   --- SOPHYA software - SkyT module ---

   R. Ansari , C. Magneville   2024 

   SOPHYA project - (C) Univ. Paris-Saclay, IJCLab CNRS/IN2P3, CEA-Irfu

   (C) UPS+LAL IN2P3/CNRS     (C) IRFU/CEA 
   (C) Univ. Paris-Saclay     (C) IJCLab CNRS/IN2P3 
   --------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//--- include SOPHYA 
#include "machdefs.h"
#include "pexceptions.h"
#include "tarray.h"
#include "angle.h"
#include "longlat.h"
//-- liste d'objets sur cartes spheriques 
#include "listosphmap.h"

#include "progbar.h"

namespace SOPHYA {

 
/*! 
  \class AssRecord
  \ingroup SkyT
  \brief A simple class to hold source association information between objects from two celestial object lists 

  The object identifier in AssRecord (idm or getId() ) correspond to the object index in the std::vector<LongitudeLatitude> 
  containing the primary list of celestial source positions. Each AssRecord, when filled, contains then the list 
  of secondary objects within the angular association distance, with the following information: 
  object index in the secondary list, angular distance and redshift difference if applicable.
  The AssRecord object holds also the object identifier with the minimum angular distance. 

  \sa AstroMatch
*/
class AssRecord {
public:
  AssRecord(size_t idm=0)  : idm_(idm), idass_(0), dstass_(9.e19), delz_(0.), maxassdist_(0.) 
    {  }
  AssRecord(AssRecord const& a)
    : idm_(a.idm_), idass_(a.idass_), dstass_(a.dstass_), delz_(a.delz_),
      vass_(a.vass_), maxassdist_(a.maxassdist_)
    {  }
  AssRecord& operator = (AssRecord const& a) {
    idm_=a.idm_; 
    idass_=a.idass_;   dstass_=a.dstass_; delz_=a.delz_;
    vass_=a.vass_;
    maxassdist_=a.maxassdist_;
    return *this; 
  }
  //! return the identifier of the object in the primary list (index in the vector) 
  inline size_t& getId() { return idm_; }
  /*! 
    \brief return the identifier of the associated object from the secondary list (index in the vector) 

    \warning the returned identifier (vector index) is valid only if the association count is larger than zero 
  */
  inline size_t& getAssId() { return  idass_; }
  //! return the association  count (number of objects with the association distance) 
  inline size_t getAssCount() const { return vass_.size(); }
  //! return true if association count > 0 
  inline bool okAss() const { return vass_.size()>0; }
  //! return the association distance (angle in radian) between the object in the primary and the best angular match in the secondary list 
  inline double& getAssDistance() { return dstass_; }
  //! return the redshift distance between the object in the primary list and the best matched object in the secondary list  
  inline double& getRedshiftDiff() { return delz_; }
  //! add an object to the list of associated objects (i.e. within the angular distance and redshift difference if applicable)
  inline void addAssObject(size_t idass, double dass, double dz=0.)
  { vass_.push_back(pair<size_t, pair<double,double> >(idass, pair<double,double>(dass, dz)) ); }
  //! return the n-th associated object - NO bound checking performed (0<=n<AssCcount) 
  inline pair<size_t, pair<double,double> > getAssObject(size_t n) const { return vass_[n]; }
  //! return the n-th associated object Id - NO bound checking performed (0<=n<AssCcount) 
  inline size_t getAssObjectId(size_t n) const { return vass_[n].first; }
  //! return the n-th associated object angular distance - NO bound checking performed (0<=n<AssCcount) 
  inline double getAssObjectDistance(size_t n) const { return vass_[n].second.first; }
  //! return the n-th associated object redshift difference - NO bound checking performed (0<=n<AssCcount) 
  inline double getAssObjectRedshiftDiff(size_t n) const { return vass_[n].second.second; }
  //! return the maximum association distance (angle in radian) for this object 
  inline double& getMaxxAssDistance() { return maxassdist_; }
  //! return the maximum association distance (angle in radian) for this object (const version)
  inline double getMaxxAssDistance() const { return maxassdist_; }
protected:
  size_t idm_;
  size_t idass_;
  double dstass_, delz_;
  std::vector< pair<size_t, pair<double,double> > >  vass_;  // list of all associated objects within the angular distance and redshift distance
  double maxassdist_;
};

/*!
  \class AsM_MaxAssDistInterface
  \ingroup SkyT
  \brief An abstract class defining the interface classes used to specify a per object maximum association distance. 

  If it is necessary to specify a maximum association distance depending on the object in the primary 
  or main object list, one has to define a class inheriting from AsM_MaxAssDistInterface, implementing 
  getMaxAssDistance() and pass an object of this class as the argument of the appropriate 
  AstroMatch::doMatch() method.
 */
class AsM_MaxAssDistInterface {
public:
  AsM_MaxAssDistInterface() { }
  virtual ~AsM_MaxAssDistInterface() { }
  /*! \brief return maximum association distance for object \b numom in the main or primary object list 

    \param numom : object index in the primary object list 
    \param z : redshift for the object (if applicable)
    \param dz : returned value, maximum redshift difference acceptable for the association (if applicable)

    The return value should be the maximum association angular distance in radian.
  */
  virtual double getMaxAssDistance(size_t numom, double z, double & dz) const = 0;
  //! return true if getMaxAssDistance() provides also the maximum redshift difference dz
  virtual bool reqRedshiftMatch() const = 0;
};

  
//--------------------------------------------------------------------------
//! Class to perform source matching based on spherical coordinates  
class AstroMatch {
 public:
  //! constructor - associates objects from list lls to object in llm, optionaly with a list of redshifts 
  AstroMatch(std::vector<LongitudeLatitude>& llm, std::vector<LongitudeLatitude>& lls, Angle skylistresol,
	     std::vector<double>* z_reds_m=nullptr, std::vector<double>* z_reds_s=nullptr);
  //! constructor - associates objects from list lls to object in llm, and a region of sky, optionaly with a list of redshifts 
  AstroMatch(std::vector<LongitudeLatitude>& llm, std::vector<LongitudeLatitude>& lls,
	     LongitudeLatitude const& llmin, LongitudeLatitude const& llmax, Angle skylistresol,
	     std::vector<double>* z_reds_m=nullptr, std::vector<double>* z_reds_s=nullptr);

  //! associate objects and fills the association record AssRecord for all objects in the primary list - No cut/selection on redshift differences 
  size_t doMatch(Angle assmaxdist, std::vector<AssRecord>& vassr)
  {
    return doMatchP(assmaxdist, -1., vassr);
  }
  /*! \brief associate objects and fills the association record AssRecord for all objects in the primary list
    Only objects with redshift difference |delta-z|<=dz are matched ( dz > 1.e-9)    */
  size_t doMatch(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr)
  {
    return doMatchP(assmaxdist, dz, vassr);
  }    
  /*! \brief associate objects and fills the association record AssRecord for all objects in the primary list with
  a per object maximum association distance defined by \b am_maxdist (and redshift difference if applicable) */
  size_t doMatch(AsM_MaxAssDistInterface const & am_maxdist, std::vector<AssRecord>& vassr)
  {
    return doMatchP(Angle(skylistresol_.ToRadian()*0.15), 0.1, vassr, &am_maxdist);
  }
  
  //! return a spherical map with the source count per pixel (valid if full sky) 
  SphereThetaPhi<int_4> getCountMapS() const
    { return maplls_.getCountMap(); }
  //! return a spherical ECP map with the source count per pixel (valid if partial sky) 
  SphereECP<int_4> getCountMapS_ECP() const
    { return ecpmaplls_.getCountMap(); }
  //! diplay of progress bar : ProgBarM_None , ProgBarM_Percent, ProgBarM_Time
  inline void setProgressBarMode(ProgressBarMode pgm=ProgBarM_None) { pgm_mode_=pgm; } 
 protected:
  size_t doMatchP(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr,
		  AsM_MaxAssDistInterface const * am_maxdist = nullptr);

  //!  perform the match for a full sky coverage, i.e. using  SphThetaPhiList<size_t> maplls_
  size_t doMatch_TP(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr, size_t& loopcnt,
		    AsM_MaxAssDistInterface const * am_maxdist_p = nullptr);
  //!  perform the match for a secondary list covering fraction of the sky, i.e. using  SphThetaPhiList<size_t> maplls_SphECPList<size_t> ecpmaplls_
  size_t doMatch_ECP(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr, size_t& loopcnt,
		     AsM_MaxAssDistInterface const * am_maxdist_p = nullptr);

  std::vector<LongitudeLatitude>& llm_;  // premiere liste de longitude-latitude , liste principale 
  std::vector<LongitudeLatitude>& lls_;  // deuxieme liste de longitude-latitude
  Angle skylistresol_;
  std::vector<double>* z_redshift_m_;
  std::vector<double>* z_redshift_s_;
  bool fg_with_redshift_; 
  SphThetaPhiList<size_t> maplls_;
  bool fgecpmap_; 
  SphECPList<size_t> ecpmaplls_;
  std::vector<int_8> vmidx_; 
  ProgressBarMode pgm_mode_; 
};
 
} // Fin du namespace
