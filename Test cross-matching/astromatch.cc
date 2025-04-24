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

#include "astromatch.h"

using namespace std;
namespace SOPHYA {

/*! 
  \class AstroMatch
  \ingroup SkyT
  
  This class can be used to cross match two list of celestical objects, based on the angular distance 
  on the sky. In addition, it is possible to impose a maximum difference along the resdhift axis. 
  For each object in the first or main list (\b llm ), defined by their angular coordinates, the 
  AstroMatch::doMatch() method finds all objects in the second list (\b lls ), with an angular
  distance less than the limit \b assmaxdist .  

  In order to perform the association, the objects from the secondary list are distributed 
  over pixels on a spherical map, with a user specified spherical pixel angular resolution. 
  This method speeds up significantly the association process. If the first list has N object,
  and the second list M objects, the association process without prior distribution over 
  spherical pixels will take a time proportional to NxM, which becomes easily prohibitive 
  with N,M >= 10^6. Assuming a roughly unform distribution, and a angular sky map with P pixels, 
  the association time is decreased by a factor P, NxM/P. 

  The AstroMatch class constructor has an angular resolution parameter (\b skylistresol ) which defines 
  the angular size of spehrical map pixels over which the objects from the second list are distributed. 
  These pixel size should be larger than the association angular limit (\b assmaxdist ) by a factor 5.
  A factor 5-10 seems a reasonable lower limit for the ratio skylistresol/assmaxdist > 4. 
  Keep also in mind that fine sky pixels leads to large memory footprint for the process and possibly 
  to decreased performance. 
  For a typical few arcsec association distance, few arcminute pixels should perform well. 

  \sa SOPHYA::SphThetaPhiList<T> 
  \sa SOPHYA::SphECPList<T> 
  \sa SOPHYA::AssRecord
*/
  
/* --Methode-- */
/*!
  Constructor, assuming that objects are distributed over the full celestial sphere 

  \param llm : vector of sky angular position for the primary celestial object list 
  \param lls : vector of sky angular position for the secondary object list 
  \param skylistresol : angular resolution of object list, organised as a function of position on sky (SphericalMapList)
  \param z_reds_m : pointer to a vector of redshift value for the objects in the primary list (if not a nullptr) 
  \param z_reds_s : pointer to a vector of redshift value for the objects in the secondary list (if not a nullptr) 
*/
AstroMatch::AstroMatch(std::vector<LongitudeLatitude>& llm, std::vector<LongitudeLatitude>& lls, Angle skylistresol,
		       std::vector<double>* z_reds_m, std::vector<double>* z_reds_s)
  : llm_(llm), lls_(lls), skylistresol_(skylistresol), maplls_(skylistresol),
    z_redshift_m_(z_reds_m), z_redshift_s_(z_reds_s), fg_with_redshift_(false),
    fgecpmap_(false), ecpmaplls_(5) 
{
  if ((llm_.size() < 1)||(lls_.size()<1))
    throw ParmError("AstroMatch::AstroMatch() ERROR  empty main source or secondary source list ");

  if (z_redshift_m_ && z_redshift_s_) {
    fg_with_redshift_=true;
    if (z_redshift_m_->size() != llm_.size() )
      throw SzMismatchError("AstroMatch::AstroMatch() ERROR  provided redshift vector size != llm.size (main/first list) ");
    if (z_redshift_s_->size() != lls_.size() )
      throw SzMismatchError("AstroMatch::AstroMatch() ERROR  provided redshift vector size != lls.size (second list) ");
    
  }
  setProgressBarMode();
  for(size_t i=0; i<lls_.size(); i++) {
    maplls_.addObject(lls_[i],i);
  }
  cout <<"AstroMatch()/Info FullSky , source list sizes: main="<<llm_.size()<<" second="<<lls_.size()<<endl;
  if (fg_with_redshift_) cout<<" [ association with redshifts ] "<<endl;
}

/* --Methode-- */
/*!
  Constructor to be used when objects are distributed in a limited region of sky, defined by \b llmin-llmax 
  longitude-latitude limits.   

  \param llm : vector of sky angular position for the primary celestial object list 
  \param lls : vector of sky angular position for the secondary object list 
  \param llmin : lower limit (longitude-latitude or alpha, delta) of the sky region 
  \param llmax : upper limit (longitude-latitude or alpha, delta) of the sky region 
  \param skylistresol : angular resolution of object list, organised as a function of position on sky  (SphericalMapList)
  \param z_reds_m : pointer to a vector of redshift value for the objects in the primary list (if not a nullptr) 
  \param z_reds_s : pointer to a vector of redshift value for the objects in the secondary list (if not a nullptr) 
*/

AstroMatch::AstroMatch(std::vector<LongitudeLatitude>& llm, std::vector<LongitudeLatitude>& lls,
		       LongitudeLatitude const& llmin, LongitudeLatitude const& llmax, Angle skylistresol,
		       std::vector<double>* z_reds_m, std::vector<double>* z_reds_s)
  : llm_(llm), lls_(lls), skylistresol_(skylistresol), maplls_(5),
    z_redshift_m_(z_reds_m), z_redshift_s_(z_reds_s), fg_with_redshift_(false),
    fgecpmap_(true), ecpmaplls_(llmin,llmax,skylistresol) 
{
  if ((llm_.size() < 1)||(lls_.size()<1))
    throw ParmError("AstroMatch::AstroMatch(llLimits) ERROR  empty main source or secondary source list ");

  if (z_redshift_m_ && z_redshift_s_) {
    fg_with_redshift_=true;
    if (z_redshift_m_->size() != llm_.size() )
      throw SzMismatchError("AstroMatch::AstroMatch(llLimits) ERROR  provided redshift vector size != llm.size (main/first list) ");
    if (z_redshift_s_->size() != lls_.size() )
      throw SzMismatchError("AstroMatch::AstroMatchllLimits() ERROR  provided redshift vector size != lls.size (second list) ");
    
  }
  setProgressBarMode();
  vmidx_.resize(llm_.size(),false);
  size_t nokm=0;
  //  cout<<"*DBG*  Checking llm_ values in ECP map ... "<<endl; 
  for(size_t i=0; i<llm_.size(); i++) {
    vmidx_[i]=ecpmaplls_.getIndexNoExc(llm_[i]);
    if (vmidx_[i] >= 0)  nokm++; 
  }
  size_t noks=0;
  //  cout<<"*DBG*  Checking/filling lls_ values in ECP map ... ecpmaplls_.size()="<<ecpmaplls_.size()<<endl; 
  for(size_t i=0; i<lls_.size(); i++) {
    if (ecpmaplls_.getIndexNoExc(lls_[i])<0)  continue;
    ecpmaplls_.addObject(lls_[i],i);
    noks++; 
  }
  cout <<"AstroMatch()/Info PartialSky , source count InRegion/Total: main="<<nokm<<"/"<<llm_.size()
       <<" second="<<noks<<"/"<<lls_.size()<<endl;
  if (fg_with_redshift_) cout<<" [ association with redshifts ] "<<endl;
}  
/* Fonction pour verifier que la valeur de l'index est presente ou pas ds le vecteur - 
    renvoie true si nouvel index et l'ajoute ds le vecteur */
static inline bool check_new_index(vector<int_8>& vidx, int_8 nidx)
{
  for(size_t i=0; i<vidx.size(); i++)
    if (nidx == vidx[i])  return false;
  vidx.push_back(nidx);
  return true; 
}

/* Fonction pour decaler l'objet longlat selon theta,phi */
static LongitudeLatitude shift_longlat(LongitudeLatitude const& ll, double shiftrad, int jj, int ii)
{
  double theta = ll.Theta()+(double)jj*shiftrad;
  double cst=sin(theta);
  cst=(cst>1.e-3)?1./cst:1e3;
  double phi = ll.Phi()+(double)ii*shiftrad*cst;
  if (theta<0.) theta=0.;
  if (theta>Angle::OnePiCst()) theta=Angle::OnePiCst()-1.e-19;
  while (phi<0.) phi += Angle::TwoPiCst();
  while (phi>Angle::TwoPiCst()) phi -= Angle::TwoPiCst();
  return LongitudeLatitude(theta,phi);
}

/* --Methode-- */
/*!
  \param assmaxdist : object association maximum angular distance 
  \param dz : maximum redshift difference for associated objects  
  \param vassr : output vector, with one AssRecord for each object in the primary list, 
  containing the list of associated object from the secondary list, within 
  the angular distance \b assmaxdist.
  \param am_maxdist : if a valid pointer (not null) provided, it will be used to define the association maximum 
  angular distance and redshift difference, if applicable 

  \sa AssRecord

  \warning throws an exception if assmaxdist >  0.2*SourceListMapResolution
*/
size_t AstroMatch::doMatchP(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr, AsM_MaxAssDistInterface const * am_maxdist)
{
  
  if (!am_maxdist && (assmaxdist.ToRadian()>0.20*maplls_.getAngularResolution().ToRadian()) )
    throw ParmError("AstroMatch::doMatchP() ERROR assmaxdist > 0.2*SourceListMapResolution");
  
  double shiftrad = maplls_.getAngularResolution().ToRadian()/4.;   // 1/4 de la taille de pixel
  double amaxrad=assmaxdist.ToRadian();

  bool fgmatchredshift = fg_with_redshift_ && (dz>1.e-9);
  if (am_maxdist) fgmatchredshift = fg_with_redshift_ && am_maxdist->reqRedshiftMatch();
  
  cout<<"AstroMatch::doMatchP()/Info "<<(fgecpmap_?"Partial Sky/ECP":"Full Sky/ThetaPhi")<<" shift="<<Angle(shiftrad).ToDegree()
      << " AssDist="<<Angle(amaxrad).ToArcSec()<<" fgmatchredshift="<<(fgmatchredshift?"TRUE":"FALSE")<<endl;
  if (am_maxdist)
    cout << " ... using an AsM_MaxAssDistInterface object for max distance association , am_maxdist="
	 << hex << am_maxdist << dec << endl; 
  size_t loopcnt=0;
  size_t totcnt=0; 
  if (fgecpmap_) {
    totcnt=doMatch_ECP(assmaxdist, dz, vassr, loopcnt, am_maxdist);
  }
  else {
    totcnt=doMatch_TP(assmaxdist, dz, vassr, loopcnt, am_maxdist);
  }
  cout << "AstroMatch::doMatch()  main list source count="<<llm_.size()<<" TotCount (associated objects)="<<totcnt
       <<" TotLoopCount="<<loopcnt<<endl;
  double moydass=0.;
  size_t nwass=0;
  for(size_t m=0; m<vassr.size(); m++) {
    if (vassr[m].okAss()) {
      moydass+=vassr[m].getAssDistance();
      nwass++;
    }
  }
  if (nwass>0) moydass/=(double)nwass;
  cout<<" AstroMatch::doMatch()/Info : NbSrc w/ Ass="<<nwass<<" Mean-D-ass= "<<Angle(moydass).ToArcMin()<<endl;
  return totcnt; 
}
  
/* --Methode-- */
size_t AstroMatch::doMatch_TP(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr, size_t& loopcnt,
			      AsM_MaxAssDistInterface const * am_maxdist_p)
{
  if (!am_maxdist_p && (assmaxdist.ToRadian()>0.20*maplls_.getAngularResolution().ToRadian()) )
    throw ParmError("AstroMatch::doMatch_TP() ERROR assmaxdist > 0.2*SourceListMapResolution");

  double shiftrad = maplls_.getAngularResolution().ToRadian()/4.;   // 1/4 de la taille de pixel
  double amaxrad=assmaxdist.ToRadian();

  bool fgmatchredshift = fg_with_redshift_ && (dz>1.e-9);
  if (am_maxdist_p) fgmatchredshift = fg_with_redshift_ && am_maxdist_p->reqRedshiftMatch();

  vassr.resize(llm_.size());

  ProgressBar pgbar(llm_.size(), pgm_mode_);

  loopcnt=0; 
  size_t totcnt=0; 
  for(size_t m=0; m<llm_.size(); m++) {  // loop over all sources in the main list
    pgbar.update(m);
    //  if (m%20 == 0)  cout << " DEBUG*doMatch_TP() m="<<m<<" llm_[m]="<<llm_[m]<<endl; 
    vector<int_8> vidx;
    vidx.push_back(maplls_.getIndex(llm_[m])); // pixel de la carte ds lequel tombe la source de depart
    for(int jj=-1; jj<=1; jj++) {
      for(int ii=-1; ii<=1; ii++) {
	if ((jj==0)&&(ii==0))  continue;
	int_8 nidx=maplls_.getIndex(shift_longlat(llm_[m], shiftrad, jj, ii));
	if (maplls_.getCount(nidx)<1) continue;
	check_new_index(vidx, nidx);
      }
    }

    // on met a jour la distance maximum d'association si l'objet AsM_MaxAssDistInterface am_maxdist fourni
    if (am_maxdist_p) {
      amaxrad=am_maxdist_p->getMaxAssDistance(m, z_redshift_m_->at(m), dz);
    }
    
    // On a donc tous les pixels ou il faut boucler sur les sources ds lls
    AssRecord arec(m);
    arec.getMaxxAssDistance()=amaxrad;
    
    double dmin=9.e9;

    for(size_t i=0; i<vidx.size(); i++) {
      vector<size_t>& vsi= maplls_.getList(vidx[i]);
      for(size_t k=0; k<vsi.size(); k++) {
	loopcnt++;
	size_t ns=vsi[k];
	double sep=llm_[m].SepAngle(lls_[ns]);
	if (sep > amaxrad) continue;
	// Selection sur redshift si applicable
	double delz=0;
	if (fgmatchredshift) {
	  delz=z_redshift_m_->at(m)-z_redshift_s_->at(ns);
	  if (fabs(delz)>dz)  continue;
	}
	arec.addAssObject(vsi[k],sep,delz);  totcnt++;
	if (sep<dmin) {
	  dmin=sep;  arec.getAssId()=vsi[k];
	  arec.getAssDistance()=dmin;
	  arec.getRedshiftDiff()=delz;
	}
      }
    }
    vassr[m]=arec;
  }

  return totcnt; 
}

/* --Methode-- */
size_t AstroMatch::doMatch_ECP(Angle assmaxdist, double dz, std::vector<AssRecord>& vassr, size_t& loopcnt,
			       AsM_MaxAssDistInterface const * am_maxdist_p)
{
  if (!am_maxdist_p && (assmaxdist.ToRadian()>0.20*maplls_.getAngularResolution().ToRadian()) )
    throw ParmError("AstroMatch::doMatch_ECP() ERROR assmaxdist > 0.2*SourceListMapResolution");

  double shiftrad = maplls_.getAngularResolution().ToRadian()/4.;   // 1/4 de la taille de pixel
  double amaxrad=assmaxdist.ToRadian();

  bool fgmatchredshift = fg_with_redshift_ && (dz>1.e-9);
  if (am_maxdist_p) fgmatchredshift = fg_with_redshift_ && am_maxdist_p->reqRedshiftMatch();

  vassr.resize(llm_.size());

  ProgressBar pgbar(llm_.size(), pgm_mode_);

  loopcnt=0; 
  size_t totcnt=0; 
  for(size_t m=0; m<llm_.size(); m++) {  // loop over all sources in the main list
    pgbar.update(m);
    if (vmidx_[m]<0)  {
      vassr[m]=AssRecord(m);
      continue;
    }
    vector<int_8> vidx;
    vidx.push_back(vmidx_[m]); // pixel de la carte ds lequel tombe la source de depart
    for(int jj=-1; jj<=1; jj++) {
      for(int ii=-1; ii<=1; ii++) {
	if ((jj==0)&&(ii==0))  continue;
	int_8 nidx=ecpmaplls_.getIndexNoExc(shift_longlat(llm_[m], shiftrad, jj, ii));
	if ( (nidx<0)||(ecpmaplls_.getCount(nidx)<1) ) continue;
	check_new_index(vidx, nidx);
      }
    }
    // on met a jour la distance maximum d'association si l'objet AsM_MaxAssDistInterface am_maxdist fourni
    if (am_maxdist_p) {
      amaxrad=am_maxdist_p->getMaxAssDistance(m, z_redshift_m_->at(m), dz);
    }
    
    // On a donc tous les pixels ou il faut boucler sur les sources ds lls
    AssRecord arec(m);
    arec.getMaxxAssDistance()=amaxrad;

    double dmin=9.e9;

    for(size_t i=0; i<vidx.size(); i++) {
      vector<size_t>& vsi= ecpmaplls_.getList(vidx[i]);
      for(size_t k=0; k<vsi.size(); k++) {
	loopcnt++;
	size_t ns=vsi[k];
	double sep=llm_[m].SepAngle(lls_[ns]);
	if (sep > amaxrad) continue;
	// Selection sur redshift si applicable
	double delz=0;
	if (fgmatchredshift) {
	  delz=z_redshift_m_->at(m)-z_redshift_s_->at(ns);
	  if (fabs(delz)>dz)  continue;
	}
	arec.addAssObject(vsi[k],sep,delz);  totcnt++;
	if (sep<dmin) {
	  dmin=sep;  arec.getAssId()=vsi[k];
	  arec.getAssDistance()=dmin;
	  arec.getRedshiftDiff()=delz;
	}
      }
    }
    vassr[m]=arec;
  }

  return totcnt; 
}

  /* --- DEBUG 
     if (m%150 == 0) {
     size_t asstcnt=0;
     cout<<"====*DBG*["<<m<<"] vidx.size()="<<vidx.size()<<"  vidx:"<<endl; 
     for(size_t nn=0; nn<vidx.size(); nn++) {
     cout<<vidx[nn]<< "-> "<<maplls_.getList(vidx[nn]).size()<<" , ";
     if (nn+1%3==0) cout<<endl;
     asstcnt += maplls_.getList(vidx[nn]).size();
     }
      cout<< "\n  ==== TotAssCnt="<<asstcnt<<endl; 
      }
  */
} // Fin du namespace
  
