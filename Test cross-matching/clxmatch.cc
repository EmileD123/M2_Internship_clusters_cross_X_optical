/* ---------------------------------------------------------------   
   SOPHYA project - (C) Univ. Paris-Saclay, IJCLab CNRS/IN2P3, CEA-Irfu
   Programme d'association d'objets entre deux catalogues d'amas  
   
   R. Ansari ,  Juin  2024  (stage Jeanine Smoyan) 
                Avril 2025  (Stage E. Dosso) 


   (C) UPS+LAL IN2P3/CNRS     (C) IRFU/CEA 
   (C) Univ. Paris-Saclay     (C) IJCLab CNRS/IN2P3 
   --------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>

//--- include SOPHYA 
#include "array.h"
#include "histats.h"
#include "fitsioserver.h"
#include "skymap.h"
#include "randfmt.h"

#include "slininterp.h"
#include "astromatch.h"
#include "luc.h"
#include "interpcosmoc.h"


using namespace std; 
using namespace SOPHYA; 


//--- Declaration des fonctions de ce fichier

using namespace std; 
using namespace SOPHYA; 

//---------- Parametres pour l'association de catalogues d'objets 
// resolution de la carte du ciel en degre, pour repartir les sources 
static double skyresol_deg = 5.;   // resolution des listes d'objets repartis sur le ciel (pour accelerer la recherche et l'association 
static double asslimite_arcmin = 15.;   // Valeur fixe de la limite angulaire d'association 
static double assd0_arcmin = 5.;        // Valeur plancher pour une limite d'association dependant du redshift 
static double ray_clus = 1;      // Parametre d'association dependant du redshift - distance transverse in Mpc
static double ass_dz0 = 0.15;    // Parametre d'association dependant du redshift - Ecart maximum en z (redshift) 
//--------------------------------------------------------------

//--- Classe permettant de 
class CluxMaxAssDist : public AsM_MaxAssDistInterface  {
public:
  // clus_rad : cluster radius in Mpc
  // d0ass : min association radius in arcmin
  CluxMaxAssDist(SimpleUniverse & su, double clus_rad=0.35, Angle d0ass=Angle(1,Angle::ArcMin), double dz0=0.1)
    : icosmo_(su), clus_rad_(clus_rad), d0ass_rad_(d0ass.ToRadian()), dz0_(dz0)
  {
    maxdass_rad_=Angle(skyresol_deg, Angle::Degree).ToRadian()*0.5;
  } 
  virtual double getMaxAssDistance(size_t numom, double z, double & dz) const 
  {
    if (z<0.)  return d0ass_rad_;
    double ascl=clus_rad_/icosmo_.AngularDiameterDistanceMpc(z);
    dz=dz0_;
    double dass=ascl+d0ass_rad_;
    if (dass>maxdass_rad_)  return maxdass_rad_;
    return dass; 
  }
  virtual bool reqRedshiftMatch() const
  {
    return true;
  }
  InterpCosmoCalc icosmo_;
  double clus_rad_;
  double d0ass_rad_, maxdass_rad_;  // in radian
  double dz0_;   

};

DataTable create_assinfo_dt(DataTable & catm, DataTable & cats, std::vector<AssRecord>& asrecs);
DataTable create_assinfo_all_dt(DataTable & catm, DataTable & cats, std::vector<AssRecord>& asrecs);
int update_mcat_assinfo(DataTable & catm, DataTable & cats, std::vector<AssRecord>& asrecs);

//---------  Main program -------------
//--------------------------------------------------------------

int main (int narg, char* arg[])
{
  int rc = 0;
  if (narg<4) {
    cout << " usage: clxmatch Catalog1-Fits-Name Catalog2-Fits-Name OutFileName" <<endl;
    cout << " OutFileName : PPF or FITS  (.ppf , .fits) " <<endl;
    return 1;
  }
  string cat1name=arg[1];
  string cat2name=arg[2];
  string outfname=arg[3];
  bool fgfits=false;
  if (outfname.length()>=5) {
    if (outfname.substr(outfname.length()-5)==".fits") fgfits=true;
  }
  try {
    SophyaInit();
    Timer tm("clxmatch.cc");

    DataTable cat1, cat2; 
    cout << "[1] Reading input catalog 1 (main) from input FITS file"<< cat1name<<endl;
    {
      FitsInOutFile fis(cat1name, FitsInOutFile::Fits_RO);
      fis.MoveAbsToHDU(2);
      fis >> cat1;
    }
    cat1.SetShowMinMaxFlag(true);
    cout << cat1;

    vector<double> vram;
    vector<double> vdecm;
    vector<double> vredzm;
    cat1.GetColumn(cat1.IndexNom("RAJ2000"),vram);
    cat1.GetColumn(cat1.IndexNom("DEJ2000"),vdecm);
    {
      vector<float> vredzf;
      cat1.GetColumn(cat1.IndexNom("zBest"),vredzf);
      vredzm.resize(vredzf.size());
      for(size_t ii=0; ii<vredzf.size(); ii++) vredzm[ii]=vredzf[ii];
    }
    vector<LongitudeLatitude> vllm(vram.size());
    for(size_t i=0; i<vram.size(); i++)  {
      vllm[i]=LongitudeLatitude(Angle(vram[i],Angle::Degree), Angle(vdecm[i],Angle::Degree), true);
    }
    
    cout << "[2] Reading input catalog 2  from input FITS file"<< cat2name<<endl;
    {
      FitsInOutFile fis(cat2name, FitsInOutFile::Fits_RO);
      fis.MoveAbsToHDU(2);
      fis >> cat2;
    }
    cat2.SetShowMinMaxFlag(true);
    cout << cat2;

    vector<double> vra;
    vector<double> vdec;
    vector<double> vredz;
    cat2.GetColumn(cat2.IndexNom("RAJ2000"),vra);
    cat2.GetColumn(cat2.IndexNom("DEJ2000"),vdec);
    {
      vector<float> vredzf;
      cat2.GetColumn(cat2.IndexNom("zph"),vredzf);
      vredz.resize(vredzf.size());
      for(size_t ii=0; ii<vredzf.size(); ii++) vredz[ii]=vredzf[ii];
    }
    vector<LongitudeLatitude> vll(vra.size());
    
    for(size_t i=0; i<vra.size(); i++)  {
      LongitudeLatitude ll(Angle(vra[i],Angle::Degree), Angle(vdec[i],Angle::Degree), true);
      vll[i]=ll;
    }
    

    tm.Split("After-read");

    //----- On fait le croisement des catalogues ici 
    Angle skyresol(skyresol_deg,Angle::Degree);
    //----- Classe pour faire l'association 
    AstroMatch asma(vllm, vll, skyresol, &vredzm, &vredz);   

    //---- Tolerance d'association dependant du redshift 
    SimpleUniverse su;
    CluxMaxAssDist cmaxassdist(su, ray_clus, Angle(assd0_arcmin, Angle::ArcMin), ass_dz0);
    std::vector<AssRecord> asrecs;

    tm.Split("Before-match");    
    asma.doMatch(cmaxassdist, asrecs);
    cout << "[3]  Source match done " << endl;
    tm.Split("After-match");

    cout << "[4.a]  Creating association DataTable " << endl;
    DataTable dtass=create_assinfo_dt(cat1, cat2, asrecs);
    DataTable dtassall=create_assinfo_all_dt(cat1, cat2, asrecs);

    cout << "[4.b]  Updating first (main) object catalog with association information " << endl;
    update_mcat_assinfo(cat1, cat2, asrecs);
			
    SphereThetaPhi<int_4> cntmap=asma.getCountMapS();
    cout << "[4.c]  Got cat2 count map " << endl;

    dtass.SetShowMinMaxFlag(true);
    cout<<dtass;
    dtassall.SetShowMinMaxFlag(true);
    cout<<dtassall;
    cat1.SetShowMinMaxFlag(true);
    cout<<cat1;


    if (fgfits) {
      cout << "[5]  saving cntmap and catass, dtass, dtassall DataTable to output FITS file: "<<outfname<< endl;
      FitsInOutFile fos(outfname, FitsInOutFile::Fits_Create);
      fos << FitsNameTag("catass") << cat1;
      fos << FitsNameTag("dtass") << dtass; 
      fos << FitsNameTag("dtassall") << dtassall;
    }
    else {
    cout << "[5]  saving cntmap and catass, dtass, dtassall DataTable to output PPF file: "<<outfname<< endl;
      POutPersist pof(outfname);
      pof << PPFNameTag("cntmap") << cntmap;
      pof << PPFNameTag("catass") << cat1;
      pof << PPFNameTag("dtass") << dtass; 
      pof << PPFNameTag("dtassall") << dtassall;
    }
  }
     // End of try bloc 
  catch (PThrowable & exc) {
    cerr << " clxmatch.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << " - Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {
    cerr << " clxmatch.cc: Catched std::exception "  
         << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {
    cerr << " clxmatch.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ------------ End of clxmatch.cc --------------- Rc=" << rc << endl;
  return rc;
};

/* --Fonction-- */
DataTable create_assinfo_dt(DataTable & catm, DataTable & cats, std::vector<AssRecord>& asrecs)
{
  DataTable dtm(256);

  dtm.AddIntegerColumn("idm");
  dtm.AddIntegerColumn("asscnt");
  dtm.AddIntegerColumn("idass");
  dtm.AddFloatColumn("dstass");
  dtm.AddFloatColumn("delz");
  dtm.AddFloatColumn("maxassdist");

  DataTableRow row=dtm.EmptyRow();
  for(size_t nn=0; nn<asrecs.size(); nn++) {
    row("idm")=(int_8)nn;
    row("asscnt")=(int_8)asrecs[nn].getAssCount();
    row("idass")=(int_8)asrecs[nn].getAssId();
    row("dstass")=Angle(asrecs[nn].getAssDistance()).ToArcMin();
    row("delz")=asrecs[nn].getRedshiftDiff();
    row("maxassdist")=Angle(asrecs[nn].getMaxxAssDistance()).ToArcMin();
    dtm.AddRow(row);
  }
  return dtm;
}

/* --Fonction-- */
DataTable create_assinfo_all_dt(DataTable & catm, DataTable & cats, std::vector<AssRecord>& asrecs)
{
  DataTable dtm(256);

  dtm.AddIntegerColumn("idm");
  dtm.AddIntegerColumn("asscnt");
  dtm.AddIntegerColumn("numass");
  dtm.AddIntegerColumn("idass");
  dtm.AddFloatColumn("dstass");
  dtm.AddFloatColumn("delz");
  dtm.AddFloatColumn("maxassdist");
  dtm.AddIntegerColumn("isbest");

  DataTableRow row=dtm.EmptyRow();
  for(size_t nn=0; nn<asrecs.size(); nn++) {
    row("idm")=(uint_8)nn;
    row("asscnt")=(uint_8)asrecs[nn].getAssCount();
    for(size_t j=0; j<asrecs[nn].getAssCount(); j++) {
      row("numass")=(uint_8)j; 
      row("idass")=(uint_8)asrecs[nn].getAssObjectId(j);
      row("dstass")=Angle(asrecs[nn].getAssObjectDistance(j)).ToArcMin();
      row("delz")=asrecs[nn].getAssObjectRedshiftDiff(j);
      row("maxassdist")=Angle(asrecs[nn].getMaxxAssDistance()).ToArcMin();
      if (asrecs[nn].getAssObjectId(j)==asrecs[nn].getAssId())  row("isbest")=(uint_8)j;
      else row("isbest")=-1;
      dtm.AddRow(row);
    }
  }
  return dtm;
}

/* --Fonction-- */
int update_mcat_assinfo(DataTable & catm, DataTable & cats, std::vector<AssRecord>& asrecs)
{

  catm.AddIntegerColumn("idm",1,true);
  catm.AddIntegerColumn("asscnt",1,true);
  catm.AddIntegerColumn("idass",1,true);
  catm.AddFloatColumn("dstass",1,true);
  catm.AddFloatColumn("delz",1,true);
  catm.AddFloatColumn("maxassdist",1,true);
  catm.AddFloatColumn("zph",1,true);
  catm.AddFloatColumn("rmag",1,true);
  catm.AddFloatColumn("r200",1,true);
  catm.AddIntegerColumn("N200",1,true);

  if (asrecs.size() != catm.NRows())
    throw SzMismatchError("update_mcat_assinfo()/ERROR asrecs.size() != catm.NRows() ");
  vector<int_4> vidm(asrecs.size()), vacnt(asrecs.size()), vidass(asrecs.size());
  vector<int_4> vidass2(asrecs.size());
  vector<r_4> vdstass(asrecs.size()), vmaxdstass(asrecs.size());
  vector<r_4> vdelz(asrecs.size()), vzph(asrecs.size()), vrmag(asrecs.size());
  vector<r_4> vr200(asrecs.size());
  vector<int_4> vN200(asrecs.size());

  DataTableRow rows=cats.EmptyRow();

  for(size_t nn=0; nn<asrecs.size(); nn++) {
    vidm[nn]=nn;  vacnt[nn]=asrecs[nn].getAssCount();
    vmaxdstass[nn]=Angle(asrecs[nn].getMaxxAssDistance()).ToArcMin();
    if (asrecs[nn].okAss()) {
      vidass[nn]=asrecs[nn].getAssId();
      vdstass[nn]=Angle(asrecs[nn].getAssDistance()).ToArcMin();
      vdelz[nn]=asrecs[nn].getRedshiftDiff();
      cats.GetRow(asrecs[nn].getAssId(), rows);
      vzph[nn]=rows("zph");
      vrmag[nn]=rows("rmag");
      vr200[nn]=rows("r200");
      vN200[nn]=rows("N200");
    }
    else {
      vidass[nn]=-1; vdstass[nn]=-1.; vzph[nn]=vrmag[nn]=vr200[nn]=vN200[nn]=-1.;
    }
  }
  catm.FillColumn(catm.IndexNom("idm"), vidm);
  catm.FillColumn(catm.IndexNom("asscnt"), vacnt);
  catm.FillColumn(catm.IndexNom("idass"), vidass);
  catm.FillColumn(catm.IndexNom("dstass"), vdstass);
  catm.FillColumn(catm.IndexNom("delz"), vdelz);
  catm.FillColumn(catm.IndexNom("maxassdist"), vmaxdstass);
  catm.FillColumn(catm.IndexNom("zph"), vzph);
  catm.FillColumn(catm.IndexNom("rmag"), vrmag);
  catm.FillColumn(catm.IndexNom("r200"), vr200);
  catm.FillColumn(catm.IndexNom("N200"), vN200);

  return 0;
}

