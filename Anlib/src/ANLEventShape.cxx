//*************************************************************************
//* =======================
//*  ANLEventShape Classes 
//* =======================
//*
//* (Description)
//*    Jet finder classes for JLC analyses.
//* (Requires)
//* 	class TLorentzVector
//* 	class Lockable
//* 	class ANL4DVector
//* (Provides)
//* 	class ANLEventShape
//* (Usage)
//* (Update Recored)
//*    1999/08/11  K.Fujii	Original version derived from the LCD 
//*				version converted from java routines 
//*				written by Gary Bower.
//*				No essential changes except for the
//*				input event format.
//*    1999/09/05  K.Ikematsu	Replaced LockableLVector with ANL4DVector.
//*    1999/09/28  K.Fujii	Plugged memory leaks as recommended by
//*				M.Iwasaki.
//*
//*************************************************************************
//
#include "ANLEventShape.h"

//_____________________________________________________________________
//  -------------------
//  ANLEventShape Class
//  -------------------
//

ClassImp(ANLEventShape)

Int_t ANLEventShape::m_maxpart = 1000;

//______________________________________________________________

ANLEventShape::ANLEventShape() 
	: m_dSphMomPower(2.0), m_dDeltaThPower(0), m_iFast(4),
	  m_dConv(0.0001), m_iGood(2),
	  m_ThrustAxis(0), m_MajorAxis(0),
	  m_MinorAxis(0), m_Thrust(0)
{
  m_dAxes.ResizeTo(4,4);
}
//______________________________________________________________

ANLEventShape::~ANLEventShape() 
{
  if (m_ThrustAxis) 	delete m_ThrustAxis;
  if (m_MajorAxis) 	delete m_MajorAxis;
  if (m_MinorAxis) 	delete m_MinorAxis;
  if (m_Thrust) 	delete m_Thrust;
}
//______________________________________________________________


void ANLEventShape::Initialize(const TObjArray& parts) 
{
  setEvent((TObjArray *)&parts);
}
//______________________________________________________________

/** 
* Call this to input a new event to the event shape routines.
* @param e An Enumeration of ANL4DVectors or their variants.
*/
// The zero element in each array in m_dAxes is ignored. Elements
// 1,2 and 3 are the x, y, and z direction cosines of the axis.
// Also the zeroth axis and thrust are ignored.
 
void ANLEventShape::setEvent(TObjArray* e)
{
  //To make this look like normal physics notation the
  //zeroth element of each array, mom[i][0], will be ignored
  //and operations will be on elements 1,2,3...
  TMatrix mom(m_maxpart,6);
  Double_t tmax = 0;
  Double_t phi = 0.;
  Double_t the = 0.;
  Double_t sgn;
  TMatrix fast(m_iFast + 1,6);
  TMatrix work(11,6);
  Double_t tdi[4] = {0.,0.,0.,0.};
  Double_t tds;
  Double_t tpr[4] = {0.,0.,0.,0.};
  Double_t thp;
  Double_t thps;
  TMatrix temp(3,5);

  // -----------------------------------
  // Store particles' 3-momenta in "mom"
  // -----------------------------------

  Int_t np = 0;
  TIter nextinev(e);
  TObject* o;
  while ((o = (TObject *)nextinev())) {
      if (np >= m_maxpart) { 
	  cerr << "Too many particles input to ANLEventShape" << endl;
	  return;
      }
      ANL4DVector *qp;
      if (o->IsA()->InheritsFrom("ANL4DVector")) {
	  qp = (ANL4DVector *)o;
      } else {
         cerr << "ANLEventShape::setEvent input is not a ANL4DVector" << endl;
	 continue;
      }
      ANL4DVector &q = *qp;
      if (q.IsLocked()) continue; 	// Skip if locked.
      mom(np,1) = q(1);
      mom(np,2) = q(2);
      mom(np,3) = q(3);
      mom(np,4) = q.GetMag();
      if ( TMath::Abs( m_dDeltaThPower ) <= 0.001 ) {
	 mom(np,5) = 1.0;
      } else {
	 mom(np,5) = TMath::Power(mom(np,4),m_dDeltaThPower);
      }
      tmax = tmax + mom(np,4)*mom(np,5);
      np++;
  }

  if ( np < 2 ) {			// Skip if only a single track there
      m_dThrust[1] = -1.0;
      m_dOblateness = -1.0;
      return;
  }
  
  // -----------------------------------
  // Find thrust/major axes
  // -----------------------------------
  // for pass = 1: find thrust axis.
  // for pass = 2: find major axis.

  for ( Int_t pass=1; pass < 3; pass++) {
      if ( pass == 2 ) {
	  phi = ulAngle(m_dAxes(1,1), m_dAxes(1,2));
	  ludbrb( &mom, 0, -phi, 0., 0., 0. );
	  for ( Int_t i = 0; i < 3; i++ ) {
	      for ( Int_t j = 1; j < 4; j++ ) temp(i,j) = m_dAxes(i+1,j);
	      temp(i,4) = 0;
	  }
	  ludbrb(&temp,0.,-phi,0.,0.,0.);
	  for ( Int_t ib = 0; ib < 3; ib++ ) {
	      for ( Int_t j = 1; j < 4; j++ ) m_dAxes(ib+1,j) = temp(ib,j);
	  }
	  the = ulAngle( m_dAxes(1,3), m_dAxes(1,1) );
	  ludbrb( &mom, -the, 0., 0., 0., 0. );
	  for ( Int_t ic = 0; ic < 3; ic++ ) {
	      for ( Int_t j = 1; j < 4; j++ ) temp(ic,j) = m_dAxes(ic+1,j);
	      temp(ic,4) = 0;
	  }
	  ludbrb(&temp,-the,0.,0.,0.,0.);
	  for ( Int_t id = 0; id < 3; id++ ) {	
	      for ( Int_t j = 1; j < 4; j++ ) m_dAxes(id+1,j) = temp(id,j);
	  }
      }

      for ( Int_t ifas = 0; ifas < m_iFast + 1 ; ifas++ ) fast(ifas,4) = 0.;

      // Find the m_iFast highest momentum particles and
      // put the highest in fast[0], next in fast[1],....fast[m_iFast-1].
      // fast[m_iFast] is just a workspace.

      for ( Int_t i = 0; i < np; i++ ) {
	  if ( pass == 2 ) {
	      mom(i,4) = TMath::Sqrt( mom(i,1)*mom(i,1) + mom(i,2)*mom(i,2) ); 
	  }
	  for ( Int_t ifas = m_iFast - 1; ifas > -1; ifas-- ) {
	      if ( mom(i,4) > fast(ifas,4) ) {
		  for ( Int_t j = 1; j < 6; j++ ) {
		      fast(ifas+1,j) = fast(ifas,j);
		      if ( ifas == 0 ) fast(ifas,j) = mom(i,j);
		  }
	      } else {
		  for ( Int_t j = 1; j < 6; j++ ) fast(ifas+1,j) = mom(i,j);
		  break;
	      }
	  }
      }

      // Find axis with highest thrust (case 1)/ highest major (case 2).

      for ( Int_t ie = 0; ie < work.GetNrows(); ie++ ) work(ie,4) = 0.;

      Int_t p = TMath::Min( m_iFast, np ) - 1;

      // Don't trust Math.pow to give right answer always.
      // Want nc = 2**p.

      Int_t nc = iPow(2,p); 
      for ( Int_t n = 0; n < nc; n++ ) {
	  for ( Int_t j = 1; j < 4; j++ ) tdi[j] = 0.;
	  for ( Int_t i = 0; i < TMath::Min(m_iFast,n); i++ ) {
	      sgn = fast(i,5);
	      if (iPow(2,(i+1))*((n+iPow(2,i))/iPow(2,(i+1))) >= i+1) {
	      	 sgn = -sgn;
	      }
	      for ( Int_t j = 1; j < 5-pass; j++ ) { 
	         tdi[j] = tdi[j] + sgn*fast(i,j);
	      }
	  }
	  tds = tdi[1]*tdi[1] + tdi[2]*tdi[2] + tdi[3]*tdi[3];
	  for ( Int_t iw = TMath::Min(n,9); iw > -1; iw--) {
	      if( tds > work(iw,4) ) {
		  for ( Int_t j = 1; j < 5; j++ ) {
		      work(iw+1,j) = work(iw,j);
		      if ( iw == 0 ) {
			  if ( j < 4 ) work(iw,j) = tdi[j];
			  else work(iw,j) = tds;
		      }
		  }
	      } else {
		  for ( Int_t j = 1; j < 4; j++ ) work(iw+1,j) = tdi[j];
		  work(iw+1,4) = tds;
	      }
	  }
      }
      
      // Iterate direction of axis until stable maximum.
      
      m_dThrust[pass] = 0;
      thp = -99999.;
      Int_t nagree = 0;
      for ( Int_t iw = 0; iw < TMath::Min(nc,10) && nagree < m_iGood; iw++ ) {
	  thp = 0.;
	  thps = -99999.;
	  while ( thp > thps + m_dConv ) {
	      thps = thp;
	      for ( Int_t j = 1; j < 4; j++ ) {
		  if ( thp <= 1E-10 ) {
		      tdi[j] = work(iw,j);
		  } else {
		      tdi[j] = tpr[j];
		      tpr[j] = 0;
		  }
	      }
	      for ( Int_t i = 0; i < np; i++ ) {
		  sgn = sign(mom(i,5), tdi[1]*mom(i,1)
			               + tdi[2]*mom(i,2)
			               + tdi[3]*mom(i,3));
		  for ( Int_t j = 1; j < 5 - pass; j++ ) {
		      tpr[j] = tpr[j] + sgn*mom(i,j);
		  }
	      }
	      thp = TMath::Sqrt(tpr[1]*tpr[1] 
				+ tpr[2]*tpr[2] 
				+ tpr[3]*tpr[3])/tmax;
	  }

	  // Save good axis. Try new initial axis until enough
	  // tries agree.

	  if ( thp < m_dThrust[pass] - m_dConv ) break;

	  if ( thp > m_dThrust[pass] + m_dConv ) {
	      nagree = 0;
	      sgn = iPow( -1, (Int_t)TMath::Nint(m_random.Rndm()) );
	      for ( Int_t j = 1; j < 4; j++ ) {
		  m_dAxes(pass,j) = sgn*tpr[j]/(tmax*thp);
	      }
	      m_dThrust[pass] = thp;
	  }
	  nagree = nagree + 1;
      }
  }
  
  // -------------------------------------------
  // Find minor axis and value by orthogonality.
  // -------------------------------------------

  sgn = iPow( -1, (Int_t)TMath::Nint(m_random.Rndm()));
  m_dAxes(3,1) = -sgn*m_dAxes(2,2);
  m_dAxes(3,2) = sgn*m_dAxes(2,1);
  m_dAxes(3,3) = 0.;
  thp = 0.;

  for ( Int_t i = 0; i < np; i++ ) {
      thp += mom(i,5)*TMath::Abs(m_dAxes(3,1)*mom(i,1) + 
				m_dAxes(3,2)*mom(i,2));
  }
  m_dThrust[3] = thp/tmax;

  // -------------------------------------------
  // Rotate back to original coordinate system.
  // -------------------------------------------

  for ( Int_t i6 = 0; i6 < 3; i6++ ) {
      for ( Int_t j = 1; j < 4; j++ ) temp(i6,j) = m_dAxes(i6+1,j);
      temp(i6,4) = 0;
  }
  ludbrb(&temp,the,phi,0.,0.,0.);
  for ( Int_t i7 = 0; i7 < 3; i7++ ) {
      for ( Int_t j = 1; j < 4; j++ ) m_dAxes(i7+1,j) = temp(i7,j);
  }
  m_dOblateness = m_dThrust[2] - m_dThrust[3];

}
//______________________________________________________________
	
// Setting and getting parameters.
	
void ANLEventShape::setThMomPower(Double_t tp)
{
  // Error if sp not positive.
  if ( tp > 0. ) m_dDeltaThPower = tp - 1.0;
  return;
}
//______________________________________________________________

Double_t ANLEventShape::getThMomPower()
{
  return 1.0 + m_dDeltaThPower;
}
//______________________________________________________________

void ANLEventShape::setFast(Int_t nf)
{
  // Error if sp not positive.
  if ( nf > 3 ) m_iFast = nf;
  return;
}
//______________________________________________________________

Int_t ANLEventShape::getFast()
{
  return m_iFast;
}
//______________________________________________________________

// Returning results

TVector3* ANLEventShape::thrustAxis()
{
  if (m_ThrustAxis) {
    m_ThrustAxis->SetXYZ(m_dAxes(1,1),m_dAxes(1,2),m_dAxes(1,3));
  } else {
    m_ThrustAxis = new TVector3(m_dAxes(1,1),m_dAxes(1,2),m_dAxes(1,3));
  }
  return m_ThrustAxis;
}
//______________________________________________________________

TVector3* ANLEventShape::majorAxis()
{
  if (m_MajorAxis) {
    m_MajorAxis->SetXYZ(m_dAxes(2,1),m_dAxes(2,2),m_dAxes(2,3));
  } else {
    m_MajorAxis = new TVector3(m_dAxes(2,1),m_dAxes(2,2),m_dAxes(2,3));
  }
  return m_MajorAxis;
}
//______________________________________________________________

TVector3* ANLEventShape::minorAxis()
{
  if (m_MinorAxis) {
    m_MinorAxis->SetXYZ(m_dAxes(3,1),m_dAxes(3,2),m_dAxes(3,3));
  } else {
    m_MinorAxis = new TVector3(m_dAxes(3,1),m_dAxes(3,2),m_dAxes(3,3));
  }
  return m_MinorAxis;
}
//______________________________________________________________

Double_t ANLEventShape::GetThrust() const
{
  return m_dThrust[1];
}
//______________________________________________________________

TVector3* ANLEventShape::thrust()
{
  if (m_Thrust) {
    m_Thrust->SetXYZ(m_dThrust[1],m_dThrust[2],m_dThrust[3]);
  } else {
    m_Thrust = new TVector3(m_dThrust[1],m_dThrust[2],m_dThrust[3]);
  }
  return m_Thrust;
}
//______________________________________________________________

Double_t ANLEventShape::oblateness()
{
  return m_dOblateness;
}
//______________________________________________________________
//
// utilities(from Jetset):

Double_t ANLEventShape::ulAngle(Double_t x, Double_t y)
{
  Double_t ulangl = 0;
  Double_t r = TMath::Sqrt(x*x + y*y);
  if ( r < 1.0E-20 ) return ulangl;
  if ( TMath::Abs(x)/r < 0.8 ) {
      ulangl = sign(TMath::ACos(x/r),y);
  } else {
      ulangl = TMath::ASin(y/r);
      if ( x < 0. && ulangl >= 0. ) {
	  ulangl = TMath::Pi() - ulangl;
      } else if ( x < 0. ) {
	  ulangl = -TMath::Pi() - ulangl;
      }
  }
  return ulangl;
}
//______________________________________________________________

Double_t ANLEventShape::sign(Double_t a, Double_t b)
{
  if ( b < 0 ) 	return -TMath::Abs(a);
  else 		return TMath::Abs(a);
}
//______________________________________________________________

void ANLEventShape::ludbrb(TMatrix* mom, 
			Double_t the, 
			Double_t phi, 
			Double_t bx, 
			Double_t by,
			Double_t bz)
{
  // Ignore "zeroth" elements in rot,pr,dp.
  // Trying to use physics-like notation.

  TMatrix rot(4,4);
  Double_t pr[4];
  Double_t dp[5];
  Int_t np = mom->GetNrows();

  if ( the*the + phi*phi > 1.0E-20 ) {
      rot(1,1) = TMath::Cos(the)*TMath::Cos(phi);
      rot(1,2) = -TMath::Sin(phi);
      rot(1,3) = TMath::Sin(the)*TMath::Cos(phi);
      rot(2,1) = TMath::Cos(the)*TMath::Sin(phi);
      rot(2,2) = TMath::Cos(phi);
      rot(2,3) = TMath::Sin(the)*TMath::Sin(phi);
      rot(3,1) = -TMath::Sin(the);
      rot(3,2) = 0.0;
      rot(3,3) = TMath::Cos(the);

      for ( Int_t i = 0; i < np; i++ ) {
	  for ( Int_t j = 1; j < 4; j++ ) {
	      pr[j] = (*mom)(i,j);
	      (*mom)(i,j) = 0;
	  }
	  for ( Int_t jb = 1; jb < 4; jb++) {
	      for ( Int_t k = 1; k < 4; k++) {
		  (*mom)(i,jb) = (*mom)(i,jb) + rot(jb,k)*pr[k];
	      }
	  }
      }
      Double_t beta = TMath::Sqrt( bx*bx + by*by + bz*bz );
      if ( beta*beta > 1.0E-20 ) {
	  if ( beta >  0.99999999 ) {
			 //send message: boost too large, resetting to <~1.0.
	      bx = bx*(0.99999999/beta);
	      by = by*(0.99999999/beta);
	      bz = bz*(0.99999999/beta);
	      beta =   0.99999999;
	  }
	  Double_t gamma = 1.0/TMath::Sqrt(1.0 - beta*beta);
	  for ( Int_t i = 0; i < np; i++ ) {
	      for ( Int_t j = 1; j < 5; j++ ) dp[j] = (*mom)(i,j);
	      Double_t bp = bx*dp[1] + by*dp[2] + bz*dp[3];
	      Double_t gbp = gamma*(gamma*bp/(1.0 + gamma) + dp[4]);
	      (*mom)(i,1) = dp[1] + gbp*bx;
	      (*mom)(i,2) = dp[2] + gbp*by;
	      (*mom)(i,3) = dp[3] + gbp*bz;
	      (*mom)(i,4) = gamma*(dp[4] + bp);
	  }
      }
  }
  return;
}
//______________________________________________________________

Int_t ANLEventShape::iPow(Int_t man, Int_t exp)
{
  Int_t ans = 1;
  for( Int_t k = 0; k < exp; k++) ans = ans*man;
  return ans;
}


