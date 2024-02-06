#include <iostream>
//#include <hash_set.h>
//#include <hash_set>
#include <set>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>


#include "hull3D.h"

using namespace std;


/* copyright 2016 Dr David Sinclair
   david@newtonapples.net 
   all rights reserved.
 

   version of 10-feb-2016.


   this code is released under GPL3,
   a copy of the license can be found at 
   http://www.gnu.org/licenses/gpl-3.0.html

   you can purchase a un-restricted license from
   http://newtonapples.net

   where algorithm details are explained.

   If you do choose to purchase a license I will do my best 
   to post you a free bag of Newton Apple Chocolates!, 
   the cleverest chocolates on the Internet.



 */



/*  
    read an ascii file of (r,c) or (r,c,z) point pairs.

    the first line of the points file should contain
    "<NUMP>  2 points" or maybe     "<NUMP>  3 points"  

    if it does not have the word points in it the first line is 
    interpretted as a point pair or triplod.

 */

int read_R3(std::vector<R3> &pts, char * fname){
  char s0[513];
  int nump =0;
  float p1,p2, p3;

  R3 pt;

  std::string line;
  std::string points_str("points");

  std::ifstream myfile;
  myfile.open(fname);

  if (myfile.is_open()){
    
    getline (myfile,line);
    //int numc = line.length();
    
    // check string for the string "points"
    int n = (int) line.find( points_str);
    if( n < line.length() ){ 
      while ( myfile.good() ){
	getline (myfile,line);
	if( line.length() <= 512){
	  copy( line.begin(), line.end(), s0);
	  s0[line.length()] = 0;
	  int v = sscanf( s0, "%g %g %g", &p1,&p2, &p3);
	  if( v == 3 ){
	    pt.id = nump; 
	    nump++;
	    pt.r = p1;
	    pt.c = p2;
	    pt.z = p3;
	    
	    pts.push_back(pt);
	  }
	  
	  else{
	    v = sscanf( s0, "%g %g", &p1,&p2);
	    if( v == 2 ){
	      pt.id = nump; 
	      nump++;
	      pt.r = p1;
	      pt.c = p2;
	      pt.z = p1*p1 + p2*p2;
	      
	      pts.push_back(pt);
	    }
	  }
	}   
      }
    }
    else{   // assume all number pairs on a line are points
      if( line.length() <= 512){
	copy( line.begin(), line.end(), s0);
	s0[line.length()] = 0;
	int v = sscanf( s0, "%g %g %g", &p1,&p2, &p3);
	if( v == 3 ){
	  pt.id = nump; 
	  nump++;
	  pt.r = p1;
	  pt.c = p2;
	  pt.z = p3;
	  
	  pts.push_back(pt);
	}

	else{
	  v = sscanf( s0, "%g %g", &p1,&p2);
	  if( v == 2 ){
	    pt.id = nump; 
	    nump++;
	    pt.r = p1;
	    pt.c = p2;
	    pt.z = p1*p1 + p2*p2;

	    pts.push_back(pt);
	  }
	}


      }   

      while ( myfile.good() ){
	getline (myfile,line);
	if( line.length() <= 512){
	  copy( line.begin(), line.end(), s0);
	  s0[line.length()] = 0;
	  int v = sscanf( s0, "%g %g %g", &p1,&p2, &p3);
	  if( v == 3 ){
	    pt.id = nump; 
	    nump++;
	    pt.r = p1;
	    pt.c = p2;
	    pt.z = p3;
	    
	    pts.push_back(pt);
	  }
	  
	  else{
	    v = sscanf( s0, "%g %g", &p1,&p2);
	    if( v == 2 ){
	      pt.id = nump; 
	      nump++;
	      pt.r = p1;
	      pt.c = p2;
	      pt.z = p1*p1 + p2*p2;
	      
	      pts.push_back(pt);
	    }
	  }
	}
      }
    }
    myfile.close();
  }

  nump = (int) pts.size();

  return(nump);
};

/*
	write out a set of points to disk


*/

void write_R3(std::vector<R3> &pts, char * fname){
   std::ofstream out(fname, ios::out);
   
   int nr = (int) pts.size();
   out << nr << " 3 points" << endl;
   
   for (int r = 0; r < nr; r++){
     out << pts[r].r << ' ' << pts[r].c <<  ' ' << pts[r].z << endl;
   }
   out.close();
   
   return;
};



/*
 write out triangle ids to be compatible with matlab/octave array numbering.

 */
void write_Tris(std::vector<Tri> &ts, char * fname){
   std::ofstream out(fname, ios::out);
   
   int nr = (int) ts.size();
   out << nr << " 6   point-ids (1,2,3)  adjacent triangle-ids ( limbs ab  ac  bc )" << endl;
   
   for (int r = 0; r < nr; r++){
     out << ts[r].a+1 << ' ' << ts[r].b+1 <<' ' << ts[r].c+1 <<' ' 
	 << ts[r].ab+1 <<' ' << ts[r].ac+1 <<' ' << ts[r].bc+1 << endl; //" " << ts[r].ro <<  endl;
   }
   out.close();
   
   return;
};



/*  To create a Delaunay triangulation from a 2D point set:
    make the 'z' coordinate = (x*x + y*y),
   find the convex hull in R3 of a point cloud.
   using a sweep-hull algorithm called the Newton Apple Wrapper.
   discard the facets that are not downward facing.

 */

int NewtonApple_Delaunay( std::vector<R3> &pts, std::vector<Tri> &hulk)
{

  int nump = (int) pts.size();


  if( nump <= 4 ){
    cerr << "less than 4 points, aborting " << endl;
    return(-1);
  }


  sort( pts.begin(), pts.end() );

  std::vector<Tri> hull;

  int num = init_hull3D(pts, hull);

  //   return(0); // exit here is you do not need to write the triangles to disk.

  // just pick out the hull triangles and renumber.
  int numh = hull.size();
  int *taken = new int [numh];
  //int taken[numh];
  memset(taken, -1, numh*sizeof(int));

  int cnt = 0;
  for(int t=0; t<numh; t++){  // create an index from old tri-id to new tri-id.
    if( hull[t].keep > 0 ){   // point index remains unchanged.
      taken[t] = cnt;
      cnt ++;
    }
  }

  for(int t=0; t<numh; t++){  // create an index from old tri-id to new tri-id.
    if( hull[t].keep > 0 ){   // point index remains unchanged.
      Tri T = hull[t];
      T.id = taken[t];
      if( taken[T.ab] < 0 ){
	cerr << "broken hull" << endl;
	delete [] taken;
	exit(0);
      }
      T.ab = taken[T.ab];


      if( taken[T.bc] < 0 ){
	cerr << "broken hull" << endl;
	delete [] taken;
	exit(0);
      }
      T.bc = taken[T.bc];

      if( taken[T.ac] < 0 ){
	cerr << "broken hull" << endl;
	delete [] taken;
	exit(0);
      }
      T.ac = taken[T.ac];


      // look at the normal to the triangle
      if( hull[t].ez < 0 ){
	hulk.push_back(T);
      }
    }
  }

  delete [] taken;

  return(1);
}




/*  
   find the convex hull in R3 of a point cloud.
   using a sweep-hull algorithm called the Newton Apple Wrapper.

 */

int NewtonApple_hull_3D( std::vector<R3> &pts, std::vector<Tri> &hulk)
{

  int nump = (int) pts.size();


  if( nump <= 4 ){
    cerr << "less than 4 points, aborting " << endl;
    return(-1);
  }


  sort( pts.begin(), pts.end() );

  std::vector<Tri> hull;

  int num = init_hull3D(pts, hull);

  //   return(0); // exit here is you do not need to write the triangles to disk.

  // just pick out the hull triangles and renumber.
  int numh = hull.size();
  int *taken = new int [numh];
  //int taken[numh];
  memset(taken, -1, numh*sizeof(int));

  int cnt = 0;
  for(int t=0; t<numh; t++){  // create an index from old tri-id to new tri-id.
    if( hull[t].keep > 0 ){   // point index remains unchanged.
      taken[t] = cnt;
      cnt ++;
    }
  }

  for(int t=0; t<numh; t++){  // create an index from old tri-id to new tri-id.
    if( hull[t].keep > 0 ){   // point index remains unchanged.
      Tri T = hull[t];
      T.id = taken[t];
      if( taken[T.ab] < 0 ){
	cerr << "broken hull" << endl;
	delete [] taken;
	exit(0);
      }
      T.ab = taken[T.ab];


      if( taken[T.bc] < 0 ){
	cerr << "broken hull" << endl;
	delete [] taken;
	exit(0);
      }
      T.bc = taken[T.bc];

      if( taken[T.ac] < 0 ){
	cerr << "broken hull" << endl;
	delete [] taken;
	exit(0);
      }
      T.ac = taken[T.ac];

      hulk.push_back(T);
    }
  }

  delete [] taken;

  return(1);
}

// if you are de-duplicating points in R2 remember to set .z to something sensible.


int de_duplicateR3( std::vector<R3> &pts, std::vector<int> &outx,std::vector<R3> &pts2 ){

  int nump = (int) pts.size();
  std::vector<R3> dpx;
  R3 d;
  for( int k=0; k<nump; k++){
    d.r = pts[k].r;
    d.c = pts[k].c;
    d.z = pts[k].z;
    d.id = k;
    dpx.push_back(d);
  }

  sort(dpx.begin(), dpx.end());
  
  cerr << "de-duplicating ";
  pts2.clear();
  pts2.push_back(pts[dpx[0].id]);
  pts2[0].id = 0;
  int cnt = 1;
  
  for( int k=0; k<nump-1; k++){
    if( dpx[k].r == dpx[k+1].r && dpx[k].c == dpx[k+1].c && dpx[k].z == dpx[k+1].z ){
      outx.push_back( dpx[k+1].id);
    }
    else{
      pts[dpx[k+1].id].id = cnt;
      pts2.push_back(pts[dpx[k+1].id]);
      cnt++;
    }
  }

  cerr << "removed  " << outx.size() << " points " << endl;

  return(outx.size());
}


// initialise the hull to the point where there is a not zero volume 
// hull.


int init_hull3D( std::vector<R3> &pts, std::vector<Tri> &hull)
{

  int nump = (int) pts.size();

  float mr=0, mc = 0, mz = 0;
  float Mr=0, Mc = 0, Mz = 0;

  Tri T1(0,1,2);
  float r0 = pts[0].r, c0 = pts[0].c, z0 = pts[0].z;
  float r1 = pts[1].r, c1 = pts[1].c, z1 = pts[1].z;
  float r2 = pts[2].r, c2 = pts[2].c, z2 = pts[2].z;
  
  Mr = r0+r1+r2;
  Mc = c0+c1+c2;
  Mz = z0+z1+z2;

  // check for colinearity 
  float r01 = r1-r0, r02 = r2-r0;
  float c01 = c1-c0, c02 = c2-c0;  
  float z01 = z1-z0, z02 = z2-z0;

  float e0 = c01*z02 - c02*z01;
  float e1 = -r01*z02 + r02*z01;
  float e2 = r01*c02 - r02*c01;

  if( e0==0 && e1==0 && e2==0 ){ // do not add a facet.
    cerr << "stop fucking me arround and give me a valid opening facet, you tit. " << endl;
    return(-1);
  }

  T1.id = 0;
  T1.er = e0;
  T1.ec = e1;
  T1.ez = e2;

  T1.ab = 1;   // adjacent facet id number
  T1.ac = 1;
  T1.bc = 1;

  hull.push_back( T1);

  T1.id = 1;
  T1.er = -e0;
  T1.ec = -e1;
  T1.ez = -e2;
  
  T1.ab = 0;
  T1.ac = 0;
  T1.bc = 0;

  hull.push_back(T1);
  std::vector<int> xlist;
  Tri Tnew;

  int numt = 2;

  for( int p=3; p<nump; p++){ // add points until a non coplanar set of points is achieved.
    R3 &pt = pts[p];

    Mr += pt.r; mr = Mr/(p+1);
    Mc += pt.c; mc = Mc/(p+1);
    Mz += pt.z; mz = Mz/(p+1);

    // find the first visible plane.
    int numh = hull.size();
    int hvis = -1;
    float r = pt.r;
    float c = pt.c;
    float z = pt.z;
    xlist.clear();

    for( int h=numh-1; h>=0; h--){
      Tri &t= hull[h];
      float R1 = pts[t.a].r;
      float C1 = pts[t.a].c;
      float Z1 = pts[t.a].z;

      float dr = r-R1;
      float dc = c-C1;
      float dz = z-Z1;

      float d = dr*t.er + dc*t.ec + dz*t.ez;

      if( d > 0 ){
	hvis = h;
	hull[h].keep = 0;
	xlist.push_back(hvis);
	break;
      }
    }
    if( hvis < 0 ){
      add_coplanar(pts, hull, p);
    }
    if( hvis >= 0 ){
      // new triangular facets are formed from neighbouring invisible planes.
      int numh = hull.size();
      int numx = xlist.size();
      for( int x=0; x<numx; x++){
	int xid = xlist[x];
	int ab = hull[xid].ab;     // facet adjacent to line ab
	Tri &tAB= hull[ab];


	float R1 = pts[tAB.a].r;  // point on next triangle
	float C1 = pts[tAB.a].c;
	float Z1 = pts[tAB.a].z;

	float dr = r-R1;
	float dc = c-C1;
	float dz = z-Z1;

	float d = dr*tAB.er + dc*tAB.ec + dz*tAB.ez;

	if( d > 0 ){ // add to xlist.
	  if( hull[ab].keep == 1){
	    hull[ab].keep = 0;
	    xlist.push_back(ab);
	    numx ++;
	  }
	}
	else{ // spawn a new triangle.
	  Tnew.id = hull.size();
	  Tnew.keep = 2;
	  Tnew.a = p;
	  Tnew.b = hull[xid].a;  
	  Tnew.c = hull[xid].b;

	  Tnew.ab = -1;
	  Tnew.ac = -1;
	  Tnew.bc = ab;

	  // make normal vector.
	  float dr1 = pts[Tnew.a].r- pts[Tnew.b].r, dr2 = pts[Tnew.a].r- pts[Tnew.c].r;
	  float dc1 = pts[Tnew.a].c- pts[Tnew.b].c, dc2 = pts[Tnew.a].c- pts[Tnew.c].c;
	  float dz1 = pts[Tnew.a].z- pts[Tnew.b].z, dz2 = pts[Tnew.a].z- pts[Tnew.c].z;

	  float er = (dc1*dz2-dc2*dz1);
	  float ec = -(dr1*dz2-dr2*dz1);
	  float ez = (dr1*dc2-dr2*dc1);

	  dr = mr-r; // points from new facet towards [mr,mc,mz]
	  dc = mc-c;
	  dz = mz-z;
	  // make it point outwards.

	  float dromadery = dr*er +  dc*ec + dz*ez;

	  if( dromadery > 0 ){
	    Tnew.er = -er;
	    Tnew.ec = -ec;
	    Tnew.ez = -ez;
	  }
	  else{
	    Tnew.er = er;
	    Tnew.ec = ec;
	    Tnew.ez = ez;
	  }	    

	  // update the touching triangle tAB
	  int A = hull[xid].a, B = hull[xid].b;
	  if( (tAB.a == A && tAB.b == B ) ||(tAB.a == B && tAB.b == A ) ){
	    tAB.ab = hull.size();
	  }
	  else if( (tAB.a == A && tAB.c == B ) ||(tAB.a == B && tAB.c == A ) ){
	    tAB.ac = hull.size();
	  }
	  else if( (tAB.b == A && tAB.c == B ) ||(tAB.b == B && tAB.c == A ) ){
	    tAB.bc = hull.size();
	  }
	  else{
	    cerr << "Oh crap, the di-lithium crystals are fucked!" << endl;
	    cerr << "numeric stability issue, exiting now." << endl;
	    exit(9);
	  }


	  hull.push_back(Tnew);

	}
      

	// second side of the struck out triangle

	int ac = hull[xid].ac;     // facet adjacent to line ac
	Tri &tAC = hull[ac];

	R1 = pts[tAC.a].r;  // point on next triangle
	C1 = pts[tAC.a].c;
	Z1 = pts[tAC.a].z;

	dr = r-R1;
	dc = c-C1;
	dz = z-Z1;

	d = dr*tAC.er + dc*tAC.ec + dz*tAC.ez;

	if( d > 0 ){ // add to xlist.
	  if( hull[ac].keep == 1){
	    hull[ac].keep = 0;
	    xlist.push_back(ac);
	    numx ++;
	  }
	}
	else{ // spawn a new triangle.
	  Tnew.id = hull.size();
	  Tnew.keep = 2;
	  Tnew.a = p;
	  Tnew.b = hull[xid].a;  
	  Tnew.c = hull[xid].c;

	  Tnew.ab = -1;
	  Tnew.ac = -1;
	  Tnew.bc = ac;

	  // make normal vector.
	  float dr1 = pts[Tnew.a].r- pts[Tnew.b].r, dr2 = pts[Tnew.a].r- pts[Tnew.c].r;
	  float dc1 = pts[Tnew.a].c- pts[Tnew.b].c, dc2 = pts[Tnew.a].c- pts[Tnew.c].c;
	  float dz1 = pts[Tnew.a].z- pts[Tnew.b].z, dz2 = pts[Tnew.a].z- pts[Tnew.c].z;

	  float er = (dc1*dz2-dc2*dz1);
	  float ec = -(dr1*dz2-dr2*dz1);
	  float ez = (dr1*dc2-dr2*dc1);

	  dr = mr-r; // points from new facet towards [mr,mc,mz]
	  dc = mc-c;
	  dz = mz-z;
	  // make it point outwards.

	  float dromadery = dr*er +  dc*ec + dz*ez;

	  if( dromadery > 0 ){
	    Tnew.er = -er;
	    Tnew.ec = -ec;
	    Tnew.ez = -ez;
	  }
	  else{
	    Tnew.er = er;
	    Tnew.ec = ec;
	    Tnew.ez = ez;
	  }	    
	  // update the touching triangle tAC
	  int A = hull[xid].a, C = hull[xid].c;
	  if( (tAC.a == A && tAC.b == C ) ||(tAC.a == C && tAC.b == A ) ){
	    tAC.ab = hull.size();
	  }
	  else if( (tAC.a == A && tAC.c == C ) ||(tAC.a == C && tAC.c == A ) ){
	    tAC.ac = hull.size();
	  }
	  else if( (tAC.b == A && tAC.c == C ) ||(tAC.b == C && tAC.c == A ) ){
	    tAC.bc = hull.size();
	  }
	  else{
	    cerr << "Oh crap, warp drive failure, dude!" << endl;
	    cerr << "numeric stability issue, exiting now." << endl;
	    exit(9);
	  }


	  hull.push_back(Tnew);
	}
      

	// third side of the struck out triangle


	int bc = hull[xid].bc;     // facet adjacent to line ac
	Tri &tBC = hull[bc];


	R1 = pts[tBC.a].r;  // point on next triangle
	C1 = pts[tBC.a].c;
	Z1 = pts[tBC.a].z;

	dr = r-R1;
	dc = c-C1;
	dz = z-Z1;

	d = dr*tBC.er + dc*tBC.ec + dz*tBC.ez;

	if( d > 0 ){ // add to xlist.
	  if( hull[bc].keep == 1){
	    hull[bc].keep = 0;
	    xlist.push_back(bc);
	    numx ++;
	  }
	}
	else{ // spawn a new triangle.
	  Tnew.id = hull.size();
	  Tnew.keep = 2;
	  Tnew.a = p;
	  Tnew.b = hull[xid].b;  
	  Tnew.c = hull[xid].c;

	  Tnew.ab = -1;
	  Tnew.ac = -1;
	  Tnew.bc = bc;

	  // make normal vector.
	  float dr1 = pts[Tnew.a].r- pts[Tnew.b].r, dr2 = pts[Tnew.a].r- pts[Tnew.c].r;
	  float dc1 = pts[Tnew.a].c- pts[Tnew.b].c, dc2 = pts[Tnew.a].c- pts[Tnew.c].c;
	  float dz1 = pts[Tnew.a].z- pts[Tnew.b].z, dz2 = pts[Tnew.a].z- pts[Tnew.c].z;

	  float er = (dc1*dz2-dc2*dz1);
	  float ec = -(dr1*dz2-dr2*dz1);
	  float ez = (dr1*dc2-dr2*dc1);

	  dr = mr-r; // points from new facet towards [mr,mc,mz]
	  dc = mc-c;
	  dz = mz-z;
	  // make it point outwards.

	  float dromadery = dr*er +  dc*ec + dz*ez;

	  if( dromadery > 0 ){
	    Tnew.er = -er;
	    Tnew.ec = -ec;
	    Tnew.ez = -ez;
	  }
	  else{
	    Tnew.er = er;
	    Tnew.ec = ec;
	    Tnew.ez = ez;
	  }	    

	  // update the touching triangle tBC
	  int B = hull[xid].b, C = hull[xid].c;
	  if( (tBC.a == B && tBC.b == C ) ||(tBC.a == C && tBC.b == B ) ){
	    tBC.ab = hull.size();
	  }
	  else if( (tBC.a == B && tBC.c == C ) ||(tBC.a == C && tBC.c == B ) ){
	    tBC.ac = hull.size();
	  }
	  else if( (tBC.b == B && tBC.c == C ) ||(tBC.b == C && tBC.c == B ) ){
	    tBC.bc = hull.size();
	  }
	  else{
	    cerr << "Oh crap, rocket engine failure" << endl;
	    cerr << "numeric stability issue, exiting now." << endl;
	    exit(9);
	  }

	  hull.push_back(Tnew);
	}

      }

      // patch up the new triangles in hull.

      int numN = hull.size();
      std::vector<Snork> norts;
      Snork snort;
      for( int q = numN-1; q>= numh; q--){
	if( hull[q].keep > 1){
	  snort.id = q;
	  snort.a = hull[q].b;
	  snort.b = 1;
	  norts.push_back(snort);

	  snort.a = hull[q].c;
	  snort.b = 0;
	  norts.push_back(snort);

	  hull[q].keep = 1;
	  
	}
      }

      sort( norts.begin(), norts.end());
      int nums = norts.size();

      if( nums >= 2 ){
	for( int s=0; s<nums-1; s++){
	  if( norts[s].a == norts[s+1].a ){
	    // link triangle sides.
	    if( norts[s].b == 1){
	      hull[norts[s].id].ab = norts[s+1].id;
	    }
	    else{
	      hull[norts[s].id].ac = norts[s+1].id;
	    }

	    if( norts[s+1].b == 1){
	      hull[norts[s+1].id].ab = norts[s].id;
	    }
	    else{
	      hull[norts[s+1].id].ac = norts[s].id;
	    }
	  }
	}
      }

    }
    /*else{
      cerr << "still in the coplanar state you fucking baboon..." << endl;
      // rather complicated and need to add points to the 2D-hull as two faced triangles.
      exit(0);
      }*/

    
  }


  return(0);
}

// add a point coplanar to the existing planar hull in 3D
// this is not an efficient routine and should only be used 
// to add in duff (coplanar) pts at the start.
// it should not be called in doing a Delaunay triangulation 
// of a 2D set raised into 3D.

void add_coplanar( std::vector<R3> &pts, std::vector<Tri> &hull, int id)
{

  int numh = hull.size();
  float er, ec, ez;
  for( int k=0; k<numh; k++){
    //find vizible edges. from external edges.

    if( hull[k].c == hull[hull[k].ab].c ){ // ->  ab is an external edge.
      // test this edge for visibility from new point pts[id].
      int A = hull[k].a;
      int B = hull[k].b;
      int C = hull[k].c;
      
      int zot = cross_test( pts, A, B, C, id, er,ec,ez); 

      if( zot < 0 ){ // visible edge facet, create 2 new hull plates.
	Tri up, down;
	up.keep = 2;
	up.id = hull.size();
	up.a = id;
	up.b = A;
	up.c = B;

	up.er = er; up.ec = ec; up.ez = ez;	
	up.ab = -1; up.ac = -1;

	down.keep = 2;
	down.id = hull.size()+1;
	down.a = id;
	down.b = A;
	down.c = B;
	
	down.ab = -1; down.ac = -1;
	down.er = -er; down.ec = -ec; down.ez = -ez;


	float xx = hull[k].er*er + hull[k].ec*ec + hull[k].ez*ez;
	if( xx > 0 ){
	  up.bc = k; 
	  down.bc = hull[k].ab; 

	  hull[k].ab = up.id;
	  hull[down.bc].ab = down.id;
	}
	else{
	  down.bc = k; 
	  up.bc = hull[k].ab; 

	  hull[k].ab = down.id;
	  hull[up.bc].ab = up.id;
	}

	hull.push_back(up);
	hull.push_back(down);
      }
    }



    if( hull[k].a == hull[hull[k].bc].a ){   // bc is an external edge.
      // test this edge for visibility from new point pts[id].
      int A = hull[k].b;
      int B = hull[k].c;
      int C = hull[k].a;
      
      int zot = cross_test( pts, A, B, C, id, er,er,ez); 

      if( zot < 0 ){ // visible edge facet, create 2 new hull plates.
	Tri up, down;
	up.keep = 2;
	up.id = hull.size();
	up.a = id;
	up.b = A;
	up.c = B;

	up.er = er; up.ec = ec; up.ez = ez;	
	up.ab = -1; up.ac = -1;

	down.keep = 2;
	down.id = hull.size()+1;
	down.a = id;
	down.b = A;
	down.c = B;
	
	down.ab = -1; down.ac = -1;
	down.er = -er; down.ec = -ec; down.ez = -ez;


	float xx = hull[k].er*er + hull[k].ec*ec + hull[k].ez*ez;
	if( xx > 0 ){
	  up.bc = k; 
	  down.bc = hull[k].bc; 

	  hull[k].bc = up.id;
	  hull[down.bc].bc = down.id;
	}
	else{
	  down.bc = k; 
	  up.bc = hull[k].bc; 

	  hull[k].bc = down.id;
	  hull[up.bc].bc = up.id;
	}

	hull.push_back(up);
	hull.push_back(down);
      }
    }




    if( hull[k].b == hull[hull[k].ac].b ){   // ac is an external edge.
      // test this edge for visibility from new point pts[id].
      int A = hull[k].a;
      int B = hull[k].c;
      int C = hull[k].b;
      
      int zot = cross_test( pts, A, B, C, id, er,er,ez); 

      if( zot < 0 ){ // visible edge facet, create 2 new hull plates.
	Tri up, down;
	up.keep = 2;
	up.id = hull.size();
	up.a = id;
	up.b = A;
	up.c = B;

	up.er = er; up.ec = ec; up.ez = ez;	
	up.ab = -1; up.ac = -1;

	down.keep = 2;
	down.id = hull.size()+1;
	down.a = id;
	down.b = A;
	down.c = B;
	
	down.ab = -1; down.ac = -1;
	down.er = -er; down.ec = -ec; down.ez = -ez;


	float xx = hull[k].er*er + hull[k].ec*ec + hull[k].ez*ez;
	if( xx > 0 ){
	  up.bc = k; 
	  down.bc = hull[k].ac; 

	  hull[k].ac = up.id;
	  hull[down.bc].ac = down.id;
	}
	else{
	  down.bc = k; 
	  up.bc = hull[k].ac; 

	  hull[k].ac = down.id;
	  hull[up.bc].ac = up.id;
	}

	hull.push_back(up);
	hull.push_back(down);
      }
    }

  }


  // fix up the non asigned hull adjecencies (correctly).
  
  int numN = hull.size();
  std::vector<Snork> norts;
  Snork snort;
  for( int q = numN-1; q>= numh; q--){
    if( hull[q].keep > 1){
      snort.id = q;
      snort.a = hull[q].b;
      snort.b = 1;
      norts.push_back(snort);
      
      snort.a = hull[q].c;
      snort.b = 0;
      norts.push_back(snort);
      
      hull[q].keep = 1;
      
    }
  }
   
  sort( norts.begin(), norts.end());
  int nums = norts.size();
  Snork snor;
  snor.id = -1; snor.a=-1; snor.b=-1;
  norts.push_back(snor);
  snor.id = -2; snor.a=-2; snor.b=-2;
  norts.push_back(snor);

  if( nums >= 2 ){
    for( int s=0; s<nums-1; s++){
      if( norts[s].a == norts[s+1].a ){
	// link triangle sides.
	if( norts[s].a != norts[s+2].a ){ // edge of figure case
	  if( norts[s].b == 1){
	    hull[norts[s].id].ab = norts[s+1].id;
	  }
	  else{
	    hull[norts[s].id].ac = norts[s+1].id;
	  }
	  
	  if( norts[s+1].b == 1){
	    hull[norts[s+1].id].ab = norts[s].id;
	  }
	  else{
	    hull[norts[s+1].id].ac = norts[s].id;
	  }
	  s++;
	}
	else{ // internal figure boundary 4 junction case.
	  int s1 = s+1, s2 = s+2, s3 = s+3;
	  int id = norts[s].id;
	  int id1 = norts[s1].id;
	  int id2 = norts[s2].id;
	  int id3 = norts[s3].id;

	  // check normal directions of id and id1..3
	  float barf = hull[id].er*hull[id1].er + hull[id].ec*hull[id1].ec + hull[id].ez*hull[id1].ez;
	  if( barf > 0){
	  }
	  else{
	    barf = hull[id].er*hull[id2].er + hull[id].ec*hull[id2].ec + hull[id].ez*hull[id2].ez;
	    if( barf > 0 ){
	      int tmp = id2; id2 = id1; id1 = tmp;
	      tmp = s2; s2 = s1; s1 = tmp;
	    }
	    else{
	      barf = hull[id].er*hull[id3].er + hull[id].ec*hull[id3].ec + hull[id].ez*hull[id3].ez;
	      if( barf > 0 ){
		int tmp = id3; id3 = id1; id1 = tmp;
		tmp = s3; s3 = s1; s1 = tmp;
	      }
	    }
	  }
	
	  if( norts[s].b == 1){
	    hull[norts[s].id].ab = norts[s1].id;
	  }
	  else{
	    hull[norts[s].id].ac = norts[s1].id;
	  }
	  
	  if( norts[s1].b == 1){
	    hull[norts[s1].id].ab = norts[s].id;
	  }
	  else{
	    hull[norts[s1].id].ac = norts[s].id;
	  }
	      

	  // use s2 and s3

	  if( norts[s2].b == 1){
	    hull[norts[s2].id].ab = norts[s3].id;
	  }
	  else{
	    hull[norts[s2].id].ac = norts[s3].id;
	  }
	  
	  if( norts[s3].b == 1){
	    hull[norts[s3].id].ab = norts[s2].id;
	  }
	  else{
	    hull[norts[s3].id].ac = norts[s2].id;
	  }

	  s+=3;
	}
      }
    }
  }
  

 
  return;
}

// cross product relative sign test.
// remmebers the cross product of (ab x cx)

int cross_test( std::vector<R3> &pts, int A, int B, int C, int X, 
		float &er, float &ec, float &ez)
{ 
  
  float Ar = pts[A].r;
  float Ac = pts[A].c;
  float Az = pts[A].z;

  float Br = pts[B].r;
  float Bc = pts[B].c;
  float Bz = pts[B].z;

  float Cr = pts[C].r;
  float Cc = pts[C].c;
  float Cz = pts[C].z;

  float Xr = pts[X].r;
  float Xc = pts[X].c;
  float Xz = pts[X].z;

  float ABr = Br-Ar;
  float ABc = Bc-Ac;
  float ABz = Bz-Az;

  float ACr = Cr-Ar;
  float ACc = Cc-Ac;
  float ACz = Cz-Az;

  float AXr = Xr-Ar;
  float AXc = Xc-Ac;
  float AXz = Xz-Az;

  er = (ABc*AXz-ABz*AXc);
  ec = -(ABr*AXz-ABz*AXr);
  ez = (ABr*AXc-ABc*AXr);
  
  float kr = (ABc*ACz-ABz*ACc);
  float kc = -(ABr*ACz-ABz*ACr);
  float kz = (ABr*ACc-ABc*ACr);

  //  look at sign of (ab x ac).(ab x ax)

  float globit =  kr * er +  kc * ec + kz * ez;

  if( globit > 0 ) return(1);
  if( globit == 0 ) return(0);
  
  return(-1);
  
}
