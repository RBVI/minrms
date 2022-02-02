//
// Copyright (c) 2002 The Regents of the University of California.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions, and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions, and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//   3. Redistributions must acknowledge that this software was
//      originally developed by the UCSF Computer Graphics Laboratory
//      under support by the NIH National Center for Research Resources,
//      grant P41-RR01081.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _DISTANCE_METRIC_H
#define _DISTANCE_METRIC_H

#include <cassert>
using namespace std;

#include <vect_3d.h>
#include <biopolymer.h>

namespace minrms {

namespace DistanceMetric
{

  extern Real max_allowed_dist_metric; //Please treat this variable as a private
                                        //member. Do not access it directly.




  //SetMaxPhysicalDistance() is used to enable or dissable the
  //"distance-clipping" feature.  If this feature is enabled,
  //then whenever the actual physical distance between the two residues
  //exceeds the amount specified by parameter "r", any functions
  //that compute distance described in this file will instead
  //return an unattainably high number (essentially infinity).
  //  (This is to prevent alignments between residues farther apart than
  //   a certain distance.)
  //To dissable this feature, pass: "NO_MAX_PHYS_DISTANCE" for r.
  void SetMaxPhysicalDistance(Real r);
  extern const Real NO_MAX_PHYS_DISTANCE;
                               //To dissable this behavior, call
                               //max_distance_allowed = NO_MAX_DISTANCE
                               //(this is the default setting)

  //GetMaxPhysicalDistance() retreives the maximum physical distance according
  //to the parameter you set in a previous call to SetMaxPhysicalDistance().
  Real GetMaxPhysicalDistance();

  //Set and GetMaxMetricDistance() do exactly the same thing, only they use
  //the units of whatever metric is being used, instead of physical distance.
  //(An alternate way to dissable the "distance-clipping" is to call
  // SetMaxMetricDistance(UPPER_BOUND).)
  void SetMaxMetricDistance(Real r);
  extern const Real UPPER_BOUND;//The square of
                                 //1 kilometer
                                 //(as measured in angstroms).
                                 //The distance between any
                                 //two residues should be less
                                 //than 1km  ("Well,"
                                 //"...one would hope, anyway"
                                 // -Eric says.).
  Real GetMaxMetricDistance();

  #ifdef COMPILER_BUG1
  //Since the "DistanceMetric" namespace is actually nested inside the "minrms"
  //namespace, the following line shouldn't be necessary.  However
  //current Digital compiler has a bug which requires this line:
  using namespace minrms;
  #endif //#ifdef COMPILER_BUG1

  //The next two functions below offer a measure of the distance
  //between two residues from two molecules.
  Real
    SquaredDistanceBetweenIandJ(Biopolymer::const_iterator i,
                                Biopolymer::const_iterator j);

  Real
    DistanceBetweenIandJ(Biopolymer::const_iterator i,
                         Biopolymer::const_iterator j);


  //For slightly faster access, use the inlined versions.
  inline Real
    InlinedSquaredDistanceBetweenIandJ(Biopolymer::const_iterator i,
                                       Biopolymer::const_iterator j)
  {
    Real output = 0.0f;

    Biopolymer::Residue::const_iterator pa[2];
    for(int q = 0; 
        q < Biopolymer::Residue::NumBackboneAtoms();
        ++q)
    {
      pa[0] = (*i).GetBackboneAtom(q);
      pa[1] = (*j).GetBackboneAtom(q);
      assert(pa[0] != i->end());
      assert(pa[1] != j->end());
      Vect3 displacement_between_ij;
      SubtractVect3((*pa[0]).second.xyz,
                    (*pa[1]).second.xyz,
                    displacement_between_ij);
      output += DotProduct3(displacement_between_ij, displacement_between_ij);
    }
    //if (output > DistanceMetric::max_allowed_dist_metric)
    //return DistanceMetric::UPPER_BOUND;
    if (output > max_allowed_dist_metric)
      return UPPER_BOUND;
    return output;
  }



  inline Real
    InlinedDistanceBetweenIandJ(Biopolymer::const_iterator i,
                                Biopolymer::const_iterator j)
  {
    Real output = 0.0f;

    Biopolymer::Residue::const_iterator pa[2];
    for(int q = 0; 
        q < Biopolymer::Residue::NumBackboneAtoms();
        ++q)
    {
      pa[0] = (*i).GetBackboneAtom(q);
      pa[1] = (*j).GetBackboneAtom(q);
      assert(pa[0] != i->end());
      assert(pa[1] != j->end());
      Vect3 displacement_between_ij;
      SubtractVect3((*pa[0]).second.xyz,
                    (*pa[1]).second.xyz,
                    displacement_between_ij);
      output += sqrt(DotProduct3(displacement_between_ij, displacement_between_ij));
    }
    //if (output > DistanceMetric::max_allowed_dist_metric)
    //return DistanceMetric::UPPER_BOUND;
    if (output > max_allowed_dist_metric)
      return UPPER_BOUND;
    return output;
  }


  // The next function is only occasionally used.  It's very simple, it does
  // what the name says.  Instead of computing RMSD between two sets of
  // of points, many functions return their answer in terms of sum-quared-
  // distance between the two sets of points.  This function is used to
  // compute RMSD from sum-squared-distance.
  Real
    DivideByNthenSqrt(Real x, int n);

  } //namespace DistanceMetric


} //namespace minrms

#endif //#ifndef _DISTANCE_METRIC_H





