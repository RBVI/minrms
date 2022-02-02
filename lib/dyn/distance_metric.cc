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


#include <cmath>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_utils.h>
#include "distance_metric.h"

//Hopefully, the following line includes an RCS-version-string in the binary.
static char rcs_id[] = "@(#)$Header: /usr/local/src/minrms/lib/dyn/RCS/distance_metric.cc,v 3.4 2002/10/15 21:41:14 conrad Exp $";

//extern global variables
//extern Real g_cost_threshold; //defined in dyn.cc

using namespace minrms;

using namespace DistanceMetric;

const Real DistanceMetric::UPPER_BOUND = 1.0e26f;//The square of
                                           //1 kilometer
                                           //(as measured in angstroms).
                                           //The distance between any
                                           //two residues should be less
                                           //than 1km  ("Well,"
                                           //"...one would hope, anyway"
                                           // -Eric says.).
const Real DistanceMetric::NO_MAX_PHYS_DISTANCE = -1.0f;
/* DEK: commented out static because it's declared extern in the class */
/* static */ Real DistanceMetric::max_allowed_dist_metric = UPPER_BOUND;
                              //"max_allowed_dist_metric" is a cutoff for the
                              //maximum value of whatever distance metric
                              //you are using.  If the distance between two
                              //residues according to the metric you are
                              //using exceeds this value, then the residues
                              //will not be matched.
                              //To dissable this behavior, set
                              //max_distance_allowed = UPPER_BOUND
                              //(this is the default setting)

void  DistanceMetric::SetMaxMetricDistance(Real r)
{  max_allowed_dist_metric = r; }

Real DistanceMetric::GetMaxMetricDistance()
{  return max_allowed_dist_metric; }

void
DistanceMetric::SetMaxPhysicalDistance(Real r)
{
  if (r == NO_MAX_PHYS_DISTANCE)
    max_allowed_dist_metric = UPPER_BOUND;
  else
  {
    assert(Biopolymer::Residue::NumBackboneAtoms() != 0);
    max_allowed_dist_metric =
      SQR(r) * Biopolymer::Residue::NumBackboneAtoms();
  }
}

Real
DistanceMetric::GetMaxPhysicalDistance()
{
  if (max_allowed_dist_metric == UPPER_BOUND)
    return NO_MAX_PHYS_DISTANCE;
  else
  {
    assert(Biopolymer::Residue::NumBackboneAtoms() != 0);
    return sqrt(max_allowed_dist_metric
                /
                Biopolymer::Residue::NumBackboneAtoms());
  }
}


Real
DistanceMetric::DivideByNthenSqrt(Real total_cost, int n)
{
  if (n < 1)
    ERR("DivBy0. something bad happened: "
        << __FILE__ << ":" << __LINE__);
  if (total_cost < 0.0)
    ERR("sqrt(" << total_cost << "). something bad happened: "
        << __FILE__ << ":" << __LINE__);
  return sqrt(total_cost / n);
}

Real
DistanceMetric::SquaredDistanceBetweenIandJ(Biopolymer::const_iterator i,
                                            Biopolymer::const_iterator j)
{
  return InlinedSquaredDistanceBetweenIandJ(i,j);
}



Real
DistanceMetric::DistanceBetweenIandJ(Biopolymer::const_iterator i,
                                     Biopolymer::const_iterator j)
{
  return InlinedDistanceBetweenIandJ(i,j);
}

