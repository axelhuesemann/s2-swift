//
//  S2EdgeUtility.swift
//  s2-swift
//

import Foundation


// SimpleCrossing reports whether edge AB crosses CD at a point that is interior
// to both edges. Properties:
//
//  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d)
//  (2) SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
func simpleCrossing(_ a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
	// We compute the equivalent of Sign for triangles ACB, CBD, BDA,
	// and DAC. All of these triangles need to have the same orientation
	// (CW or CCW) for an intersection to exist.

	let ab = a.v.cross(b.v)
	let acb = -(ab.dot(c.v))
	let bda = ab.dot(d.v)
	if acb * bda <= 0 {
		return false
	}

	let cd = c.v.cross(d.v)
	let cbd = -(cd.dot(b.v))
	let dac = cd.dot(a.v)
	return (acb * cbd > 0) && (acb * dac > 0)
}

// VertexCrossing reports whether two edges "cross" in such a way that point-in-polygon
// containment tests can be implemented by counting the number of edge crossings.
//
// Given two edges AB and CD where at least two vertices are identical
// (i.e. CrossingSign(a,b,c,d) == 0), the basic rule is that a "crossing"
// occurs if AB is encountered after CD during a CCW sweep around the shared
// vertex starting from a fixed reference point.
//
// Note that according to this rule, if AB crosses CD then in general CD
// does not cross AB.  However, this leads to the correct result when
// counting polygon edge crossings.  For example, suppose that A,B,C are
// three consecutive vertices of a CCW polygon.  If we now consider the edge
// crossings of a segment BP as P sweeps around B, the crossing number
// changes parity exactly when BP crosses BA or BC.
//
// Useful properties of VertexCrossing (VC):
//
//  (1) VC(a,a,c,d) == VC(a,b,c,c) == false
//  (2) VC(a,b,a,b) == VC(a,b,b,a) == true
//  (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
//  (3) If exactly one of a,b equals one of c,d, then exactly one of
//      VC(a,b,c,d) and VC(c,d,a,b) is true
//
// It is an error to call this method with 4 distinct vertices.
func vertexCrossing(_ a: S2Point, _ b: S2Point, _ c: S2Point, _ d: S2Point) -> Bool {
	// If A == B or C == D there is no intersection. We need to check this
	// case first in case 3 or more input points are identical.
	if a.approxEquals(b) || c.approxEquals(d) {
		return false
	}

	// If any other pair of vertices is equal, there is a crossing if and only
	// if OrderedCCW indicates that the edge AB is further CCW around the
	// shared vertex O (either A or B) than the edge CD, starting from an
	// arbitrary fixed reference point.
	if a.approxEquals(d) {
    return S2Point.orderedCCW(S2Point(raw: a.v.ortho()), c, b, a)
	} else if b.approxEquals(c) {
		return S2Point.orderedCCW(S2Point(raw: b.v.ortho()), d, a, b)
	} else if a.approxEquals(c) {
		return S2Point.orderedCCW(S2Point(raw: a.v.ortho()), d, b, a)
	} else if  b.approxEquals(d){
		return S2Point.orderedCCW(S2Point(raw: b.v.ortho()), c, a, b)
	}

	return false
}

// Interpolate returns the point X along the line segment AB whose distance from A
// is the given fraction "t" of the distance AB. Does NOT require that "t" be
// between 0 and 1. Note that all distances are measured on the surface of
// the sphere, so this is more complicated than just computing (1-t)*a + t*b
// and normalizing the result.
func interpolate(t: Double, a: S2Point, b: S2Point) -> S2Point {
	if t == 0.0 {
		return a
	}
	if t == 1.0 {
		return b
	}
	let ab = a.angle(b)
	return interpolateAtDistance(t * ab, a: a, b: b)
}

// InterpolateAtDistance returns the point X along the line segment AB whose
// distance from A is the angle ax.
func interpolateAtDistance(_ ax: Double, a: S2Point, b: S2Point) -> S2Point {
	let aRad = ax

	// Use PointCross to compute the tangent vector at A towards B. The
	// result is always perpendicular to A, even if A=B or A=-B, but it is not
	// necessarily unit length. (We effectively normalize it below.)
	let normal = a.pointCross(b)
	let tangent = normal.v.cross(a.v)

	// Now compute the appropriate linear combination of A and "tangent". With
	// infinite precision the result would always be unit length, but we
	// normalize it anyway to ensure that the error is within acceptable bounds.
	// (Otherwise errors can build up when the result of one interpolation is
	// fed into another interpolation.)
  let v = a.v.mul(cos(aRad)).add(tangent.mul(sin(aRad) / tangent.norm))
  return S2Point(raw: v)
}

// RectBounder is used to compute a bounding rectangle that contains all edges
// defined by a vertex chain (v0, v1, v2, ...). All vertices must be unit length.
// Note that the bounding rectangle of an edge can be larger than the bounding
// rectangle of its endpoints, e.g. consider an edge that passes through the North Pole.
//
// The bounds are calculated conservatively to account for numerical errors
// when points are converted to LatLngs. More precisely, this function
// guarantees the following:
// Let L be a closed edge chain (Loop) such that the interior of the loop does
// not contain either pole. Now if P is any point such that L.ContainsPoint(P),
// then RectBound(L).ContainsPoint(LatLngFromPoint(P)).
struct RectBounder {

  // The previous vertex in the chain.
  var a: S2Point
	// The previous vertex latitude longitude.
  var aLL: LatLng
  var bound: S2Rect

  //
  init() {
    a = S2Point.origin
    aLL = LatLng(lat: 0.0, lng: 0.0)
  	bound = S2Rect.empty
  }

  // AddPoint adds the given point to the chain. The Point must be unit length.
  mutating func add(point b: S2Point) {
    let bLL = LatLng(point: b)

    if bound.isEmpty {
      a = b
      aLL = bLL
      bound = bound.add(bLL)
      return
    }

    // First compute the cross product N = A x B robustly. This is the normal
    // to the great circle through A and B. We don't use RobustSign
    // since that method returns an arbitrary vector orthogonal to A if the two
    // vectors are proportional, and we want the zero vector in that case.
    let n = a.v.sub(b.v).cross(a.v.add(b.v)) // N = 2 * (A x B)

    // The relative error in N gets large as its norm gets very small (i.e.,
    // when the two points are nearly identical or antipodal). We handle this
    // by choosing a maximum allowable error, and if the error is greater than
    // this we fall back to a different technique. Since it turns out that
    // the other sources of error add up to at most 1.16 * dblEpsilon, and it
    // is desirable to have the total error be a multiple of dblEpsilon, we
    // have chosen the maximum error threshold here to be 3.84 * dblEpsilon.
    // It is possible to show that the error is less than this when
    //
    // n.norm >= 8 * sqrt(3) / (3.84 - 0.5 - sqrt(3)) * dblEpsilon
    //          = 1.91346e-15 (about 8.618 * dblEpsilon)
    let nNorm = n.norm
    if nNorm < 1.91346e-15 {
      // A and B are either nearly identical or nearly antipodal (to within
      // 4.309 * dblEpsilon, or about 6 nanometers on the earth's surface).
      if a.v.dot(b.v) < 0 {
        // The two points are nearly antipodal. The easiest solution is to
        // assume that the edge between A and B could go in any direction
        // around the sphere.
        bound = S2Rect.full
      } else {
        // The two points are nearly identical (to within 4.309 * dblEpsilon).
        // In this case we can just use the bounding rectangle of the points,
        // since after the expansion done by GetBound this Rect is
        // guaranteed to include the (lat,lng) values of all points along AB.
        bound = bound.union(S2Rect(latLng: aLL).add(bLL))
      }
      a = b
      aLL = bLL
      return
    }

    // Compute the longitude range spanned by AB.
    var lngAB = S1Interval.empty.add(aLL.lng).add(bLL.lng)
    if lngAB.length >= .pi - 2 * Cell.dblEpsilon {
      // The points lie on nearly opposite lines of longitude to within the
      // maximum error of the calculation. The easiest solution is to assume
      // that AB could go on either side of the pole.
      lngAB = S1Interval.full
    }

    // Next we compute the latitude range spanned by the edge AB. We start
    // with the range spanning the two endpoints of the edge:
    var latAB = R1Interval(point: aLL.lat).add(bLL.lat)

    // This is the desired range unless the edge AB crosses the plane
    // through N and the Z-axis (which is where the great circle through A
    // and B attains its minimum and maximum latitudes). To test whether AB
    // crosses this plane, we compute a vector M perpendicular to this
    // plane and then project A and B onto it.
    let m = n.cross(R3Vector(x: 0, y: 0, z: 1))
    let mA = m.dot(a.v)
    let mB = m.dot(b.v)

    // We want to test the signs of "mA" and "mB", so we need to bound
    // the error in these calculations. It is possible to show that the
    // total error is bounded by
    //
    // (1 + sqrt(3)) * dblEpsilon * nNorm + 8 * sqrt(3) * (dblEpsilon**2)
    //   = 6.06638e-16 * nNorm + 6.83174e-31

    let mError  = 6.06638e-16 * nNorm + 6.83174e-31
    if mA * mB < 0 || abs(mA) <= mError || abs(mB) <= mError {
      // Minimum/maximum latitude *may* occur in the edge interior.
      //
      // The maximum latitude is 90 degrees minus the latitude of N. We
      // compute this directly using atan2 in order to get maximum accuracy
      // near the poles.
      //
      // Our goal is compute a bound that contains the computed latitudes of
      // all S2Points P that pass the point-in-polygon containment test.
      // There are three sources of error we need to consider:
      // - the directional error in N (at most 3.84 * dblEpsilon)
      // - converting N to a maximum latitude
      // - computing the latitude of the test point P
      // The latter two sources of error are at most 0.955 * dblEpsilon
      // individually, but it is possible to show by a more complex analysis
      // that together they can add up to at most 1.16 * dblEpsilon, for a
      // total error of 5 * dblEpsilon.
      //
      // We add 3 * dblEpsilon to the bound here, and GetBound() will pad
      // the bound by another 2 * dblEpsilon.
      let maxLat = min(atan2(sqrt(n.x*n.x + n.y*n.y), abs(n.z)) + 3 * Cell.dblEpsilon, .pi/2)

      // In order to get tight bounds when the two points are close together,
      // we also bound the min/max latitude relative to the latitudes of the
      // endpoints A and B. First we compute the distance between A and B,
      // and then we compute the maximum change in latitude between any two
      // points along the great circle that are separated by this distance.
      // This gives us a latitude change "budget". Some of this budget must
      // be spent getting from A to B; the remainder bounds the round-trip
      // distance (in latitude) from A or B to the min or max latitude
      // attained along the edge AB.
      let latBudget = 2 * asin(0.5 * (a.v.sub(b.v)).norm * sin(maxLat))
      let maxDelta = 0.5 * (latBudget - latAB.length) + Cell.dblEpsilon

      // Test whether AB passes through the point of maximum latitude or
      // minimum latitude. If the dot product(s) are small enough then the
      // result may be ambiguous.
      if mA <= mError && mB >= -mError {
        latAB = R1Interval(lo: latAB.lo, hi: min(maxLat, latAB.hi+maxDelta))
      }
      if mB <= mError && mA >= -mError {
        latAB = R1Interval(lo: max(-maxLat, latAB.lo-maxDelta), hi:latAB.hi)
      }
    }
    a = b
    aLL = bLL
    bound = bound.union(S2Rect(lat: latAB, lng: lngAB))
  }

  // RectBound returns the bounding rectangle of the edge chain that connects the
  // vertices defined so far. This bound satisfies the guarantee made
  // above, i.e. if the edge chain defines a Loop, then the bound contains
  // the LatLng coordinates of all Points contained by the loop.
  func rectBound() -> S2Rect {
    return bound.expanded(LatLng(lat: 2 * Cell.dblEpsilon, lng: 0)).polarClosure()
  }

}

// ExpandForSubregions expands a bounding Rect so that it is guaranteed to
// contain the bounds of any subregion whose bounds are computed using
// ComputeRectBound. For example, consider a loop L that defines a square.
// GetBound ensures that if a point P is contained by this square, then
// LatLngFromPoint(P) is contained by the bound. But now consider a diamond
// shaped loop S contained by L. It is possible that GetBound returns a
// *larger* bound for S than it does for L, due to rounding errors. This
// method expands the bound for L so that it is guaranteed to contain the
// bounds of any subregion S.
//
// More precisely, if L is a loop that does not contain either pole, and S
// is a loop such that L.Contains(S), then
//
//   ExpandForSubregions(L.RectBound).Contains(S.RectBound).
//
extension S2Rect {
 
  func expandForSubregions() -> S2Rect {
    // Empty bounds don't need expansion.
    if isEmpty {
      return self
    }
    // First we need to check whether the bound B contains any nearly-antipodal
    // points (to within 4.309 * dblEpsilon). If so then we need to return
    // FullRect, since the subregion might have an edge between two
    // such points, and AddPoint returns Full for such edges. Note that
    // this can happen even if B is not Full for example, consider a loop
    // that defines a 10km strip straddling the equator extending from
    // longitudes -100 to +100 degrees.
    //
    // It is easy to check whether B contains any antipodal points, but checking
    // for nearly-antipodal points is trickier. Essentially we consider the
    // original bound B and its reflection through the origin B', and then test
    // whether the minimum distance between B and B' is less than 4.309 * dblEpsilon.
    // lngGap is a lower bound on the longitudinal distance between B and its
    // reflection B'. (2.5 * dblEpsilon is the maximum combined error of the
    // endpoint longitude calculations and the Length call.)
    let lngGap = max(0, .pi - lng.length - 2.5 * Cell.dblEpsilon)
    // minAbsLat is the minimum distance from B to the equator (if zero or
    // negative, then B straddles the equator).
    let minAbsLat = max(lat.lo, -lat.hi)
    // latGapSouth and latGapNorth measure the minimum distance from B to the
    // south and north poles respectively.
    let latGapSouth = .pi/2 + lat.lo
    let latGapNorth = .pi/2 - lat.hi
    if minAbsLat >= 0 {
      // The bound B does not straddle the equator. In this case the minimum
      // distance is between one endpoint of the latitude edge in B closest to
      // the equator and the other endpoint of that edge in B'. The latitude
      // distance between these two points is 2*minAbsLat, and the longitude
      // distance is lngGap. We could compute the distance exactly using the
      // Haversine formula, but then we would need to bound the errors in that
      // calculation. Since we only need accuracy when the distance is very
      // small (close to 4.309 * dblEpsilon), we substitute the Euclidean
      // distance instead. This gives us a right triangle XYZ with two edges of
      // length x = 2*minAbsLat and y ~= lngGap. The desired distance is the
      // length of the third edge z, and we have
      //
      //         z  ~=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
      //
      // Therefore the region may contain nearly antipodal points only if
      //
      //  2*minAbsLat + lngGap  <  sqrt(2) * 4.309 * dblEpsilon
      //                        ~= 1.354e-15
      //
      // Note that because the given bound B is conservative, minAbsLat and
      // lngGap are both lower bounds on their true values so we do not need
      // to make any adjustments for their errors.
      if 2 * minAbsLat + lngGap < 1.354e-15 {
        return S2Rect.full
      }
    } else if lngGap >= .pi/2 {
      // B spans at most Pi/2 in longitude. The minimum distance is always
      // between one corner of B and the diagonally opposite corner of B'. We
      // use the same distance approximation that we used above; in this case
      // we have an obtuse triangle XYZ with two edges of length x = latGapSouth
      // and y = latGapNorth, and angle Z >= Pi/2 between them. We then have
      //
      //         z  >=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
      //
      // Unlike the case above, latGapSouth and latGapNorth are not lower bounds
      // (because of the extra addition operation, and because math.Pi/2 is not
      // exactly equal to Pi/2); they can exceed their true values by up to
      // 0.75 * dblEpsilon. Putting this all together, the region may contain
      // nearly antipodal points only if
      //
      //   latGapSouth + latGapNorth  <  (sqrt(2) * 4.309 + 1.5) * dblEpsilon
      //                              ~= 1.687e-15
      if latGapSouth+latGapNorth < 1.687e-15 {
        return S2Rect.full
      }
    } else {
      // Otherwise we know that (1) the bound straddles the equator and (2) its
      // width in longitude is at least Pi/2. In this case the minimum
      // distance can occur either between a corner of B and the diagonally
      // opposite corner of B' (as in the case above), or between a corner of B
      // and the opposite longitudinal edge reflected in B'. It is sufficient
      // to only consider the corner-edge case, since this distance is also a
      // lower bound on the corner-corner distance when that case applies.
      //
      // Consider the spherical triangle XYZ where X is a corner of B with
      // minimum absolute latitude, Y is the closest pole to X, and Z is the
      // point closest to X on the opposite longitudinal edge of B'. This is a
      // right triangle (Z = Pi/2), and from the spherical law of sines we have
      //
      //     sin(z) / sin(Z)  =  sin(y) / sin(Y)
      //     sin(maxLatGap) / 1  =  sin(dMin) / sin(lngGap)
      //     sin(dMin)  =  sin(maxLatGap) * sin(lngGap)
      //
      // where "maxLatGap" = max(latGapSouth, latGapNorth) and "dMin" is the
      // desired minimum distance. Now using the facts that sin(t) >= (2/Pi)*t
      // for 0 <= t <= Pi/2, that we only need an accurate approximation when
      // at least one of "maxLatGap" or lngGap is extremely small (in which
      // case sin(t) ~= t), and recalling that "maxLatGap" has an error of up
      // to 0.75 * dblEpsilon, we want to test whether
      //
      //   maxLatGap * lngGap  <  (4.309 + 0.75) * (Pi/2) * dblEpsilon
      //                       ~= 1.765e-15
      if max(latGapSouth, latGapNorth) * lngGap < 1.765e-15 {
        return S2Rect.full
      }
    }
    // Next we need to check whether the subregion might contain any edges that
    // span (math.Pi - 2 * dblEpsilon) radians or more in longitude, since AddPoint
    // sets the longitude bound to Full in that case. This corresponds to
    // testing whether (lngGap <= 0) in lngExpansion below.
    // Otherwise, the maximum latitude error in AddPoint is 4.8 * dblEpsilon.
    // In the worst case, the errors when computing the latitude bound for a
    // subregion could go in the opposite direction as the errors when computing
    // the bound for the original region, so we need to double this value.
    // (More analysis shows that it's okay to round down to a multiple of
    // dblEpsilon.)
    //
    // For longitude, we rely on the fact that atan2 is correctly rounded and
    // therefore no additional bounds expansion is necessary.
    let latExpansion = 9 * Cell.dblEpsilon
    var lngExpansion = 0.0
    if lngGap <= 0 {
      lngExpansion = .pi
    }
    return expanded(LatLng(lat: latExpansion, lng: lngExpansion)).polarClosure()
  }
}
