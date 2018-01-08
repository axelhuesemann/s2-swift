//
//  S2EdgeUtility.swift
//  s2-swift
//

import Foundation


// edgeClipErrorUVCoord is the maximum error in a u- or v-coordinate
// compared to the exact result, assuming that the points A and B are in
// the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
let edgeClipErrorUVCoord = 2.25 * Cell.dblEpsilon

// edgeClipErrorUVDist is the maximum distance from a clipped point to
// the corresponding exact result. It is equal to the error in a single
// coordinate because at most one coordinate is subject to error.
let edgeClipErrorUVDist = 2.25 * Cell.dblEpsilon

// faceClipErrorRadians is the maximum angle between a returned vertex
// and the nearest point on the exact edge AB. It is equal to the
// maximum directional error in PointCross, plus the error when
// projecting points onto a cube face.
let faceClipErrorRadians = 3 * Cell.dblEpsilon

// faceClipErrorDist is the same angle expressed as a maximum distance
// in (u,v)-space. In other words, a returned vertex is at most this far
// from the exact edge AB projected into (u,v)-space.
let faceClipErrorUVDist = 9 * Cell.dblEpsilon

// faceClipErrorUVCoord is the maximum angle between a returned vertex
// and the nearest point on the exact edge AB expressed as the maximum error
// in an individual u- or v-coordinate. In other words, for each
// returned vertex there is a point on the exact edge AB whose u- and
// v-coordinates differ from the vertex by at most this amount.
let faceClipErrorUVCoord = 9.0 * (1.0 / sqrt(2.0)) * Cell.dblEpsilon

// intersectsRectErrorUVDist is the maximum error when computing if a point
// intersects with a given Rect. If some point of AB is inside the
// rectangle by at least this distance, the result is guaranteed to be true;
// if all points of AB are outside the rectangle by at least this distance,
// the result is guaranteed to be false. This bound assumes that rect is
// a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it
// (e.g., by 1e-10 or less).
let intersectsRectErrorUVDist = 3 * sqrt(2.0) * Cell.dblEpsilon

// intersectionError can be set somewhat arbitrarily, because the algorithm
// uses more precision if necessary in order to achieve the specified error.
// The only strict requirement is that intersectionError >= dblEpsilon
// radians. However, using a larger error tolerance makes the algorithm more
// efficient because it reduces the number of cases where exact arithmetic is
// needed.
let intersectionError = 4 * Cell.dblEpsilon

// intersectionMergeRadius is used to ensure that intersection points that
// are supposed to be coincident are merged back together into a single
// vertex. This is required in order for various polygon operations (union,
// intersection, etc) to work correctly. It is twice the intersection error
// because two coincident intersection points might have errors in
// opposite directions.
let intersectionMergeRadius = 2 * intersectionError

// cellPadding defines the total error when clipping an edge which comes
// from two sources:
// (1) Clipping the original spherical edge to a cube face (the face edge).
//     The maximum error in this step is faceClipErrorUVCoord.
// (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
//     The maximum error in this step is edgeClipErrorUVCoord.
// Finally, since we encounter the same errors when clipping query edges, we
// double the total error so that we only need to pad edges during indexing
// and not at query time.
let cellPadding = 2.0 * (faceClipErrorUVCoord + edgeClipErrorUVCoord)

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

// A Crossing indicates how edges cross.
enum Crossing: Int {
  // Cross means the edges cross.
  case cross = 0
  // MaybeCross means two vertices from different edges are the same.
  case maybeCross = 1
  // DoNotCross means the edges do not cross.
  case doNotCross = 2
}


// EdgeCrosser allows edges to be efficiently tested for intersection with a
// given fixed edge AB. It is especially efficient when testing for
// intersection with an edge chain connecting vertices v0, v1, v2, ...
struct EdgeCrosser {
  let a: S2Point
  let b: S2Point
  let aXb: S2Point

	// To reduce the number of calls to expensiveSign, we compute an
	// outward-facing tangent at A and B if necessary. If the plane
	// perpendicular to one of these tangents separates AB from CD (i.e., one
	// edge on each side) then there is no intersection.
	let aTangent: S2Point // Outward-facing tangent at A.
  let bTangent: S2Point // Outward-facing tangent at B.

	// The fields below are updated for each vertex in the chain.
	var c: S2Point     // Previous vertex in the vertex chain.
  var acb: Direction // The orientation of triangle ACB.
  
  // EdgeCrosser with the fixed edge AB
  init(a: S2Point, b: S2Point, c: S2Point, acb: Direction) {
    let norm = a.pointCross(b).v
    self.a = a
    self.b = b
    aXb = a.pointCross(b)
    aTangent = S2Point(raw: a.v.cross(norm))
    bTangent = S2Point(raw: norm.cross(b.v))
    //
    self.c = c
    self.acb = acb
  }
  
}

extension EdgeCrosser {

  // NewChainEdgeCrosser is a convenience constructor that uses AB as the fixed edge,
  // and C as the first vertex of the vertex chain (equivalent to calling RestartAt(c)).
  //
  // You don't need to use this or any of the chain functions unless you're trying to
  // squeeze out every last drop of performance. Essentially all you are saving is a test
  // whether the first vertex of the current edge is the same as the second vertex of the
  // previous edge.
  init(a: S2Point, b: S2Point, c: S2Point) {
    let acb = -S2Point.triageSign(a, b, c)
    self.init(a: a, b: b, c: c, acb: acb)
  }
  
  // CrossingSign reports whether the edge AB intersects the edge CD.
  // If any two vertices from different edges are the same, returns MaybeCross.
  // If either edge is degenerate (A == B or C == D), returns DoNotCross or MaybeCross.
  //
  // Properties of CrossingSign:
  //
  //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
  //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
  //  (3) CrossingSign(a,b,c,d) == MaybeCross if a==c, a==d, b==c, b==d
  //  (3) CrossingSign(a,b,c,d) == DoNotCross or MaybeCross if a==b or c==d
  //
  // Note that if you want to check an edge against a chain of other edges,
  // it is slightly more efficient to use the single-argument version
  // ChainCrossingSign below.
  mutating func crossingSign(_ c: S2Point, d: S2Point) -> Crossing {
    if c != self.c {
      restartAt(c)
    }
    return chainCrossingSign(d)
  }

  // EdgeOrVertexCrossing reports whether if CrossingSign(c, d) > 0, or AB and
  // CD share a vertex and VertexCrossing(a, b, c, d) is true.
  //
  // This method extends the concept of a "crossing" to the case where AB
  // and CD have a vertex in common. The two edges may or may not cross,
  // according to the rules defined in VertexCrossing above. The rules
  // are designed so that point containment tests can be implemented simply
  // by counting edge crossings. Similarly, determining whether one edge
  // chain crosses another edge chain can be implemented by counting.
  mutating func edgeOrVertexCrossing(_ c: S2Point, d: S2Point) -> Bool {
    if c != self.c {
      restartAt(c)
    }
    return edgeOrVertexChainCrossing(d)
  }

  // RestartAt sets the current point of the edge crosser to be c.
  // Call this method when your chain 'jumps' to a new place.
  // The argument must point to a value that persists until the next call.
  mutating func restartAt(_ c: S2Point) {
    self.c = c
    acb = -S2Point.triageSign(a, b, self.c)
  }

  // ChainCrossingSign is like CrossingSign, but uses the last vertex passed to one of
  // the crossing methods (or RestartAt) as the first vertex of the current edge.
  mutating func chainCrossingSign(_ d: S2Point) -> Crossing {
    // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
    // all be oriented the same way (CW or CCW). We keep the orientation of ACB
    // as part of our state. When each new point D arrives, we compute the
    // orientation of BDA and check whether it matches ACB. This checks whether
    // the points C and D are on opposite sides of the great circle through AB.

    // Recall that triageSign is invariant with respect to rotating its
    // arguments, i.e. ABD has the same orientation as BDA.
    let bda = S2Point.triageSign(a, b, d)
    if acb == -bda && bda != .indeterminate {
      // The most common case -- triangles have opposite orientations. Save the
      // current vertex D as the next vertex C, and also save the orientation of
      // the new triangle ACB (which is opposite to the current triangle BDA).
      c = d
      acb = -bda
      return .doNotCross
    }
    return crossingSign(d, bda: bda)
  }

  // EdgeOrVertexChainCrossing is like EdgeOrVertexCrossing, but uses the last vertex
  // passed to one of the crossing methods (or RestartAt) as the first vertex of the current edge.
  mutating func edgeOrVertexChainCrossing(_ d: S2Point) -> Bool {
    // We need to copy e.c since it is clobbered by ChainCrossingSign.
    let c = self.c
    switch chainCrossingSign(d) {
    case .doNotCross: return false
    case .cross: return true
    default: break
    }
    return vertexCrossing(a, b, c, d)
  }

  // crossingSign handle the slow path of CrossingSign.
  mutating func crossingSign(_ d: S2Point, bda: Direction) -> Crossing {
    // Compute the actual result, and then save the current vertex D as the next
    // vertex C, and save the orientation of the next triangle ACB (which is
    // opposite to the current triangle BDA).
    defer {
      c = d
      acb = -bda
    }

    // RobustSign is very expensive, so we avoid calling it if at all possible.
    // First eliminate the cases where two vertices are equal.
    if a == c || a == d || b == c || b == d {
      return .maybeCross
    }

    // At this point, a very common situation is that A,B,C,D are four points on
    // a line such that AB does not overlap CD. (For example, this happens when
    // a line or curve is sampled finely, or when geometry is constructed by
    // computing the union of S2CellIds.) Most of the time, we can determine
    // that AB and CD do not intersect using the two outward-facing
    // tangents at A and B (parallel to AB) and testing whether AB and CD are on
    // opposite sides of the plane perpendicular to one of these tangents. This
    // is moderately expensive but still much cheaper than expensiveSign.

    // The error in RobustCrossProd is insignificant. The maximum error in
    // the call to CrossProd (i.e., the maximum norm of the error vector) is
    // (0.5 + 1/sqrt(3)) * dblEpsilon. The maximum error in each call to
    // DotProd below is dblEpsilon. (There is also a small relative error
    // term that is insignificant because we are comparing the result against a
    // constant that is very close to zero.)
    let maxError = (1.5 + 1/sqrt(3)) * Cell.dblEpsilon
    if (c.v.dot(aTangent.v) > maxError && d.v.dot(aTangent.v) > maxError) || (c.v.dot(bTangent.v) > maxError && d.v.dot(bTangent.v) > maxError) {
      return .doNotCross
    }

    // Otherwise it's time to break out the big guns.
    if acb == .indeterminate {
      acb = -S2Point.expensiveSign(a, b, c)
    }
    var bda = bda
    if bda == .indeterminate {
      bda = S2Point.expensiveSign(a, b, d)
    }
    if bda != acb {
      return .doNotCross
    }
    let cbd = -S2Point.robustSign(c, d, b)
    if cbd != acb {
      return .doNotCross
    }
    let dac = S2Point.robustSign(c, d, a)
    if dac == acb {
      return .cross
    }
    return .doNotCross
  }

}

// pointUVW represents a Point in (u,v,w) coordinate space of a cube face.
struct PointUVW {
  
  //
  let p: S2Point

  // intersectsFace reports whether a given directed line L intersects the cube face F.
  // The line L is defined by its normal N in the (u,v,w) coordinates of F.
  func intersectsFace() -> Bool {
    // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot
    // products of N with the four corner vertices (-1,-1,1), (1,-1,1), (1,1,1),
    // and (-1,1,1) do not all have the same sign. This is true exactly when
    // |Nu| + |Nv| >= |Nw|. The code below evaluates this expression exactly.
    let u = abs(p.x)
    let v = abs(p.y)
    let w = abs(p.z)
    // We only need to consider the cases where u or v is the smallest value,
    // since if w is the smallest then both expressions below will have a
    // positive LHS and a negative RHS.
    return (v >= w - u) && (u >= w - v)
  }

  // intersectsOppositeEdges reports whether a directed line L intersects two
  // opposite edges of a cube face F. This includs the case where L passes
  // exactly through a corner vertex of F. The directed line L is defined
  // by its normal N in the (u,v,w) coordinates of F.
  func intersectsOppositeEdges() -> Bool {
    // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if
    // and only exactly two of the corner vertices lie on each side of L. This
    // is true exactly when ||Nu| - |Nv|| >= |Nw|. The code below evaluates this
    // expression exactly.
    let u = abs(p.x)
    let v = abs(p.y)
    let w = abs(p.z)
    // If w is the smallest, the following line returns an exact result.
    if abs(u - v) != w {
      return abs(u - v) >= w
    }
    // Otherwise u - v = w exactly, or w is not the smallest value. In either
    // case the following returns the correct result.
    if u >= v {
      return u - w >= v
    }
    return v - w >= u
  }

  // axis represents the possible results of exitAxis.
  enum UVAxis: Int {
    case u = 0
    case v = 1
  }

  // exitAxis reports which axis the directed line L exits the cube face F on.
  // The directed line L is represented by its CCW normal N in the (u,v,w) coordinates
  // of F. It returns axisU if L exits through the u=-1 or u=+1 edge, and axisV if L exits
  // through the v=-1 or v=+1 edge. Either result is acceptable if L exits exactly
  // through a corner vertex of the cube face.
  func exitAxis() -> UVAxis {
    if intersectsOppositeEdges() {
      // The line passes through through opposite edges of the face.
      // It exits through the v=+1 or v=-1 edge if the u-component of N has a
      // larger absolute magnitude than the v-component.
      if abs(p.x) >= abs(p.y) {
        return .v
      }
      return .u
    }
    // The line passes through through two adjacent edges of the face.
    // It exits the v=+1 or v=-1 edge if an even number of the components of N
    // are negative. We test this using signbit() rather than multiplication
    // to avoid the possibility of underflow.
    let x = p.x.sign == .plus ? 1 : 0
    let y = p.y.sign == .plus ? 1 : 0
    let z = p.z.sign == .plus ? 1 : 0
    if x ^ y ^ z == 0 {
      return .v
    }
    return .u
  }

  // exitPoint returns the UV coordinates of the point where a directed line L (represented
  // by the CCW normal of this point), exits the cube face this point is derived from along
  // the given axis.
  func exitPoint(axis a: UVAxis) -> R2Point {
    switch a  {
    case .u:
      var u = -1.0
      if p.y > 0 {
        u = 1.0
      }
      return R2Point(x: u, y: (-u*p.x - p.z) / p.y)
    case .v:
      var v = -1.0
      if p.x < 0 {
        v = 1.0
      }
      return R2Point(x: (-v * p.y - p.z) / p.x, y: v)
    }
  }

}
