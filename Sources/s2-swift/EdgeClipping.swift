//
//  EdgeClipping.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/12/18.
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

// ClipToFace returns the (u,v) coordinates for the portion of the edge AB that
// intersects the given face, or false if the edge AB does not intersect.
// This method guarantees that the clipped vertices lie within the [-1,1]x[-1,1]
// cube face rectangle and are within faceClipErrorUVDist of the line AB, but
// the results may differ from those produced by FaceSegments.
func clipToFace(a: S2Point, b: S2Point, face: Int) -> (aUV: R2Point, bUV: R2Point, intersects: Bool) {
  return clipToPaddedFace(a: a, b: b, f: face, padding: 0.0)
}

// ClipToPaddedFace returns the (u,v) coordinates for the portion of the edge AB that
// intersects the given face, but rather than clipping to the square [-1,1]x[-1,1]
// in (u,v) space, this method clips to [-R,R]x[-R,R] where R=(1+padding).
// Padding must be non-negative.
func clipToPaddedFace(a: S2Point, b: S2Point, f: Int, padding: Double) -> (aUV: R2Point, bUV: R2Point, intersects: Bool) {
  // Fast path: both endpoints are on the given face.
  if S2Cube.face(point: a) == f && S2Cube.face(point: b) == f {
    let (au, av) = S2Cube.validFaceXYZToUV(face: f, point: a)
    let (bu, bv) = S2Cube.validFaceXYZToUV(face: f, point: b)
    return (R2Point(x: au, y: av), R2Point(x: bu, y: bv), true)
  }
  // Convert everything into the (u,v,w) coordinates of the given face. Note
  // that the cross product *must* be computed in the original (x,y,z)
  // coordinate system because PointCross (unlike the mathematical cross
  // product) can produce different results in different coordinate systems
  // when one argument is a linear multiple of the other, due to the use of
  // symbolic perturbations.
  var normUVW = PointUVW(p: S2Cube.faceXYZtoUVW(face: f, point: a.pointCross(b)))
  let aUVW = PointUVW(p: S2Cube.faceXYZtoUVW(face: f, point: a))
  let bUVW = PointUVW(p: S2Cube.faceXYZtoUVW(face: f, point: b))
  // Padding is handled by scaling the u- and v-components of the normal.
  // Letting R=1+padding, this means that when we compute the dot product of
  // the normal with a cube face vertex (such as (-1,-1,1)), we will actually
  // compute the dot product with the scaled vertex (-R,-R,1). This allows
  // methods such as intersectsFace, exitAxis, etc, to handle padding
  // with no further modifications.
  let scaleUV = 1 + padding
  let p = S2Point(raw: R3Vector(x: scaleUV * normUVW.p.x, y: scaleUV * normUVW.p.y, z: normUVW.p.z).normalized())
  let scaledN = PointUVW(p: p)
  if !scaledN.intersectsFace() {
    let aUV = R2Point(x: 0, y: 0)
    let bUV = R2Point(x: 0, y: 0)
    return (aUV, bUV, false)
  }
  // TODO(roberts): This is a workaround for extremely small vectors where some
  // loss of precision can occur in Normalize causing underflow. When PointCross
  // is updated to work around this, this can be removed.
  if max(abs(normUVW.p.x), max(abs(normUVW.p.y), abs(normUVW.p.z))) < ldexp(1, -511) {
    normUVW = PointUVW(p: S2Point(raw: normUVW.p.mul(ldexp(1, 563)).normalized()))
  }
  normUVW = PointUVW(p: normUVW.p)
  let aTan = PointUVW(p: S2Point(raw: normUVW.p.cross(aUVW.p).normalized()))
  let bTan = PointUVW(p: S2Point(raw: bUVW.p.cross(normUVW.p).normalized()))
  // As described in clipDestination, if the sum of the scores from clipping the two
  // endpoints is 3 or more, then the segment does not intersect this face.
  let scaledM = PointUVW(p: S2Point(raw: scaledN.p.mul(-1)))
  let (aUV, aScore) = clipDestination(a: bUVW, b: aUVW, scaledN: scaledM, aTan: bTan, bTan: aTan, scaleUV: scaleUV)
  let (bUV, bScore) = clipDestination(a: aUVW, b: bUVW, scaledN: scaledN, aTan: aTan, bTan: bTan, scaleUV: scaleUV)
  return (aUV, bUV, aScore + bScore < 3)
}

// ClipEdge returns the portion of the edge defined by AB that is contained by the
// given rectangle. If there is no intersection, false is returned and aClip and bClip
// are undefined.
func clipEdge(a: R2Point, b: R2Point, clip: R2Rect) -> (aClip: R2Point, bClip: R2Point)? {
  // Compute the bounding rectangle of AB, clip it, and then extract the new
  // endpoints from the clipped bound.
  let bound2 = R2Rect(p0: a, p1: b)
  let (bound, intersects) = clipEdgeBound(a: a, b: b, clip: clip, bound: bound2)
  if !intersects {
    return nil
  }
  let ai = (a.x > b.x) ? 1 : 0
  let aj = (a.y > b.y) ? 1 : 0
  return (bound.vertex(i: ai, j: aj), bound.vertex(i: 1 - ai, j: 1 - aj))
}

// The three functions below (sumEqual, intersectsFace, intersectsOppositeEdges)
// all compare a sum (u + v) to a third value w. They are implemented in such a
// way that they produce an exact result even though all calculations are done
// with ordinary floating-point operations. Here are the principles on which these
// functions are based:
//
// A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
//
// B. If u + v < w in exact arithmetic, then at least one of the following
//    expressions is true in floating-point:
//       u + v < w
//       u < w - v
//       v < w - u
//
// Proof: By rearranging terms and substituting ">" for "<", we can assume
// that all values are non-negative.  Now clearly "w" is not the smallest
// value, so assume WLOG that "u" is the smallest.  We want to show that
// u < w - v in floating-point.  If v >= w/2, the calculation of w - v is
// exact since the result is smaller in magnitude than either input value,
// so the result holds.  Otherwise we have u <= v < w/2 and w - v >= w/2
// (even in floating point), so the result also holds.

// sumEqual reports whether u + v == w exactly.
func sumEqual(u: Double, v: Double, w: Double) -> Bool {
  return (u + v == w) && (u == w - v) && (v == w - u)
}

// axis represents the possible results of exitAxis.
enum UVWAxis: Int {
  case u = 0, v, w
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
  
  // exitAxis reports which axis the directed line L exits the cube face F on.
  // The directed line L is represented by its CCW normal N in the (u,v,w) coordinates
  // of F. It returns axisU if L exits through the u=-1 or u=+1 edge, and axisV if L exits
  // through the v=-1 or v=+1 edge. Either result is acceptable if L exits exactly
  // through a corner vertex of the cube face.
  func exitAxis() -> UVWAxis {
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
  func exitPoint(axis a: UVWAxis) -> R2Point {
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
    default: fatalError("???")
    }
  }
  
}

// clipDestination returns a score which is used to indicate if the clipped edge AB
// on the given face intersects the face at all. This function returns the score for
// the given endpoint, which is an integer ranging from 0 to 3. If the sum of the scores
// from both of the endpoints is 3 or more, then edge AB does not intersect this face.
//
// First, it clips the line segment AB to find the clipped destination B' on a given
// face. (The face is specified implicitly by expressing *all arguments* in the (u,v,w)
// coordinates of that face.) Second, it partially computes whether the segment AB
// intersects this face at all. The actual condition is fairly complicated, but it
// turns out that it can be expressed as a "score" that can be computed independently
// when clipping the two endpoints A and B.
func clipDestination(a: PointUVW, b: PointUVW, scaledN: PointUVW, aTan: PointUVW, bTan: PointUVW, scaleUV: Double) -> (R2Point, Int) {
  var uv: R2Point
  // Optimization: if B is within the safe region of the face, use it.
  let maxSafeUVCoord = 1 - faceClipErrorUVCoord
  if b.p.z > 0 {
    uv = R2Point(x: b.p.x / b.p.z, y: b.p.y / b.p.z)
    if max(abs(uv.x), abs(uv.y)) <= maxSafeUVCoord {
      return (uv, 0)
    }
  }
  // Otherwise find the point B' where the line AB exits the face.
  uv = scaledN.exitPoint(axis: scaledN.exitAxis()).mul(scaleUV)
  let p = PointUVW(p: S2Point(raw: R3Vector(x: uv.x, y: uv.y, z: 1.0)))
  // Determine if the exit point B' is contained within the segment. We do this
  // by computing the dot products with two inward-facing tangent vectors at A
  // and B. If either dot product is negative, we say that B' is on the "wrong
  // side" of that point. As the point B' moves around the great circle AB past
  // the segment endpoint B, it is initially on the wrong side of B only; as it
  // moves further it is on the wrong side of both endpoints; and then it is on
  // the wrong side of A only. If the exit point B' is on the wrong side of
  // either endpoint, we can't use it; instead the segment is clipped at the
  // original endpoint B.
  //
  // We reject the segment if the sum of the scores of the two endpoints is 3
  // or more. Here is what that rule encodes:
  //  - If B' is on the wrong side of A, then the other clipped endpoint A'
  //    must be in the interior of AB (otherwise AB' would go the wrong way
  //    around the circle). There is a similar rule for A'.
  //  - If B' is on the wrong side of either endpoint (and therefore we must
  //    use the original endpoint B instead), then it must be possible to
  //    project B onto this face (i.e., its w-coordinate must be positive).
  //    This rule is only necessary to handle certain zero-length edges (A=B).
  var score = 0
  if p.p.sub(a.p).dot(aTan.p.v) < 0 {
    score = 2 // B' is on wrong side of A.
  } else if p.p.sub(b.p).dot(bTan.p.v) < 0 {
    score = 1 // B' is on wrong side of B.
  }
  if score > 0 { // B' is not in the interior of AB.
    if b.p.z <= 0 {
      score = 3 // B cannot be projected onto this face.
    } else {
      uv = R2Point(x: b.p.x / b.p.z, y: b.p.y / b.p.z)
    }
  }
  return (uv, score)
}

/// Returns the interval with the specified endpoint updated to
/// the given value. If the value lies beyond the opposite endpoint, nothing is
/// changed and false is returned.
func updateEndpoint(bound: R1Interval, highEndpoint: Bool, value: Double) -> R1Interval? {
  if highEndpoint {
    if bound.lo > value {
      return nil
    }
    if bound.hi > value {
      return R1Interval(lo: bound.lo, hi: value)
    }
  } else {
    if bound.hi < value {
      return nil
    }
    if bound.lo < value {
      return R1Interval(lo: value, hi: bound.hi)
    }
  }
  return bound
}

// clipBoundAxis returns the clipped versions of the bounding intervals for the given
// axes for the line segment from (a0,a1) to (b0,b1) so that neither extends beyond the
// given clip interval. negSlope is a precomputed helper variable that indicates which
// diagonal of the bounding box is spanned by AB; it is false if AB has positive slope,
// and true if AB has negative slope. If the clipping interval doesn't overlap the bounds,
// false is returned.
func clipBoundAxis(a0: Double, b0: Double, bound0: R1Interval, a1: Double, b1: Double, bound1: R1Interval, negSlope: Bool, clip: R1Interval) -> (bound0c: R1Interval, bound1c: R1Interval, updated: Bool) {
  var bound0 = bound0
  var bound1 = bound1
  if bound0.lo < clip.lo {
    // If the upper bound is below the clips lower bound, there is nothing to do.
    if bound0.hi < clip.lo {
      return (bound0, bound1, false)
    }
    // narrow the intervals lower bound to the clip bound.
    bound0 = R1Interval(lo: clip.lo, hi: bound0.hi)
    let value = interpolateFloat64(x: clip.lo, a: a0, b: b0, a1: a1, b1: b1)
    bound1 = updateEndpoint(bound: bound1, highEndpoint: negSlope, value: value) ?? bound1
  }
  if bound0.hi > clip.hi {
    // If the lower bound is above the clips upper bound, there is nothing to do.
    if bound0.lo > clip.hi {
      return (bound0, bound1, false)
    }
    // narrow the intervals upper bound to the clip bound.
    bound0 = R1Interval(lo: bound0.lo, hi: clip.hi)
    let value = interpolateFloat64(x: clip.hi, a: a0, b: b0, a1: a1, b1: b1)
    bound1 = updateEndpoint(bound: bound1, highEndpoint: !negSlope, value: value) ?? bound1
  }
  return (bound0, bound1, true)
}

// edgeIntersectsRect reports whether the edge defined by AB intersects the
// given closed rectangle to within the error bound.
func edgeIntersectsRect(a: R2Point, b: R2Point, r: R2Rect) -> Bool {
  // First check whether the bounds of a Rect around AB intersects the given rect.
  if !r.intersects(R2Rect(p0: a, p1: b)) {
    return false
  }
  // Otherwise AB intersects the rect if and only if all four vertices of rect
  // do not lie on the same side of the extended line AB. We test this by finding
  // the two vertices of rect with minimum and maximum projections onto the normal
  // of AB, and computing their dot products with the edge normal.
  let n = b.sub(a).ortho()
  var i = 0
  if n.x >= 0 {
    i = 1
  }
  var j = 0
  if n.y >= 0 {
    j = 1
  }
  let max = n.dot(r.vertex(i: i, j: j).sub(a))
  let min = n.dot(r.vertex(i: 1 - i, j: 1 - j).sub(a))
  return (max >= 0) && (min <= 0)
}

// clippedEdgeBound returns the bounding rectangle of the portion of the edge defined
// by AB intersected by clip. The resulting bound may be empty. This is a convenience
// function built on top of clipEdgeBound.
func clippedEdgeBound(a: R2Point, b: R2Point, clip: R2Rect) -> R2Rect {
  let bound = R2Rect(p0: a, p1: b)
  let (b1, intersects) = clipEdgeBound(a: a, b: b, clip: clip, bound: bound)
  if intersects {
    return b1
  }
  return R2Rect.empty
}

// clipEdgeBound clips an edge AB to sequence of rectangles efficiently.
// It represents the clipped edges by their bounding boxes rather than as a pair of
// endpoints. Specifically, let A'B' be some portion of an edge AB, and let bound be
// a tight bound of A'B'. This function returns the bound that is a tight bound
// of A'B' intersected with a given rectangle. If A'B' does not intersect clip,
// it returns false and the original bound.
func clipEdgeBound(a: R2Point, b: R2Point, clip: R2Rect, bound: R2Rect) -> (R2Rect, Bool) {
  // negSlope indicates which diagonal of the bounding box is spanned by AB: it
  // is false if AB has positive slope, and true if AB has negative slope. This is
  // used to determine which interval endpoints need to be updated each time
  // the edge is clipped.
  let negSlope = (a.x > b.x) != (a.y > b.y)
  let (b0x, b0y, up1) = clipBoundAxis(a0: a.x, b0: b.x, bound0: bound.x, a1: a.y, b1: b.y, bound1: bound.y, negSlope: negSlope, clip: clip.x)
  if !up1 {
    return (bound, false)
  }
  let (b1y, b1x, up2) = clipBoundAxis(a0: a.y, b0: b.y, bound0: b0y, a1: a.x, b1: b.x, bound1: b0x, negSlope: negSlope, clip: clip.y)
  if !up2 {
    return (R2Rect(x: b0x, y: b0y), false)
  }
  return (R2Rect(x: b1x, y: b1y), true)
}

// interpolateFloat64 returns a value with the same combination of a1 and b1 as the
// given value x is of a and b. This function makes the following guarantees:
//  - If x == a, then x1 = a1 (exactly).
//  - If x == b, then x1 = b1 (exactly).
//  - If a <= x <= b, then a1 <= x1 <= b1 (even if a1 == b1).
// This requires a != b.
func interpolateFloat64(x: Double, a: Double, b: Double, a1: Double, b1: Double) -> Double {
  // To get results that are accurate near both A and B, we interpolate
  // starting from the closer of the two points.
  if abs(a - x) <= abs(b - x) {
    return a1 + (b1 - a1) * (x - a) / (b - a)
  }
  return b1 + (a1 - b1) * (x - b) / (a - b)
}

// FaceSegment represents an edge AB clipped to an S2 cube face. It is
// represented by a face index and a pair of (u,v) coordinates.
struct FaceSegment {
  let face: Int
  let a: R2Point
  let b: R2Point
}

// FaceSegments subdivides the given edge AB at every point where it crosses the
// boundary between two S2 cube faces and returns the corresponding FaceSegments.
// The segments are returned in order from A toward B. The input points must be
// unit length.
//
// This function guarantees that the returned segments form a continuous path
// from A to B, and that all vertices are within faceClipErrorUVDist of the
// line AB. All vertices lie within the [-1,1]x[-1,1] cube face rectangles.
// The results are consistent with Sign, i.e. the edge is well-defined even its
// endpoints are antipodal.
// TODO(roberts): Extend the implementation of PointCross so that this is true.
func faceSegments(a: S2Point, b: S2Point) -> [FaceSegment] {
  // Fast path: both endpoints are on the same face.
  let aFaceUV = S2Cube(point: a)
  let bFaceUV = S2Cube(point: b)
  let sa1 = R2Point(x: aFaceUV.u, y: aFaceUV.v)
  let sb1 = R2Point(x: bFaceUV.u, y: bFaceUV.v)
  if aFaceUV.face == bFaceUV.face {
    let segment = FaceSegment(face: aFaceUV.face, a: sa1, b: sb1)
    return [segment]
  }
  // Starting at A, we follow AB from face to face until we reach the face
  // containing B. The following code is designed to ensure that we always
  // reach B, even in the presence of numerical errors.
  //
  // First we compute the normal to the plane containing A and B. This normal
  // becomes the ultimate definition of the line AB; it is used to resolve all
  // questions regarding where exactly the line goes. Unfortunately due to
  // numerical errors, the line may not quite intersect the faces containing
  // the original endpoints. We handle this by moving A and/or B slightly if
  // necessary so that they are on faces intersected by the line AB.
  let ab = a.pointCross(b)
  let (aFace2, sa2) = moveOriginToValidFace(face: aFaceUV.face, a: a, ab: ab, aUV: sa1)
  let (bFace2, sb2) = moveOriginToValidFace(face: bFaceUV.face, a: b, ab: S2Point(raw: ab.mul(-1)), aUV: sb1)
  // Now we simply follow AB from face to face until we reach B.
  var segments: [FaceSegment] = []
  var sFace = aFace2
  var sa = sa2
  var sb = sb2
  let bSaved = sb2
  var face = aFace2
  while face != bFace2 {
    // Complete the current segment by finding the point where AB
    // exits the current face.
    let z = S2Cube.faceXYZtoUVW(face: face, point: ab)
    let n = PointUVW(p: z)
    let exitAxis = n.exitAxis()
    sb = n.exitPoint(axis: exitAxis)
    segments.append(FaceSegment(face: sFace, a: sa, b: sb))
    // Compute the next face intersected by AB, and translate the exit
    // point of the current segment into the (u,v) coordinates of the
    // next face. This becomes the first point of the next segment.
    let exitXyz = S2Cube.faceUVToXYZ(face: face, u: sb.x, v: sb.y)
    face = nextFace(face: face, exit: sb, axis: exitAxis, n: n, targetFace: bFace2)
    let exitUvw = S2Cube.faceXYZtoUVW(face: face, point: S2Point(raw: exitXyz))
    sFace = face
    sa = R2Point(x: exitUvw.x, y: exitUvw.y)
  }
  // Finish the last segment.
  segments.append(FaceSegment(face: sFace, a: sa, b: bSaved))
  return segments
}

// moveOriginToValidFace updates the origin point to a valid face if necessary.
// Given a line segment AB whose origin A has been projected onto a given cube
// face, determine whether it is necessary to project A onto a different face
// instead. This can happen because the normal of the line AB is not computed
// exactly, so that the line AB (defined as the set of points perpendicular to
// the normal) may not intersect the cube face containing A. Even if it does
// intersect the face, the exit point of the line from that face may be on
// the wrong side of A (i.e., in the direction away from B). If this happens,
// we reproject A onto the adjacent face where the line AB approaches A most
// closely. This moves the origin by a small amount, but never more than the
// error tolerances.
func moveOriginToValidFace(face: Int, a: S2Point, ab: S2Point, aUV: R2Point) -> (Int, R2Point) {
  // Fast path: if the origin is sufficiently far inside the face, it is
  // always safe to use it.
  let maxSafeUVCoord = 1 - faceClipErrorUVCoord
  if max(abs((aUV).x), abs((aUV).y)) <= maxSafeUVCoord {
    return (face, aUV)
  }
  
  // Otherwise check whether the normal AB even intersects this face.
  let z = S2Cube.faceXYZtoUVW(face: face, point: ab)
  let n = PointUVW(p: z)
  if n.intersectsFace() {
    // Check whether the point where the line AB exits this face is on the
    // wrong side of A (by more than the acceptable error tolerance).
    let uv = n.exitPoint(axis: n.exitAxis())
    let exit = S2Cube.faceUVToXYZ(face: face, u: uv.x, v: uv.y)
    let aTangent = ab.cross(a)
    // We can use the given face.
    if exit.sub(a.v).dot(aTangent) >= -faceClipErrorRadians {
      return (face, aUV)
    }
  }
  // Otherwise we reproject A to the nearest adjacent face. (If line AB does
  // not pass through a given face, it must pass through all adjacent faces.)
  var dir = 0
  var face2 = face
  if abs(aUV.x) >= abs(aUV.y) {
    // U-axis
    if aUV.x > 0 {
      dir = 1
    }
    face2 = S2Cube.uvwFace(face: face, axis: 0, direction: dir)
  } else {
    // V-axis
    if aUV.y > 0 {
      dir = 1
    }
    face2 = S2Cube.uvwFace(face: face, axis: 1, direction: dir)
  }
  let (aUVx, aUVy) = S2Cube.validFaceXYZToUV(face: face2, point: a)
  let x = max(-1.0, min(1.0, aUVx))
  let y = max(-1.0, min(1.0, aUVy))
  return (face2, R2Point(x: x, y: y))
}

// nextFace returns the next face that should be visited by FaceSegments, given that
// we have just visited face and we are following the line AB (represented
// by its normal N in the (u,v,w) coordinates of that face). The other
// arguments include the point where AB exits face, the corresponding
// exit axis, and the target face containing the destination point B.
func nextFace(face: Int, exit: R2Point, axis: UVWAxis, n: PointUVW, targetFace: Int) -> Int {
  // this bit is to work around C++ cleverly casting bools to ints for you.
  var exitA = exit.x
  var exit1MinusA = exit.y
  if axis == .v {
    exitA = exit.y
    exit1MinusA = exit.x
  }
  var exitAPos = 0
  if exitA > 0 {
    exitAPos = 1
  }
  var exit1MinusAPos = 0
  if exit1MinusA > 0 {
    exit1MinusAPos = 1
  }
  // We return the face that is adjacent to the exit point along the given
  // axis. If line AB exits *exactly* through a corner of the face, there are
  // two possible next faces. If one is the target face containing B, then
  // we guarantee that we advance to that face directly.
  //
  // The three conditions below check that (1) AB exits approximately through
  // a corner, (2) the adjacent face along the non-exit axis is the target
  // face, and (3) AB exits *exactly* through the corner. (The sumEqual
  // code checks whether the dot product of (u,v,1) and n is exactly zero.)
  if abs(exit1MinusA) == 1 {
    if S2Cube.uvwFace(face: face, axis: 1 - axis.rawValue, direction: exit1MinusAPos) == targetFace {
      if sumEqual(u: exit.x * n.p.x, v: exit.y * n.p.y, w: -n.p.z) {
        return targetFace
      }
    }
  }
  // Otherwise return the face that is adjacent to the exit point in the
  // direction of the exit axis.
  return S2Cube.uvwFace(face: face, axis: axis.rawValue, direction: exitAPos)
}
