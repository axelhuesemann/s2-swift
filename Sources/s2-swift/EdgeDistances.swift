//
//  EdgeDistances.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/19/18.
//

import Foundation

/// DistanceFromSegment returns the distance of point X from line segment AB.
/// The points are expected to be normalized. The result is very accurate for small
/// distances but may have some numerical error if the distance is large
/// (approximately pi/2 or greater). The case A == B is handled correctly.
func distanceFromSegment(x: S2Point, a: S2Point, b: S2Point) -> S1Angle {
  var minDist = S1ChordAngle.zero
  (minDist, _) = updateMinDistance(x: x, a: a, b: b, minDist: minDist, alwaysUpdate: true)
  return minDist.angle()
}

/// IsDistanceLess reports whether the distance from X to the edge AB is less
/// than limit. This method is faster than DistanceFromSegment(). If you want to
/// compare against a fixed s1.Angle, you should convert it to an s1.ChordAngle
/// once and save the value, since this conversion is relatively expensive.
func isDistanceLess(x: S2Point, a: S2Point, b: S2Point, limit: S1ChordAngle) -> Bool {
  let (_, less) = updateMinDistance(x: x, a: a, b: b, minDist: limit)
  return less
}

/// UpdateMinDistance checks if the distance from X to the edge AB is less
/// then minDist, and if so, returns the updated value and true.
/// The case A == B is handled correctly.
//
/// Use this method when you want to compute many distances and keep track of
/// the minimum. It is significantly faster than using DistanceFromSegment
/// because (1) using s1.ChordAngle is much faster than s1.Angle, and (2) it
/// can save a lot of work by not actually computing the distance when it is
/// obviously larger than the current minimum.
func updateMinDistance(x: S2Point, a: S2Point, b: S2Point, minDist: S1ChordAngle) -> (S1ChordAngle, Bool) {
  return updateMinDistance(x: x, a: a, b: b, minDist: minDist, alwaysUpdate: false)
}

/// IsInteriorDistanceLess reports whether the minimum distance from X to the
/// edge AB is attained at an interior point of AB (i.e., not an endpoint), and
/// that distance is less than limit.
func isInteriorDistanceLess(x: S2Point, a: S2Point, b: S2Point, limit: S1ChordAngle) -> Bool {
  let (_, less) = updateMinInteriorDistance(x: x, a: a, b: b, minDist: limit)
  return less
}

/// UpdateMinInteriorDistance reports whether the minimum distance from X to AB
/// is attained at an interior point of AB (i.e., not an endpoint), and that distance
/// is less than minDist. If so, the value of minDist is updated and true is returned.
/// Otherwise it is unchanged and returns false.
func updateMinInteriorDistance(x: S2Point, a: S2Point, b: S2Point, minDist: S1ChordAngle) -> (S1ChordAngle, Bool) {
  return interiorDist(x: x, a: a, b: b, minDist: minDist, alwaysUpdate: false)
}

/// Project returns the point along the edge AB that is closest to the point X.
/// The fractional distance of this point along the edge AB can be obtained
/// using DistanceFraction.
//
/// This requires that all points are unit length.
func project(x: S2Point, a: S2Point, b: S2Point) -> S2Point {
  let aXb = a.pointCross(b)
  // Find the closest point to X along the great circle through AB.
  let p = x.v.sub(aXb.mul(x.dot(aXb) / aXb.v.norm2))
  // If this point is on the edge AB, then it's the closest point.
  if S2Point.sign(aXb, b: a, c: S2Point(raw: p)) && S2Point.sign(S2Point(raw: p), b: b, c: aXb) {
    return S2Point(raw: p)
  }
  // Otherwise, the closest point is either A or B.
  if x.sub(a).norm2 <= x.sub(b).norm2 {
    return a
  }
  return b
}

/// DistanceFraction returns the distance ratio of the point X along an edge AB.
/// If X is on the line segment AB, this is the fraction T such
/// that X == Interpolate(T, A, B).
//
/// This requires that A and B are distinct.
func distanceFraction(x: S2Point, a: S2Point, b: S2Point) -> Double {
  let d0 = x.angle(a)
  let d1 = x.angle(b)
  return d0 / (d0 + d1)
}

/// Returns the point X along the line segment AB whose distance from A
/// is the given fraction "t" of the distance AB. Does NOT require that "t" be
/// between 0 and 1. Note that all distances are measured on the surface of
/// the sphere, so this is more complicated than just computing (1-t)*a + t*b
/// and normalizing the result.
func interpolate(t: Double, a: S2Point, b: S2Point) -> S2Point {
  if t == 0.0 {
    return a
  }
  if t == 1.0 {
    return b
  }
  let ab = a.angle(b)
  return interpolateAtDistance(ax: t * ab, a: a, b: b)
}

/// Returns the point X along the line segment AB whose
/// distance from A is the angle ax.
func interpolateAtDistance(ax: S1Angle, a: S2Point, b: S2Point) -> S2Point {
  let aRad = ax
  // Use PointCross to compute the tangent vector at A towards B. The
  // result is always perpendicular to A, even if A=B or A=-B, but it is not
  // necessarily unit length. (We effectively normalize it below.)
  let normal = a.pointCross(b)
  let tangent = normal.cross(a)
  // Now compute the appropriate linear combination of A and "tangent". With
  // infinite precision the result would always be unit length, but we
  // normalize it anyway to ensure that the error is within acceptable bounds.
  // (Otherwise errors can build up when the result of one interpolation is
  // fed into another interpolation.)
  return S2Point(raw: a.mul(cos(aRad)).add(tangent.mul(sin(aRad) / tangent.norm)))
}

/// minUpdateDistanceMaxError returns the maximum error in the result of
/// UpdateMinDistance (and the associated functions such as
/// UpdateMinInteriorDistance, IsDistanceLess, etc), assuming that all
/// input points are normalized to within the bounds guaranteed by r3.Vector's
/// Normalize. The error can be added or subtracted from an s1.ChordAngle
/// using its Expanded method.
func minUpdateDistanceMaxError(dist: S1ChordAngle) -> Double {
  // There are two cases for the maximum error in UpdateMinDistance(),
  // depending on whether the closest point is interior to the edge.
  return max(minUpdateInteriorDistanceMaxError(dist: dist), dist.maxPointError())
}

/// minUpdateInteriorDistanceMaxError returns the maximum error in the result of
/// UpdateMinInteriorDistance, assuming that all input points are normalized
/// to within the bounds guaranteed by Point's Normalize. The error can be added
/// or subtracted from an s1.ChordAngle using its Expanded method.
func minUpdateInteriorDistanceMaxError(dist: S1ChordAngle) -> Double {
  // This bound includes all source of error, assuming that the input points
  // are normalized. a and b are components of chord length that are
  // perpendicular and parallel to a plane containing the edge respectively.
  let b = 0.5 * dist.value * dist.value
  let a = dist.value * sqrt(1 - 0.5 * b)
  let d0 = (2.5 + 2 * sqrt(3) + 8.5 * a) * a
  let d1 = (2 + 2 * sqrt(3) / 3 + 6.5 * (1 - b)) * b
  let d2 = (23 + 16 / sqrt(3)) * S1Interval.dblEpsilon
  return (d0 + d1 + d2) * S1Interval.dblEpsilon
}

/// updateMinDistance computes the distance from a point X to a line segment AB,
/// and if either the distance was less than the given minDist, or alwaysUpdate is
/// true, the value and whether it was updated are returned.
func updateMinDistance(x: S2Point, a: S2Point, b: S2Point, minDist: S1ChordAngle, alwaysUpdate: Bool) -> (S1ChordAngle, Bool) {
  let (d, ok) = interiorDist(x: x, a: a, b: b, minDist: minDist, alwaysUpdate: alwaysUpdate)
  if ok {
    // Minimum distance is attained along the edge interior.
    return (d, true)
  }
  // Otherwise the minimum distance is to one of the endpoints.
  let (xa2, xb2) = (x.sub(a).norm2, x.sub(b).norm2)
  let dist = S1ChordAngle(value: min(xa2, xb2))
  if !alwaysUpdate && dist >= minDist {
    return (minDist, false)
  }
  return (dist, true)
}

/// interiorDist returns the shortest distance from point x to edge ab, assuming
/// that the closest point to X is interior to AB. If the closest point is not
/// interior to AB, interiorDist returns (minDist, false). If alwaysUpdate is set to
/// false, the distance is only updated when the value exceeds certain the given minDist.
func interiorDist(x: S2Point, a: S2Point, b: S2Point, minDist: S1ChordAngle, alwaysUpdate: Bool) -> (S1ChordAngle, Bool) {
  // Chord distance of x to both end points a and b.
  let (xa2, xb2) = (x.sub(a).norm2, x.sub(b).norm2)
  // The closest point on AB could either be one of the two vertices (the
  // vertex case) or in the interior (the interior case). Let C = A x B.
  // If X is in the spherical wedge extending from A to B around the axis
  // through C, then we are in the interior case. Otherwise we are in the
  // vertex case.
  //
  // Check whether we might be in the interior case. For this to be true, XAB
  // and XBA must both be acute angles. Checking this condition exactly is
  // expensive, so instead we consider the planar triangle ABX (which passes
  // through the sphere's interior). The planar angles XAB and XBA are always
  // less than the corresponding spherical angles, so if we are in the
  // interior case then both of these angles must be acute.
  //
  // We check this by computing the squared edge lengths of the planar
  // triangle ABX, and testing acuteness using the law of cosines:
  //
  //   max(XA^2, XB^2) < min(XA^2, XB^2) + AB^2
  if max(xa2, xb2) >= min(xa2, xb2) + (a.sub(b)).norm2 {
    return (minDist, false)
  }
  // The minimum distance might be to a point on the edge interior. Let R
  // be closest point to X that lies on the great circle through AB. Rather
  // than computing the geodesic distance along the surface of the sphere,
  // instead we compute the "chord length" through the sphere's interior.
  //
  // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
  // is the point X projected onto the plane through the great circle AB.
  // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
  // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
  // is faster and the corresponding distance on the Earth's surface is
  // accurate to within 1% for distances up to about 1800km.
  let c = a.pointCross(b)
  let c2 = c.norm2
  let xDotC = x.dot(c)
  let xDotC2 = xDotC * xDotC
  if !alwaysUpdate && xDotC2 >= c2 * minDist.value {
    // The closest point on the great circle AB is too far away.
    return (minDist, false)
  }
  // Otherwise we do the exact, more expensive test for the interior case.
  // This test is very likely to succeed because of the conservative planar
  // test we did initially.
  let cx = c.cross(x)
  if a.v.dot(cx) >= 0 || b.v.dot(cx) <= 0 {
    return (minDist, false)
  }
  // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
  // This calculation has good accuracy for all chord lengths since it
  // is based on both the dot product and cross product (rather than
  // deriving one from the other). However, note that the chord length
  // representation itself loses accuracy as the angle approaches π.
  let qr = 1 - sqrt(cx.norm2 / c2)
  let dist = S1ChordAngle(value: (xDotC2 / c2) + (qr * qr))
  if !alwaysUpdate && dist >= minDist {
    return (minDist, false)
  }
  return (dist, true)
}

func updateEdgePairMinDistance(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point, minDist: S1ChordAngle) -> (S1ChordAngle, Bool) {
  if minDist == .zero {
    return (.zero, false)
  }
  var minDist = minDist
  if EdgeCrosser.crossingSign(a: a0, b: a1, c: b0, d: b1) == .cross {
    minDist = .zero
    return (minDist, true)
  }
  // Otherwise, the minimum distance is achieved at an endpoint of at least
  // one of the two edges. We ensure that all four possibilities are always checked.
  //
  // The calculation below computes each of the six vertex-vertex distances
  // twice (this could be optimized).
  var ok1: Bool
  var ok2: Bool
  var ok3: Bool
  var ok4: Bool
  (minDist, ok1) = updateMinDistance(x: a0, a: b0, b: b1, minDist: minDist)
  (minDist, ok2) = updateMinDistance(x: a1, a: b0, b: b1, minDist: minDist)
  (minDist, ok3) = updateMinDistance(x: b0, a: a0, b: a1, minDist: minDist)
  (minDist, ok4) = updateMinDistance(x: b1, a: a0, b: a1, minDist: minDist)
  return (minDist, ok1 || ok2 || ok3 || ok4)
}

/// EdgePairClosestPoints returns the pair of points (a, b) that achieves the
/// minimum distance between edges a0a1 and b0b1, where a is a point on a0a1 and
/// b is a point on b0b1. If the two edges intersect, a and b are both equal to
/// the intersection point. Handles a0 == a1 and b0 == b1 correctly.
func edgePairClosestPoints(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> (S2Point, S2Point) {
  if EdgeCrosser.crossingSign(a: a0, b: a1, c: b0, d: b1) == .cross {
    let x = EdgeCrosser.intersection(a0: a0, a1: a1, b0: b0, b1: b1)
    return (x, x)
  }
  // We save some work by first determining which vertex/edge pair achieves
  // the minimum distance, and then computing the closest point on that edge.
  var minDist = S1ChordAngle.zero
  var ok: Bool
  (minDist, ok) = updateMinDistance(x: a0, a: b0, b: b1, minDist: minDist, alwaysUpdate: true)
  var closestVertex = 0
  (minDist, ok) = updateMinDistance(x: a1, a: b0, b: b1, minDist: minDist)
  if ok {
    closestVertex = 1
  }
  (minDist, ok) = updateMinDistance(x: b0, a: a0, b: a1, minDist: minDist)
  if ok {
    closestVertex = 2
  }
  (minDist, ok) = updateMinDistance(x: b1, a: a0, b: a1, minDist: minDist)
  if ok {
    closestVertex = 3
  }
  switch closestVertex {
  case 0:
    return (a0, project(x: a0, a: b0, b: b1))
  case 1:
    return (a1, project(x: a1, a: b0, b: b1))
  case 2:
    return (project(x: b0, a: a0, b: a1), b0)
  case 3:
    return (project(x: b1, a: a0, b: a1), b1)
  default:
    fatalError("illegal case reached")
  }
  
}
