//
//  EdgeCrosser.swift
//  s2-swiftPackageDescription
//
//  Created by Axel Huesemann on 1/11/18.
//

import Foundation


// A Crossing indicates how edges cross.
enum Crossing: Int {
  // Cross means the edges cross.
  case cross = 0
  // MaybeCross means two vertices from different edges are the same.
  case maybeCross = 1
  // DoNotCross means the edges do not cross.
  case doNotCross = 2
}

/// Allows edges to be efficiently tested for intersection with a
/// given fixed edge AB. It is especially efficient when testing for
/// intersection with an edge chain connecting vertices v0, v1, v2, ...
///
/// Example usage:
///  func countIntersections(a: S2Point, b: S2Point, edges: [Edge]) -> Int {
///    var count = 0
///    var crosser = EdgeCrosser(a, b)
///    for edge in edges {
///      if crosser.crossingSign(edge.v0, edge.v1) != .doNotCross {
///        count += 1
///      }
///    }
///    return count
///  }
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
  var c: S2Point?     // Previous vertex in the vertex chain.
  var acb: Direction? // The orientation of triangle ACB.

  // Contructs an EdgeCrosser with the fixed edge AB.
  init(a: S2Point, b: S2Point) {
    self.a = a
    self.b = b
    aXb = S2Point(raw: a.cross(b))
    let norm = a.pointCross(b)
    aTangent = S2Point(raw: a.cross(norm))
    bTangent = S2Point(raw: norm.cross(b))
    c = nil
    acb = nil
  }

  /// Reports whether the edge AB intersects the edge CD. If any two
  /// vertices from different edges are the same, returns MaybeCross. If either edge
  /// is degenerate (A == B or C == D), returns either DoNotCross or MaybeCross.
  /// Properties of CrossingSign:
  ///  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
  ///  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
  ///  (3) CrossingSign(a,b,c,d) == MaybeCross if a==c, a==d, b==c, b==d
  ///  (3) CrossingSign(a,b,c,d) == DoNotCross or MaybeCross if a==b or c==d
  /// Note that if you want to check an edge against a chain of other edges,
  /// it is slightly more efficient to use the single-argument version
  /// ChainCrossingSign below.
  mutating func crossingSign(c: S2Point, d: S2Point) -> Crossing {
    if c != self.c {
      restart(at: c)
    }
    return chainCrossingSign(d: d)
  }
  
  /// Reports whether if CrossingSign(c, d) > 0, or AB and
  /// CD share a vertex and VertexCrossing(a, b, c, d) is true.
  //
  /// This method extends the concept of a "crossing" to the case where AB
  /// and CD have a vertex in common. The two edges may or may not cross,
  /// according to the rules defined in VertexCrossing above. The rules
  /// are designed so that point containment tests can be implemented simply
  /// by counting edge crossings. Similarly, determining whether one edge
  /// chain crosses another edge chain can be implemented by counting.
  mutating func isEdgeOrVertexCrossing(c: S2Point, d: S2Point) -> Bool {
    if c != self.c {
      restart(at: c)
    }
    return isEdgeOrVertexChainCrossing(d: d)
  }
  
  /// A convenience constructor that uses AB as the fixed edge,
  /// and C as the first vertex of the vertex chain (equivalent to calling restart(at: c)).
  /// rerturn a ChainEdgeCrosser which is an EdgeCrosser
  /// You don't need to use this or any of the chain functions unless you're trying to
  /// squeeze out every last drop of performance. Essentially all you are saving is a test
  /// whether the first vertex of the current edge is the same as the second vertex of the
  /// previous edge.
  init(a: S2Point, b: S2Point, c: S2Point) {
    self.init(a: a, b: b)
    restart(at: c)
  }
  
  /// Sets the current point of the edge crosser to be c.
  /// Call this method when your chain 'jumps' to a new place.
  /// The argument must point to a value that persists until the next call.
  mutating func restart(at c: S2Point) {
    self.c = c
    acb = -S2Point.triageSign(a, b, c)
  }
  
  /// Like crossingSign, but uses the last vertex passed to one of
  /// the crossing methods (or RestartAt) as the first vertex of the current edge.
  mutating func chainCrossingSign(d: S2Point) -> Crossing {
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
    return crossingSign(d: d, bda: bda)
  }
  
  /// EdgeOrVertexChainCrossing is like EdgeOrVertexCrossing, but uses the last vertex
  /// passed to one of the crossing methods (or RestartAt) as the first vertex of the current edge.
  mutating func isEdgeOrVertexChainCrossing(d: S2Point) -> Bool {
    guard let c = c else { return false }
    // We need to copy e.c since it is clobbered by ChainCrossingSign.???
    // c := e.c
    switch chainCrossingSign(d: d) {
    case .doNotCross:
      return false
    case .cross:
      return true
    default:
      return EdgeCrosser.vertexCrossing(a: a, b: b, c: c, d: d)
    }
  }
  
  /// Handles the slow path of CrossingSign.
  mutating func crossingSign(d: S2Point, bda: Direction) -> Crossing {
    guard let c = c else { return .maybeCross }
    // Compute the actual result, and then save the current vertex D as the next
    // vertex C, and save the orientation of the next triangle ACB (which is
    // opposite to the current triangle BDA).
    defer {
      self.c = d
      acb = -bda
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
    let maxError = (1.5 + 1 / sqrt(3)) * S1Interval.dblEpsilon
    if c.dot(aTangent) > maxError && d.dot(aTangent) > maxError {
      return .doNotCross
    }
    if c.dot(bTangent) > maxError && d.dot(bTangent) > maxError {
      return .doNotCross
    }
    // Otherwise, eliminate the cases where two vertices from different edges are
    // equal. (These cases could be handled in the code below, but we would rather
    // avoid calling ExpensiveSign if possible.)
    if a == c || a == d || b == c || b == d {
      return .maybeCross
    }
    // Eliminate the cases where an input edge is degenerate. (Note that in
    // most cases, if CD is degenerate then this method is not even called
    // because acb and bda have different signs.)
    if a == b || c == d {
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
    if dac != acb {
      return .doNotCross
    }
    return .cross
  }
  
  // intersectionError can be set somewhat arbitrarily, because the algorithm
  // uses more precision if necessary in order to achieve the specified error.
  // The only strict requirement is that intersectionError >= dblEpsilon
  // radians. However, using a larger error tolerance makes the algorithm more
  // efficient because it reduces the number of cases where exact arithmetic is
  // needed.
  static let intersectionError = S1Angle(8 * Cell.dblEpsilon)

  // intersectionMergeRadius is used to ensure that intersection points that
  // are supposed to be coincident are merged back together into a single
  // vertex. This is required in order for various polygon operations (union,
  // intersection, etc) to work correctly. It is twice the intersection error
  // because two coincident intersection points might have errors in
  // opposite directions.
  static let intersectionMergeRadius = 2 * intersectionError


  // CrossingSign reports whether the edge AB intersects the edge CD.
  // If AB crosses CD at a point that is interior to both edges, Cross is returned.
  // If any two vertices from different edges are the same it returns MaybeCross.
  // Otherwise it returns DoNotCross.
  // If either edge is degenerate (A == B or C == D), the return value is MaybeCross
  // if two vertices from different edges are the same and DoNotCross otherwise.
  //
  // Properties of CrossingSign:
  //
  //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
  //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
  //  (3) CrossingSign(a,b,c,d) == MaybeCross if a==c, a==d, b==c, b==d
  //  (3) CrossingSign(a,b,c,d) == DoNotCross or MaybeCross if a==b or c==d
  //
  // This method implements an exact, consistent perturbation model such
  // that no three points are ever considered to be collinear. This means
  // that even if you have 4 points A, B, C, D that lie exactly in a line
  // (say, around the equator), C and D will be treated as being slightly to
  // one side or the other of AB. This is done in a way such that the
  // results are always consistent (see RobustSign).
  static func crossingSign(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Crossing {
    var crosser = EdgeCrosser(a: a, b: b, c: c)
    return crosser.chainCrossingSign(d: d)
  }

  // SimpleCrossing reports whether edge AB crosses CD at a point that is interior
  // to both edges. Properties:
  //
  //  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d)
  //  (2) SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
  static func simpleCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
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
  
  // Reports whether two edges "cross" in such a way that point-in-polygon
  // containment tests can be implemented by counting the number of edge crossings.
  //
  // Given two edges AB and CD where at least two vertices are identical
  // (i.e. CrossingSign(a,b,c,d) == 0), the basic rule is that a "crossing"
  // occurs if AB is encountered after CD during a CCW sweep around the shared
  // vertex starting from a fixed reference point.
  //
  // Note that according to this rule, if AB crosses CD then in general CD
  // does not cross AB. However, this leads to the correct result when
  // counting polygon edge crossings. For example, suppose that A,B,C are
  // three consecutive vertices of a CCW polygon. If we now consider the edge
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
  static func vertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
    // If A == B or C == D there is no intersection. We need to check this
    // case first in case 3 or more input points are identical.
    if a == b || c == d {
      return false
    }
    // If any other pair of vertices is equal, there is a crossing if and only
    // if OrderedCCW indicates that the edge AB is further CCW around the
    // shared vertex O (either A or B) than the edge CD, starting from an
    // arbitrary fixed reference point.
    if a == d {
      return S2Point.orderedCCW(S2Point(raw: a.ortho()), c, b, a)
    } else if b == c {
      return S2Point.orderedCCW(S2Point(raw: b.ortho()), d, a, b)
    } else if a == c {
      return S2Point.orderedCCW(S2Point(raw: a.ortho()), d, b, a)
    } else if b == d {
      return S2Point.orderedCCW(S2Point(raw: b.ortho()), c, a, b)
    }
    return false
  }

  // EdgeOrVertexCrossing is a convenience function that calls CrossingSign to
  // handle cases where all four vertices are distinct, and VertexCrossing to
  // handle cases where two or more vertices are the same. This defines a crossing
  // function such that point-in-polygon containment tests can be implemented
  // by simply counting edge crossings.
  static func edgeOrVertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
    switch crossingSign(a: a, b: b, c: c, d: d) {
    case .doNotCross:
      return false
    case .cross:
      return true
    default:
      return
        vertexCrossing(a: a, b: b, c: c, d: d)
    }
  }

  // Intersection returns the intersection point of two edges AB and CD that cross
  // (CrossingSign(a,b,c,d) == Crossing).
  //
  // Useful properties of Intersection:
  //
  //  (1) Intersection(b,a,c,d) == Intersection(a,b,d,c) == Intersection(a,b,c,d)
  //  (2) Intersection(c,d,a,b) == Intersection(a,b,c,d)
  //
  // The returned intersection point X is guaranteed to be very close to the
  // true intersection point of AB and CD, even if the edges intersect at a
  // very small angle.
  static func intersection(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> S2Point {
    // It is difficult to compute the intersection point of two edges accurately
    // when the angle between the edges is very small. Previously we handled
    // this by only guaranteeing that the returned intersection point is within
    // intersectionError of each edge. However, this means that when the edges
    // cross at a very small angle, the computed result may be very far from the
    // true intersection point.
    //
    // Instead this function now guarantees that the result is always within
    // intersectionError of the true intersection. This requires using more
    // sophisticated techniques and in some cases extended precision.
    //
    //  - intersectionStable computes the intersection point using
    //    projection and interpolation, taking care to minimize cancellation
    //    error.
    //
    //  - intersectionExact computes the intersection point using precision
    //    arithmetic and converts the final result back to an Point.
    let pt = intersectionStable(a0: a0, a1: a1, b0: b0, b1: b1) ?? intersectionExact(a0: a0, a1: a1, b0: b0, b1: b1)
    // Make sure the intersection point is on the correct side of the sphere.
    // Since all vertices are unit length, and edges are less than 180 degrees,
    // (a0 + a1) and (b0 + b1) both have positive dot product with the
    // intersection point.  We use the sum of all vertices to make sure that the
    // result is unchanged when the edges are swapped or reversed.
    let sumA = a0 + a1
    let sumB = b0 + b1
    if pt.v.dot(sumA + sumB) < 0 {
      return pt.inverse()
    }
    return pt
  }

  // Computes the cross product of two vectors, normalized to be unit length.
  // Also returns the length of the cross
  // product before normalization, which is useful for estimating the amount of
  // error in the result.  For numerical stability, the vectors should both be
  // approximately unit length.
  static func robustNormalWithLength(x: R3Vector, y: R3Vector) -> (R3Vector, Double) {
    // This computes 2 * (x.Cross(y)), but has much better numerical
    // stability when x and y are unit length.
    let tmp = (x - y).cross(x + y)
    let length = tmp.norm
    guard length != 0 else {
      return (R3Vector(x: 0, y: 0, z: 0), 0) // Since tmp == 2 * (x.Cross(y))
    }
    return (tmp * (1 / length), length * 0.5)
  }

  /*
   // intersectionSimple is not used by the C++ so it is skipped here.
   */

  // projection returns the projection of aNorm onto X (x.Dot(aNorm)), and a bound
  // on the error in the result. aNorm is not necessarily unit length.
  //
  // The remaining parameters (the length of aNorm (aNormLen) and the edge endpoints
  // a0 and a1) allow this dot product to be computed more accurately and efficiently.
  static func projection(x: R3Vector, aNorm: R3Vector, aNormLen: Double, a0: S2Point, a1: S2Point) -> (proj: Double, bound: Double) {
    // The error in the dot product is proportional to the lengths of the input
    // vectors, so rather than using x itself (a unit-length vector) we use
    // the vectors from x to the closer of the two edge endpoints. This
    // typically reduces the error by a huge factor.
    let x0 = x - a0.v
    let x1 = x - a1.v
    let x0Dist2 = x0.norm2
    let x1Dist2 = x1.norm2
    // If both distances are the same, we need to be careful to choose one
    // endpoint deterministically so that the result does not change if the
    // order of the endpoints is reversed.
    var dist: Double
    var proj: Double
    if x0Dist2 < x1Dist2 || (x0Dist2 == x1Dist2 && x0 < x1) {
      dist = sqrt(x0Dist2)
      proj = x0.dot(aNorm)
    } else {
      dist = sqrt(x1Dist2)
      proj = x1.dot(aNorm)
    }
    // This calculation bounds the error from all sources: the computation of
    // the normal, the subtraction of one endpoint, and the dot product itself.
    // dblEpsilon appears because the input points are assumed to be
    // normalized in double precision.
    //
    // For reference, the bounds that went into this calculation are:
    // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * dblEpsilon) * epsilon
    // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * epsilon
    // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * epsilon
    let bound = (((3.5 + 2 * sqrt(3)) * aNormLen + 32 * sqrt(3) * Cell.dblEpsilon) * dist + 1.5 * abs(proj)) * R1Interval.epsilon
    return (proj, bound)
  }

  // compareEdges reports whether (a0,a1) is less than (b0,b1) with respect to a total
  // ordering on edges that is invariant under edge reversals.
  static func compareEdges(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> Bool {
    let (a0, a1) = (a0 >= a1) ? (a1, a0) : (a0, a1)
    let (b0, b1) = (b0 >= b1) ? (b1, b0) : (b0, b1)
    return a0 < b0 || (a0 == b0 && a1 < b1)
  }

  // intersectionStable returns the intersection point of the edges (a0,a1) and
  // (b0,b1) if it can be computed to within an error of at most intersectionError
  // by this function.
  //
  // The intersection point is not guaranteed to have the correct sign because we
  // choose to use the longest of the two edges first. The sign is corrected by
  // Intersection.
  static func intersectionStable(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> S2Point? {
    // Sort the two edges so that (a0,a1) is longer, breaking ties in a
    // deterministic way that does not depend on the ordering of the endpoints.
    // This is desirable for two reasons:
    //  - So that the result doesn't change when edges are swapped or reversed.
    //  - It reduces error, since the first edge is used to compute the edge
    //    normal (where a longer edge means less error), and the second edge
    //    is used for interpolation (where a shorter edge means less error).
    let aLen2 = a1.sub(a0).norm2
    let bLen2 = b1.sub(b0).norm2
    if aLen2 < bLen2 || (aLen2 == bLen2 && compareEdges(a0: a0, a1: a1, b0: b0, b1: b1)) {
      return intersectionStableSorted(a0: b0, a1: b1, b0: a0, b1: a1)
    }
    return intersectionStableSorted(a0: a0, a1: a1, b0: b0, b1: b1)
  }

  // intersectionStableSorted is a helper function for intersectionStable.
  // It expects that the edges (a0,a1) and (b0,b1) have been sorted so that
  // the first edge passed in is longer.
  static func intersectionStableSorted(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> S2Point? {
    // Compute the normal of the plane through (a0, a1) in a stable way.
    let aNorm = a0.sub(a1).cross(a0.add(a1))
    let aNormLen = aNorm.norm
    let bLen = b1.sub(b0).norm
    // Compute the projection (i.e., signed distance) of b0 and b1 onto the
    // plane through (a0, a1).  Distances are scaled by the length of aNorm.
    let (b0Dist, b0Error) = projection(x: b0.v, aNorm: aNorm, aNormLen: aNormLen, a0: a0, a1: a1)
    let (b1Dist, b1Error) = projection(x: b1.v, aNorm: aNorm, aNormLen: aNormLen, a0: a0, a1: a1)
    // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is
    // |b0Dist - b1Dist|.  Note that b0Dist and b1Dist generally have
    // opposite signs because b0 and b1 are on opposite sides of (a0, a1).  The
    // code below finds the intersection point by interpolating along the edge
    // (b0, b1) to a fractional distance of b0Dist / (b0Dist - b1Dist).
    //
    // It can be shown that the maximum error in the interpolation fraction is
    //
    //   (b0Dist * b1Error - b1Dist * b0Error) / (distSum * (distSum - errorSum))
    //
    // We save ourselves some work by scaling the result and the error bound by
    // "distSum", since the result is normalized to be unit length anyway.
    let distSum = abs(b0Dist - b1Dist)
    let errorSum = b0Error + b1Error
    if distSum <= errorSum {
      return nil // Error is unbounded in this case.
    }
    let x = b1.mul(b0Dist).sub(b0.mul(b1Dist))
    let err = bLen * abs(b0Dist * b1Error - b1Dist * b0Error) / (distSum - errorSum) + 2 * distSum * R1Interval.epsilon
    // Finally we normalize the result, compute the corresponding error, and
    // check whether the total error is acceptable.
    let xLen = x.norm
    let maxError = intersectionError
    if err > (Double(maxError) - R1Interval.epsilon) * xLen {
      return nil
    }
    return S2Point(raw: x.mul(1 / xLen))
  }

  // intersectionExact returns the intersection point of (a0, a1) and (b0, b1)
  // using precise arithmetic. Note that the result is not exact because it is
  // rounded down to double precision at the end. Also, the intersection point
  // is not guaranteed to have the correct sign (i.e., the return value may need
  // to be negated).
  static func intersectionExact(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point) -> S2Point {
    // Since we are using presice arithmetic, we don't need to worry about
    // numerical stability.
    let a0P = PreciseVector(v: a0.v)
    let a1P = PreciseVector(v: a1.v)
    let b0P = PreciseVector(v: b0.v)
    let b1P = PreciseVector(v: b1.v)
    let aNormP = a0P.cross(a1P)
    let bNormP = b0P.cross(b1P)
    let xP = aNormP.cross(bNormP)
    // The final Normalize() call is done in double precision, which creates a
    // directional error of up to 2*dblEpsilon. (Precise conversion and Normalize()
    // each contribute up to dblEpsilon of directional error.)
    var x = xP.vector()
    if x == R3Vector(x: 0, y: 0, z: 0) {
      // The two edges are exactly collinear, but we still consider them to be
      // "crossing" because of simulation of simplicity. Out of the four
      // endpoints, exactly two lie in the interior of the other edge. Of
      // those two we return the one that is lexicographically smallest.
      x = R3Vector(x: 10, y: 10, z: 10) // Greater than any valid S2Point
      let aNorm = S2Point(raw: aNormP.vector())
      let bNorm = S2Point(raw: bNormP.vector())
      if S2Point.orderedCCW(b0, a0, b1, bNorm) && a0.v < x {
        return a0
      }
      if S2Point.orderedCCW(b0, a1, b1, bNorm) && a1.v < x {
        return a1
      }
      if S2Point.orderedCCW(a0, b0, a1, aNorm) && b0.v < x {
        return b0
      }
      if S2Point.orderedCCW(a0, b1, a1, aNorm) && b1.v < x {
        return b1
      }
    }
    return S2Point(raw: x.normalized())
  }

}
