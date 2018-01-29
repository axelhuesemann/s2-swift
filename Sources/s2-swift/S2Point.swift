//
//  S2Point.swift
//  s2-swift
//

import Foundation

// maxDeterminantError is the maximum error in computing (AxB).C where all vectors
// are unit length. Using standard inequalities, it can be shown that
//
//  fl(AxB) = AxB + D where |D| <= (|AxB| + (2/sqrt(3))*|A|*|B|) * e
//
// where "fl()" denotes a calculation done in floating-point arithmetic,
// |x| denotes either absolute value or the L2-norm as appropriate, and
// e is a reasonably small value near the noise level of floating point
// number accuracy. Similarly,
//
//  fl(B.C) = B.C + d where |d| <= (|B.C| + 2*|B|*|C|) * e .
//
// Applying these bounds to the unit-length vectors A,B,C and neglecting
// relative error (which does not affect the sign of the result), we get
//
//  fl((AxB).C) = (AxB).C + d where |d| <= (3 + 2/sqrt(3)) * e
let maxDeterminantError = 4.6125e-16

// detErrorMultiplier is the factor to scale the magnitudes by when checking
// for the sign of set of points with certainty. Using a similar technique to
// the one used for maxDeterminantError, the error is at most:
//
//   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
//
// If the determinant magnitude is larger than this value then we know its sign with certainty.
let detErrorMultiplier = 7.1767e-16

/// Represents a point on the unit sphere as a normalized 3D vector.
/// Points are guaranteed to be close to normalized.
/// Fields should be treated as read-only. Use one of the factory methods for creation.
/// Represents a point in RxRxR.
public struct S2Point {
  
  let x: Double
  let y: Double
  let z: Double
  
  //
  static let epsilon = 1e-14
  
  // MARK: inits / factory
  
  init(x: Double, y: Double, z: Double, normalize: Bool) {
    if normalize {
      (self.x, self.y, self.z) = S2Point.normalize(x: x, y: y, z: z)
    } else {
      (self.x, self.y, self.z) = (x, y, z)
    }
  }
  
  init(x: Double, y: Double, z: Double) {
    // normalize explicitly to prevent needing recursive construction
    (self.x, self.y, self.z) = S2Point.normalize(x: x, y: y, z: z)
  }
  
  /// Creates a new normalized point from coordinates.
  /// This always returns a valid point. If the given coordinates can not be normalized
  /// the origin point will be returned.
  init(raw: R3Vector) {
    self.init(x: raw.x, y: raw.y, z: raw.z)
  }
  
  static func normalize(x: Double, y: Double, z: Double) -> (Double, Double, Double) {
    let norm2 = x * x + y * y + z * z
    guard norm2 != 0.0 else {
      // this is different from R3Vector, instead of 0,0,0 we pick a point that's on the unit sphere
      // but not quite on an edge so this remains stable(?)
      return (0.00456762077230, 0.99947476613078, 0.03208315302933)
    }
    let norm = sqrt(norm2)
    return (x / norm, y / norm, z / norm)
  }
  
  /// A unique "origin" on the sphere for operations that need a fixed
  /// reference point. In particular, this is the "point at infinity" used for
  /// point-in-polygon testing (by counting the number of edge crossings).
  ///
  /// It should *not* be a point that is commonly used in edge tests in order
  /// to avoid triggering code to handle degenerate cases (this rules out the
  /// north and south poles). It should also not be on the boundary of any
  /// low-level S2Cell for the same reason.
  static let origin = S2Point(x: 0, y: 0, z: 0)
  
  // MARK: tests
  
  /// Returns whether this S2Point is of approximately unit length.
  var isUnit: Bool {
    return abs(v.norm2 - 1.0) <= S2Point.epsilon
  }
  
  // MARK: computed members
  
  var v: R3Vector {
    return R3Vector(x: x, y: y, z:z)
  }
  
  /// Returns the S2Point's norm.
  var norm: Double {
    return sqrt(dot(self))
  }
  
  /// Returns the square of the norm.
  var norm2: Double {
    return dot(self)
  }
  
  /// Returns a unit S2Point in the same direction as
  func normalized() -> R3Vector {
    if x == 0.0 && y == 0.0 && z == 0.0 {
      return mul(1.0) //self
    }
    return mul(1.0 / norm)
  }
  
  /// Returns the S2Point with nonnegative components.
  func absolute() -> S2Point {
    return S2Point(x: abs(x), y: abs(y), z: abs(z))
  }
  
  func largestAbsComponent() -> Int {
    let temp = absolute()
    if temp.x > temp.y {
      if temp.x > temp.z { return 0 }
      return 2
    }
    if temp.y > temp.z { return 1 }
    return 2
  }

  /// Returns a unit S2Point that is orthogonal to
  /// Ortho(-v) = -Ortho(v) for all
  func ortho() -> R3Vector {
    // Grow a component other than the largest in v, to guarantee that they aren't
    // parallel (which would make the cross product zero).
    let other: S2Point
    if abs(x) > abs(y) {
      other = S2Point(x: 0.012, y: 1.0, z: 0.00457)
    } else {
      other = S2Point(x: 1.0, y: 0.0053, z: 0.00457)
    }
    return cross(other).normalized()
  }
  
  // MARK: arithmetic
  // Some of the below methods are convenience methods so the R3Vector does not have to be exposed all the time.
  // Those return R3Vector because normalization into S2 is expensive and often not wanted at this point.
  
  /// The inverse vector
  func inverse() -> S2Point {
    return S2Point(x: -x, y: -y, z: -z)
  }
  
  /// Returns the standard S2Point difference of v and other.
  func add(_ point: S2Point) -> R3Vector {
    return R3Vector(x: x + point.x, y: y + point.y, z: z + point.z)
  }
  
  static func +(lhs: S2Point, rhs: S2Point) -> R3Vector {
    return R3Vector(x: lhs.x + rhs.x, y: lhs.y + rhs.y, z: lhs.z + rhs.z)
  }
  
  /// Returns the standard S2Point difference of v and other.
  func sub(_ point: S2Point) -> R3Vector {
    return R3Vector(x: x - point.x, y: y - point.y, z: z - point.z)
  }
  
  static func -(lhs: S2Point, rhs: S2Point) -> R3Vector {
    return R3Vector(x: lhs.x - rhs.x, y: lhs.y - rhs.y, z: lhs.z - rhs.z)
  }
  
  /// Returns the standard scalar product of v and m.
  func mul(_ m: Double) -> R3Vector {
    return R3Vector(x: m * x, y: m * y, z: m * z)
  }
  
  static func *(lhs: S2Point, rhs: Double) -> R3Vector {
    return R3Vector(x: lhs.x * rhs, y: lhs.y * rhs, z: lhs.z * rhs)
  }
  
  /// Returns the standard dot product of v and other.
  func dot(_ point: S2Point) -> Double {
    return x * point.x + y * point.y + z * point.z
  }

  static func *(lhs: S2Point, rhs: S2Point) -> Double {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z
  }
  
  /// Returns the standard cross product of v and other.
  func cross(_ point: S2Point) -> R3Vector {
    return R3Vector(x: y * point.z - z * point.y, y: z * point.x - x * point.z, z: x * point.y - y * point.x)
  }
  
  /// Returns the Euclidean distance between v and other.
  func euclideanDistance(_ point: S2Point) -> Double {
    return v.sub(point.v).norm
  }

  /// Returns the angle between two points.
  func distance(_ b: S2Point) -> Double {
    return angle(b)
  }
  
  /// Returns the angle between v and other.
  func angle(_ point: S2Point) -> Double {
    return atan2(v.cross(point.v).norm, v.dot(point.v))
  }
  
  /// Returns the exterior angle at vertex B in the triangle ABC. The
  /// return value is positive if ABC is counterclockwise and negative otherwise.
  /// If you imagine an ant walking from A to B to C, this is the angle that the
  /// ant turns at vertex B (positive = left = CCW, negative = right = CW).
  /// This quantity is also known as the "geodesic curvature" at B.
  ///
  /// Ensures that TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all distinct
  /// a,b,c. The result is undefined if (a == b || b == c), but is either
  /// -Pi or Pi if (a == c). All points should be normalized.
  static func turnAngle(a: S2Point, b: S2Point, c: S2Point) -> S1Angle {
    // We use PointCross to get good accuracy when two points are very
    // close together, and RobustSign to ensure that the sign is correct for
    // turns that are close to 180 degrees.
    let angle = a.pointCross(b).angle(b.pointCross(c))
    // Don't return RobustSign * angle because it is legal to have (a == c).
    if S2Point.robustSign(a, b, c) == .counterClockwise {
      return angle
    }
    return -angle
  }
  
  /// Rotates the given point about the given axis by the given angle. p and
  /// axis must be unit length; angle has no restrictions (e.g., it can be
  /// positive, negative, greater than 360 degrees, etc).
  func rotate(p: S2Point, axis: S2Point, angle: S1Angle) -> S2Point {
    // Let M be the plane through P that is perpendicular to axis, and let
    // center be the point where M intersects axis. We construct a
    // right-handed orthogonal frame (dx, dy, center) such that dx is the
    // vector from center to P, and dy has the same length as dx. The
    // result can then be expressed as (cos(angle)*dx + sin(angle)*dy + center).
    let center = axis.mul(p.dot(axis))
    let dx = p.v.sub(center)
    let dy = axis.cross(p)
    // Mathematically the result is unit length, but normalization is necessary
    // to ensure that numerical errors don't accumulate.
    let v = dx.mul(cos(angle)).add(dy.mul(sin(angle))).add(center)
    return S2Point(raw: v)
  }
  
  /// Returns a Point that is orthogonal to both p and op. This is similar to
  /// p.Cross(op) (the true cross product) except that it does a better job of
  /// ensuring orthogonality when the Point is nearly parallel to op, it returns
  /// a non-zero result even when p == op or p == -op and the result is a Point,
  /// so it will have norm 1.
  ///
  /// It satisfies the following properties (f == PointCross):
  ///
  ///   (1) f(p, op) != 0 for all p, op
  ///   (2) f(op,p) == -f(p,op) unless p == op or p == -op
  ///   (3) f(-p,op) == -f(p,op) unless p == op or p == -op
  ///   (4) f(p,-op) == -f(p,op) unless p == op or p == -op
  func pointCross(_ op: S2Point) -> S2Point {
    // NOTE(dnadasi): In the C++ API the equivalent method here was known as "RobustCrossProd",
    // but PointCross more accurately describes how this method is used.
    let x = v.add(op.v).cross(op.v.sub(v))
    if x.approxEquals(R3Vector(x: 0, y: 0, z: 0)) {
      // The only result that makes sense mathematically is to return zero, but
      // we find it more convenient to return an arbitrary orthogonal vector.
      return S2Point(raw: v.ortho())
    }
    return S2Point(raw: x)
  }

  /// Returns true if the points A, B, C are strictly counterclockwise,
  /// and returns false if the points are clockwise or collinear (i.e. if they are all
  /// contained on some great circle).
  ///
  /// Due to numerical errors, situations may arise that are mathematically
  /// impossible, e.g. ABC may be considered strictly CCW while BCA is not.
  /// However, the implementation guarantees the following:
  ///
  ///   If Sign(a,b,c), then !Sign(c,b,a) for all a,b,c.
  static func sign(_ a: S2Point, b: S2Point, c: S2Point) -> Bool {
    // NOTE(dnadasi): In the C++ API the equivalent method here was known as "SimpleSign".
    // We compute the signed volume of the parallelepiped ABC. The usual
    // formula for this is (A ⨯ B) · C, but we compute it here using (C ⨯ A) · B
    // in order to ensure that ABC and CBA are not both CCW. This follows
    // from the following identities (which are true numerically, not just
    // mathematically):
    //
    //     (1) x ⨯ y == -(y ⨯ x)
    //     (2) -x · y == -(x · y)
    return c.v.cross(a.v).dot(b.v) > 0
  }

  /// Returns a Direction representing the ordering of the points.
  /// CounterClockwise is returned if the points are in counter-clockwise order,
  /// Clockwise for clockwise, and Indeterminate if any two points are the same (collinear),
  /// or the sign could not completely be determined.
  ///
  /// This function has additional logic to make sure that the above properties hold even
  /// when the three points are coplanar, and to deal with the limitations of
  /// floating-point arithmetic.
  ///
  /// RobustSign satisfies the following conditions:
  ///
  ///  (1) RobustSign(a,b,c) == Indeterminate if and only if a == b, b == c, or c == a
  ///  (2) RobustSign(b,c,a) == RobustSign(a,b,c) for all a,b,c
  ///  (3) RobustSign(c,b,a) == -RobustSign(a,b,c) for all a,b,c
  ///
  /// In other words:
  ///
  ///  (1) The result is Indeterminate if and only if two points are the same.
  ///  (2) Rotating the order of the arguments does not affect the result.
  ///  (3) Exchanging any two arguments inverts the result.
  ///
  /// On the other hand, note that it is not true in general that
  /// RobustSign(-a,b,c) == -RobustSign(a,b,c), or any similar identities
  /// involving antipodal points.
  static func robustSign(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> S2Direction {
    let sign = triageSign(a, b, c)
    if sign == .indeterminate {
      return expensiveSign(a, b, c)
    }
    return sign
  }

  /// Returns the direction sign of the points. It returns Indeterminate if two
  /// points are identical or the result is uncertain. Uncertain cases can be resolved, if
  /// desired, by calling expensiveSign.
  ///
  /// The purpose of this method is to allow additional cheap tests to be done without
  /// calling expensiveSign.
  static func triageSign(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> S2Direction {
    let det = c.v.cross(a.v).dot(b.v)
    if det > maxDeterminantError {
      return .counterClockwise
    }
    if det < -maxDeterminantError {
      return .clockwise
    }
    return .indeterminate
  }

  /// Reports the direction sign of the points. It returns Indeterminate
  /// if two of the input points are the same. It uses multiple-precision arithmetic
  /// to ensure that its results are always self-consistent.
  static func expensiveSign(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> S2Direction {
    // Return Indeterminate if and only if two points are the same.
    // This ensures RobustSign(a,b,c) == Indeterminate if and only if a == b, b == c, or c == a.
    // ie. Property 1 of RobustSign.
    if a == b || b == c || c == a {
      return .indeterminate
    }
    // Next we try recomputing the determinant still using floating-point
    // arithmetic but in a more precise way. This is more expensive than the
    // simple calculation done by triageSign, but it is still *much* cheaper
    // than using arbitrary-precision arithmetic. This optimization is able to
    // compute the correct determinant sign in virtually all cases except when
    // the three points are truly collinear (e.g., three points on the equator).
    let detSign = stableSign(a, b, c)
    if detSign != .indeterminate {
      return detSign
    }
    // Otherwise fall back to exact arithmetic and symbolic permutations.
    return exactSign(a, b, c)
  }

  /// Reports the direction sign of the points in a numerically stable way.
  /// Unlike triageSign, this method can usually compute the correct determinant sign even when all
  /// three points are as collinear as possible. For example if three points are
  /// spaced 1km apart along a random line on the Earth's surface using the
  /// nearest representable points, there is only a 0.4% chance that this method
  /// will not be able to find the determinant sign. The probability of failure
  /// decreases as the points get closer together; if the collinear points are
  /// 1 meter apart, the failure rate drops to 0.0004%.
  ///
  /// This method could be extended to also handle nearly-antipodal points (and
  /// in fact an earlier version of this code did exactly that), but antipodal
  /// points are rare in practice so it seems better to simply fall back to
  /// exact arithmetic in that case.
  static func stableSign(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> S2Direction {
    let ab = a.v.sub(b.v)
    let ab2 = ab.norm2
    let bc = b.v.sub(c.v)
    let bc2 = bc.norm2
    let ca = c.v.sub(a.v)
    let ca2 = ca.norm2
    // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
    // cyclically permuted if necessary so that AB is the longest edge. (This
    // minimizes the magnitude of cross product.)  At the same time we also
    // compute the maximum error in the determinant.
    // The two shortest edges, pointing away from their common point.
    let e1: R3Vector
    let e2: R3Vector
    let op: S2Point
    if ab2 >= bc2 && ab2 >= ca2 {
      // AB is the longest edge.
      e1 = ca
      e2 = bc
      op = c
    } else if bc2 >= ca2 {
      // BC is the longest edge.
      e1 = ab
      e2 = ca
      op = a
    } else {
      // CA is the longest edge.
      e1 = bc
      e2 = ab
      op = b
    }
    let det = e1.cross(e2).dot(op.v)
    let maxErr = detErrorMultiplier * sqrt(e1.norm2 * e2.norm2)
    // If the determinant isn't zero, within maxErr, we know definitively the point ordering.
    if det > maxErr {
      return .counterClockwise
    }
    if det < -maxErr {
      return .clockwise
    }
    return .indeterminate
  }

  /// Reports the direction sign of the points using exact precision arithmetic.
  static func exactSign(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> S2Direction {
    // In the C++ version, the final computation is performed using OpenSSL's
    // Bignum exact precision math library. The existence of an equivalent
    // library in Go is indeterminate. In C++, using the exact precision library
    // to solve this stage is ~300x slower than the above checks.
    // TODO: Select and incorporate an appropriate Go exact precision
    // floating point library for the remaining calculations.
    return .indeterminate
  }

  /// Returns true if the edges OA, OB, and OC are encountered in that
  /// order while sweeping CCW around the point O.
  ///
  /// You can think of this as testing whether A <= B <= C with respect to the
  /// CCW ordering around O that starts at A, or equivalently, whether B is
  /// contained in the range of angles (inclusive) that starts at A and extends
  /// CCW to C. Properties:
  ///
  ///  (1) If OrderedCCW(a,b,c,o) && OrderedCCW(b,a,c,o), then a == b
  ///  (2) If OrderedCCW(a,b,c,o) && OrderedCCW(a,c,b,o), then b == c
  ///  (3) If OrderedCCW(a,b,c,o) && OrderedCCW(c,b,a,o), then a == b == c
  ///  (4) If a == b or b == c, then OrderedCCW(a,b,c,o) is true
  ///  (5) Otherwise if a == c, then OrderedCCW(a,b,c,o) is false
  static func orderedCCW(_ a: S2Point, _ b: S2Point, _ c: S2Point, _ o: S2Point) -> Bool {
    var sum = 0
    if robustSign(b, o, a) != .clockwise {
      sum += 1
    }
    if robustSign(c, o, b) != .clockwise {
      sum += 1
    }
    if robustSign(a, o, c) == .counterClockwise {
      sum += 1
    }
    return sum >= 2
  }

  static func signedArea(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> Double {
    return pointArea(a, b, c) * Double(robustSign(a, b, c).rawValue)
  }
  
  /// Returns the area on the unit sphere for the triangle defined by the
  /// given points.
  ///
  /// This method is based on l'Huilier's theorem,
  ///
  ///   tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
  ///
  /// where E is the spherical excess of the triangle (i.e. its area),
  ///       a, b, c are the side lengths, and
  ///       s is the semiperimeter (a + b + c) / 2.
  ///
  /// The only significant source of error using l'Huilier's method is the
  /// cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
  /// *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
  /// to a relative error of about 1e-15 / E using Girard's formula, where E is
  /// the true area of the triangle. Girard's formula can be even worse than
  /// this for very small triangles, e.g. a triangle with a true area of 1e-30
  /// might evaluate to 1e-5.
  ///
  /// So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
  /// dmin = min(s-a, s-b, s-c). This basically includes all triangles
  /// except for extremely long and skinny ones.
  ///
  /// Since we don't know E, we would like a conservative upper bound on
  /// the triangle area in terms of s and dmin. It's possible to show that
  /// E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
  /// Using this, it's easy to show that we should always use l'Huilier's
  /// method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
  /// if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
  /// k3 is about 0.1. Since the best case error using Girard's formula
  /// is about 1e-15, this means that we shouldn't even consider it unless
  /// s >= 3e-4 or so.
  static func pointArea(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> Double {
    let sa = b.angle(c)
    let sb = c.angle(a)
    let sc = a.angle(b)
    let s = 0.5 * (sa + sb + sc)
    if s >= 3e-4 {
      // Consider whether Girard's formula might be more accurate.
      let dmin = s - max(sa, max(sb, sc))
      if dmin < 1e-2*s*s*s*s*s {
        // This triangle is skinny enough to use Girard's formula.
        let ab = a.pointCross(b)
        let bc = b.pointCross(c)
        let ac = a.pointCross(c)
        let area = max(0.0, ab.angle(ac)-ab.angle(bc)+bc.angle(ac))
        if dmin < s * 0.1 * area {
          return area
        }
      }
    }
    // Use l'Huilier's formula.
    return 4.0 * atan(sqrt(max(0.0, tan(0.5*s) * tan(0.5*(s-sa)) * tan(0.5*(s-sb)) * tan(0.5*(s-sc)))))
  }

  /// Returns the true centroid of the spherical triangle ABC multiplied by the
  /// signed area of spherical triangle ABC. The result is not normalized.
  /// The reasons for multiplying by the signed area are (1) this is the quantity
  /// that needs to be summed to compute the centroid of a union or difference of triangles,
  /// and (2) it's actually easier to calculate this way. All points must have unit length.
  ///
  /// The true centroid (mass centroid) is defined as the surface integral
  /// over the spherical triangle of (x,y,z) divided by the triangle area.
  /// This is the point that the triangle would rotate around if it was
  /// spinning in empty space.
  ///
  /// The best centroid for most purposes is the true centroid. Unlike the
  /// planar and surface centroids, the true centroid behaves linearly as
  /// regions are added or subtracted. That is, if you split a triangle into
  /// pieces and compute the average of their centroids (weighted by triangle
  /// area), the result equals the centroid of the original triangle. This is
  /// not true of the other centroids.
  static func trueCentroid(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> R3Vector {
    var ra = 1.0
    let sa = b.distance(c)
    if sa != 0 {
      ra = sa / sin(sa)
    }
    var rb = 1.0
    let sb = c.distance(a)
    if sb != 0 {
      rb = sb / sin(sb)
    }
    var rc = 1.0
    let sc = a.distance(b)
    if sc != 0 {
      rc = sc / sin(sc)
    }
    // Now compute a point M such that:
    //
    //  [Ax Ay Az] [Mx]                       [ra]
    //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
    //  [Cx Cy Cz] [Mz]                       [rc]
    //
    // To improve the numerical stability we subtract the first row (A) from the
    // other two rows; this reduces the cancellation error when A, B, and C are
    // very close together. Then we solve it using Cramer's rule.
    //
    // This code still isn't as numerically stable as it could be.
    // The biggest potential improvement is to compute B-A and C-A more
    // accurately so that (B-A)x(C-A) is always inside triangle ABC.
    let x = R3Vector(x: a.x, y: b.x - a.x, z: c.x - a.x)
    let y = R3Vector(x: a.y, y: b.y - a.y, z: c.y - a.y)
    let z = R3Vector(x: a.z, y: b.z - a.z, z: c.z - a.z)
    let r = R3Vector(x: ra, y: rb - ra, z: rc - ra)
    return R3Vector(x: y.cross(z).dot(r), y: z.cross(x).dot(r), z: x.cross(y).dot(r)).mul(0.5)
  }

  /// Returns the centroid of the planar triangle ABC, which is not normalized.
  /// It can be normalized to unit length to obtain the "surface centroid" of the corresponding
  /// spherical triangle, i.e. the intersection of the three medians. However,
  /// note that for large spherical triangles the surface centroid may be
  // nowhere near the intuitive "center" (see example in TrueCentroid comments).
  ///
  /// Note that the surface centroid may be nowhere near the intuitive
  /// "center" of a spherical triangle. For example, consider the triangle
  /// with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
  /// The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
  /// within a distance of 2*eps of the vertex B. Note that the median from A
  /// (the segment connecting A to the midpoint of BC) passes through S, since
  /// this is the shortest path connecting the two endpoints. On the other
  /// hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
  /// the surface is a much more reasonable interpretation of the "center" of
  /// this triangle.
  static func planarCentroid(_ a: S2Point, _ b: S2Point, _ c: S2Point) -> S2Point {
    return S2Point(raw: a.v.add(b.v).add(c.v).mul(1.0 / 3.0))
  }

  /// Constructs a ChordAngle corresponding to the distance
  /// between the two given points. The points must be unit length.
  func chordAngleBetweenPoints(x: S2Point, y: S2Point) -> S1ChordAngle {
    return S1ChordAngle(value: min(4.0, (x - y).norm2))
  }
  
  /// Generates a slice of points shaped as a regular polygon with
  /// the numVertices vertices, all located on a circle of the specified angular radius
  /// around the center. The radius is the actual distance from center to each vertex.
  static func regularPoints(center: S2Point, radius: S1Angle, numVertices: Int) -> [S2Point] {
    return S2Point.regularPointsForFrame(frame: S2Point.getFrame(center), radius: radius, numVertices: numVertices)
  }
  
  /// Generates a slice of points shaped as a regular polygon
  /// with numVertices vertices, all on a circle of the specified angular radius around
  /// the center. The radius is the actual distance from the center to each vertex.
  static func regularPointsForFrame(frame: R3Matrix, radius: S1Angle, numVertices: Int) -> [S2Point] {
    // We construct the loop in the given frame coordinates, with the center at
    // (0, 0, 1). For a loop of radius r, the loop vertices have the form
    // (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r). The distance on the
    // sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
    let z = cos(radius)
    let r = sin(radius)
    let radianStep = 2 * .pi / Double(numVertices)
    var vertices: [S2Point] = []
    for i in 0..<numVertices {
      let angle = Double(i) * radianStep
      let p = S2Point(x: r * cos(angle), y: r * sin(angle), z: z)
      vertices.append(S2Point.fromFrame(frame, point: p))
    }
    return vertices
  }
  
  // MARK: lat lng
  
  /// Latitude
  var lat: Double {
    return atan2(z, sqrt(x * x + y * y))
  }
  
  /// Longitude
  var lng: Double {
    return atan2(y, x)
  }
  
  /// LatLng representation
  var latLng: LatLng {
    return LatLng(lat: lat, lng: lng)
  }
  
  // MARK: frames
  
  /// Returns the orthonormal frame for the given point on the unit sphere.
  static func getFrame(_ point: S2Point) -> R3Matrix {
    // Given the point p on the unit sphere, extend this into a right-handed
    // coordinate frame of unit-length column vectors m = (x,y,z).  Note that
    // the vectors (x,y) are an orthonormal frame for the tangent space at point p,
    // while p itself is an orthonormal frame for the normal space at p.
    let o = point.v.ortho()
    let p1 = S2Point(raw: o.cross(point.v))
    let p2 = S2Point(raw: o)
    let p3 = point
    return R3Matrix(
      r1: [p1.x, p2.x, p3.x],
      r2: [p1.y, p2.y, p3.y],
      r3: [p1.z, p2.z, p3.z])
  }
  
  /// Returns the coordinates of the given point with respect to its orthonormal basis m.
  /// The resulting point q satisfies the identity (m * q == p).
  static func toFrame(_ matrix: R3Matrix, point: S2Point) -> S2Point {
    // The inverse of an orthonormal matrix is its transpose.
    return matrix.transpose().mul(point: point)
  }
  
  /// Returns the coordinates of the given point in standard axis-aligned basis
  /// from its orthonormal basis m.
  /// The resulting point p satisfies the identity (p == m * q).
  static func fromFrame(_ matrix: R3Matrix, point: S2Point) -> S2Point {
    return matrix.mul(point: point)
  }
  
}

extension S2Point: Equatable, CustomStringConvertible, Approximatable, Hashable, Comparable {
  
  public static func ==(lhs: S2Point, rhs: S2Point) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z
  }
  
  public var description: String {
    return "(\(x), \(y), \(z))"
  }
  
  /// Reports whether self and point are equal within a small epsilon.
  public func approxEquals(_ point: S2Point) -> Bool {
    return angle(point) <= S2Point.epsilon
    // as opposed to the vector implementation
    // return abs(x-point.x) < epsilon && abs(y-point.y) < epsilon && abs(z-point.z) < epsilon
  }
  
  public var hashValue: Int {
    return x.hashValue ^ y.hashValue ^ z.hashValue
  }
  
  /// This is a lexicographic comparison. Kinda a hack, makes them sortable
  public static func <(lhs: S2Point, rhs: S2Point) -> Bool {
    if lhs.x != rhs.x {
      return lhs.x < rhs.x
    }
    if lhs.y != rhs.y {
      return lhs.y < rhs.y
    }
    return lhs.z < rhs.z
  }
  
}

