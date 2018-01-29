//
//  R3Vector.swift
//  s2-swift
//

import Foundation


/// Enumerates the 3 axes of ℝ³.
public enum Axis {
  case x, y, z
}

/// Represents a point in RxRxR.
public struct R3Vector {
  
  //
  let x: Double
  let y: Double
  let z: Double
    
  // MARK: inits

  public init(x: Double, y: Double, z: Double) {
    self.x = x
    self.y = y
    self.z = z
  }

  // MARK: tests

  /// Returns whether this vector is of approximately unit length.
  var isUnit: Bool {
    return abs(norm2 - 1) <= R1Interval.epsilon
  }
  
  // MARK: computed members
  
  /// Returns the vector's norm.
  var norm: Double {
    return sqrt(dot(self))
  }

  /// Returns the square of the norm.
  var norm2: Double {
    return dot(self)
  }

  /// Returns a unit vector in the same direction as
  func normalized() -> R3Vector {
    if x == 0.0 && y == 0.0 && z == 0.0 {
      return self
    }
    return mul(1.0 / norm)
  }

  /// Returns the vector with nonnegative components.
  func absolute() -> R3Vector {
    return R3Vector(x: abs(x), y: abs(y), z: abs(z))
  }

  /// Returns a unit vector that is orthogonal to
  /// Ortho(-v) = -Ortho(v) for all
  func ortho() -> R3Vector {
    // Grow a component other than the largest in v, to guarantee that they aren't
    // parallel (which would make the cross product zero).
    let vector: R3Vector
    if abs(x) > abs(y) {
      vector = R3Vector(x: 0.012, y: 1.0, z: 0.00457)
    } else {
      vector = R3Vector(x: 1.0, y: 0.0053, z: 0.00457)
    }
    return cross(vector).normalized()
  }
  
  var s2: S2Point {
    return S2Point(raw: self)
  }
  
  // MARK: arithmetic
  
  /// Returns the standard vector sum of self and vector.
  func add(_ vector: R3Vector) -> R3Vector {
    return R3Vector(x: x + vector.x, y: y + vector.y, z: z + vector.z)
  }
  
  static func +(lhs: R3Vector, rhs: R3Vector) -> R3Vector {
    return R3Vector(x: lhs.x + rhs.x, y: lhs.y + rhs.y, z: lhs.z + rhs.z)
  }
  
  /// Returns the standard vector difference of self and vector.
  func sub(_ vector: R3Vector) -> R3Vector {
    return R3Vector(x: x - vector.x, y: y - vector.y, z: z - vector.z)
  }

  static func -(lhs: R3Vector, rhs: R3Vector) -> R3Vector {
    return R3Vector(x: lhs.x - rhs.x, y: lhs.y - rhs.y, z: lhs.z - rhs.z)
  }
  
  /// Returns the standard scalar product of self and m.
  func mul(_ m: Double) -> R3Vector {
    return R3Vector(x: m * x, y: m * y, z: m * z)
  }

  static func *(lhs: R3Vector, rhs: Double) -> R3Vector {
    return R3Vector(x: lhs.x * rhs, y: lhs.y * rhs, z: lhs.z * rhs)
  }
  
  /// Returns the standard dot product of self and vector.
  func dot(_ vector: R3Vector) -> Double {
    return x*vector.x + y*vector.y + z*vector.z
  }

  static func *(lhs: R3Vector, rhs: R3Vector) -> Double {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z
  }
  
  /// Returns the standard cross product of self and vector.
  func cross(_ vector: R3Vector) -> R3Vector {
    return R3Vector(x: y*vector.z - z*vector.y, y: z*vector.x - x*vector.z, z: x*vector.y - y*vector.x)
  }

  /// Returns the Euclidean distance between self and vector.
  func distance(_ vector: R3Vector) -> Double {
    return sub(vector).norm
  }

  /// Returns the angle between self and vector.
  func angle(_ vector: R3Vector) -> Double {
    return atan2(cross(vector).norm, dot(vector))
  }

}

extension R3Vector: Equatable, CustomStringConvertible, Approximatable, Hashable, Comparable {
  
  public static func ==(lhs: R3Vector, rhs: R3Vector) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z
  }
  
  public var description: String {
    return "(\(x), \(y), \(z))"
  }
  
  /// Reports whether self and vector are equal within a small epsilon.
  public func approxEquals(_ vector: R3Vector) -> Bool {
    return abs(x - vector.x) < R1Interval.epsilon && abs(y - vector.y) < R1Interval.epsilon && abs(z - vector.z) < R1Interval.epsilon
  }
  
  public var hashValue: Int {
    return x.hashValue ^ y.hashValue ^ z.hashValue
  }

  /// This is a lexicographic comparison. Kinda a hack, makes vectors sortable
  public static func <(lhs: R3Vector, rhs: R3Vector) -> Bool {
    if lhs.x != rhs.x {
      return lhs.x < rhs.x
    }
    if lhs.y != rhs.y {
      return lhs.y < rhs.y
    }
    return lhs.z < rhs.z
  }

}

typealias BigFloat = Double

// prec is the number of bits of precision to use for the Float values.
// To keep things simple, we use the maximum allowable precision on big
// values. This allows us to handle all values we expect in the s2 library.
//let prec = big.MaxPrec

// define some commonly referenced values.
let precise0 = BigFloat(0)
let precise1 = BigFloat(1)

func precAdd(_ lhs: BigFloat, _ rhs: BigFloat) -> BigFloat {
  return lhs + rhs
}

func precSub(_ lhs: BigFloat, _ rhs: BigFloat) -> BigFloat {
  return lhs - rhs
}

func precMul(_ lhs: BigFloat, _ rhs: BigFloat) -> BigFloat {
  return lhs * rhs
}
// PreciseVector represents a point in ℝ³ using high-precision values.
// Note that this is NOT a complete implementation because there are some
// operations that Vector supports that are not feasible with arbitrary precision
// math. (e.g., methods that need divison like Normalize, or methods needing a
// square root operation such as Norm)
struct PreciseVector {
  let x: BigFloat
  let y: BigFloat
  let z: BigFloat
}

extension PreciseVector {

  /// Creates a high precision vector from the given Vector.
  init(v: R3Vector) {
    x = v.x
    y = v.y
    z = v.z
  }

  /// Creates a high precision vector from the given floating point values.
//  init(x: Double, y: Double, z: Double) {
//    self.x = precFloat(x)
//    self.y = precFloat(y)
//    self.z = precFloat(z)
//  }
  
  // Vector returns this precise vector converted to a Vector.
  func vector() -> R3Vector {
    // The accuracy flag is ignored on these conversions back to float64.
    let xd = Double(x)
    let yd = Double(y)
    let zd = Double(z)
    return R3Vector(x: xd, y: yd, z: zd)
  }

  /// Reports whether v and ov are equal.
  static func ==(lhs: PreciseVector, rhs: PreciseVector) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z
  }

  var description: String {
    return "(\(x), \(y), \(z)"
  }

  /// Returns the square of the norm.
  var norm2: BigFloat { return dot(self) }
  
  // IsUnit reports whether this vector is of unit length.
  var isUnit: Bool {
    return norm2 == precise1
  }
  
  // Abs returns the vector with nonnegative components.
  func abs() -> PreciseVector {
    return PreciseVector(x: Swift.abs(x), y: Swift.abs(y), z: Swift.abs(z))
  }
  
  // Add returns the standard vector sum of v and ov.
  func add(_ ov: PreciseVector) -> PreciseVector {
    return PreciseVector(x: precAdd(x, ov.x), y: precAdd(y, ov.y), z: precAdd(z, ov.z))
  }
  
  // Sub returns the standard vector difference of v and ov.
  func sub(_ ov: PreciseVector) -> PreciseVector {
    return PreciseVector(x: precSub(x, ov.x), y: precSub(y, ov.y), z: precSub(z, ov.z))
  }
  
  // Mul returns the standard scalar product of v and f.
  func mul(_ f: BigFloat) -> PreciseVector {
    return PreciseVector(x: precMul(x, f), y: precMul(y, f), z: precMul(z, f))
  }

  // MulByFloat64 returns the standard scalar product of v and f.
//  func mul(_ f: Double) -> PreciseVector {
//    return mul(precFloat(f))
//  }

  // Dot returns the standard dot product of v and ov.
  func dot(_ ov: PreciseVector) -> BigFloat {
    return precAdd(precMul(x, ov.x), precAdd(precMul(y, ov.y), precMul(z, ov.z)))
  }
  
  // Cross returns the standard cross product of v and ov.
  func cross(_ ov: PreciseVector) -> PreciseVector {
    return PreciseVector(
      x: precSub(precMul(y, ov.z), precMul(z, ov.y)),
      y: precSub(precMul(z, ov.x), precMul(x, ov.z)),
      z: precSub(precMul(x, ov.y), precMul(y, ov.x)))
  }
  
  // LargestComponent returns the axis that represents the largest component in this vector.
  func largestComponent() -> Axis {
    let t = abs()
    if t.x > t.y {
      if t.x > t.z {
        return .x
      }
      return .z
    }
    if t.y > t.z {
      return .y
    }
    return .z
  }
  
  // SmallestComponent returns the axis that represents the smallest component in this vector.
  func smallestComponent() -> Axis {
    let t = abs()
    if t.x < t.y {
      if t.x < t.z {
        return .x
      }
      return .z
    }
    if t.y < t.z {
      return .y
    }
    return .z
  }

}
