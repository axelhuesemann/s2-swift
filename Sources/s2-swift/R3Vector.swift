//
//  R3Vector.swift
//  s2-swift
//

import Foundation


/// Represents a point in RxRxR.
public struct R3Vector {
  
  static let epsilon = 1e-14

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
    return abs(norm2 - 1) <= R3Vector.epsilon
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

  /// Returns the standard vector difference of self and vector.
  func sub(_ vector: R3Vector) -> R3Vector {
    return R3Vector(x: x - vector.x, y: y - vector.y, z: z - vector.z)
  }

  /// Returns the standard scalar product of self and m.
  func mul(_ m: Double) -> R3Vector {
    return R3Vector(x: m * x, y: m * y, z: m * z)
  }

  /// Returns the standard dot product of self and vector.
  func dot(_ vector: R3Vector) -> Double {
    return x*vector.x + y*vector.y + z*vector.z
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

extension R3Vector: Equatable, CustomStringConvertible, Approximatable, Hashable {
  
  public static func ==(lhs: R3Vector, rhs: R3Vector) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z
  }
  
  public var description: String {
    return "(\(x), \(y), \(z))"
  }
  
  /// Reports whether self and vector are equal within a small epsilon.
  public func approxEquals(_ vector: R3Vector) -> Bool {
    return abs(x-vector.x) < R3Vector.epsilon && abs(y-vector.y) < R3Vector.epsilon && abs(z-vector.z) < R3Vector.epsilon
  }
  
  public var hashValue: Int {
    return x.hashValue ^ y.hashValue ^ z.hashValue
  }
  
}
