//
//  S1ChrodAngle.swift
//  s2-swift
//

import Foundation

/// ChordAngle represents the angle subtended by a chord (i.e., the straight
/// line segment connecting two points on the sphere). Its representation
/// makes it very efficient for computing and comparing distances, but unlike
/// Angle it is only capable of representing angles between 0 and π radians.
/// Generally, ChordAngle should only be used in loops where many angles need
/// to be calculated and compared. Otherwise it is simpler to use Angle.
/// ChordAngles are represented by the squared chord length, which can
/// range from 0 to 4. Positive infinity represents an infinite squared length.
struct ChordAngle {
  
  let value: Double
  
  /// Zero
  public static let zero = ChordAngle(value: 0)
  
  /// Represents a chord angle smaller than the zero angle.
  /// The only valid operations on a NegativeChordAngle are comparisons and
  /// Angle conversions.
  public static let negative = ChordAngle(value: -1)
  
  /// Represents a chord angle of 90 degrees (a "right angle").
  public static let right = ChordAngle(value: 2)
  
  /// Represents a chord angle of 180 degrees (a "straight angle").
  /// This is the maximum finite chord angle.
  public static let straight = ChordAngle(value: 4)
  
  /// Represents a chord angle larger than any finite chord angle.
  /// The only valid operations on an InfChordAngle are comparisons and Angle conversions.
  public static let infinite = ChordAngle(value: Double.greatestFiniteMagnitude)
  
  /// The square of the maximum length allowed in a ChordAngle.
  static let maxLength2 = 4.0
  
  init(value: Double) {
    self.value = value
  }
  
  /// Returns a ChordAngle from the given Angle.
  public init(angle: S1Angle) {
    if angle < 0 {
      self = .negative
      return
    }
    if angle.isInfinite {
      self = .infinite
      return
    }
    let l = 2 * Darwin.sin(0.5 * min(.pi, angle))
    value = l * l
  }
  
  /// Returns a ChordAngle from the squared chord length.
  /// Note that the argument is automatically clamped to a maximum of 4 to
  /// handle possible roundoff errors. The argument must be non-negative.
  public init(squaredLength length2: Double) {
    if length2 > ChordAngle.maxLength2 {
      self = .straight
      return
    }
    value = length2
  }
  
  /// Returns a new ChordAngle that has been adjusted by the given error
  /// bound (which can be positive or negative). Error should be the value
  /// returned by either MaxPointError or MaxAngleError. For example:
  ///    a := ChordAngleFromPoints(x, y)
  ///    a1 := a.Expanded(a.MaxPointError())
  public func expanded(e: Double) -> ChordAngle {
    // If the angle is special, don't change it. Otherwise clamp it to the valid range.
    if isSpecial() {
      return self
    }
    return ChordAngle(value: max(0.0, min(ChordAngle.maxLength2, value + e)))
  }
  
  /// Converts this ChordAngle to an Angle.
  public func angle() -> S1Angle {
    if value < 0 {
      return -toRadians
    }
    if isInfinite() {
      return Double.greatestFiniteMagnitude
    }
    return 2 * asin(0.5 * sqrt(value))
  }
  
  /// Reports whether this ChordAngle is infinite.
  public func isInfinite() -> Bool {
    return value == Double.greatestFiniteMagnitude
  }
  
  /// Reports whether this ChordAngle is one of the special cases.
  public func isSpecial() -> Bool {
    return value < 0.0 || isInfinite()
  }
  
  /// Reports whether this ChordAngle is valid or not.
  public func isValid() -> Bool {
    return (value >= 0 && value <= ChordAngle.maxLength2) || isSpecial()
  }
  
  /// Returns the smallest representable ChordAngle larger than this one.
  /// This can be used to convert a "<" comparison to a "<=" comparison.
  ///
  /// Note the following special cases:
  ///   NegativeChordAngle.Successor == 0
  ///   StraightChordAngle.Successor == InfChordAngle
  ///   InfChordAngle.Successor == InfChordAngle
  public func successor() -> ChordAngle {
    if value >= ChordAngle.maxLength2 {
      return .infinite
    }
    if value < 0 {
      return .zero
    }
    return ChordAngle(value: nextafter(value, 10.0))
  }
  
  /// Returns the largest representable ChordAngle less than this one.
  ///
  /// Note the following special cases:
  ///   InfChordAngle.Predecessor == StraightChordAngle
  ///   ChordAngle(0).Predecessor == NegativeChordAngle
  ///   NegativeChordAngle.Predecessor == NegativeChordAngle
  public func predecessor() -> ChordAngle {
    if value <= 0 {
      return .negative
    }
    if value > ChordAngle.maxLength2 {
      return .straight
    }
    return ChordAngle(value: nextafter(value, -10.0))
  }
  
  /// Returns the maximum error size for a ChordAngle constructed
  /// from 2 Points x and y, assuming that x and y are normalized to within the
  /// bounds guaranteed by s2.Point.Normalize. The error is defined with respect to
  /// the true distance after the points are projected to lie exactly on the sphere.
  public func maxPointError() -> Double {
    // There is a relative error of (2.5*dblEpsilon) when computing the squared
    // distance, plus an absolute error of (16 * dblEpsilon**2) because the
    // lengths of the input points may differ from 1 by up to (2*dblEpsilon) each.
    return 2.5 * S1Interval.dblEpsilon * value + 16 * S1Interval.dblEpsilon * S1Interval.dblEpsilon
  }
  
  /// Returns the maximum error for a ChordAngle constructed
  /// as an Angle distance.
  public func maxAngleError() -> Double {
    return S1Interval.dblEpsilon * value
  }
  
  /// Adds the other ChordAngle to this one and returns the resulting value.
  /// This method assumes the ChordAngles are not special.
  public func add(other: ChordAngle) -> ChordAngle {
    // Note that this method (and Sub) is much more efficient than converting
    // the ChordAngle to an Angle and adding those and converting back. It
    // requires only one square root plus a few additions and multiplications.
    // Optimization for the common case where b is an error tolerance
    // parameter that happens to be set to zero.
    if other.value == 0 {
      return self
    }
    // Clamp the angle sum to at most 180 degrees.
    if value + other.value >= ChordAngle.maxLength2 {
      return .straight
    }
    // Let a and b be the (non-squared) chord lengths, and let c = a+b.
    // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
    // Then the formula below can be derived from c = 2 * sin(A+B) and the
    // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
    //                 cos(X) = sqrt(1 - sin^2(X))
    let x = value * (1 - 0.25 * other.value)
    let y = other.value * (1 - 0.25 * value)
    return ChordAngle(value: min(ChordAngle.maxLength2, x + y + 2 * sqrt(x * y)))
  }
  
  /// Subtracts the other ChordAngle from this one and returns the resulting
  /// value. This method assumes the ChordAngles are not special.
  public func sub(other: ChordAngle) -> ChordAngle {
    if other.value == 0 {
      return self
    }
    if value <= other.value {
      return .zero
    }
    let x = value * (1 - 0.25 * other.value)
    let y = other.value * (1 - 0.25 * value)
    return ChordAngle(value: max(0.0, x + y - 2 * sqrt(x * y)))
  }
  
  // Returns the sine of this chord angle. This method is more efficient
  /// than converting to Angle and performing the computation.
  public func sin() -> Double {
    return sqrt(sin2())
  }
  
  /// Returns the square of the sine of this chord angle.
  /// It is more efficient than Sin.
  public func sin2() -> Double {
    // Let a be the (non-squared) chord length, and let A be the corresponding
    // half-angle (a = 2*sin(A)).  The formula below can be derived from:
    //   sin(2*A) = 2 * sin(A) * cos(A)
    //   cos^2(A) = 1 - sin^2(A)
    // This is much faster than converting to an angle and computing its sine.
    return value * (1 - 0.25 * value)
  }
  
  /// Returns the cosine of this chord angle. This method is more efficient
  /// than converting to Angle and performing the computation.
  public func cos() -> Double {
    // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
    return 1 - 0.5 * value
  }
  
  // Returns the tangent of this chord angle.
  public func tan() -> Double {
    return sin() / cos()
  }
  
}

extension ChordAngle: Comparable {
  
  public static func ==(lhs: ChordAngle, rhs: ChordAngle) -> Bool {
    return lhs.value == rhs.value
  }
  
  public static func <(lhs: ChordAngle, rhs: ChordAngle) -> Bool {
    return lhs.value < rhs.value
  }
  
}
