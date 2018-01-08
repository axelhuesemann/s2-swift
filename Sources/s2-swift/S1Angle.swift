//
//  S1Angle.swift
//  s2-swift
//

import Foundation

/// Represents a 1D angle.
public typealias S1Angle = Double

public let toRadians = .pi / 180.0
public let toDegrees = 180.0 / .pi

extension S1Angle {

  init(degrees: Double) {
    self = degrees * toRadians
  }

  /// Returns an equivalent angle in [0, 2π).
  func normalized() -> S1Angle {
    var rad = fmod(self, 2.0 * .pi)
    if rad < 0.0 {
      rad += 2.0 * .pi
    }
    return rad
  }

}


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
  
  /// Represents a chord angle smaller than the zero angle.
  /// The only valid operations on a NegativeChordAngle are comparisons and
  /// Angle conversions.
  static let negativeChordAngle = ChordAngle(value: -1)
  
  /// Represents a chord angle of 180 degrees (a "straight angle").
  /// This is the maximum finite chord angle.
  static let straightChordAngle = ChordAngle(value: 4)
  
  /// Represents a chord angle larger than any finite chord angle.
  /// The only valid operations on an InfChordAngle are comparisons and Angle conversions.
  static let infiniteChordAngle = ChordAngle(value: Double.greatestFiniteMagnitude)

  init(value: Double) {
    self.value = value
  }
  
  /// Reports whether this ChordAngle is infinite.
  func isInfinite() -> Bool {
    return value == Double.greatestFiniteMagnitude
  }

  /// Reports whether this ChordAngle is one of the special cases.
  func isSpecial() -> Bool {
    return value < 0.0 || isInfinite()
  }

}
