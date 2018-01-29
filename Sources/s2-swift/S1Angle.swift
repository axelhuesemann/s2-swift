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
