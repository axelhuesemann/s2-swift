//
//  R2Point.swift
//  s2-swift
//

import Foundation


/// Represents a point in ℝ²
public struct R2Point {
  
  let x: Double
  let y: Double
  
  public init(x: Double, y: Double) {
    self.x = x
    self.y = y
  }
  
}

extension R2Point: Equatable, CustomStringConvertible, Approximatable {
  
  public static func ==(lhs: R2Point, rhs: R2Point) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y
  }
  
  public var description: String {
    return String(format: "[%f, %f]", x, y)
  }
  
  /// Returns true if the x and y of the two points are
  /// the same up to the given tolerance.
  public func approxEquals(_ point: R2Point) -> Bool {
    return x.approxEquals(point.x) && y.approxEquals(point.y)
  }
  
}
