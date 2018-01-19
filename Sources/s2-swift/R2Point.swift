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
  
  func add(_ p: R2Point) -> R2Point {
    return R2Point(x: x + p.x, y: y + p.y)
  }
  
  func sub(_ p: R2Point) -> R2Point {
    return R2Point(x: x - p.x, y: y - p.y)
  }
  
  func mul(_ value: Double) -> R2Point {
    return R2Point(x: x * value, y: y * value)
  }
  
  func ortho() -> R2Point {
    return R2Point(x: -y, y: x)
  }
  
  func dot(_ p: R2Point) -> Double {
    return x * p.x + y * p.y
  }
  
  func cross(_ p: R2Point) -> Double {
    return x * p.x - y * p.y
  }
  
  func norm() -> Double {
    return hypot(x, y)
  }
  
  func normalized() -> R2Point {
    if x == 0 && y == 0 {
      return self
    }
    return mul(1 / norm())
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
