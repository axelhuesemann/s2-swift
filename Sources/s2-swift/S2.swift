//
//  S2.swift
//  s2-swift
//

import Foundation

public protocol Approximatable {
  func approxEquals(_ other: Self) -> Bool
}

extension Double: Approximatable {
  
  public func approxEquals(_ other: Double) -> Bool {
    return abs(self - other) <= 1e-15
  }
  
  public func approxEquals(_ other: Double, _ epsilon: Double) -> Bool {
    return abs(self - other) <= epsilon
  }
  
}

public protocol IntervalType {
  associatedtype Member
  associatedtype Distance
  init(lo: Member, hi: Member)
  init(p1: Member, p2: Member)
  // empty, full
  static var empty: Self { get }
  static var full: Self { get }
  var isEmpty: Bool { get }
  var isFull: Bool { get }
  // contains, intersects
  func contains(_: Member) -> Bool
  func interiorContains(_: Member) -> Bool
  func contains(_: Self) -> Bool
  func interiorContains(_: Self) -> Bool
  func intersects(_: Self) -> Bool
  func interiorIntersects(_: Self) -> Bool
  // computed
  var center: Member { get }
  var length: Distance { get }
  // arithmetic
  func union(_: Self) -> Self
  func intersection(_: Self) -> Self
  func add(_: Member) -> Self
  func expanded(_: Distance) -> Self
  //  func clamp(_: Member) -> Member
}

/// Represents a two-dimensional region on the unit sphere.
/// The purpose of this interface is to allow complex regions to be
/// approximated as simpler regions. The interface is restricted to methods
/// that are useful for computing approximations.
public protocol S2RegionType  {
  /// Returns a bounding spherical cap. This is not guaranteed to be exact.
  func capBound() -> S2Cap
  /// Returns a bounding latitude-longitude rectangle that contains
  /// the region. The bounds are not guaranteed to be tight.
  func rectBound() -> S2Rect
  /// Reports whether the region completely contains the given region.
  /// It returns false if containment could not be determined.
  func contains(_ cell: Cell) -> Bool
  /// Reports whether the region intersects the given cell or
  /// if intersection could not be determined. It returns false if the region
  /// does not intersect.
  func intersects(_ cell: Cell) -> Bool
}

/// Defines an interface for any S2 type that needs to be indexable.
public protocol S2ShapeType {
  /// Returns the number of edges in this shape.
  func numEdges() -> Int
  /// Returns endpoints for the given edge index.
  func edge(_ i: Int) -> (S2Point, S2Point)
  /// Returns true if this shape has an interior.
  /// i.e. the Shape consists of one or more closed non-intersecting loops.
  func hasInterior() -> Bool
  /// Returns true if this shape contains s2.Origin.
  /// Shapes that do not have an interior will return false.
  func containsOrigin() -> Bool
}

/// Ordering of a set of points
public enum Direction: Int {
  case clockwise = -1
  case indeterminate = 0
  case counterClockwise = 1
}

extension Direction {
  static prefix func -(d: Direction) -> Direction {
    switch d {
    case .clockwise: return .counterClockwise
    case .indeterminate: return .indeterminate
    case .counterClockwise: return .clockwise
    }
  }
}

/// Returns number closest to x within the range min..max.
func clamp(_ x: Int, min: Int, max: Int) -> Int {
  if x < min {
    return min
  }
  if x > max {
    return max
  }
  return x
}
