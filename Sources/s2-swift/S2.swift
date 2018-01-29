//
//  S2.swift
//  s2-swift
//

import Foundation

public protocol Approximatable {
  func approxEquals(_ other: Self) -> Bool
}

public protocol IntervalType {
  associatedtype Member
  associatedtype Distance
  init(lo: Member, hi: Member)
  init(p0: Member, p1: Member)
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
public protocol S2Region {
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

/// Ordering of a set of points
public enum S2Direction: Int {
  case clockwise = -1
  case indeterminate = 0
  case counterClockwise = 1
}

extension S2Direction {
  static prefix func -(d: S2Direction) -> S2Direction {
    switch d {
    case .clockwise: return .counterClockwise
    case .indeterminate: return .indeterminate
    case .counterClockwise: return .clockwise
    }
  }
}

/// Returns element closest to x within the range min..max.
func clamp<T: Comparable>(_ x: T, min: T, max: T) -> T {
  if x < min {
    return min
  }
  if x > max {
    return max
  }
  return x
}
