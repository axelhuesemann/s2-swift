//
//  R1Interval.swift
//  s2-swift
//

import Foundation

/// Represents a closed interval on ℝ.
/// Zero-length intervals (where lo == hi) represent single points.
/// If lo > hi then the interval is empty.
public struct R1Interval: IntervalType {
  
  public typealias Member = Double
  public typealias Distance = Double

  //
  let lo: Double
  let hi: Double
  
  // MARK: inits / factory
  
  public init(lo: Double, hi: Double) {
    self.lo = lo
    self.hi = hi
  }
  
  // Convenience method to construct an interval containing a single point.
  public init(point: Double)  {
    self.lo = point
    self.hi = point
  }
  
  // Convenience method to construct the minimal interval containing the two
  // given points. This is equivalent to starting with an empty interval and
  // calling AddPoint() twice, but it is more efficient.
  public init(p0: Double, p1: Double) {
    if p0 <= p1 {
      self.init(lo: p0, hi: p1)
    } else {
      self.init(lo: p1, hi: p0)
    }
  }
  
  public static let empty = R1Interval(lo: 1.0, hi: 0.0)
  public static let full = R1Interval(lo: -Double.greatestFiniteMagnitude, hi: Double.greatestFiniteMagnitude)

  // MARK: tests
  
  /// Reports whether the interval is empty.
  public var isEmpty: Bool {
    return lo > hi
  }
  
  /// Reports whether the interval is full.
  public var isFull: Bool {
    return lo == -Double.greatestFiniteMagnitude && hi == Double.greatestFiniteMagnitude
  }
  
  // MARK: contains/ intersects

  /// Returns true iff the interval contains p.
  public func contains(_ point: Double) -> Bool {
    return lo <= point && point <= hi
  }
  
  /// Returns true iff the interval contains other.
  public func contains(_ interval: R1Interval) -> Bool {
    if interval.isEmpty {
      return true
    }
    return lo <= interval.lo && interval.hi <= hi
  }
  
  /// Returns true iff the the interval strictly contains p.
  public func interiorContains(_ point: Double) -> Bool {
    return lo < point && point < hi
  }
  
  /// Returns true iff the interval strictly contains other.
  public func interiorContains(_ interval: R1Interval) -> Bool {
    if interval.isEmpty {
      return true
    }
    return lo < interval.lo && interval.hi < hi
  }
  
  /// Returns true iff the interval contains any points in common with other.
  public func intersects(_ interval: R1Interval) -> Bool {
    if lo <= interval.lo {
      return interval.lo <= hi && interval.lo <= interval.hi // interval.lo ∈ i and interval is not empty
    }
    return lo <= interval.hi && lo <= hi // lo ∈ interval and i is not empty
  }
  
  /// Returns true iff the interior of the interval contains any points in common with other, including the latter's boundary.
  public func interiorIntersects(_ interval: R1Interval) -> Bool {
    return interval.lo < hi && lo < interval.hi && lo < hi && interval.lo <= hi
  }
  
  // MARK: computed members
  
  /// Returns the midpoint of the interval.
  // It is undefined for empty intervals.
  public var center: Double {
    return 0.5 * (lo + hi)
  }
  
  /// Returns the length of the interval.
  /// The length of an empty interval is negative.
  public var length: Double {
    return hi - lo
  }
  
  // MARK: arithmetic
  
  /// Returns the interval containing all points common to i and j.
  public func intersection(_ interval: R1Interval) -> R1Interval {
    // Empty intervals do not need to be special-cased.
    return R1Interval(lo: max(lo, interval.lo), hi: min(hi, interval.hi))
  }
  
  /// Returns the smallest interval that contains this interval and the given interval.
  public func union(_ interval: R1Interval) -> R1Interval {
    if isEmpty {
      return interval
    }
    if interval.isEmpty {
      return self
    }
    return R1Interval(lo: min(lo, interval.lo), hi: max(hi, interval.hi))
  }
  
  /// Returns the interval expanded so that it contains the given point.
  public func add(_ point: Double) -> R1Interval {
    if isEmpty {
      return R1Interval(lo: point, hi: point)
    } else if point < lo {
      return R1Interval(lo: point, hi: hi)
    } else if point > hi {
      return R1Interval(lo: lo, hi: point)
    }
    return self
  }
  
  /// Returns an interval that has been expanded on each side by margin.
  /// If margin is negative, then the function shrinks the interval on
  /// each side by margin instead. The resulting interval may be empty. Any
  /// expansion of an empty interval remains empty.
  public func expanded(_ margin: Double) -> R1Interval {
    if isEmpty {
      return self
    }
    return R1Interval(lo: lo - margin, hi: hi + margin)
  }
  
  /// Returns the closest point in the interval to the given point.
  /// The interval must be non-empty.
  func clamp(_ point: Double) -> Double {
    return max(lo, min(hi, point))
  }
  
}

extension R1Interval: Equatable, CustomStringConvertible, Approximatable {

  /// A small number that represents a reasonable level of noise between two
  /// values that can be considered to be equal.
  static let epsilon = 1e-14
  
  public static func ==(lhs: R1Interval, rhs: R1Interval) -> Bool {
    return lhs.lo == rhs.lo && lhs.hi == rhs.hi || (lhs.isEmpty && rhs.isEmpty)
  }
  
  public var description: String {
    let l = String(format: "%.7f", lo * toDegrees)
    let h = String(format: "%.7f", hi * toDegrees)
    return "[\(l), \(h)]"
  }
  
  /// Reports whether the interval can be transformed into the
  /// given interval by moving each endpoint a small distance.
  /// The empty interval is considered to be positioned arbitrarily on the
  /// real line, so any interval with a small enough length will match
  /// the empty interval.
  public func approxEquals(_ other: R1Interval) -> Bool {
    if isEmpty {
      return other.length <= 2 * R1Interval.epsilon
    }
    if other.isEmpty {
      return length <= 2 * R1Interval.epsilon
    }
    return abs(other.lo - lo) <= R1Interval.epsilon && abs(other.hi - hi) <= R1Interval.epsilon
  }
  
}
