//
//  S1Interval.swift
//  s2-swift
//

import Foundation

/// Represents a closed interval on a unit circle.
/// If lo > hi then the interval is "inverted".
/// Zero-length intervals, where lo == hi, represent single points.
/// The point at (-1, 0) on the unit circle has two valid representations:
///   - [π,π]
///   - [-π,-π], this is normalized to the former
/// There are two special intervals that take advantage of that:
///   - the full interval, [-π,π], and
///   - the empty interval, [π,-π].
public struct S1Interval: IntervalType {
  
  public typealias Member = S1Angle
  public typealias Distance = S1Angle
  
  // MARK: constants
  
  static let epsilon = 1e-14

  // MARK: state 
  
  let lo: S1Angle
  let hi: S1Angle
  
  // MARK: inits / factory
  
  /// Constructs a new interval from without range checks.
  public init(lo: S1Angle, hi: S1Angle) {
    self.lo = lo
    self.hi = hi
  }
  
  /// Constructs a new interval from endpoints.
  /// This function allows inverted intervals to be created.
  /// Both arguments must be in the range [-π,π].
  public init(p1: S1Angle, p2: S1Angle) {
    self.lo = (p1 == -.pi && p2 != .pi) ? .pi : p1
    self.hi = (p2 == -.pi && p1 != .pi) ? .pi : p2
  }
  
  /// Returns the empty interval.
  public static let empty = S1Interval(lo: .pi, hi: -.pi)

  /// Returns the full interval.
  public static let full = S1Interval(lo: -.pi, hi: .pi)

  // MARK: tests
  
  /// Reports whether the interval is valid.
  var isValid: Bool {
    return abs(lo) <= .pi && abs(hi) <= .pi && !(lo == -.pi && hi != .pi) && !(hi == -.pi && lo != .pi)
  }

  /// Reports whether the interval is full.
  public var isFull: Bool {
    return hi - lo == 2 * .pi
  }

  /// Reports whether the interval is empty.
  public var isEmpty: Bool {
    return lo - hi == 2 * .pi
  }

  /// Reports whether the interval is inverted, that is, whether lo > hi.
  var isInverted: Bool {
    return lo > hi
  }

  // MARK: contains / intersects
  
  /// Fast version of contains(point:)
  /// Assumes p ∈ (-π,π].
  func fastContains(_ point: S1Angle) -> Bool {
    if isInverted {
      return (point >= lo || point <= hi) && !isEmpty
    }
    return point >= lo && point <= hi
  }
  
  /// Returns true iff the interval contains p.
  /// Assumes p ∈ [-π,π].
  public func contains(_ point: S1Angle) -> Bool {
    if point == -.pi {
      return fastContains(.pi)
    }
    return fastContains(point)
  }
  
  /// Returns true iff the interior of the interval contains p.
  /// Assumes p ∈ [-π,π].
  public func interiorContains(_ point: S1Angle) -> Bool {
    var point = point
    if point == -.pi {
      point = .pi
    }
    if isInverted {
      return point > lo || point < hi
    }
    return (point > lo && point < hi) || isFull
  }
  
  /// Returns true iff the interval contains other.
  public func contains(_ interval: S1Interval) -> Bool {
    if isInverted {
      if interval.isInverted {
        return interval.lo >= lo && interval.hi <= hi
      }
      return (interval.lo >= lo || interval.hi <= hi) && !isEmpty
    }
    if interval.isInverted {
      return isFull || interval.isEmpty
    }
    return interval.lo >= lo && interval.hi <= hi
  }
  
  /// Returns true iff the interior of the interval contains other.
  public func interiorContains(_ interval: S1Interval) -> Bool {
    if isInverted {
      if interval.isInverted {
        return (interval.lo > lo && interval.hi < hi) || interval.isEmpty
      }
      return interval.lo > lo || interval.hi < hi
    }
    if interval.isInverted {
      return isFull || interval.isEmpty
    }
    return (interval.lo > lo && interval.hi < hi) || isFull
  }
  
  /// Returns true iff the interval contains any points in common with interval.
  public func intersects(_ interval: S1Interval) -> Bool {
    if isEmpty || interval.isEmpty {
      return false
    }
    if isInverted {
      return interval.isInverted || interval.lo <= hi || interval.hi >= lo
    }
    if interval.isInverted {
      return interval.lo <= hi || interval.hi >= lo
    }
    return interval.lo <= hi && interval.hi >= lo
  }
  
  /// Returns true iff the interior of the interval contains any points in common with other, including the latter's boundary.
  public func interiorIntersects(_ interval: S1Interval) -> Bool {
    if isEmpty || interval.isEmpty || lo == hi {
      return false
    }
    if isInverted {
      return interval.isInverted || interval.lo < hi || interval.hi > lo
    }
    if interval.isInverted {
      return interval.lo < hi || interval.hi > lo
    }
    return (interval.lo < hi && interval.hi > lo) || isFull
  }
  
  // MARK: computed members
  
  /// Returns the length of the interval.
  /// The length of an empty interval is negative.
  public var length: S1Angle {
    var l = hi - lo
    if l >= 0.0 {
      return l
    }
    l += 2 * .pi
    if l > 0 {
      return l
    }
    return -1
  }

  /// Returns the midpoint of the interval.
  /// It is undefined for full and empty intervals.
  public var center: S1Angle {
    let c = 0.5 * (lo + hi)
    if !isInverted {
      return c
    }
    if c <= 0 {
      return c + .pi
    }
    return c - .pi
  }
  
  // MARK: arithmetic
  
  /// Computes distance from a to b in [0,2π], in a numerically stable way.
  static func positiveDistance(_ a: S1Angle, _ b: S1Angle) -> S1Angle {
    let d = b - a
    if d >= 0 {
      return d
    }
    return (b + .pi) - (a - .pi)
  }

  /// Returns the smallest interval that contains both the interval and other.
  public func union(_ interval: S1Interval) -> S1Interval {
    if interval.isEmpty {
      return self
    }
    if fastContains(interval.lo) {
      if fastContains(interval.hi) {
        // Either interval ⊂ i, or i ∪ interval is the full interval.
        if contains(interval) {
          return self
        }
        return S1Interval.full
      }
      return S1Interval(lo: lo, hi: interval.hi)
    }
    if fastContains(interval.hi) {
      return S1Interval(lo: interval.lo, hi: hi)
    }
    // Neither endpoint of interval is in  Either i ⊂ other, or i and other are disjoint.
    if isEmpty || interval.fastContains(lo) {
      return interval
    }
    // This is the only hard case where we need to find the closest pair of endpoints.
    if S1Interval.positiveDistance(interval.hi, lo) < S1Interval.positiveDistance(hi, interval.lo) {
      return S1Interval(lo: interval.lo, hi: hi)
    }
    return S1Interval(lo: lo, hi: interval.hi)
  }

  /// Returns the smallest interval that contains the intersection of the interval and other.
  public func intersection(_ interval: S1Interval) -> S1Interval {
    if interval.isEmpty {
      return S1Interval.empty
    }
    if fastContains(interval.lo) {
      if fastContains(interval.hi) {
        // Either other ⊂ i, or i and other intersect twice. Neither are empty.
        // In the first case we want to return i (which is shorter than other).
        // In the second case one of them is inverted, and the smallest interval
        // that covers the two disjoint pieces is the shorter of i and other.
        // We thus want to pick the shorter of i and other in both cases.
        if interval.length < length {
          return interval
        }
        return self
      }
      return S1Interval(lo: interval.lo, hi: hi)
    }
    if fastContains(interval.hi) {
      return S1Interval(lo: lo, hi: interval.hi)
    }
    // Neither endpoint of other is in  Either i ⊂ other, or i and other are disjoint.
    if interval.fastContains(lo) {
      return self
    }
    return S1Interval.empty
  }

  /// Returns the interval expanded by the minimum amount necessary such
  /// that it contains the given point "p" (an angle in the range [-Pi, Pi]).
  public func add(_ point: S1Angle) -> S1Interval {
    var point = point
    if abs(point) > .pi {
      return self
    }
    if point == -.pi {
      point = .pi
    }
    if fastContains(point) {
      return self
    }
    if isEmpty {
      return S1Interval(lo: point, hi: point)
    }
    if S1Interval.positiveDistance(point, lo) < S1Interval.positiveDistance(hi, point) {
      return S1Interval(lo: point, hi: hi)
    }
    return S1Interval(lo: lo, hi: point)
  }

  //Returns an interval that has been expanded on each side by margin.
  /// If margin is negative, then the function shrinks the interval on each side by margin instead.
  /// The resulting interval may be empty or full.
  /// Any expansion (positive or negative) of a full interval remains full,
  /// and any expansion of an empty interval remains empty.
  public func expanded(_ margin: S1Angle) -> S1Interval {
    if margin >= 0 {
      if isEmpty {
        return self
      }
      // Check whether this interval will be full after expansion, allowing
      // for a 1-bit rounding error when computing each endpoint.
      if length + 2 * margin + 2 * S1Interval.epsilon >= 2.0 * .pi {
        return S1Interval.full
      }
    } else {
      if isFull {
        return self
      }
      // Check whether this interval will be empty after expansion, allowing
      // for a 1-bit rounding error when computing each endpoint.
      if length + 2 * margin - 2 * S1Interval.epsilon <= 0 {
        return S1Interval.empty
      }
    }
    let l = (lo - margin).truncatingRemainder(dividingBy: 2.0 * .pi)
    let h = (hi + margin).truncatingRemainder(dividingBy: 2.0 * .pi)
    if l <= -.pi {
      return S1Interval(p1: .pi, p2: h)
    }
    return S1Interval(p1: l, p2: h)
  }

}

extension S1Interval: CustomStringConvertible {
  
  public var description: String {
    let l = String(format: "%.7f", lo * toDegrees)
    let h = String(format: "%.7f", hi * toDegrees)
    return "[\(l), \(h)]"
  }
  
}
