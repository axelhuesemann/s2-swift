//
//  R2Rect.swift
//  s2-swift
//

import Foundation


/// Represents a closed axis-aligned rectangle in the (x,y) plane
public struct R2Rect {
  
  let x: R1Interval
  let y: R1Interval

  // MARK: inits / factory
  
  public init(x: R1Interval, y: R1Interval) {
    self.x = x
    self.y = y
  }
  
  public init(p0: R2Point, p1: R2Point) {
    let x = R1Interval(lo: p0.x, hi: p1.x)
    let y = R1Interval(lo: p0.y, hi: p1.y)
    self.init(x: x, y: y)
  }
  
  /// Constructs a rect that contains the given points.
  init(points: [R2Point]) {
    // Because the default value on interval is 0,0, we need to manually
    // define the interval from the first point passed in as our starting
    // interval, otherwise we end up with the case of passing in
    // R2Point{0.2, 0.3} and getting the starting Rect of {0, 0.2}, {0, 0.3}
    // instead of the Rect {0.2, 0.2}, {0.3, 0.3} which is not correct.
    guard let p0 = points.first else {
      self.init(x: R1Interval.empty, y: R1Interval.empty)
      return
    }
    var x = R1Interval(lo: p0.x, hi: p0.x)
    var y = R1Interval(lo: p0.y, hi: p0.y)
    for p in points.dropFirst() {
      x = x.add(p.x)
      y = y.add(p.y)
    }
    self.init(x: x, y: y)
  }

  /// Constructs a rectangle with the given center and size.
  /// Both dimensions of size must be non-negative.
  public init(center: R2Point, size: R2Point) {
    let x = R1Interval(lo: center.x - size.x/2, hi: center.x + size.x/2)
    let y = R1Interval(lo: center.y - size.y/2, hi: center.y + size.y/2)
    self.init(x: x, y: y)
  }

  /// Constructs the canonical empty rectangle. Use isEmpty to test
  /// for empty rectangles, since they have more than one representation. A Rect
  /// is not the same as the EmptyRect.
  public static let empty = R2Rect(x: R1Interval.empty, y: R1Interval.empty)

}

extension R2Rect {
  
  // MARK: tests
  
  /// Reports whether the rectangle is valid.
  // This requires the width to be empty iff the height is empty.
  var isValid: Bool {
    return x.isEmpty == y.isEmpty
  }

  /// Reports whether the rectangle is empty.
  var isEmpty: Bool {
    return x.isEmpty
  }

  // MARK: computed members
  
  /// Returns all four vertices of the rectangle. Vertices are returned in
  /// CCW direction starting with the lower left corner.
  var vertices: [R2Point] {
    return [
      R2Point(x: x.lo, y: y.lo),
      R2Point(x: x.hi, y: y.lo),
      R2Point(x: x.hi, y: y.hi),
      R2Point(x: x.lo, y: y.hi)]
  }

  /// Returns the center of the rectangle in (x,y)-space
  var center: R2Point {
    return R2Point(x: x.center, y: y.center)
  }

  /// Returns the width and height of this rectangle in (x,y)-space. Empty
  /// rectangles have a negative width and height.
  var size: R2Point {
    return R2Point(x: x.length, y: y.length)
  }

  /// Returns the vertex in direction i along the X-axis (0=left, 1=right) and
  /// direction j along the Y-axis (0=down, 1=up).
  func vertex(i: Int, j: Int) -> R2Point {
    let px = (i == 1) ? x.lo : x.hi
    let py = (j == 1) ? y.lo : y.hi
    return R2Point(x: px, y: py)
  }
  
  /// Returns the low corner of the rect.
  var lo: R2Point {
    return R2Point(x: x.lo, y: y.lo)
  }
  
  /// Returns the high corner of the rect.
  var i: R2Point {
    return R2Point(x: x.hi, y: y.hi)
  }
  
  // MARK: contains/ intersects
  
  /// Reports whether the rectangle contains the given point.
  /// Rectangles are closed regions, i.e. they contain their boundary.
  func contains(_ point: R2Point) -> Bool {
    return x.contains(point.x) && y.contains(point.y)
  }

  /// Returns true iff the given point is contained in the interior
  /// of the region (i.e. the region excluding its boundary).
  func interiorContains(_ point: R2Point) -> Bool {
    return x.interiorContains(point.x) && y.interiorContains(point.y)
  }

  /// Reports whether the rectangle contains the given rectangle.
  func contains(_ rect: R2Rect) -> Bool {
    return x.contains(rect.x) && y.contains(rect.y)
  }

  /// Reports whether the interior of this rectangle contains all of the
  /// points of the given other rectangle (including its boundary).
  func interiorContains(_ rect: R2Rect) -> Bool {
    return x.interiorContains(rect.x) && y.interiorContains(rect.y)
  }

  /// Reports whether this rectangle and the other rectangle have any points in common.
  func intersects(_ rect: R2Rect) -> Bool {
    return x.intersects(rect.x) && y.intersects(rect.y)
  }

  /// Reports whether the interior of this rectangle intersects
  /// any point (including the boundary) of the given other rectangle.
  func interiorIntersects(_ rect: R2Rect) -> Bool {
    return x.interiorIntersects(rect.x) && y.interiorIntersects(rect.y)
  }

  // MARK: arithmetic
  
  /// Expands the rectangle to include the given point. The rectangle is
  // expanded by the minimum amount possible.
  func add(_ point: R2Point) -> R2Rect {
    return R2Rect(x: x.add(point.x), y: y.add(point.y))
  }

  /// Expands the rectangle to include the given rectangle. This is the
  /// same as replacing the rectangle by the union of the two rectangles, but
  /// is more efficient.
  func add(_ rect: R2Rect) -> R2Rect {
    return R2Rect(x: x.union(rect.x), y: y.union(rect.y))
  }

  /// Returns the closest point in the rectangle to the given point.
  /// The rectangle must be non-empty.
  func clamp(_ point: R2Point) -> R2Point {
    return R2Point(x: x.clamp(point.x), y: y.clamp(point.y))
  }

  /// Returns a rectangle that has been expanded in the x-direction
  /// by margin.x, and in y-direction by margin.y. If either margin is empty,
  /// then shrink the interval on the corresponding sides instead. The resulting
  /// rectangle may be empty. Any expansion of an empty rectangle remains empty.
  func expanded(_ margin: R2Point) -> R2Rect {
    let xx = x.expanded(margin.x)
    let yy = y.expanded(margin.y)
    if xx.isEmpty || yy.isEmpty {
      return R2Rect.empty
    }
    return R2Rect(x: xx, y: yy)
  }

  // ExpandedByMargin returns a R2Rect that has been expanded by the amount on all sides.
  func expanded(_ margin: Double) -> R2Rect {
    return expanded(R2Point(x: margin, y: margin))
  }

  // Union returns the smallest rectangle containing the union of this rectangle and
  // the given rectangle.
  func union(_ rect: R2Rect) -> R2Rect {
    return R2Rect(x: x.union(rect.x), y: y.union(rect.y))
  }

  // Intersection returns the smallest rectangle containing the intersection of this
  // rectangle and the given rectangle.
  func intersection(_ rect: R2Rect) -> R2Rect {
    let xx = x.intersection(rect.x)
    let yy = y.intersection(rect.y)
    if xx.isEmpty || yy.isEmpty {
      return R2Rect.empty
    }
    return R2Rect(x: xx, y: yy)
  }

}

extension R2Rect: Equatable, CustomStringConvertible, Approximatable {
  
  public static func ==(lhs: R2Rect, rhs: R2Rect) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y
  }
  
  public var description: String {
    return "X \(x), Y \(y)"
//    return "[\(lo), \(hi)]"
  }
  
  // ApproxEquals returns true if the x- and y-intervals of the two rectangles are
  // the same up to the given tolerance.
  public func approxEquals(_ rect: R2Rect) -> Bool {
    return x.approxEquals(rect.x) && y.approxEquals(rect.y)
  }
  
}
