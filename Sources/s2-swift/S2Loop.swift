//
//  S2Loop.swift
//  s2-swift
//

import Foundation

// package s2
// import math, r1, r3, s1

// Loop represents a simple spherical polygon. It consists of a sequence
// of vertices where the first vertex is implicitly connected to the
// last. All loops are defined to have a CCW orientation, i.e. the interior of
// the loop is on the left side of the edges. This implies that a clockwise
// loop enclosing a small area is interpreted to be a CCW loop enclosing a
// very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not), and non-adjacent edges are not allowed to intersect. Loops must have
// at least 3 vertices (except for the "empty" and "full" loops discussed
// below).
//
// There are two special loops: the "empty" loop contains no points and the
// "full" loop contains all points. These loops do not have any edges, but to
// preserve the invariant that every loop can be represented as a vertex
// chain, they are defined as having exactly one vertex each (see EmptyLoop
// and FullLoop).
// The major differences from the C++ version is pretty much everything.
public struct S2Loop: S2ShapeType, S2RegionType {
  
  //
  let vertices: [S2Point]
  
  // originInside keeps a precomputed value whether this loop contains the origin
  // versus computing from the set of vertices every time.
  var originInside: Bool
  
  // bound is a conservative bound on all points contained by this loop.
  // If l.ContainsPoint(P), then l.bound.ContainsPoint(P).
  var bound: S2Rect
  
  // Since "bound" is not exact, it is possible that a loop A contains
  // another loop B whose bounds are slightly larger. subregionBound
  // has been expanded sufficiently to account for this error, i.e.
  // if A.Contains(B), then A.subregionBound.Contains(B.bound).
  var subregionBound: S2Rect
  
  // MARK: inits / factory
  
  init(vertices: [S2Point], originInside: Bool, bound: S2Rect, subregionBound: S2Rect) {
    self.vertices = vertices
    self.originInside = originInside
    self.bound = bound
    self.subregionBound = subregionBound
  }
  
//  init(vertices: [S2Point], originInside: Bool) {
//    self.vertices = vertices
//    self.originInside = originInside
//    self.bound = S2Rect.empty
//    self.subregionBound = S2Rect.empty
//  }

  // constructs a loop from the given points
  init(points: [S2Point]) {
    var vertices = points
    // create preliminary loop object with empty bounds
//    var l =
    // figure out origin
    let originInside: Bool
    switch vertices.count {
    case 1:
      // This is the special empty or full loop, so the origin depends on if
      // the vertex is in the southern hemisphere or not.
      originInside = vertices[0].z < 0
    case 0, 2:
      // these are incomplete loops
      self.init(vertices: points, originInside: false, bound: S2Rect.empty, subregionBound: S2Rect.empty)
      return
    default:
      // Point containment testing is done by counting edge crossings starting
      // at a fixed point on the sphere (OriginPoint). We need to know whether
      // the reference point (OriginPoint) is inside or outside the loop before
      // we can construct the S2ShapeIndex. We do this by first guessing that
      // it is outside, and then seeing whether we get the correct containment
      // result for vertex 1. If the result is incorrect, the origin must be
      // inside the loop.
      //
      // A loop with consecutive vertices A,B,C contains vertex B if and only if
      // the fixed vector R = B.Ortho is contained by the wedge ABC. The
      // wedge is closed at A and open at C, i.e. the point B is inside the loop
      // if A = R but not if C = R. This convention is required for compatibility
      // with VertexCrossing. (Note that we can't use OriginPoint
      // as the fixed vector because of the possibility that B == OriginPoint.)
      let ortho = S2Point(raw: vertices[1].v.ortho())
      let v1Inside = S2Point.orderedCCW(ortho, vertices[0], vertices[2], vertices[1])
      originInside = v1Inside != S2Loop.contains(vertices: vertices, originInside: false, point: vertices[1])
    }
//    l = S2Loop(vertices: points, originInside: originInside, bound: S2Rect.empty, subregionBound: S2Rect.empty)
    
    // We *must* call initBound before initIndex, because initBound calls
    // ContainsPoint(s2.Point), and ContainsPoint(s2.Point) does a bounds check whenever the
    // index is not fresh (i.e., the loop has been added to the index but the
    // index has not been updated yet).
    // Check for the special "empty" and "full" loops.
    if vertices.count == 1 { // empty of full
      let bound = originInside ? S2Rect.full : S2Rect.empty
      self.init(vertices: vertices, originInside: originInside, bound: bound, subregionBound: bound)
      return
    }
    
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices. First, the maximal latitude may be
    // attained along the interior of an edge. Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe). Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    var bounder = RectBounder()
    for p in vertices {
      bounder.add(point: p)
    }
    bounder.add(point: vertices[0])
    var b = bounder.rectBound()
    
    if S2Loop.contains(vertices: vertices, originInside: originInside, point: S2Point(x: 0, y: 0, z: 1)) {
      b = S2Rect(lat: R1Interval(lo: b.lat.lo, hi: .pi / 2.0), lng: S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.Lng.isFull due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.Lng.isFull.
    if b.lng.isFull && S2Loop.contains(vertices: vertices, originInside: originInside, point: S2Point(x: 0, y: 0, z: -1)) {
      b = S2Rect(lat: R1Interval(lo: -.pi / 2.0, hi: b.lat.hi), lng: S1Interval.full)
    }
    self.init(vertices: vertices, originInside: originInside, bound: b, subregionBound: b.expandForSubregions())
    // TODO: Depends on s2shapeindex being implemented.
    // l.initIndex()
  }

  // returns true if the loop contains the point
  static func contains(vertices: [S2Point], originInside: Bool, point: S2Point) -> Bool {
    // TODO: Move to bruteForceContains and update with ShapeIndex when available.
    // Empty and full loops don't need a special case, but invalid loops with
    // zero vertices do, so we might as well handle them all at once.
    if vertices.count < 3 {
      return originInside
    }
    let origin = S2Point.origin
    var inside = originInside
    var crosser = EdgeCrosser(a: origin, b: point, c: vertices[0])
    for i in 1..<vertices.count {
      inside = inside != crosser.edgeOrVertexChainCrossing(vertices[i])
    }
    // Test the closing edge of the loop too.
    inside = inside != crosser.edgeOrVertexChainCrossing(vertices[0])
    return inside
  }
  
  // constructs a loop from the given points
  static func loopFromPoints(_ points: [S2Point]) -> S2Loop {
    var vertices = points
    // create preliminary loop object with empty bounds
    var l = S2Loop(vertices: points, originInside: false, bound: S2Rect.empty, subregionBound: S2Rect.empty)
    // figure out origin
    let originInside: Bool
    switch l.vertices.count {
    case 1:
      // This is the special empty or full loop, so the origin depends on if
      // the vertex is in the southern hemisphere or not.
      originInside = vertices[0].z < 0
    case 0, 2:
      // these are incomplete loops
      return l
    default:
        // Point containment testing is done by counting edge crossings starting
      // at a fixed point on the sphere (OriginPoint). We need to know whether
      // the reference point (OriginPoint) is inside or outside the loop before
      // we can construct the S2ShapeIndex. We do this by first guessing that
      // it is outside, and then seeing whether we get the correct containment
      // result for vertex 1. If the result is incorrect, the origin must be
      // inside the loop.
      //
      // A loop with consecutive vertices A,B,C contains vertex B if and only if
      // the fixed vector R = B.Ortho is contained by the wedge ABC. The
      // wedge is closed at A and open at C, i.e. the point B is inside the loop
      // if A = R but not if C = R. This convention is required for compatibility
      // with VertexCrossing. (Note that we can't use OriginPoint
      // as the fixed vector because of the possibility that B == OriginPoint.)
      let ortho = S2Point(raw: vertices[1].v.ortho())
      let v1Inside = S2Point.orderedCCW(ortho, vertices[0], vertices[2], vertices[1])
      originInside = v1Inside != l.contains(vertices[1])
    }
    l = S2Loop(vertices: points, originInside: originInside, bound: S2Rect.empty, subregionBound: S2Rect.empty)
    
    // We *must* call initBound before initIndex, because initBound calls
    // ContainsPoint(s2.Point), and ContainsPoint(s2.Point) does a bounds check whenever the
    // index is not fresh (i.e., the loop has been added to the index but the
    // index has not been updated yet).
    // Check for the special "empty" and "full" loops.
    if l.isEmptyOrFull() {
      let bound = l.isEmpty ? S2Rect.empty : S2Rect.full
      return S2Loop(vertices: vertices, originInside: originInside, bound: bound, subregionBound: bound)
    }
    
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices. First, the maximal latitude may be
    // attained along the interior of an edge. Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe). Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    var bounder = RectBounder()
    for p in vertices {
      bounder.add(point: p)
    }
    bounder.add(point: vertices[0])
    var b = bounder.rectBound()
    
    if l.contains(S2Point(x: 0, y: 0, z: 1)) {
      b = S2Rect(lat: R1Interval(lo: b.lat.lo, hi: .pi / 2.0), lng: S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.Lng.isFull due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.Lng.isFull.
    if b.lng.isFull && l.contains(S2Point(x: 0, y: 0, z: -1)) {
      b = S2Rect(lat: R1Interval(lo: -.pi / 2.0, hi: b.lat.hi), lng: S1Interval.full)
    }
    return S2Loop(vertices: vertices, originInside: originInside, bound: b, subregionBound: b.expandForSubregions())
    // TODO: Depends on s2shapeindex being implemented.
    // l.initIndex()
  }

  // EmptyLoop returns a special "empty" loop.
  static let empty = S2Loop.loopFromPoints([S2Point(x: 0, y: 0, z: 1)])

  // FullLoop returns a special "full" loop.
  static let full = S2Loop.loopFromPoints([S2Point(x: 0, y: 0, z: -1)])

  // MARK: tests
  
  // ContainsOrigin reports true if this loop contains s2.OriginPoint().
  public func containsOrigin() -> Bool {
    return originInside
  }

  // HasInterior returns true because all loops have an interior.
  public func hasInterior() -> Bool {
    return true
  }

  // IsEmpty reports true if this is the special "empty" loop that contains no points.
  var isEmpty: Bool {
    return isEmptyOrFull() && !containsOrigin()
  }
  
  // IsFull reports true if this is the special "full" loop that contains all points.
  var isFull: Bool {
    return isEmptyOrFull() && containsOrigin()
  }
  
  // isEmptyOrFull reports true if this loop is either the "empty" or "full" special loops.
  func isEmptyOrFull() -> Bool {
    return vertices.count == 1
  }
  
  // MARK: computed members
  
  // NumEdges returns the number of edges in this shape.
  public func numEdges() -> Int {
    if isEmptyOrFull() {
      return 0
    }
    return vertices.count
  }

  // Edge returns the endpoints for the given edge index.
  public func edge(_ i: Int) -> (S2Point, S2Point) {
    return (vertex(i), vertex(i + 1))
  }

  // RectBound returns a tight bounding rectangle. If the loop contains the point,
  // the bound also contains it.
  public func rectBound() -> S2Rect {
    return bound
  }

  // CapBound returns a bounding cap that may have more padding than the corresponding
  // RectBound. The bound is conservative such that if the loop contains a point P,
  // the bound also contains it.
  public func capBound() -> S2Cap {
    return bound.capBound()
  }

  // Vertex returns the vertex for the given index. For convenience, the vertex indices
  // wrap automatically for methods that do index math such as Edge.
  // i.e., Vertex(NumEdges() + n) is the same as Vertex(n).
  func vertex(_ i: Int) -> S2Point {
    return vertices[i % vertices.count]
  }

  // returns true if the loop contains the point
  func contains(_ point: S2Point) -> Bool {
    // TODO: Move to bruteForceContains and update with ShapeIndex when available.
    // Empty and full loops don't need a special case, but invalid loops with
    // zero vertices do, so we might as well handle them all at once.
    if vertices.count < 3 {
      return originInside
    }
    let origin = S2Point.origin
    var inside = originInside
    var crosser = EdgeCrosser(a: origin, b: point, c: vertices[0])
    for i in 1..<vertices.count {
      inside = inside != crosser.edgeOrVertexChainCrossing(vertices[i])
    }
    // Test the closing edge of the loop too.
    inside = inside != crosser.edgeOrVertexChainCrossing(vertices[0])
    return inside
  }

  // reports whether the loop contains the given cell
  public func contains(_ cell: Cell) -> Bool {
    let cellVertices = (0..<4).map { cell.vertex($0) }
    // if the loop does not contain all cell vertices, return false
    for k in 0..<4 {
      if !contains(cellVertices[k]) {
        return false
      }
    }
    // if there are any edge crossing, it is not containing
    for j in 0..<4 {
      var crosser = EdgeCrosser(a: cellVertices[j], b: cellVertices[(j+1)&3], c: vertices[0])
      for i in 1..<vertices.count {
        if crosser.chainCrossingSign(vertices[i]) != .doNotCross {
          // There is a proper crossing, or two vertices were the same.
          return false
        }
      }
    }
    return true
  }
  
  // reports whether the loop intersects the cell
  public func intersects(_ cell: Cell) -> Bool {
    let cellVertices = (0..<4).map { cell.vertex($0) }
    // intersects if the loop contains any cell vertex
    for k in 0..<4 {
      let vertex = cell.vertex(k)
      if contains(vertex) {
        return true
      }
    }
    // intersects if any loop edge crosses with any cell edge
    for j in 0..<4 {
      var crosser = EdgeCrosser(a: cellVertices[j], b: cellVertices[(j+1)&3], c: vertices[0])
      for i in 1..<vertices.count {
        if crosser.chainCrossingSign(vertices[i]) != .doNotCross {
          // There is a proper crossing, or two vertices were the same.
          return true
        }
      }
    }
    // intersects if tightest rect is enclosed in cell
    var bounder = RectBounder()
    for p in vertices {
      bounder.add(point: p)
    }
    bounder.add(point: vertices[0])
    let rect = bounder.rectBound()
    if cell.rectBound().contains(rect) {
      return true
    }
    // otherwise fals
    return false
  }
  
}
