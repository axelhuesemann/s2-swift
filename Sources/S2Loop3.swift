//
//  S2Loop3.swift
//  Sphere2Go
//
//  Created by Axel Huesemann on 3/21/17.
//  Copyright © 2017 Axel Huesemann. All rights reserved.
//

import Foundation

struct S2AreaCentroid {
  let area: Double
  let centroid: S2Point?
}

/**
 *
 * An S2Loop represents a simple spherical polygon. It consists of a single
 * chain of vertices where the first vertex is implicitly connected to the last.
 * All loops are defined to have a CCW orientation, i.e. the interior of the
 * polygon is on the left side of the edges. This implies that a clockwise loop
 * enclosing a small area is interpreted to be a CCW loop enclosing a very large
 * area.
 *
 *  Loops are not allowed to have any duplicate vertices (whether adjacent or
 * not), and non-adjacent edges are not allowed to intersect. Loops must have at
 * least 3 vertices. Although these restrictions are not enforced in optimized
 * code, you may get unexpected results if they are violated.
 *
 *  Point containment is defined such that if the sphere is subdivided into
 * faces (loops), every point is contained by exactly one face. This implies
 * that loops do not necessarily contain all (or any) of their vertices An
 * S2Rect represents a latitude-longitude rectangle. It is capable of
 * representing the empty and full rectangles as well as single points.
 *
 */

extension S2Loop: S2Region, Comparable, CustomStringConvertible {

  var description: String {
    let verticess = vertices.map { $0.description }
    return "S1Loop \(vertices.count) points. [\(verticess.joined(separator: " "))]"
  }
  
  static func ==(lhs: S2Loop, rhs: S2Loop) -> Bool { return false }
  static func <(lhs: S2Loop, rhs: S2Loop) -> Bool { return false }
  
}

class S2Loop {
  /**
   * Max angle that intersections can be off by and yet still be considered
   * colinear.
   */
  let MaxIntersectionError = 1e-15
  
  /**
   * Edge index used for performance-critical operations. For example,
   * contains() can determine whether a point is inside a loop in nearly
   * constant time, whereas without an edge index it is forced to compare the
   * query point against every edge in the loop.
   */
  var numVertices: Int
  var vertices: [S2Point]
  
  /**
   * The index (into "vertices") of the vertex that comes first in the total
   * ordering of all vertices in this loop.
   */
  
  var bound: S2Rect
  var depth: Int
  var firstLogicalVertex: Int
  var originInside: Bool
  var index: S2EdgeIndex?
  var vertexToIndex: [S2Point: Int]?
  
  // MARK: inits
  
  /**
   * Initialize a loop connecting the given vertices. The last vertex is
   * implicitly connected to the first. All points should be unit length. Loops
   * must have at least 3 vertices.
   *
   * @param vertices
   */
  init(vertices: [S2Point]) {
    self.vertices = vertices
    numVertices = vertices.count
    bound = S2Rect.full
    depth = 0
    // if (debugMode) {
    //  assert (isValid(vertices, DEFAULT_MAX_ADJACENT))
    // }
    // initOrigin() must be called before InitBound() because the latter
    // function expects contains() to work properly.
    InitOrigin()
    InitBound()
    InitFirstLogicalVertex()
  }
  
  /**
   * Initialize a loop corresponding to the given cell.
   */
  convenience init(cell: Cell) {
    self.init(cell: cell, bound: cell.rectBound())
  }
  
  /**
   * Like the constructor above, but assumes that the cell's bounding rectangle
   * has been precomputed.
   *
   * @param cell
   * @param bound
   */
  init(cell: Cell, bound: S2Rect) {
    self.bound = bound
    numVertices = 4
    vertices = (0..<4).map { cell.vertex($0) }
    vertexToIndex = nil
    index = nil
    depth = 0
    InitOrigin()
    InitFirstLogicalVertex()
  }
  
  // Copy constructor.
  init(src: S2Loop) {
    numVertices = src.numVertices
    vertices = src.vertices
    vertexToIndex = src.vertexToIndex
    index = src.index
    firstLogicalVertex = src.firstLogicalVertex
    bound = src.rectBound()
    originInside = src.originInside
    depth = src.depth
  }
  
  // MARK: tests
  
  // Return true if this loop represents a hole in its containing polygon.
  var IsHole: Bool { return (depth & 1) != 0 }
  
  // The sign of a loop is -1 if the loop represents a hole in its containing polygon, and +1 otherwise.
  var Sign: Int { return IsHole ? -1 : 1 }
  
  var IsNormalized: Bool {
    // We allow a bit of error so that exact hemispheres are considered normalized.
    return Area <= 2 * .pi + 1e-14
  }
  
  var AreaAndCentroid: S2AreaCentroid { return GetAreaCentroid(doCentroid: true) }
  
  // Return the area of the polygon interior, i.e. the region on the left side of an odd number of loops. The return value is between 0 and 4*Pi.
  var Area: Double { return GetAreaCentroid(doCentroid: false).area }
  
  // Return the true centroid of the polygon multiplied by the area of the
  // polygon (see {@link S2} for details on centroids). Note that the centroid
  //may not be contained by the polygon.
  var Centroid: S2Point? { return GetAreaCentroid(doCentroid: true).centroid }

  var IsValid: Bool {
    if (numVertices < 3) {
      NSLog("Degenerate loop")
      return false
    }
    // All vertices must be unit length.
    for i in 0..<numVertices {
      if !vertex(i).isUnit {
        NSLog("Vertex \(i) is not unit length")
        return false
      }
    }
    // Loops are not allowed to have any duplicate vertices.
    var vmap = [S2Point: Int]()
    for i in 0..<numVertices {
      var key = vertex(i)
      if let prevIndex = vmap[key] {
        NSLog("Duplicate vertices: \(prevIndex) and \(i)")
        // update always
        vmap[key] = i
        return false
      }
      vmap[key] = i
    }
    // Non-adjacent edges are not allowed to intersect.
    // var crosses = false
    var it = GetEdgeIterator(expectedQueries: numVertices)
    for a1 in 0..<numVertices {
      var a2 = (a1 + 1) % numVertices
      var crosser = EdgeCrosser(a: vertex(a1), b: vertex(a2), c: vertex(0))
      var previousIndex = -2
      it.GetCandidates(a: vertex(a1), b: vertex(a2))
      for b1 in it { // it.GetCandidates(vertex(a1), vertex(a2)); it.HasNext; it.Next()) {
        //var b1 = it.Index
        var b2 = (b1 + 1)%numVertices
        // If either 'a' index equals either 'b' index, then these two edges
        // share a vertex. If a1==b1 then it must be the case that a2==b2, e.g.
        // the two edges are the same. In that case, we skip the test, since we
        // don't want to test an edge against itself. If a1==b2 or b1==a2 then
        // we have one edge ending at the start of the other, or in other words,
        // the edges share a vertex -- and in S2 space, where edges are always
        // great circle segments on a sphere, edges can only intersect at most
        // once, so we don't need to do further checks in that case either.
        if (a1 != b2 && a2 != b1 && a1 != b1) {
          // WORKAROUND(shakusa, ericv): S2.robustCCW() currently
          // requires arbitrary-precision arithmetic to be truly robust. That
          // means it can give the wrong answers in cases where we are trying
          // to determine edge intersections. The workaround is to ignore
          // intersections between edge pairs where all four points are
          // nearly colinear.
          var abc = S2.Angle(Vertex(a1), Vertex(a2), Vertex(b1))
          var abcNearlyLinear = abc.approxEquals(0.0, S2Loop.maxIntersectionError) ||
            abc.pproxEquals(.pi, MaxIntersectionError)
          var abd = S2.Angle(Vertex(a1), Vertex(a2), Vertex(b2))
          var abdNearlyLinear = abd.approxEquals(0.0, S2Loop.maxIntersectionError) ||
            abd.approxEquals(.pi, MaxIntersectionError)
          if (abcNearlyLinear && abdNearlyLinear) {
            continue
          }
          if (previousIndex != b1) {
            crosser.RestartAt(Vertex(b1))
          }
          // Beware, this may return the loop is valid if there is a
          // "vertex crossing".
          // TODO(user): Fix that.
          var crosses = crosser.RobustCrossing(Vertex(b2)) > 0
          previousIndex = b2
          if (crosses) {
            NSLog("Edges " + a1 + " and " + b1 + " cross")
            NSLog("Edge locations in degrees: " + "{0}-{1} and {2}-{3}",
                  S2LatLng(Vertex(a1)).ToStringDegrees(),
                  S2LatLng(Vertex(a2)).ToStringDegrees(),
                  S2LatLng(Vertex(b1)).ToStringDegrees(),
                  S2LatLng(Vertex(b2)).ToStringDegrees())
            return false
          }
        }
      }
    }
    return true
  }

  static let maxIntersectionError = 1e-15

  func CompareTo(other: S2Loop) -> Int {
    if numVertices != other.numVertices {
      return numVertices - other.numVertices
    }
    // Compare the two loops' vertices, starting with each loop's
    // firstLogicalVertex. This allows us to always catch cases where logically
    // identical loops have different vertex orderings (e.g. ABCD and BCDA).
    var maxVertices = numVertices
    var iThis = firstLogicalVertex
    var iOther = other.firstLogicalVertex
    for i in 0..<maxVertices {
      iThis += 1
      iOther += 1
      var compare = vertex(iThis).CompareTo(other.Vertex(iOther))
      if (compare != 0) {
        return compare
      }
    }
    return 0
  }

  func capBound() -> S2Cap { return bound.capBound() }

  /** Return a bounding latitude-longitude rectangle. */
  
  func rectBound() -> S2Rect { return bound }

  /**
   * If this method returns true, the region completely contains the given cell.
   * Otherwise, either the region does not contain the cell or the containment
   * relationship could not be determined.
   */
  func contains(_ cell: Cell) -> Bool {
    // It is faster to construct a bounding rectangle for an S2Cell than for
    // a general polygon. A future optimization could also take advantage of
    // the fact than an S2Cell is convex.
    let cellBound = cell.rectBound()
    if !bound.contains(cellBound) {
      return false
    }
    let cellLoop = S2Loop(cell: cell, bound: cellBound)
    return contains(cellLoop)
  }
  
  /**
   * If this method returns false, the region does not intersect the given cell.
   * Otherwise, either region intersects the cell, or the intersection
   * relationship could not be determined.
   */
  
  func intersects(_ cell: Cell) -> Bool {
    // It is faster to construct a bounding rectangle for an S2Cell than for
    // a general polygon. A future optimization could also take advantage of
    // the fact than an S2Cell is convex.
    let cellBound = cell.rectBound()
    if !bound.intersects(cellBound) {
      return false
    }
    return S2Loop(cell: cell, bound: cellBound).intersects(self)
  }
  
  /**
   * The depth of a loop is defined as its nesting level within its containing
   * polygon. "Outer shell" loops have depth 0, holes within those loops have
   * depth 1, shells within those holes have depth 2, etc. This field is only
   * used by the S2Polygon implementation.
   *
   * @param depth
   */
  
  /**
   * For convenience, we make two entire copies of the vertex list available:
   * vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == numVertices().
   */
  func vertex(_ i: Int) -> S2Point {
    let index = i >= vertices.count ? i - vertices.count : i
    return vertices[index]
  }
  
  /**
   * Comparator (needed by Comparable interface)
   */
  
  /**
   * Calculates firstLogicalVertex, the vertex in this loop that comes first in
   * a total ordering of all vertices (by way of S2Point's compareTo function).
   */
  
  func InitFirstLogicalVertex() {
    var first = 0
    for i in 1..<numVertices {
      if vertex(i).CompareTo(vertex(first)) < 0 {
        first = i
      }
    }
    firstLogicalVertex = first
  }
  
  /**
   * Return true if the loop area is at most 2*Pi.
   */
  /**
   * Invert the loop if necessary so that the area enclosed by the loop is at
   * most 2*Pi.
   */
  func Normalize() {
    if (!IsNormalized) {
      Invert()
    }
  }
  
  /**
   * Reverse the order of the loop vertices, effectively complementing the
   * region represented by the loop.
   */
  func Invert() {
    let last = numVertices - 1
    for j in 0..<(last - 1) / 2 {
      let i = (last - 1) / 2 - j
      let t = vertices[i]
      vertices[i] = vertices[last - i]
      vertices[last - i] = t
    }
    vertexToIndex = nil
    index = nil
    originInside = !originInside
    if bound.lat.lo > -.pi / 2 && bound.lat.hi < .pi / 2 {
      // The complement of this loop contains both poles.
      bound = S2Rect.full
    } else {
      InitBound()
    }
    InitFirstLogicalVertex()
  }
  
  /**
   * Helper method to get area and optionally centroid.
   */
  func GetAreaCentroid(doCentroid: Bool) -> S2AreaCentroid {
    // Don't crash even if loop is not well-defined.
    if numVertices < 3 {
      return S2AreaCentroid(area: 0.0, centroid: nil)
    }
    // The triangle area calculation becomes numerically unstable as the length
    // of any edge approaches 180 degrees. However, a loop may contain vertices
    // that are 180 degrees apart and still be valid, e.g. a loop that defines
    // the northern hemisphere using four points. We handle this case by using
    // triangles centered around an origin that is slightly displaced from the
    // first vertex. The amount of displacement is enough to get plenty of
    // accuracy for antipodal points, but small enough so that we still get
    // accurate areas for very tiny triangles.
    //
    // Of course, if the loop contains a point that is exactly antipodal from
    // our slightly displaced vertex, the area will still be unstable, but we
    // expect this case to be very unlikely (i.e. a polygon with two vertices on
    // opposite sides of the Earth with one of them displaced by about 2mm in
    // exactly the right direction). Note that the approximate point resolution
    // using the E7 or S2CellId representation is only about 1cm.
    var origin = vertex(0)
    var axis = (origin.largestAbsComponent + 1) % 3
    var slightlyDisplaced = origin[axis] + .e * 1e-10
    origin = S2Point((axis == 0) ? slightlyDisplaced : origin.X, (axis == 1) ? slightlyDisplaced : origin.Y, (axis == 2) ? slightlyDisplaced : origin.Z)
//    origin = origin.normalized()
    var areaSum = 0.0
    var centroidSum = S2Point(x: 0, y: 0, z: 0)
    for i in 1...numVertices {
      areaSum += S2Point.signedArea(origin, vertex(i - 1), vertex(i))
      if (doCentroid) {
        // The true centroid is already premultiplied by the triangle area.
        var trueCentroid = S2Point.trueCentroid(origin, vertex(i - 1), vertex(i))
        centroidSum = centroidSum.v.add(trueCentroid).s2
      }
    }
    // The calculated area at this point should be between -4*Pi and 4*Pi,
    // although it may be slightly larger or smaller than this due to
    // numerical errors.
    // assert (Math.abs(areaSum) <= 4 * S2.M_PI + 1e-12)
    if areaSum < 0 {
      // If the area is negative, we have computed the area to the right of the
      // loop. The area to the left is 4*Pi - (-area). Amazingly, the centroid
      // does not need to be changed, since it is the negative of the integral
      // of position over the region to the right of the loop. This is the same
      // as the integral of position over the region to the left of the loop,
      // since the integral of position over the entire sphere is (0, 0, 0).
      areaSum += 4 * .pi
    }
    // The loop's sign() does not affect the return result and should be taken
    // into account by the caller.
    var centroid: S2Point?
    if (doCentroid) {
      centroid = centroidSum
    }
    return S2AreaCentroid(area: areaSum, centroid: centroid)
  }
  
  /**
   * Return the area of the loop interior, i.e. the region on the left side of
   * the loop. The return value is between 0 and 4*Pi and the true centroid of
   * the loop multiplied by the area of the loop (see S2.java for details on
   * centroids). Note that the centroid may not be contained by the loop.
   */
  
  // The following are the possible relationships between two loops A and B:
  //
  // (1) A and B do not intersect.
  // (2) A contains B.
  // (3) B contains A.
  // (4) The boundaries of A and B cross (i.e. the boundary of A
  // intersects the interior and exterior of B and vice versa).
  // (5) (A union B) is the entire sphere (i.e. A contains the
  // complement of B and vice versa).
  //
  // More than one of these may be true at the same time, for example if
  // A == B or A == Complement(B).
  
  /**
   * Return true if the region contained by this loop is a superset of the
   * region contained by the given other loop.
   */
  func contains(_ b: S2Loop) -> Bool {
    // For this loop A to contains the given loop B, all of the following must
    // be true:
    //
    // (1) There are no edge crossings between A and B except at vertices.
    //
    // (2) At every vertex that is shared between A and B, the local edge
    // ordering implies that A contains B.
    //
    // (3) If there are no shared vertices, then A must contain a vertex of B
    // and B must not contain a vertex of A. (An arbitrary vertex may be
    // chosen in each case.)
    //
    // The second part of (3) is necessary to detect the case of two loops whose
    // union is the entire sphere, i.e. two loops that contains each other's
    // boundaries but not each other's interiors.
    
    if !bound.contains(b.rectBound()) {
      return false
    }
    
    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B. Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick rejection test.
    if !contains(b.vertex(0)) && FindVertex(p: b.vertex(0)) < 0 {
      return false
    }
    
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    if CheckEdgeCrossings(b: b, relation: WedgeContains()) <= 0 {
      return false
    }
    
    // At this point we know that the boundaries of A and B do not intersect,
    // and that A contains a vertex of B. However we still need to check for
    // the case mentioned above, where (A union B) is the entire sphere.
    // Normally this check is very cheap due to the bounding box precondition.
    if bound.union(b.rectBound()).isFull {
      if b.contains(vertex(0)) && b.FindVertex(p: vertex(0)) < 0 {
        return false
      }
    }
    return true
  }
  
  /**
   * Return true if the region contained by this loop intersects the region
   * contained by the given other loop.
   */
  func intersects(_ b: S2Loop) -> Bool {
    // a->intersects(b) if and only if !a->Complement()->contains(b).
    // This code is similar to contains(), but is optimized for the case
    // where both loops enclose less than half of the sphere.
    if !bound.intersects(b.rectBound()) {
      return false
    }
    // Normalize the arguments so that B has a smaller longitude span than A.
    // This makes intersection tests much more efficient in the case where
    // longitude pruning is used (see CheckEdgeCrossings).
    if b.rectBound().lng.length > bound.lng.length {
      return b.intersects(self)
    }
    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B. Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick acceptance test.
    if contains(b.vertex(0)) && FindVertex(p: b.vertex(0)) < 0 {
      return true
    }
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    if CheckEdgeCrossings(b: b, relation: WedgeIntersects()) < 0 {
      return true
    }
    // We know that A does not contain a vertex of B, and that there are no edge
    // crossings. Therefore the only way that A can intersect B is if B
    // entirely contains A. We can check this by testing whether B contains an
    // arbitrary non-shared vertex of A. Note that this check is cheap because
    // of the bounding box precondition and the fact that we normalized the
    // arguments so that A's longitude span is at least as long as B's.
    if b.rectBound().contains(bound) {
      if b.contains(vertex(0)) && b.FindVertex(p: vertex(0)) < 0 {
        return true
      }
    }
    return false
  }
  
  /**
   * Given two loops of a polygon, return true if A contains B. This version of
   * contains() is much cheaper since it does not need to check whether the
   * boundaries of the two loops cross.
   */
  func ContainsNested(b: S2Loop) -> Bool {
    if !bound.contains(b.rectBound()) {
      return false
    }
    // We are given that A and B do not share any edges, and that either one
    // loop contains the other or they do not intersect.
    let m = FindVertex(p: b.vertex(1))
    if (m < 0) {
      // Since b->vertex(1) is not shared, we can check whether A contains it.
      return contains(b.vertex(1))
    }
    // Check whether the edge order around b->vertex(1) is compatible with
    // A containin B.
    return WedgeContains().test(a0: vertex(m - 1), ab1: vertex(m), a2: vertex(m + 1), b0: b.vertex(0), b2: b.vertex(2)) > 0
  }
  
  /**
   * Return +1 if A contains B (i.e. the interior of B is a subset of the
   * interior of A), -1 if the boundaries of A and B cross, and 0 otherwise.
   * Requires that A does not properly contain the complement of B, i.e. A and B
   * do not contain each other's boundaries. This method is used for testing
   * whether multi-loop polygons contain each other.
   */
  func ContainsOrCrosses(b: S2Loop) -> Int {
    // There can be containment or crossing only if the bounds intersect.
    if !bound.intersects(b.rectBound()) {
      return 0
    }
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices. Note that unlike contains() or
    // intersects(), we can't do a point containment test as a shortcut because
    // we need to detect whether there are any edge crossings.
    let result = CheckEdgeCrossings(b: b, relation: WedgeContainsOrCrosses())
    // If there was an edge crossing or a shared vertex, we know the result
    // already. (This is true even if the result is 1, but since we don't
    // bother keeping track of whether a shared vertex was seen, we handle this
    // case below.)
    if (result <= 0) {
      return result
    }
    // At this point we know that the boundaries do not intersect, and we are
    // given that (A union B) is a proper subset of the sphere. Furthermore
    // either A contains B, or there are no shared vertices (due to the check
    // above). So now we just need to distinguish the case where A contains B
    // from the case where B contains A or the two loops are disjoint.
    if !bound.contains(b.rectBound()) {
      return 0
    }
    if !contains(b.vertex(0)) && FindVertex(p: b.vertex(0)) < 0 {
      return 0
    }
    return 1
  }
  
  /**
   * Returns true if two loops have the same boundary except for vertex
   * perturbations. More precisely, the vertices in the two loops must be in the
   * same cyclic order, and corresponding vertex pairs must be separated by no
   * more than maxError. Note: This method mostly useful only for testing
   * purposes.
   */
  
  func BoundaryApproxEquals(b: S2Loop, maxError: Double) -> Bool {
    if numVertices != b.numVertices {
      return false
    }
    var maxVertices = numVertices
    var iThis = firstLogicalVertex
    var iOther = b.firstLogicalVertex
    for i in 0..<maxVertices {
      iThis += 1
      iOther += 1
      if !vertex(iThis).approxEquals(b.Vertex(iOther)) { //, S2Loop.maxError) {
        return false
      }
    }
    return true
  }
  
  // S2Region interface (see {@code S2Region} for details):
  
  /** Return a bounding spherical cap. */
  
  // The point 'p' does not need to be normalized.
  func contains(_ p: S2Point) -> Bool {
    if (!bound.contains(p)) {
      return false
    }
    var inside = originInside
    var origin = S2Point.origin
    var crosser = EdgeCrosser(a: origin, b: p, c: vertices[numVertices - 1])
    // The s2edgeindex library is not optimized yet for long edges,
    // so the tradeoff to using it comes with larger loops.
    if (numVertices < 2000) {
      for i in 0..<numVertices {
        inside ^= crosser.EdgeOrVertexCrossing(vertices[i])
      }
    }
    else {
      var previousIndex = -2
      var it = GetEdgeIterator(expectedQueries: numVertices)
      it.GetCandidates(a: origin, b: p)
      for ai in it { // it.GetCandidates(origin, p); it.HasNext; it.Next()) {
        //var ai = it.Index
        if (previousIndex != ai - 1) {
          crosser.RestartAt(vertices[ai])
        }
        previousIndex = ai
        inside ^= crosser.EdgeOrVertexCrossing(vertex(ai + 1))
      }
    }
    return inside
  }

  /**
   * Returns the shortest distance from a point P to this loop, given as the
   * angle formed between P, the origin and the nearest point on the loop to P.
   * This angle in radians is equivalent to the arclength along the unit sphere.
   */
  func GetDistance(p: S2Point) -> S1Angle {
    var normalized = p.normalized() // not necessary?
    
    // The furthest point from p on the sphere is its antipode, which is an
    // angle of PI radians. This is an upper bound on the angle.
    var minDistance = Double.pi
    for i in 0..<numVertices {
      minDistance = min(minDistance, S2EdgeUtil.GetDistance(normalized, vertex(i), vertex(i + 1)))
    }
    return minDistance
  }
  
  /**
   * Creates an edge index over the vertices, which by itself takes no time.
   * Then the expected number of queries is used to determine whether brute
   * force lookups are likely to be slower than really creating an index, and if
   * so, we do so. Finally an iterator is returned that can be used to perform
   * edge lookups.
   */
  func GetEdgeIterator(expectedQueries: Int) -> DataEdgeIterator {
    if (index == nil) {
      index = AnonS2EdgeIndex(loop: self)
    }
    guard let index = index else { return nil }
    index.PredictAdditionalCalls(n: expectedQueries)
    return DataEdgeIterator(edgeIndex: index)
  }
  
  /** Return true if this loop is valid. */
  
  /**
   * Static version of isValid(), to be used only when an S2Loop instance is not
   * available, but validity of the points must be checked.
   *
   * @return true if the given loop is valid. Creates an instance of S2Loop and
   *         defers this call to {@link #isValid()}.
   */
  static func IsValidLoop(vertices: [S2Point]) -> Bool {
    return S2Loop(vertices: vertices).IsValid
  }
  
  func InitOrigin() {
    // The bounding box does not need to be correct before calling this
    // function, but it must at least contain vertex(1) since we need to
    // do a contains() test on this point below.
    assert(bound.contains(vertex(1)))
    
    // To ensure that every point is contained in exactly one face of a
    // subdivision of the sphere, all containment tests are done by counting the
    // edge crossings starting at a fixed point on the sphere (S2::Origin()).
    // We need to know whether this point is inside or outside of the loop.
    // We do this by first guessing that it is outside, and then seeing whether
    // we get the correct containment result for vertex 1. If the result is
    // incorrect, the origin must be inside the loop.
    //
    // A loop with consecutive vertices A,B,C contains vertex B if and only if
    // the fixed vector R = S2::Ortho(B) is on the left side of the wedge ABC.
    // The test below is written so that B is inside if C=R but not if A=R.
    
    originInside = false; // Initialize before calling contains().
    var v1Inside = S2Point.orderedCCW(vertex(1).ortho().s2, vertex(0), vertex(2), vertex(1))
    if (v1Inside != contains(vertex(1))) {
      originInside = true
    }
  }
  
  func InitBound() {
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices. First, the loop may wrap entirely
    // around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe). Second, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    
    var bounder = RectBounder()
    for i in 0...numVertices {
      bounder.add(point: vertex(i))
    }
    var b = bounder.bound
    // Note that we need to initialize bound with a temporary value since
    // contains() does a bounding rectangle check before doing anything else.
    bound = S2Rect.full
    if contains(S2Point(x: 0, y: 0, z: 1)) {
      let lat = R1Interval(lo: b.lat.lo, hi: .pi / 2.0)
      b = S2Rect(lat: lat, lng: S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.lng().isFull() due to the test above.
    
    if (b.lng.isFull && contains(S2Point(x: 0, y: 0, z: -1))) {
      let lat = R1Interval(lo: -.pi / 2.0, hi: b.lat.hi)
      b = S2Rect(lat: lat, lng: b.lng)
    }
    bound = b
  }
  
  /**
   * Return the index of a vertex at point "p", or -1 if not found. The return
   * value is in the range 1..num_vertices_ if found.
   */
  func FindVertex(p: S2Point) -> Int {
    if vertexToIndex == nil {
      vertexToIndex = [S2Point: Int]()
      for i in 1...numVertices {
        vertexToIndex?[vertex(i)] = i
      }
    }
    return vertexToIndex![p] ?? -1
  }
  
  /**
   * This method encapsulates the common code for loop containment and
   * intersection tests. It is used in three slightly different variations to
   * implement contains(), intersects(), and containsOrCrosses().
   *
   *  In a nutshell, this method checks all the edges of this loop (A) for
   * intersection with all the edges of B. It returns -1 immediately if any edge
   * intersections are found. Otherwise, if there are any shared vertices, it
   * returns the minimum value of the given WedgeRelation for all such vertices
   * (returning immediately if any wedge returns -1). Returns +1 if there are no
   * intersections and no shared vertices.
   */
  func CheckEdgeCrossings(b: S2Loop, relation: WedgeRelation) -> Int {
    var it = GetEdgeIterator(expectedQueries: b.numVertices)
    var result = 1
    // since 'this' usually has many more vertices than 'b', use the index on
    // 'this' and loop over 'b'
    for j in 0..<b.numVertices {
      var crosser = EdgeCrosser(a: b.vertex(j), b: b.vertex(j + 1), c: vertex(0))
      var previousIndex = -2
      it.GetCandidates(a: b.vertex(j), b: b.vertex(j + 1))
      for i in it { 
        //    var i = it.Index
        if (previousIndex != i - 1) {
          crosser.RestartAt(vertex(i))
        }
        previousIndex = i
        var crossing = crosser.RobustCrossing(vertex(i + 1))
        if (crossing < 0) {
          continue
        }
        if (crossing > 0) {
          return -1; // There is a proper edge crossing.
        }
        if (vertex(i + 1).Equals(b.vertex(j + 1))) {
          result = min(result, relation.test(vertex(i), vertex(i + 1), vertex(i + 2), b.vertex(j), b.vertex(j + 2)))
          if (result < 0) {
            return result
          }
        }
      }
    }
    return result
  }

}

protocol WedgeRelation {
  func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int
}

class WedgeContainsOrCrosses: WedgeRelation {
  /**
   * Given two edge chains (see WedgeRelation above), this function returns +1
   * if A contains B, 0 if B contains A or the two wedges do not intersect,
   * and -1 if the edge chains A and B cross each other (i.e. if A intersects
   * both the interior and exterior of the region to the left of B). In
   * degenerate cases where more than one of these conditions is satisfied,
   * the maximum possible result is returned. For example, if A == B then the
   * result is +1.
   */
  func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
    // There are 6 possible edge orderings at a shared vertex (all
    // of these orderings are circular, i.e. abcd == bcda):
    //
    // (1) a2 b2 b0 a0: A contains B
    // (2) a2 a0 b0 b2: B contains A
    // (3) a2 a0 b2 b0: A and B are disjoint
    // (4) a2 b0 a0 b2: A and B intersect in one wedge
    // (5) a2 b2 a0 b0: A and B intersect in one wedge
    // (6) a2 b0 b2 a0: A and B intersect in two wedges
    //
    // In cases (4-6), the boundaries of A and B cross (i.e. the boundary
    // of A intersects the interior and exterior of B and vice versa).
    // Thus we want to distinguish cases (1), (2-3), and (4-6).
    //
    // Note that the vertices may satisfy more than one of the edge
    // orderings above if two or more vertices are the same. The tests
    // below are written so that we take the most favorable
    // interpretation, i.e. preferring (1) over (2-3) over (4-6). In
    // particular note that if orderedCCW(a,b,c,o) returns true, it may be
    // possible that orderedCCW(c,b,a,o) is also true (if a == b or b == c).
    if S2Point.orderedCCW(a0, a2, b2, ab1) {
      // The cases with this vertex ordering are 1, 5, and 6,
      // although case 2 is also possible if a2 == b2.
      if S2Point.orderedCCW(b2, b0, a0, ab1) {
        return 1 // Case 1 (A contains B)
      }
      // We are in case 5 or 6, or case 2 if a2 == b2.
      return a2 == b2 ? 0 : -1 // Case 2 vs. 5,6.
    }
    // We are in case 2, 3, or 4.
    return S2Point.orderedCCW(a0, b0, a2, ab1) ? 0 : -1; // Case 2,3 vs. 4.
  }
}

class WedgeContainsOrIntersects: WedgeRelation {
  /**
   * Given two edge chains (see WedgeRelation above), this function returns +1
   * if A contains B, 0 if A and B are disjoint, and -1 if A intersects but
   * does not contain B.
   */
  func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
    // This is similar to WedgeContainsOrCrosses, except that we want to
    // distinguish cases (1) [A contains B], (3) [A and B are disjoint],
    // and (2,4,5,6) [A intersects but does not contain B].
    if S2Point.orderedCCW(a0, a2, b2, ab1) {
      // We are in case 1, 5, or 6, or case 2 if a2 == b2.
      return S2Point.orderedCCW(b2, b0, a0, ab1) ? 1 : -1 // Case 1 vs. 2,5,6.
    }
    // We are in cases 2, 3, or 4.
    if !S2Point.orderedCCW(a2, b0, b2, ab1) {
      return 0 // Case 3.
    }
    // We are in case 2 or 4, or case 3 if a2 == b0.
    return a2 == b0 ? 0 : -1 // Case 3 vs. 2,4.
  }
}

class WedgeIntersects: WedgeRelation {
  /**
   * Given two edge chains (see WedgeRelation above), this function returns -1
   * if the region to the left of A intersects the region to the left of B,
   * and 0 otherwise. Note that regions are defined such that points along a
   * boundary are contained by one side or the other, not both. So for
   * example, if A,B,C are distinct points ordered CCW around a vertex O, then
   * the wedges BOA, AOC, and COB do not intersect.
   */
  func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
  // For A not to intersect B (where each loop interior is defined to be
  // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
  // Note that it's important to write these conditions as negatives
  // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
  // results when two vertices are the same.
  return S2Point.orderedCCW(a0, b2, b0, ab1) && S2Point.orderedCCW(b0, a2, a0, ab1) ? 0 : -1
  }
}

class WedgeContains: WedgeRelation {
  /**
   * Given two edge chains (see WedgeRelation above), this function returns +1
   * if the region to the left of A contains the region to the left of B, and
   * 0 otherwise.
   */
  func test(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Int {
    // For A to contain B (where each loop interior is defined to be its left
    // side), the CCW edge order around ab1 must be a2 b2 b0 a0. We split
    // this test into two parts that test three vertices each.
    return S2Point.orderedCCW(a2, b2, b0, ab1) && S2Point.orderedCCW(b0, a0, a2, ab1) ? 1 : 0
  }
}


class AnonS2EdgeIndex: S2EdgeIndex {
  let loop: S2Loop
  
  init(loop: S2Loop) {
    self.loop = loop
  }
  
  override func NumEdges() -> Int {
    return loop.numVertices
  }
  
  override func EdgeFrom(index: Int) -> S2Point? {
    return loop.vertex(index)
  }
  
  override func EdgeTo(index: Int) -> S2Point? {
    return loop.vertex(index + 1)
  }
}
