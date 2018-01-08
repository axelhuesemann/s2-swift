//
//  S2Loop2.swift
//  Sphere2Go
//
//  Created by Axel Huesemann on 3/18/17.
//  Copyright © 2017 Axel Huesemann. All rights reserved.
//

import Foundation

//  #include "base/commandlineflags.h"
//  #include "base/logging.h"
//  #include "base/scoped_ptr.h"
//  #include "util/coding/coder.h"
//  #include "s2cap.h"
//  #include "s2cell.h"
//  #include "s2edgeindex.h"

//#include "base/logging.h"
//#include "base/macros.h"
//#include "s2edgeindex.h"
//#include "s2region.h"
//#include "s2latlngrect.h"
//#include "s2edgeutil.h"

// TODO find these
class S2 {
  func turnAngle(_ a: Double, _ b: Double, _ c: Double) -> Double { return 0.0 }
}

let kCurrentEncodingVersionNumber = 1 // uint8

// Indexing structure to efficiently compute intersections.
class S2LoopIndex: S2EdgeIndex {
  let loop: S2Loop
  init(loop: S2Loop) {
    self.loop = loop
  }
  
  // There is no need to overwrite Reset(), as the only data that's
  // used to implement this class is all contained in the loop data.
  // void Reset()

  func edge_from(index: Int) -> S2Point {
    return loop.vertex(index)
  }
  
  func edge_to(index: Int) -> S2Point {
    return loop.vertex(index+1)
  }

  func num_edges() -> Int {
    return loop.num_vertices
  }
  
  func reset() {
    // TODO implement
  }
  
  func predictAdditionalCalls(_ n: Int) {
    // TODO
  }
  
}

// An S2Loop represents a simple spherical polygon.  It consists of a single
// chain of vertices where the first vertex is implicitly connected to the
// last.  All loops are defined to have a CCW orientation, i.e. the interior
// of the polygon is on the left side of the edges.  This implies that a
// clockwise loop enclosing a small area is interpreted to be a CCW loop
// enclosing a very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not), and non-adjacent edges are not allowed to intersect.  Loops must have
// at least 3 vertices.  Although these restrictions are not enforced in
// optimized code, you may get unexpected results if they are violated.
//
// Point containment is defined such that if the sphere is subdivided into
// faces (loops), every point is contained by exactly one face.  This implies
// that loops do not necessarily contain all (or any) of their vertices.
//
// TODO(user): When doing operations on two loops, always create the
// edgeindex for the bigger of the two.  Same for polygons.
class S2Loop: S2Region {
  
  // MARK: data
  
  // We store the vertices in an array rather than a vector because we don't
  // need any STL methods, and computing the number of vertices using size()
  // would be relatively expensive (due to division by sizeof(S2Point) == 24).
  // When DecodeWithinScope is used to initialize the loop, we do not
  // take ownership of the memory for vertices_, and the owns_verticesfield
  // is used to prevent deallocation and overwriting.
  var num_vertices: Int
  var vertices: [S2Point]
  var owns_vertices: Bool
  
  var bound: S2Rect
  var origin_inside: Bool
  // The depth of a loop is defined as its nesting level within its containing
  // polygon.  "Outer shell" loops have depth 0, holes within those loops have
  // depth 1, shells within those holes have depth 2, etc.  This field is only
  // used by the S2Polygon implementation.
  var depth: Int
  

  // Quadtree index structure of this loop's edges.
  var index: S2LoopIndex
  
  // Map for speeding up FindVertex: We will compute a map from vertex to
  // index in the vertex array as soon as there has been enough calls.
  var num_find_vertex_calls: Int
  var vertex_to_index: [S2Point: Int]
  
  // MARK: inits
  
  // Create an empty S2Loop that should be initialized by calling Init() or
  // Decode().
  init() {
    num_vertices = 0
    vertices = []
    owns_vertices = false
    bound = S2Rect.empty
    depth = 0
    index = S2LoopIndex(loop: self)
    num_find_vertex_calls = 0
  }
  
  // Convenience constructor that calls Init() with the given vertices.
  init(vertices: [S2Point]) {
    num_vertices = 0
    self.vertices = []
    owns_vertices = false
    bound = S2Rect.full
    depth = 0
    index = S2LoopIndex(loop: self)
    num_find_vertex_calls = 0
    initialize(vertices: vertices)
  }

  // Initialize a loop corresponding to the given cell.
  init(cell: Cell) {
    bound = cell.rectBound()
    index = S2LoopIndex(loop: self)
    num_find_vertex_calls = 0
    num_vertices = 4
    vertices = (0..<4).map { cell.vertex($0) }
    depth = 0
    owns_vertices = true
    initOrigin()
    initBound()
  }
  
  // Internal constructor used only by Clone() that makes a deep copy of its argument.
  init(loop: S2Loop) {
    num_vertices = loop.num_vertices
    vertices = loop.vertices
    owns_vertices = true
    bound = loop.bound
    origin_inside = loop.origin_inside
    depth = loop.depth
    index = loop.index
    num_find_vertex_calls = 0
  }
  
  // Initialize a loop connecting the given vertices.  The last vertex is
  // implicitly connected to the first.  All points should be unit length.
  // Loops must have at least 3 vertices.
  func initialize(vertices: [S2Point]) {
    resetMutableFields()
    // necessary?
    if owns_vertices { self.vertices.removeAll() }
    num_vertices = vertices.count
    if vertices.isEmpty {
      self.vertices = []
    } else {
      self.vertices = vertices
    }
    owns_vertices = true
    bound = S2Rect.full
    // initOrigin() must be called before InitBound() because the latter
    // function expects Contains() to work properly.
    initOrigin()
    initBound()
  }
  
  // When the loop is modified (the only cae being Invert() at this time),
  // the indexing structures need to be deleted as they become invalid.
  func resetMutableFields() {
    index.reset()
    num_find_vertex_calls = 0
    vertex_to_index.removeAll()
  }

  // Check whether this loop is valid.  Note that in debug mode, validity
  // is checked at loop creation time, so IsValid()
  // should always return true.
  var isValid: Bool {
    // Loops must have at least 3 vertices.
    if num_vertices < 3 {
      NSLog("Degenerate loop")
      return false
    }
    // All vertices must be unit length.
    for i in 0..<num_vertices {
      if S2.IsUnitLength(vertex(i)) {
        NSLog("Vertex \(i) is not unit length")
        return false
      }
    }
    // Loops are not allowed to have any duplicate vertices.
    var vmap = [S2Point: Int] = [:]
    for  i in 0..<num_vertices {
      if !vmap.insert(make_pair(vertex(i), i)).second {
        NSLog("Duplicate vertices: \(vmap[vertex(i)]) and \(i)")
        return false
      }
    }
    // Non-adjacent edges are not allowed to intersect.
    var crosses = false
    index.predictAdditionalCalls(num_vertices)
    let it = S2EdgeIndex.Iterator(index)
    for i in 0..<num_vertices {
      let crosser = EdgeCrosser(&vertex(i), &vertex(i+1), &vertex(0))
      var previous_index = -2
      for (it.GetCandidates(vertex(i), vertex(i+1)); !it.Done(); it.Next()) {
        int ai = it.Index()
        // There is no need to test the same thing twice.  Moreover, two edges
        // that abut at ai+1 will have been tested for equality above.
        if ai > i + 1 {
          if previous_index != ai { crosser.restartAt(&vertex(ai)) }
          // Beware, this may return the loop is valid if there is a
          // "vertex crossing".
          // TODO(user): Fix that.
          crosses = crosser.robustCrossing(&vertex(ai+1)) > 0
          previous_index = ai + 1
          if crosses {
            NSLog("Edges \(i) and \(ai) cross")
            // additional debugging information:
            NLog("Edge locations in degrees: \(S2LatLng(vertex(i)))-\(S2LatLng(vertex(i+1))) and \(S2LatLng(vertex(ai)))-\(S2LatLng(vertex(ai + 1)))")
            break
          }
        }
      }
      if crosses { break }
    }
    return !crosses
  }


  // This parameter should be removed as soon as people stop using the
  // deprecated version of IsValid.
  static let kDefaultMaxAdjacent = 0
  
  // These two versions are deprecated and ignore max_adjacent.
  // DEPRECATED.
  static func isValid(vertices: [S2Point], max_adjacent: Int) -> Bool {
    if vertices.count < 3 { return false }
    let loop = S2Loop(vertices: vertices)
    return loop.isValid
  }

  // DEPRECATED.
  func isValid(max_adjacent: Int) -> Bool {
    return isValid
  }
  
  // Return true if this loop represents a hole in its containing polygon.
  var isHole: Bool { return (depth & 1) != 0 }
  
  // The sign of a loop is -1 if the loop represents a hole in its containing
  // polygon, and +1 otherwise.
  var sign: Int { return isHole ? -1 : 1 }
  
  // For convenience, we make two entire copies of the vertex list available:
  // vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == num_vertices().
  func vertex(_ i: Int) -> S2Point {
    assert(i >= 0)
    assert(i <= 2 * num_vertices)
    let j = i - num_vertices
    return vertices[j >= 0 ? j : i]
  }
  
  // Return true if the loop area is at most 2*Pi.  Degenerate loops are
  // handled consistently with S2::RobustCCW(), i.e., if a loop can be
  // expressed as the union of degenerate or nearly-degenerate CCW triangles,
  // then it will always be considered normalized.
  var isNormalized: Bool {
    // Optimization: if the longitude span is less than 180 degrees, then the
    // loop covers less than half the sphere and is therefore normalized.
    if bound.lng.length < .pi { return true }
    // We allow some error so that hemispheres are always considered normalized.
    // TODO(user): This might not be necessary once S2Polygon is enhanced so
    // that it does not require its input loops to be normalized.
    return getTurningAngle() >= -1e-14
  }
  
  // Invert the loop if necessary so that the area enclosed by the loop is at
  // most 2*Pi.
  func normalize() {
    // TODO make this functional
    if isNormalized { return }
    assert(owns_vertices)
    invert()
    assert(isNormalized)
  }
  
  // Reverse the order of the loop vertices, effectively complementing
  // the region represented by the loop.
  func invert() { // TODO make functional
    assert(owns_vertices)
    resetMutableFields()
    reverse(vertices, vertices + num_vertices())
    origin_inside ^= true
    if (bound.lat.lo > -.pi / 2 && bound.lat.hi < .pi / 2) {
      // The complement of this loop contains both poles.
      bound = S2Rect.full
    } else {
      initBound()
    }
  }
  
  // Return the area of the loop interior, i.e. the region on the left side of
  // the loop.  The return value is between 0 and 4*Pi.  (Note that the return
  // value is not affected by whether this loop is a "hole" or a "shell".)
  func getArea() -> Double {
    let area = getSurfaceIntegral(S2.signedArea)
    // The signed area should be between approximately -4*Pi and 4*Pi.
    assert(fabs(area) <= 4 * .pi + 1e-12)
    if area < 0 {
      // We have computed the negative of the area of the loop exterior.
      area += 4 * .pi
    }
    return max(0.0, min(4 * .pi, area))
  }
  
  // Return the true centroid of the loop multiplied by the area of the loop
  // (see s2.h for details on centroids).  The result is not unit length, so
  // you may want to normalize it.  Also note that in general, the centroid
  // may not be contained by the loop.
  //
  // We prescale by the loop area for two reasons: (1) it is cheaper to
  // compute this way, and (2) it makes it easier to compute the centroid of
  // more complicated shapes (by splitting them into disjoint regions and
  // adding their centroids).
  //
  // Note that the return value is not affected by whether this loop is a
  // "hole" or a "shell".
  func getCentroid() -> S2Point {
    // GetSurfaceIntegral() returns either the integral of position over loop
    // interior, or the negative of the integral of position over the loop
    // exterior.  But these two values are the same (!), because the integral of
    // position over the entire sphere is (0, 0, 0).
    return getSurfaceIntegral(S2.trueCentroid)
  }
  
  // Return the sum of the turning angles at each vertex.  The return value is
  // positive if the loop is counter-clockwise, negative if the loop is
  // clockwise, and zero if the loop is a great circle.  Degenerate and
  // nearly-degenerate loops are handled consistently with S2::RobustCCW().
  // So for example, if a loop has zero area (i.e., it is a very small CCW
  // loop) then the turning angle will always be negative.
  //
  // This quantity is also called the "geodesic curvature" of the loop.
  func getTurningAngle() -> Double {
    // Don't crash even if the loop is not well-defined.
    if num_vertices < 3 { return 0 }
    // To ensure that we get the same result when the loop vertex order is
    // rotated, and that we get the same result with the opposite sign when the
    // vertices are reversed, we need to be careful to add up the individual
    // turn angles in a consistent order.  In general, adding up a set of
    // numbers in a different order can change the sum due to rounding errors.
    let (i0, dir) = getCanonicalFirstVertex()
    var n = num_vertices
    var i = i0
    let angle = S2.turnAngle(vertex((i + n - dir) % n), vertex(i), vertex((i + dir) % n))
    while true {
      n -= 1
      if n <= 0 { break }
      i += dir
      angle += S2.turnAngle(vertex(i - dir), vertex(i), vertex(i + dir))
    }
    return dir * angle
  }
  
  // Return true if the region contained by this loop is a superset of the
  // region contained by the given other loop.
  func contains(_ b: S2Loop) -> Bool {
    // For this loop A to contains the given loop B, all of the following must
    // be true:
    //
    //  (1) There are no edge crossings between A and B except at vertices.
    //
    //  (2) At every vertex that is shared between A and B, the local edge
    //      ordering implies that A contains B.
    //
    //  (3) If there are no shared vertices, then A must contain a vertex of B
    //      and B must not contain a vertex of A.  (An arbitrary vertex may be
    //      chosen in each case.)
    //
    // The second part of (3) is necessary to detect the case of two loops whose
    // union is the entire sphere, i.e. two loops that contains each other's
    // boundaries but not each other's interiors.
    if !bound.contains(b.bound) { return false }
    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B.  Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick rejection test.
    if !contains(b.vertex(0)) && findVertex(b.vertex(0)) < 0 { return false }
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    var p_wedge = ContainsWedgeProcessor()
    if areBoundariesCrossing(b: b, wedge_processor: p_wedge) || p_wedge.doesntContain() {
      return false
    }
    // At this point we know that the boundaries of A and B do not intersect,
    // and that A contains a vertex of B.  However we still need to check for
    // the case mentioned above, where (A union B) is the entire sphere.
    // Normally this check is very cheap due to the bounding box precondition.
    if bound.union(b.bound).isFull {
      if b.contains(vertex(0)) && b.findVertex(vertex(0)) < 0 { return false }
    }
    return true
  }
  
  // Return true if the region contained by this loop intersects the region
  // contained by the given other loop.
  func intersects(_ b: S2Loop) -> Bool {
    // a->Intersects(b) if and only if !a->Complement()->Contains(b).
    // This code is similar to Contains(), but is optimized for the case
    // where both loops enclose less than half of the sphere.
    // The largest of the two loops should be edgeindex'd.
    if b.num_vertices > num_vertices { return intersects(self) }
    if !bound.intersects(b.bound) { return false }
    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B.  Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick acceptance test.
    if contains(b.vertex(0)) && findVertex(b.vertex(0)) < 0 { return true }
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    let p_wedge = IntersectsWedgeProcessor()
    if areBoundariesCrossing(b: b, wedge_processor: p_wedge) || p_wedge.intersects { return true }
    // We know that A does not contain a vertex of B, and that there are no edge
    // crossings.  Therefore the only way that A can intersect B is if B
    // entirely contains A.  We can check this by testing whether B contains an
    // arbitrary non-shared vertex of A.  Note that this check is usually cheap
    // because of the bounding box precondition.
    if b.bound.contains(bound) {
      if b.contains(vertex(0)) && b.findVertex(vertex(0)) < 0 { return true }
    }
    return false
  }
  
  // Given two loops of a polygon (see s2polygon.h for requirements), return
  // true if A contains B.  This version of Contains() is much cheaper since
  // it does not need to check whether the boundaries of the two loops cross.
  func containsNested(loop b: S2Loop) -> Bool {
    if !bound.contains(b.bound) { return false }
    // We are given that A and B do not share any edges, and that either one
    // loop contains the other or they do not intersect.
    let m = findVertex(b.vertex(1))
    if m < 0 {
      // Since b.vertex(1) is not shared, we can check whether A contains it.
      return contains(b.vertex(1))
    }
    // Check whether the edge order around b.vertex(1) is compatible with
    // A containing B.
    return S2EdgeUtil.wedgeContains(vertex(m-1), vertex(m), vertex(m+1), b.vertex(0), b.vertex(2))
  }
  
  // Return +1 if A contains B (i.e. the interior of B is a subset of the
  // interior of A), -1 if the boundaries of A and B cross, and 0 otherwise.
  // Requires that A does not properly contain the complement of B, i.e.
  // A and B do not contain each other's boundaries.  This method is used
  // for testing whether multi-loop polygons contain each other.
  //
  // WARNING: there is a bug in this function - it does not detect loop
  // crossings in certain situations involving shared edges.  CL 2926852 works
  // around this bug for polygon intersection, but it probably effects other
  // tests.  TODO: fix ContainsOrCrosses() and revert CL 2926852.
  func containsOrCrosses(loop b: S2Loop) -> Int {
    // There can be containment or crossing only if the bounds intersect.
    if !bound.intersects(b.bound) { return 0 }
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.  Note that unlike Contains() or
    // Intersects(), we can't do a point containment test as a shortcut because
    // we need to detect whether there are any edge crossings.
    let p_wedge = ContainsOrCrossesProcessor()
    if areBoundariesCrossing(b: b, wedge_processor: p_wedge) { return -1 }
    let result = p_wedge.crossesOrMayContain()
    if result <= 0 { return result }
    // At this point we know that the boundaries do not intersect, and we are
    // given that (A union B) is a proper subset of the sphere.  Furthermore
    // either A contains B, or there are no shared vertices (due to the check
    // above).  So now we just need to distinguish the case where A contains B
    // from the case where B contains A or the two loops are disjoint.
    if !bound.contains(b.bound) { return 0 }
    if !contains(b.vertex(0)) && findVertex(b.vertex(0)) < 0 { return 0 }
    return 1
  }
  
  // Return true if two loops have the same boundary.  This is true if and
  // only if the loops have the same vertices in the same cyclic order.
  // (For testing purposes.)
  func boundaryEquals(loop b: S2Loop) -> Bool {
    if num_vertices != b.num_vertices { return false }
    for offset in 0..<num_vertices {
      if vertex(offset) == b.vertex(0) {
        // There is at most one starting offset since loop vertices are unique.
        for i in 0..<num_vertices {
          if vertex(i + offset) != b.vertex(i) { return false }
        }
        return true
      }
    }
    return false
  }
  
  // Return true if two loops have the same boundary except for vertex
  // perturbations.  More precisely, the vertices in the two loops must be in
  // the same cyclic order, and corresponding vertex pairs must be separated
  // by no more than "max_error".  (For testing purposes.)
  func boundaryApproxEquals(loop: S2Loop, max_error: Double = 1e-15) -> Bool {}
  
  // Return true if the two loop boundaries are within "max_error" of each
  // other along their entire lengths.  The two loops may have different
  // numbers of vertices.  More precisely, this method returns true if the two
  // loops have parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2 such that
  // distance(a(t), b(t)) <= max_error for all t.  You can think of this as
  // testing whether it is possible to drive two cars all the way around the
  // two loops such that no car ever goes backward and the cars are always
  // within "max_error" of each other.  (For testing purposes.)
  func boundaryNear(loop: S2Loop, max_error: Double = 1e-15) -> Bool {}
  
  // This method computes the oriented surface integral of some quantity f(x)
  // over the loop interior, given a function f_tri(A,B,C) that returns the
  // corresponding integral over the spherical triangle ABC.  Here "oriented
  // surface integral" means:
  //
  // (1) f_tri(A,B,C) must be the integral of f if ABC is counterclockwise,
  //     and the integral of -f if ABC is clockwise.
  //
  // (2) The result of this function is *either* the integral of f over the
  //     loop interior, or the integral of (-f) over the loop exterior.
  //
  // Note that there are at least two common situations where it easy to work
  // around property (2) above:
  //
  //  - If the integral of f over the entire sphere is zero, then it doesn't
  //    matter which case is returned because they are always equal.
  //
  //  - If f is non-negative, then it is easy to detect when the integral over
  //    the loop exterior has been returned, and the integral over the loop
  //    interior can be obtained by adding the integral of f over the entire
  //    unit sphere (a constant) to the result.
  //
  // Also requires that the default constructor for T must initialize the
  // value to zero.  (This is true for built-in types such as "double".)
  // TODO  huh?
//  template <class T>
//  T GetSurfaceIntegral(T f_tri(S2Point const&, S2Point const&, S2Point const&))
//  const
  
  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):
  
  // GetRectBound() is guaranteed to return exact results, while GetCapBound()
  // is conservative.
  func clone() -> S2Loop { return S2Loop(loop: self) }
  func getCapBound() -> S2Cap { return bound.capBound() }
  func getRectBound() -> S2Rect { return bound }
  
  func contains(_ cell: Cell) -> Bool {
    // A future optimization could also take advantage of the fact than an S2Cell
    // is convex.
    // It's not necessarily true that bound.Contains(cell.GetRectBound())
    // because S2Cell bounds are slightly conservative.
    if !bound.contains(cell.center) { return false }
    let cell_loop = S2Loop(cell)
    return contains(cell_loop)
  }

  func mayIntersect(_ cell: Cell) -> Bool {
    // It is faster to construct a bounding rectangle for an S2Cell than for
    // a general polygon.  A future optimization could also take advantage of
    // the fact than an S2Cell is convex.
    if !bound.intersects(cell.rectBound()) { return false }
    return S2Loop(cell: cell).intersects(self)
  }
  
  // The point 'p' does not need to be normalized.
  func contains(_ p: S2Point) -> Bool {
    if !bound.contains(p) { return false }
    var inside = origin_inside
    let origin = S2Point.origin
    let crosser = EdgeCrosser(&origin, &p, &vertex(0))
    // The s2edgeindex library is not optimized yet for long edges,
    // so the tradeoff to using it comes later.
    if num_vertices < 2000 {
      for i in 1...num_vertices {
        inside ^= crosser.EdgeOrVertexCrossing(&vertex(i))
      }
      return inside
    }
    let it = S2EdgeIndex.Iterator(index)
    var previous_index = -2
    for (it.GetCandidates(origin, p); !it.Done(); it.Next()) {
      let ai = it.index()
      if previous_index != ai - 1 { crosser.restartAt(&vertex(ai)) }
      previous_index = ai
      inside ^= crosser.EdgeOrVertexCrossing(&vertex(ai+1))
    }
    return inside
  }

  class Coder {
    func put8(_ value: UInt8) {}
    func put32(_ value: UInt32) {}
    func get8() -> UInt8 { return 0 }
    func get32() -> UInt32 { return 0 }
  }
  
  func encode(encoder: Coder) {
//    encoder.ensure(num_vertices * sizeof(*vertices_) + 20);  // sufficient
    encoder.put8(UInt8(kCurrentEncodingVersionNumber))
    encoder.put32(UInt32(num_vertices))
//    encoder.putn(vertices_, sizeof(*vertices_) * num_vertices_)
    encoder.put8(UInt8(origin_inside))
    encoder.put32(UInt32(depth))
////    assert(encoder.avail() >= 0)
//    bound.encode(encoder)
  }
  
  func decode(decoder: Coder) -> Bool {
    return decodeInternal(decoder: decoder, withinScope: false)
  }
  
  func decodeInternal(decoder: Coder, withinScope: Bool) -> Bool {
    let version = Int(decoder.get8())
    if version > kCurrentEncodingVersionNumber { return false }
    num_vertices = Int(decoder.get32())
//    if owns_vertices { vertices = [] } // necessary?
//    if withinScope {
//      vertices = const_cast<S2Point *>(reinterpret_cast<S2Point const*>(decoder.ptr()))
//      decoder.skip(num_vertices * sizeof(vertices[0]))
//      owns_vertices = false
//    } else {
//      vertices = new S2Point[num_vertices]
//      decoder.getn(vertices, num_vertices * sizeof(vertices[0]))
//      owns_vertices = true
//    }
    origin_inside = Bool(decoder.get8())
    depth = Int(decoder.get32())
//    if !bound.decode(decoder) { return false }
    assert(isValid)
    return true // decoder.avail() >= 0
  }
 
  func initOrigin() {
    // The bounding box does not need to be correct before calling this
    // function, but it must at least contain vertex(1) since we need to
    // do a Contains() test on this point below.
    assert(bound.contains(vertex(1)))
    // To ensure that every point is contained in exactly one face of a
    // subdivision of the sphere, all containment tests are done by counting the
    // edge crossings starting at a fixed point on the sphere (S2::Origin()).
    // We need to know whether this point is inside or outside of the loop.
    // We do this by first guessing that it is outside, and then seeing whether
    // we get the correct containment result for vertex 1.  If the result is
    // incorrect, the origin must be inside the loop.
    //
    // A loop with consecutive vertices A,B,C contains vertex B if and only if
    // the fixed vector R = S2::Ortho(B) is on the left side of the wedge ABC.
    // The test below is written so that B is inside if C=R but not if A=R.
    // Initialize before calling contains().
    origin_inside = false
    let v1_inside = S2Point.orderedCCW(S2.ortho(vertex(1)), vertex(0), vertex(2), vertex(1))
    if v1_inside != contains(vertex(1)) {
      origin_inside = true
    }
  }

  func initBound() {
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices.  First, the loop may wrap entirely
    // around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe).  Second, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    var bounder = RectBounder()
    for i in 0...num_vertices {
      bounder.add(point: vertex(i))
    }
    var b = bounder.rectBound()
    // Note that we need to initialize bound with a temporary value since
    // Contains() does a bounding rectangle check before doing anything else.
    bound = S2Rect.full
    if contains(S2Point(x: 0, y: 0, z: 1)) {
      b = S2Rect(R1Interval(b.lat().lo, .pi / 2), S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.lng().is_full() due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.lng().is_full().
    if b.lng().isFull && contains(S2Point(0, 0, -1)) {
      b.mutable_lat().set_lo(-.pi / 2)
    }
    bound = b
  }
  
  // Internal implementation of the Decode and DecodeWithinScope methods above.
  // If within_scope is true, memory is allocated for vertices_ and data
  // is copied from the decoder using memcpy. If it is false, vertices_
  // will point to the memory area inside the decoder, and the field
  // owns_vertices is set to false.
//  bool DecodeInternal(Decoder* const decoder, bool within_scope)
  
  // Internal implementation of the Intersects() method above.
  func intersectsInternal(b: S2Loop) -> Bool {}

  // Return an index "first" and a direction "dir" (either +1 or -1) such that
  // the vertex sequence (first, first+dir, ..., first+(n-1)*dir) does not
  // change when the loop vertex order is rotated or inverted.  This allows
  // the loop vertices to be traversed in a canonical order.  The return
  // values are chosen such that (first, ..., first+n*dir) are in the range
  // [0, 2*n-1] as expected by the vertex() method.
  func getCanonicalFirstVertex() -> (Int, Int) {
    var first = 0
    var n = num_vertices
    for i in 1..<n {
      if vertex(i) < vertex(first) { first = i }
    }
    let dir: Int
    if vertex(first + 1) < vertex(first + n - 1) {
      dir = 1
      // 0 <= first <= n-1, so (first+n*dir) <= 2*n-1.
    } else {
      dir = -1
      first += n
      // n <= first <= 2*n-1, so (first+n*dir) >= 0.
    }
    return (first, dir)
  }

  // Return the index of a vertex at point "p", or -1 if not found.
  // The return value is in the range 1..num_vertices if found.
  func findVertex(_ p: S2Point) -> Int {
    num_find_vertex_calls += 1
    if num_vertices < 10 || num_find_vertex_calls < 20 {
      // Exhaustive search
      for i in 1...num_vertices {
        if vertex(i) == p { return i }
      }
      return -1
    }
    if vertex_to_index.isEmpty {  // We haven't computed it yet.
      for i0 in 0..<num_vertices {
        let i = num_vertices - i0 + 1
        vertex_to_index[vertex(i)] = i
      }
    }
    return vertex_to_index[p] ?? -1
  }
  
  // This method checks all edges of this loop (A) for intersection
  // against all edges of B.  If there is any shared vertex , the
  // wedges centered at this vertex are sent to wedge_processor.
  //
  // Returns true only when the edges intersect in the sense of
  // S2EdgeUtil::RobustCrossing, returns false immediately when the
  // wedge_processor returns true: this means the wedge processor
  // knows the value of the property that the caller wants to compute,
  // and no further inspection is needed.  For instance, if the
  // property is "loops intersect", then a wedge intersection is all
  // it takes to return true.
  //
  // See Contains(), Intersects() and ContainsOrCrosses() for the
  // three uses of this function.
  func areBoundariesCrossing(b: S2Loop, wedge_processor: WedgeProcessor) -> Bool {
    // See the header file for a description of what this method does.
    index.predictAdditionalCalls(b.num_vertices)
    let it = S2EdgeIndex.Iterator(index)
    for j in 0..<b.num_vertices {
      let crosser = EdgeCrosser(b.vertex(j), b.vertex(j+1), b.vertex(0))
      var previous_index = -2
      for (it.GetCandidates(b.vertex(j), b.vertex(j+1)); !it.Done(); it.Next()) {
          let ai = it.Index()
        if previous_index != ai - 1 { crosser.RestartAt(&vertex(ai)) }
          previous_index = ai
          let crossing = crosser.robustCrossing(&vertex(ai + 1))
          if crossing < 0 { continue }
          if crossing > 0 { return true }
          // We only need to check each shared vertex once, so we only
          // consider the case where vertex(i+1) == b.vertex(j+1).
          if vertex(ai + 1) == b.vertex(j + 1) &&
            wedge_processor.processWedge(vertex(ai), vertex(ai + 1), vertex(ai + 2), b.vertex(j), b.vertex(j + 2)) {
            return false
          }
      }
    }
    return false
  }
  
  struct IntPair: Hashable {
    let a: Int
    let b: Int
    var hashValue: Int {
      return a.hashValue + b.hashValue * 53
    }
    static func ==(lhs: IntPair, rhs: IntPair) -> Bool {
      return lhs.a == rhs.a && lhs.b == rhs.b
    }
  }

  static func matchBoundaries(a: S2Loop, b: S2Loop, offset: Int, maxError: Double) -> Bool {
  // The state consists of a pair (i,j).  A state transition consists of
  // incrementing either "i" or "j".  "i" can be incremented only if
  // a(i+1+a_offset) is near the edge from b(j) to b(j+1), and a similar rule
  // applies to "j".  The function returns true iff we can proceed all the way
  // around both loops in this way.
  //
  // Note that when "i" and "j" can both be incremented, sometimes only one
  // choice leads to a solution.  We handle this using a stack and
  // backtracking.  We also keep track of which states have already been
  // explored to avoid duplicating work.
    var pending: [IntPair] = []
    var done = Set<IntPair>()
    pending.append(IntPair(a: 0, b:0))
    while !pending.isEmpty {
      let pair = pending.removeLast()
      let (i, j) = (pair.a, pair.b)
      if i == a.num_vertices && j == b.num_vertices {
        return true
      }
      done.insert(IntPair(a: i, b: j))
      
      // If (i == na && offset == na-1) where na == a.num_vertices(), then
      // then (i+1+offset) overflows the [0, 2*na-1] range allowed by vertex().
      // So we reduce the range if necessary.
      var io = i + offset
      if io >= a.num_vertices { io -= a.num_vertices }
      
      if (i < a.num_vertices && done.count((i+1, j)) == 0 &&
        S2EdgeUtil.getDistance(a.vertex(io + 1), b.vertex(j), b.vertex(j + 1)).radians() <= maxError) {
        pending.append(IntPair(a: i+1, b: j))
      }
      if (j < b.num_vertices && done.count((i, j+1)) == 0 &&
        S2EdgeUtil.getDistance(b.vertex(j+1), a.vertex(io), a.vertex(io + 1)).radians() <= maxError) {
        pending.append(IntPair(a: i, b: j+1))
      }
    }
    return false
  }
  
  func boundaryNear(_ b: S2Loop, maxError: Double) -> Bool {
    for offset in 0..<num_vertices {
      if S2Loop.matchBoundaries(a: self, b: b, offset: offset, maxError: maxError) { return true }
    }
    return false
  }

  // We sum "f_tri" over a collection T of oriented triangles, possibly
  // overlapping.  Let the sign of a triangle be +1 if it is CCW and -1
  // otherwise, and let the sign of a point "x" be the sum of the signs of the
  // triangles containing "x".  Then the collection of triangles T is chosen
  // such that either:
  //
  //  (1) Each point in the loop interior has sign +1, and sign 0 otherwise; or
  //  (2) Each point in the loop exterior has sign -1, and sign 0 otherwise.
  //
  // The triangles basically consist of a "fan" from vertex 0 to every loop
  // edge that does not include vertex 0.  These triangles will always satisfy
  // either (1) or (2).  However, what makes this a bit tricky is that
  // spherical edges become numerically unstable as their length approaches
  // 180 degrees.  Of course there is not much we can do if the loop itself
  // contains such edges, but we would like to make sure that all the triangle
  // edges under our control (i.e., the non-loop edges) are stable.  For
  // example, consider a loop around the equator consisting of four equally
  // spaced points.  This is a well-defined loop, but we cannot just split it
  // into two triangles by connecting vertex 0 to vertex 2.
  //
  // We handle this type of situation by moving the origin of the triangle fan
  // whenever we are about to create an unstable edge.  We choose a new
  // location for the origin such that all relevant edges are stable.  We also
  // create extra triangles with the appropriate orientation so that the sum
  // of the triangle signs is still correct at every point.
  func getSurfaceIntegral<T: Summable>(f_tri: (S2Point, S2Point, S2Point) -> T, zero: T) -> T {
    // The maximum length of an edge for it to be considered numerically stable.
    // The exact value is fairly arbitrary since it depends on the stability of
    // the "f_tri" function.  The value below is quite conservative but could be
    // reduced further if desired.
    let kMaxLength = Double.pi - 1e-5
    // The default constructor for T must initialize the value to zero.
    // (This is true for built-in types such as "double".)
    var sum = zero
    var origin = vertex(0)
    for i in 1..<num_vertices - 1 {
      // Let V_i be vertex(i), let O be the current origin, and let length(A,B)
      // be the length of edge (A,B).  At the start of each loop iteration, the
      // "leading edge" of the triangle fan is (O,V_i), and we want to extend
      // the triangle fan so that the leading edge is (O,V_i+1).
      //
      // Invariants:
      //  1. length(O,V_i) < kMaxLength for all (i > 1).
      //  2. Either O == V_0, or O is approximately perpendicular to V_0.
      //  3. "sum" is the oriented integral of f over the area defined by
      //     (O, V_0, V_1, ..., V_i).
      assert(i == 1 || origin.angle(vertex(i)) < kMaxLength)
      assert(origin == vertex(0) || fabs(origin.dot(vertex(0))) < 1e-15)
      if vertex(i+1).angle(origin) > kMaxLength {
        // We are about to create an unstable edge, so choose a new origin O'
        // for the triangle fan.
        let old_origin = origin
        if origin == vertex(0) {
          // The following point is well-separated from V_i and V_0 (and
          // therefore V_i+1 as well).
          origin = vertex(0).pointCross(vertex(i))
        } else if vertex(i).angle(vertex(0)) < kMaxLength {
          // All edges of the triangle (O, V_0, V_i) are stable, so we can
          // revert to using V_0 as the origin.
          origin = vertex(0)
        } else {
          // (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are
          // perpendicular.  Therefore V_0.CrossProd(O) is approximately
          // perpendicular to all of {O, V_0, V_i, V_i+1}, and we can choose
          // this point O' as the new origin.
          origin = S2Point(raw: vertex(0).cross(old_origin))
          // Advance the edge (V_0,O) to (V_0,O').
          sum = sum + f_tri(vertex(0), old_origin, origin)
        }
        // Advance the edge (O,V_i) to (O',V_i).
        sum = sum + f_tri(old_origin, vertex(i), origin)
      }
      // Advance the edge (O,V_i) to (O,V_i+1).
      sum = sum + f_tri(origin, vertex(i), vertex(i+1))
    }
    // If the origin is not V_0, we need to sum one more triangle.
    if (origin != vertex(0)) {
      // Advance the edge (O,V_n-1) to (O,V_0).
      sum = sum + f_tri(origin, vertex(num_vertices - 1), vertex(0))
    }
    return sum
  }
  
}

protocol Summable {
  static func +(lhs: Self, rhs: Self) -> Self
}
