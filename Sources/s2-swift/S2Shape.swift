//
//  S2Shape.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/9/18.
//

import Foundation

/// Defines the types of geometry dimensions that a S2Shape supports.
public enum ShapeDimension: Int {
  case pointGeometry
  case polylineGeometry
  case polygonGeometry
}

/// Represents a geodesic edge consisting of two vertices. Zero-length edges are
/// allowed, and can be used to represent points.
public struct Edge {
  var v0: S2Point
  var v1: S2Point
}

extension Edge: Equatable, Comparable {
  
  public static func ==(lhs: Edge, rhs: Edge) -> Bool {
    return lhs.v0 == rhs.v0 && lhs.v1 == rhs.v1
  }
  
  public static func <(lhs: Edge, rhs: Edge) -> Bool {
    if lhs.v0 == rhs.v0 { return lhs.v1 < rhs.v1 }
    return lhs.v0 < rhs.v0
  }
  
}

/// Represents a range of edge IDs corresponding to a chain of connected
/// edges, specified as a (start, length) pair. The chain is defined to consist of
/// edge IDs {start, start + 1, ..., start + length - 1}.
public struct Chain {
  let start: Int
  let length: Int
}

/// Represents the position of an edge within a given edge chain,
/// specified as a (chainID, offset) pair. Chains are numbered sequentially
/// starting from zero, and offsets are measured from the start of each chain.
public struct ChainPosition {
  let chainId: Int
  let offset: Int
}

/// Consists of a point and a boolean indicating whether the point
/// is contained by a particular shape.
public struct ReferencePoint {
  let point: S2Point
  let contained: Bool
}

extension ReferencePoint {
  
  /// Returns a ReferencePoint with the given value for
  /// contained and the origin point. It should be used when all points or no
  /// points are contained.
  init(origin: Bool, contained: Bool) {
    self.init(point: S2Point.origin, contained: contained)
  }

}

/// Represents polygonal geometry in a flexible way. It is organized as a
/// collection of edges that optionally defines an interior. All geometry
/// represented by a given Shape must have the same dimension, which means that
/// an Shape can represent either a set of points, a set of polylines, or a set
/// of polygons.
///
/// Shape is defined as an interface in order to give clients control over the
/// underlying data representation. Sometimes an Shape does not have any data of
/// its own, but instead wraps some other type.
///
/// Shape operations are typically defined on a ShapeIndex rather than
/// individual shapes. An ShapeIndex is simply a collection of Shapes,
/// possibly of different dimensions (e.g. 10 points and 3 polygons), organized
/// into a data structure for efficient edge access.
///
/// The edges of a Shape are indexed by a contiguous range of edge IDs
/// starting at 0. The edges are further subdivided into chains, where each
/// chain consists of a sequence of edges connected end-to-end (a polyline).
/// For example, a Shape representing two polylines AB and CDE would have
/// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
/// Similarly, an Shape representing 5 points would have 5 chains consisting
/// of one edge each.
///
/// Shape has methods that allow edges to be accessed either using the global
/// numbering (edge ID) or within a particular chain. The global numbering is
/// sufficient for most purposes, but the chain representation is useful for
/// certain algorithms such as intersection (see BooleanOperation).
/// Defines an interface for any S2 type that needs to be indexable.
public protocol S2Shape {
  /// Returns the number of edges in this shape.
  func numEdges() -> Int
  /// Returns endpoints for the given edge index.
  func edge(_ i: Int) -> Edge
  /// Returns true if this shape has an interior.
  /// i.e. the Shape consists of one or more closed non-intersecting loops.
  func hasInterior() -> Bool
  /// Returns true if this shape contains s2.Origin.
  /// Shapes that do not have an interior will return false.
  func containsOrigin() -> Bool
  /// Returns an arbitrary reference point for the shape. (The
  /// containment boolean value must be false for shapes that do not have an interior.)
  /// This reference point may then be used to compute the containment of other
  /// points by counting edge crossings.
  func referencePoint() -> ReferencePoint
  /// Reports the number of contiguous edge chains in the shape.
  /// For example, a shape whose edges are [AB, BC, CD, AE, EF] would consist
  /// of two chains (AB,BC,CD and AE,EF). Every chain is assigned a chain Id
  /// numbered sequentially starting from zero.
  /// Note that it is always acceptable to implement this method by returning
  /// NumEdges, i.e. every chain consists of a single edge, but this may
  /// reduce the efficiency of some algorithms.
  func numChains() -> Int
  /// Returns the range of edge IDs corresponding to the given edge chain.
  /// Edge chains must form contiguous, non-overlapping ranges that cover
  /// the entire range of edge IDs. This is spelled out more formally below:
  ///  0 <= i < NumChains()
  ///  Chain(i).length > 0, for all i
  ///  Chain(0).start == 0
  ///  Chain(i).start + Chain(i).length == Chain(i+1).start, for i < NumChains()-1
  ///  Chain(i).start + Chain(i).length == NumEdges(), for i == NumChains()-1
  func chain(_ chainId: Int) -> Chain
  /// Returns the edge at offset "offset" within edge chain "chainID".
  /// Equivalent to "shape.Edge(shape.Chain(chainID).start + offset)"
  /// but more efficient.
  func chainEdge(chainId: Int, offset: Int) -> Edge
  /// Finds the chain containing the given edge, and returns the
  /// position of that edge as a ChainPosition(chainID, offset) pair.
  ///  shape.Chain(pos.chainID).start + pos.offset == edgeID
  ///  shape.Chain(pos.chainID+1).start > edgeID
  /// where pos == shape.ChainPosition(edgeID).
  func chainPosition(_ edgeId: Int) -> ChainPosition
  /// Returns the dimension of the geometry represented by this shape.
  ///  pointGeometry: Each point is represented as a degenerate edge.
  ///  polylineGeometry:  Polyline edges may be degenerate. A shape may
  ///      represent any number of polylines. Polylines edges may intersect.
  ///  polygonGeometry:  Edges should be oriented such that the polygon
  ///      interior is always on the left. In theory the edges may be returned
  ///      in any order, but typically the edges are organized as a collection
  ///      of edge chains where each chain represents one polygon loop.
  ///      Polygons may have degeneracies (e.g., degenerate edges or sibling
  ///      pairs consisting of an edge and its corresponding reversed edge).
  /// Note that this method allows degenerate geometry of different dimensions
  /// to be distinguished, e.g. it allows a point to be distinguished from a
  /// polyline or polygon that has been simplified to a single point.
  func dimension() -> ShapeDimension
}

/// Defines different ways of reporting edge intersections.
enum CrossingType {
  // CrossingTypeInterior reports intersections that occur at a point
  // interior to both edges (i.e., not at a vertex).
  case interior
  // CrossingTypeAll reports all intersections, even those where two edges
  // intersect only because they share a common vertex.
  case all
  // CrossingTypeNonAdjacent reports all intersections except for pairs of
  // the form (AB, BC) where both edges are from the same ShapeIndex.
  case nonAdjacent
}

/// A wrapper over ShapeIndexIterator with extra methods
/// that are useful for merging the contents of two or more ShapeIndexes.
struct RangeIterator {
  var it: ShapeIndexIterator
  // The min and max leaf cell ids covered by the current cell. If done() is
  // true, these methods return a value larger than any valid cell id.
  var rangeMin: CellId
  var rangeMax: CellId
}

extension RangeIterator {

  /// Creates a new rangeIterator positioned at the first cell of the given index.
  init(index: ShapeIndex) {
    it = index.iterator()
    rangeMin = it.cellId().rangeMin()
    rangeMax = it.cellId().rangeMax()
  }

  func cellId() -> CellId {
    return it.cellId()
  }
  
  func indexCell() -> ShapeIndexCell? {
    return it.indexCell()
  }
  
  mutating func next() {
    it.next()
    refresh()
  }
  
  func done() -> Bool {
    return it.done()
  }

  /// Positions the iterator at the first cell that overlaps or follows
  /// the current range minimum of the target iterator, i.e. such that its
  /// rangeMax >= target.rangeMin.
  mutating func seek(to target: RangeIterator) {
    it.seek(target: target.rangeMin)
    // If the current cell does not overlap target, it is possible that the
    // previous cell is the one we are looking for. This can only happen when
    // the previous cell contains target but has a smaller CellID.
    if it.done() || it.cellId().rangeMin() > target.rangeMax {
      if it.prev() && it.cellId().rangeMax() < target.cellId() {
        it.next()
      }
    }
    refresh()
  }
  
  /// Positions the iterator at the first cell that follows the current
  /// range minimum of the target iterator. i.e. the first cell such that its
  /// rangeMin > target.rangeMax.
  mutating func seek(beyond target: RangeIterator) {
    it.seek(target: target.rangeMax.next())
    if !it.done() && it.cellId().rangeMin() <= target.rangeMax {
      it.next()
    }
    refresh()
  }
  
  /// Updates the iterators min and max values.
  mutating func refresh() {
    rangeMin = cellId().rangeMin()
    rangeMax = cellId().rangeMax()
  }

}

/// referencePointForShape is a helper function for implementing various Shapes
/// ReferencePoint functions.
///
/// Given a shape consisting of closed polygonal loops, the interior of the
/// shape is defined as the region to the left of all edges (which must be
/// oriented consistently). This function then chooses an arbitrary point and
/// returns true if that point is contained by the shape.
///
/// Unlike Loop and Polygon, this method allows duplicate vertices and
/// edges, which requires some extra care with definitions. The rule that we
/// apply is that an edge and its reverse edge cancel each other: the result
/// is the same as if that edge pair were not present. Therefore shapes that
/// consist only of degenerate loop(s) are either empty or full; by convention,
/// the shape is considered full if and only if it contains an empty loop (see
/// laxPolygon for details).
///
/// Determining whether a loop on the sphere contains a point is harder than
/// the corresponding problem in 2D plane geometry. It cannot be implemented
/// just by counting edge crossings because there is no such thing as a point
/// at infinity that is guaranteed to be outside the loop.
///
/// This function requires that the given Shape have an interior.
func referencePoint(shape: S2Shape) -> ReferencePoint {
  if shape.numEdges() == 0 {
    // A shape with no edges is defined to be full if and only if it
    // contains an empty loop.
    return ReferencePoint(origin: true, contained: shape.numChains() > 0)
  }
  // Define a "matched" edge as one that can be paired with a corresponding
  // reversed edge. Define a vertex as "balanced" if all of its edges are
  // matched. In order to determine containment, we must find an unbalanced
  // vertex. Often every vertex is unbalanced, so we start by trying an
  // arbitrary vertex.
  let edge = shape.edge(0)
  if let ref = referencePointAtVertex(shape: shape, vTest: edge.v0) {
    return ref
  }
  // That didn't work, so now we do some extra work to find an unbalanced
  // vertex (if any). Essentially we gather a list of edges and a list of
  // reversed edges, and then sort them. The first edge that appears in one
  // list but not the other is guaranteed to be unmatched.
  let n = shape.numEdges()
  var edges = (0..<n).map { shape.edge($0) }
  var revEdges = edges.map { Edge(v0: $0.v1, v1: $0.v0) }
  edges.sort()
  revEdges.sort()
  for i in 0..<n {
    if edges[i] < revEdges[i] { // edges[i] is unmatched
      if let ref = referencePointAtVertex(shape: shape, vTest: edges[i].v0) {
        return ref
      }
    }
    if revEdges[i] < edges[i] { // revEdges[i] is unmatched
      if let ref = referencePointAtVertex(shape: shape, vTest: revEdges[i].v0) {
        return ref
      }
    }
  }
  // All vertices are balanced, so this polygon is either empty or full. By
  // convention it is defined to be full if it contains any empty loop.
  for i in 0..<shape.numChains() {
    if shape.chain(i).length == 0 {
      return ReferencePoint(origin: true, contained: true)
    }
  }
  return ReferencePoint(origin: true, contained: false)
}

/// referencePointAtVertex reports whether the given vertex is unbalanced, and
/// returns a ReferencePoint indicating if the point is contained.
/// Otherwise returns false.
func referencePointAtVertex(shape: S2Shape, vTest: S2Point) -> ReferencePoint? {
  // Let P be an unbalanced vertex. Vertex P is defined to be inside the
  // region if the region contains a particular direction vector starting from
  // P, namely the direction p.Ortho(). This can be calculated using
  // ContainsVertexQuery.
  var containsQuery = ContainsVertexQuery(target: vTest)
  let n = shape.numEdges()
  for e in 0..<n {
    let edge = shape.edge(e)
    if edge.v0 == vTest {
      containsQuery.addEdge(v: edge.v1, direction: 1)
    }
    if edge.v1 == vTest {
      containsQuery.addEdge(v: edge.v0, direction: -1)
    }
  }
  let containsSign = containsQuery.containsVertex()
  if containsSign == 0 {
    // There are no unmatched edges incident to this vertex.
    return nil
  }
  return ReferencePoint(point: vTest, contained: containsSign > 0)
}

/// containsBruteForce reports whether the given shape contains the given point.
/// Most clients should not use this method, since its running time is linear in
/// the number of shape edges. Instead clients should create a ShapeIndex and use
/// ContainsPointQuery, since this strategy is much more efficient when many
/// points need to be tested.
///
/// Polygon boundaries are treated as being semi-open (see ContainsPointQuery
/// and VertexModel for other options).
func containsBruteForce(shape: S2Shape, point: S2Point) -> Bool {
  if !shape.hasInterior() {
    return false
  }
  let refPoint = shape.referencePoint()
  if refPoint.point == point {
    return refPoint.contained
  }
  var crosser = EdgeCrosser(a: refPoint.point, b: point)
  var inside = refPoint.contained
  for e in 0..<shape.numEdges() {
    let edge = shape.edge(e)
    inside = inside != crosser.isEdgeOrVertexCrossing(c: edge.v0, d: edge.v1)
  }
  return inside
}

/// ContainsVertexQuery is used to track the edges entering and leaving the
/// given vertex of a Polygon in order to be able to determine if the point is
/// contained by the Polygon.
//
/// Point containment is defined according to the semi-open boundary model
/// which means that if several polygons tile the region around a vertex,
/// then exactly one of those polygons contains that vertex.
struct ContainsVertexQuery {
  let target: S2Point
  var edgeMap: [S2Point: Int]
}

extension ContainsVertexQuery {

  /// Returns a new query for the given vertex whose
  /// containment will be determined.
  init(target: S2Point) {
    self.target = target
    edgeMap = [:]
  }

  /// Adds the edge between target and v with the given direction.
  /// (+1 = outgoing, -1 = incoming, 0 = degenerate).
  mutating func addEdge(v: S2Point, direction: Int) {
    let value = edgeMap[v] ?? 0
    edgeMap[v] = value + direction
  }

  /// ContainsVertex reports a +1 if the target vertex is contained, -1 if it is
  /// not contained, and 0 if the incident edges consisted of matched sibling pairs.
  func containsVertex() -> Int {
    // Find the unmatched edge that is immediately clockwise from Ortho(P).
    let referenceDir = S2Point(raw: target.ortho())
    var bestPoint = referenceDir
    var bestDir = 0
    for (k, v) in edgeMap {
      if v == 0 {
        continue // This is a "matched" edge.
      }
      if S2Point.orderedCCW(referenceDir, bestPoint, k, target) {
        bestPoint = k
        bestDir = v
      }
    }
    return bestDir
  }

}

