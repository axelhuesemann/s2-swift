//
//  S2Polyline.swift
//  s2-swift
//

import Foundation


// Polyline represents a sequence of zero or more vertices connected by
// straight edges (geodesics). Edges of length 0 and 180 degrees are not
// allowed, i.e. adjacent vertices should not be identical or antipodal.
public struct S2Polyline: Shape, S2Region {
 
  private let points: [S2Point]
  
  // MARK: inits
  
  public init(points: [S2Point]) {
    self.points = points
  }

  // PolylineFromLatLngs creates a new Polyline from the given LatLngs.
  public init(latLngs: [LatLng]) {
    let points = latLngs.map { S2Point(latLng: $0) }
    self.init(points: points)
  }

  // MARK: 
  
  // Reverse reverses the order of the Polyline vertices.
  public func reversed() -> S2Polyline {
    return S2Polyline(points: points.reversed())
  }
  
  public func vertex(_ i: Int) -> S2Point {
    return points[i]
  }
  
  // MARK: computed members
  
  /// Returns the length of this Polyline.
  var length: Double {
    var length = 0.0
    for i in 1..<points.count {
      length += points[i-1].distance(points[i])
    }
    return length
  }
  
  /// Returns the dimension of the geometry represented by this Polyline.
  public func dimension() -> ShapeDimension {
    return .polylineGeometry
  }
  
  /// Reports the number of contiguous edge chains in this Polyline.
  public func numChains() -> Int {
    if numEdges() >= 1 {
      return 1
    }
    return 0
  }
  
  /// Returns the i-th edge Chain in the Shape.
  public func chain(_ chainId: Int) -> Chain {
    return Chain(start: 0, length: numEdges())
  }
  
  /// Returns the j-th edge of the i-th edge Chain.
  public func chainEdge(chainId: Int, offset: Int) -> Edge {
    return Edge(v0: points[offset], v1: points[offset+1])
  }
  
  /// Returns a pair (i, j) such that edgeID is the j-th edge
  public func chainPosition(_ edgeId: Int) -> ChainPosition {
    return ChainPosition(chainId: 0, offset: edgeId)
  }
  
  /// Returns the true centroid of the polyline multiplied by the length of the
  /// polyline. The result is not unit length, so you may wish to normalize it.
  /// Scaling by the Polyline length makes it easy to compute the centroid
  /// of several Polylines (by simply adding up their centroids).
  func centroid() -> R3Vector {
    var centroid = R3Vector.init(x: 0.0, y: 0.0, z: 0.0)
    for i in 1..<points.count {
      // The centroid (multiplied by length) is a vector toward the midpoint
      // of the edge, whose length is twice the sin of half the angle between
      // the two vertices. Defining theta to be this angle, we have:
      let vSum = points[i-1].v.add(points[i].v)  // Length == 2*cos(theta)
      let vDiff = points[i-1].v.sub(points[i].v) // Length == 2*sin(theta)
      // Length == 2*sin(theta)
      centroid = centroid.add(vSum.mul(sqrt(vDiff.norm2 / vSum.norm2)))
    }
    return centroid
  }
  
  /// Returns the id of the first edge in the i-th edge chain in this Polyline.
  func chainStart(i: Int) -> Int {
    if i == 0 {
      return 0
    }
    return numEdges()
  }

  // MARK: region protocol
  
  /// Returns the bounding Cap for this Polyline.
  public func capBound() -> S2Cap {
    return rectBound().capBound()
  }
  
  /// Returns the bounding Rect for this Polyline.
  public func rectBound() -> S2Rect {
    var rb = RectBounder()
    for v in points {
      rb.add(point: v)
    }
    return rb.rectBound()
  }
  
  /// Computes a covering of the Polyline.
  func cellUnionBound() -> CellUnion {
    return capBound().cellUnionBound()
  }
  
  /// Returns the default reference point with negative containment because Polylines are not closed.
  public func referencePoint() -> ReferencePoint {
    return ReferencePoint(origin: true, contained: false)
  }
  
  /// Reports whether this Polyline contains the given Cell. Always returns false
  /// because "containment" is not numerically well-defined except at the Polyline vertices.
  public func contains(_ cell: Cell) -> Bool {
    return false
  }
  
  /// Reports whether this Polyline intersects the given Cell.
  public func intersects(_ cell: Cell) -> Bool {
    if points.count == 0 {
      return false
    }
    // We only need to check whether the cell contains vertex 0 for correctness,
    // but these tests are cheap compared to edge crossings so we might as well
    // check all the vertices.
    for v in points {
      if cell.contains(v) {
        return true
      }
    }
    let cellVertices = (0..<4).map { cell.vertex($0) }
    for j in 0..<4 {
      var crosser = EdgeCrosser(a: cellVertices[j], b: cellVertices[(j+1)&3], c: points[0])
      for i in 1..<points.count {
        if crosser.chainCrossingSign(d: points[i]) != .doNotCross {
          // There is a proper crossing, or two vertices were the same.
          return true
        }
      }
    }
    return false
  }
  
  /// Returns false since Polylines are not closed.
  func contains(_ point: S2Point) -> Bool {
    return false
  }
  
  var numVertices: Int {
    return points.count
  }
  
  // MARK: Shape protocol
  
  /// Returns the number of edges in this shape.
  public func numEdges() -> Int {
    if points.count == 0 {
      return 0
    }
    return points.count - 1
  }
  
  /// Returns endpoints for the given edge index.
  public func edge(_ i: Int) -> Edge {
    return Edge(v0:points[i], v1: points[i+1])
  }
  
  /// Returns false as Polylines are not closed.
  public func hasInterior() -> Bool {
    return false
  }
  
  /// Ceturns false because there is no interior to contain s2.Origin.
  public func containsOrigin() -> Bool {
    return false
  }
  
  // MARK: ?
  
  /// Reports the maximal end index such that the line segment between
  /// the start index and this one such that the line segment between these two
  /// vertices passes within the given tolerance of all interior vertices, in order.
  func findEndVertex(tolerance:  S1Angle, index: Int) -> Int {
    // The basic idea is to keep track of the "pie wedge" of angles
    // from the starting vertex such that a ray from the starting
    // vertex at that angle will pass through the discs of radius
    // tolerance centered around all vertices processed so far.
    //
    // First we define a coordinate frame for the tangent and normal
    // spaces at the starting vertex. Essentially this means picking
    // three orthonormal vectors X,Y,Z such that X and Y span the
    // tangent plane at the starting vertex, and Z is up. We use
    // the coordinate frame to define a mapping from 3D direction
    // vectors to a one-dimensional ray angle in the range (-π,
    // π]. The angle of a direction vector is computed by
    // transforming it into the X,Y,Z basis, and then calculating
    // atan2(y,x). This mapping allows us to represent a wedge of
    // angles as a 1D interval. Since the interval wraps around, we
    // represent it as an Interval, i.e. an interval on the unit
    // circle.
    let origin = points[index]
    let frame = S2Point.getFrame(origin)
    // As we go along, we keep track of the current wedge of angles
    // and the distance to the last vertex (which must be
    // non-decreasing).
    var currentWedge = S1Interval.full
    var lastDistance: S1Angle = 0
    var index = index + 1
    while index < points.count {
      let candidate = points[index]
      let distance = origin.distance(candidate)
      // We don't allow simplification to create edges longer than
      // 90 degrees, to avoid numeric instability as lengths
      // approach 180 degrees. We do need to allow for original
      // edges longer than 90 degrees, though.
      if distance > .pi / 2 && lastDistance > 0 {
        break
      }
      // Vertices must be in increasing order along the ray, except
      // for the initial disc around the origin.
      if distance < lastDistance && lastDistance > tolerance {
        break
      }
      lastDistance = distance
      // Points that are within the tolerance distance of the origin
      // do not constrain the ray direction, so we can ignore them.
      if distance <= tolerance {
        index += 1
        continue
      }
      // If the current wedge of angles does not contain the angle
      // to this vertex, then stop right now. Note that the wedge
      // of possible ray angles is not necessarily empty yet, but we
      // can't continue unless we are willing to backtrack to the
      // last vertex that was contained within the wedge (since we
      // don't create new vertices). This would be more complicated
      // and also make the worst-case running time more than linear.
      let direction = S2Point.toFrame(frame, point: candidate)
      let center = atan2(direction.y, direction.x)
      if !currentWedge.contains(center) {
        break
      }
      // To determine how this vertex constrains the possible ray
      // angles, consider the triangle ABC where A is the origin, B
      // is the candidate vertex, and C is one of the two tangent
      // points between A and the spherical cap of radius
      // tolerance centered at B. Then from the spherical law of
      // sines, sin(a)/sin(A) = sin(c)/sin(C), where a and c are
      // the lengths of the edges opposite A and C. In our case C
      // is a 90 degree angle, therefore A = asin(sin(a) / sin(c)).
      // Angle A is the half-angle of the allowable wedge.
      let halfAngle = asin(sin(tolerance) / sin(distance))
      let target = S1Interval(lo: center, hi: center).expanded(halfAngle)
      currentWedge = currentWedge.intersection(target)
      index += 1
    }
    // We break out of the loop when we reach a vertex index that
    // can't be included in the line segment, so back up by one
    // vertex.
    return index - 1
  }

  /// SubsampleVertices returns a subsequence of vertex indices such that the
  /// polyline connecting these vertices is never further than the given tolerance from
  /// the original polyline. Provided the first and last vertices are distinct,
  /// they are always preserved; if they are not, the subsequence may contain
  /// only a single index.
  /// Some useful properties of the algorithm:
  ///  - It runs in linear time.
  ///  - The output always represents a valid polyline. In particular, adjacent
  ///    output vertices are never identical or antipodal.
  ///  - The method is not optimal, but it tends to produce 2-3% fewer
  ///    vertices than the Douglas-Peucker algorithm with the same tolerance.
  ///  - The output is parametrically equivalent to the original polyline to
  ///    within the given tolerance. For example, if a polyline backtracks on
  ///    itself and then proceeds onwards, the backtracking will be preserved
  ///    (to within the given tolerance). This is different than the
  ///    Douglas-Peucker algorithm which only guarantees geometric equivalence.
  func subsampleVertices(tolerance: S1Angle) -> [Int] {
    if points.count < 1 {
      return []
    }
    var result = [0]
    let clampedTolerance = S1Angle(max(tolerance, 0))
    var index = 0
    while index < points.count - 1 {
      let nextIndex = findEndVertex(tolerance: clampedTolerance, index: index)
      // Don't create duplicate adjacent vertices.
      if points[nextIndex] != points[index] {
        result.append(nextIndex)
      }
      index = nextIndex
    }
    return result
  }
  
}

extension S2Polyline: Equatable {
  
  /// Reports whether the given Polyline is exactly the same as this one.
  public static func ==(lhs: S2Polyline, rhs: S2Polyline) -> Bool {
    return lhs.points == rhs.points
  }

}
