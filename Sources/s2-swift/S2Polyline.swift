//
//  S2Polyline.swift
//  s2-swift
//

import Foundation


// Polyline represents a sequence of zero or more vertices connected by
// straight edges (geodesics). Edges of length 0 and 180 degrees are not
// allowed, i.e. adjacent vertices should not be identical or antipodal.
public struct S2Polyline: S2ShapeType, S2RegionType {
 
  let points: [S2Point]
  
  // MARK: inits
  
  init(points: [S2Point]) {
    self.points = points
  }

  // PolylineFromLatLngs creates a new Polyline from the given LatLngs.
  init(latLngs: [LatLng]) {
    let points = latLngs.map { S2Point(latLng: $0) }
    self.init(points: points)
  }

  // MARK: 
  
  // Reverse reverses the order of the Polyline vertices.
  func reversed() -> S2Polyline {
    return S2Polyline(points: points.reversed())
  }
  
  // MARK: computed members
  
  // Length returns the length of this Polyline.
  var length: Double {
    var length = 0.0
    for i in 1..<points.count {
      length += points[i-1].distance(points[i])
    }
    return length
  }
  
  // dimension returns the dimension of the geometry represented by this Polyline.
  func dimension() -> ShapeDimension {
    return .polylineGeometry
  }
  
  // numChains reports the number of contiguous edge chains in this Polyline.
  func numChains() -> Int {
    if numEdges() >= 1 {
      return 1
    }
    return 0
  }
  
  // Centroid returns the true centroid of the polyline multiplied by the length of the
  // polyline. The result is not unit length, so you may wish to normalize it.
  //
  // Scaling by the Polyline length makes it easy to compute the centroid
  // of several Polylines (by simply adding up their centroids).
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
  
  // Equals reports whether the given Polyline is exactly the same as this one.
  func equals(_ b: S2Polyline) -> Bool {
    if points.count != b.points.count {
      return false
    }
    for i in 0..<points.count {
      if points[i] != b.points[i] {
        return false
      }
    }
    return true
  }
  
  // chainStart returns the id of the first edge in the i-th edge chain in this Polyline.
  func chainStart(i: Int) -> Int {
    if i == 0 {
      return 0
    }
    
    return numEdges()
  }

  // MARK: region protocol
  
  // CapBound returns the bounding Cap for this Polyline.
  public func capBound() -> S2Cap {
    return rectBound().capBound()
  }
  
  // RectBound returns the bounding Rect for this Polyline.
  public func rectBound() -> S2Rect {
    var rb = RectBounder()
    for v in points {
      rb.add(point: v)
    }
    return rb.rectBound()
  }
  
  // ContainsCell reports whether this Polyline contains the given Cell. Always returns false
  // because "containment" is not numerically well-defined except at the Polyline vertices.
  public func contains(_ cell: Cell) -> Bool {
    return false
  }
  
  // IntersectsCell reports whether this Polyline intersects the given Cell.
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
    let cellVertices = (0...3).map { cell.vertex($0) }
    for j in 0..<4 {
      var crosser = EdgeCrosser(a: cellVertices[j], b: cellVertices[(j+1)&3], c: points[0])
      for i in 1..<points.count {
        if crosser.chainCrossingSign(points[i]) != .doNotCross {
          // There is a proper crossing, or two vertices were the same.
          return true
        }
      }
    }
    return false
  }

  // MARK: Shape protocol
  
  // NumEdges returns the number of edges in this shape.
  public func numEdges() -> Int {
    if points.count == 0 {
      return 0
    }
    return points.count - 1
  }
  
  // Edge returns endpoints for the given edge index.
  public func edge(_ i: Int) -> (S2Point, S2Point) {
    return (points[i], points[i+1])
  }
  
  // HasInterior returns false as Polylines are not closed.
  public func hasInterior() -> Bool {
    return false
  }
  
  // ContainsOrigin returns false because there is no interior to contain s2.Origin.
  public func containsOrigin() -> Bool {
    return false
  }
}
