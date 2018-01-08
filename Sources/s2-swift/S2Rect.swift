//
//  S2Rect.swift
//  s2-swift
//

import Foundation


/// Represents a closed latitude-longitude rectangle.
public struct S2Rect: S2RegionType {

  //
  static let validRectLatRange = R1Interval(lo: -.pi / 2, hi: .pi / 2)
  static let validRectLngRange = S1Interval.full
  
  //
  let lat: R1Interval
  let lng: S1Interval

  // MARK: inits / factory
  
  init(lat: R1Interval, lng: S1Interval) {
    self.lat = lat
    self.lng = lng
  }
  
  public init(lo: LatLng, hi: LatLng) {
    lat = R1Interval(lo: lo.lat, hi: hi.lat)
    lng = S1Interval(lo: lo.lng, hi: hi.lng)
  }
  
  public init(p1: LatLng, p2: LatLng) {
    lat = R1Interval(lo: p1.lat, hi: p2.lat)
    lng = S1Interval(lo: p1.lng, hi: p2.lng)
  }
  
  /// Constructs a rectangle containing a single point p.
  init(latLng: LatLng) {
    lat = R1Interval(lo: latLng.lat, hi: latLng.lat)
    lng = S1Interval(lo: latLng.lng, hi: latLng.lng)
  }
  
  /// Constructs a rectangle with the given size and center.
  /// center needs to be normalized, but size does not. The latitude
  /// interval of the result is clamped to [-90,90] degrees, and the longitude
  /// interval of the result is FullRect() if and only if the longitude size is
  /// 360 degrees or more.
  ///
  /// Examples of clamping (in degrees):
  ///   center=(80,170),  size=(40,60)   -> lat=[60,90],   lng=[140,-160]
  ///   center=(10,40),   size=(210,400) -> lat=[-90,90],  lng=[-180,180]
  ///   center=(-90,180), size=(20,50)   -> lat=[-90,-80], lng=[155,-155]
  init(center: LatLng, size: LatLng) {
    let half = LatLng(lat: size.lat / 2, lng: size.lng / 2)
    self.init(latLng: center)
    self = self.expanded(half)
  }
  
  /// Returns the empty rectangle.
  public static let empty = S2Rect(lat: R1Interval.empty, lng: S1Interval.empty)
  
  /// Returns the full rectangle.
  public static let full = S2Rect(lat: validRectLatRange, lng: validRectLngRange)
  
  // MARK: tests
  
  /// Returns true iff the rectangle is valid.
  /// This requires Lat ⊆ [-π/2,π/2] and Lng ⊆ [-π,π], and Lat = ∅ iff Lng = ∅
  var isValid: Bool {
    return abs(lat.lo) <= .pi/2 && abs(lat.hi) <= .pi/2 && lng.isValid && lat.isEmpty == lng.isEmpty
  }

  /// Reports whether the rectangle is empty.
  public var isEmpty: Bool {
    return lat.isEmpty
  }

  /// Reports whether the rectangle is full.
  public var isFull: Bool {
    return lat == S2Rect.validRectLatRange && lng.isFull
  }

  /// Reports whether the rectangle is a single point.
  var isPoint: Bool {
    return lat.lo == lat.hi && lng.lo == lng.hi
  }

  // MARK: computed members
  
  /// Returns the i-th vertex of the rectangle (i = 0,1,2,3) in CCW order
  /// (lower left, lower right, upper right, upper left).
  func vertex(_ i: Int) -> LatLng {
    switch i {
    case 0:
      return LatLng(lat: lat.lo, lng: lng.lo)
    case 1:
      return LatLng(lat: lat.lo, lng: lng.hi)
    case 2:
      return LatLng(lat: lat.hi, lng: lng.hi)
    default:
      return LatLng(lat: lat.hi, lng: lng.lo)
    }
  }
  
  /// Returns the vertices of the rectangle (i = 0,1,2,3) in CCW order
  /// (lower left, lower right, upper right, upper left).
  func vertices() -> [LatLng] {
    return [
      LatLng(lat: lat.lo, lng: lng.lo),
      LatLng(lat: lat.lo, lng: lng.hi),
      LatLng(lat: lat.hi, lng: lng.hi),
      LatLng(lat: lat.hi, lng: lng.lo)]
  }
  
  /// Returns one corner of the rectangle.
  func lo() -> LatLng {
    return LatLng(lat: lat.lo, lng: lng.lo)
  }

  /// Returns the other corner of the rectangle.
  func hi() -> LatLng {
    return LatLng(lat: lat.hi, lng: lng.hi)
  }

  /// Returns the center of the rectangle.
  public var center: LatLng {
    return LatLng(lat: lat.center, lng: lng.center)
  }

  /// Returns the size of the Rect.
  public func size() -> LatLng {
    return LatLng(lat: lat.length, lng: lng.length)
  }
  
  /// Returns the surface area of the Rect.
  public func area() -> Double {
    if isEmpty {
      return 0
    }
    let capDiff = abs(sin(lat.hi) - sin(lat.lo))
    return lng.length * capDiff
  }

  /// Returns a cap that countains Rect.
  public func capBound() -> S2Cap {
    // We consider two possible bounding caps, one whose axis passes
    // through the center of the lat-long rectangle and one whose axis
    // is the north or south pole.  We return the smaller of the two caps.
    if isEmpty {
      return S2Cap.empty
    }
    var poleZ: Double
    var poleAngle: Double
    if lat.hi + lat.lo < 0 {
      // South pole axis yields smaller cap.
      poleZ = -1
      poleAngle = .pi/2 + lat.hi
    } else {
      poleZ = 1
      poleAngle = .pi/2 - lat.lo
    }
    let poleCap = S2Cap(center: S2Point(x: 0, y: 0, z: poleZ), angle: poleAngle)
    // For bounding rectangles that span 180 degrees or less in longitude, the
    // maximum cap size is achieved at one of the rectangle vertices.  For
    // rectangles that are larger than 180 degrees, we punt and always return a
    // bounding cap centered at one of the two poles.
    if (lng.hi-lng.lo).truncatingRemainder(dividingBy: .pi*2) >= 0 && lng.hi-lng.lo < .pi*2 {
      let midCap = S2Cap(point: center.toPoint()).add(lo().toPoint()).add(hi().toPoint())
      if midCap.height < poleCap.height {
        return midCap
      }
    }
    return poleCap
  }
  
  /// Returns itself.
  public func rectBound() -> S2Rect {
    return self
  }
  
  // MARK: arithmetic
  
  /// Increases the size of the rectangle to include the given point.
  public func add(_ latLng: LatLng) -> S2Rect {
    if !latLng.isValid {
      return self
    }
    return S2Rect(lat: lat.add(latLng.lat), lng: lng.add(latLng.lng))
  }

  /// Returns a rectangle that has been expanded by margin.Lat on each side
  /// in the latitude direction, and by margin.Lng on each side in the longitude
  /// direction. If either margin is negative, then it shrinks the rectangle on
  /// the corresponding sides instead. The resulting rectangle may be empty.
  ///
  /// The latitude-longitude space has the topology of a cylinder. Longitudes
  /// "wrap around" at +/-180 degrees, while latitudes are clamped to range [-90, 90].
  /// This means that any expansion (positive or negative) of the full longitude range
  /// remains full (since the "rectangle" is actually a continuous band around the
  /// cylinder), while expansion of the full latitude range remains full only if the
  /// margin is positive.
  ///
  /// If either the latitude or longitude interval becomes empty after
  /// expansion by a negative margin, the result is empty.
  ///
  /// Note that if an expanded rectangle contains a pole, it may not contain
  /// all possible lat/lng representations of that pole, e.g., both points [π/2,0]
  /// and [π/2,1] represent the same pole, but they might not be contained by the
  /// same Rect.
  ///
  /// If you are trying to grow a rectangle by a certain distance on the
  /// sphere (e.g. 5km), refer to the ExpandedByDistance() C++ method implementation
  /// instead.
  public func expanded(_ margin: LatLng) -> S2Rect {
    let lat = self.lat.expanded(margin.lat)
    let lng = self.lng.expanded(margin.lng)
    if lat.isEmpty || lng.isEmpty {
      return S2Rect.empty
    }
    return S2Rect(lat: lat.intersection(S2Rect.validRectLatRange), lng: lng)
  }

  /// Returns the rectangle unmodified if it does not include either pole.
  /// If it includes either pole, PolarClosure returns an expansion of the rectangle along
  /// the longitudinal range to include all possible representations of the contained poles.
  func polarClosure() -> S2Rect {
    if lat.lo == -.pi/2 || lat.hi == .pi/2 {
      return S2Rect(lat: lat, lng: S1Interval.full)
    }
    return self
  }

  /// Returns the smallest Rect containing the union of this rectangle and the given rectangle.
  public func union(_ rect: S2Rect) -> S2Rect {
    return S2Rect(lat: lat.union(rect.lat), lng: lng.union(rect.lng))
  }

  /// Returns the smallest rectangle containing the intersection of
  /// this rectangle and the given rectangle. Note that the region of intersection
  /// may consist of two disjoint rectangles, in which case a single rectangle
  /// spanning both of them is returned.
  public func intersection(_ rect: S2Rect) -> S2Rect {
    let lat = self.lat.intersection(rect.lat)
    let lng = self.lng.intersection(rect.lng)
    if lat.isEmpty || lng.isEmpty {
      return S2Rect.empty
    }
    return S2Rect(lat: lat, lng: lng)
  }

  // MARK: contains / intersects
  
  /// Reports whether this rectangle and the other have any points in common.
  public func intersects(_ rect: S2Rect) -> Bool {
    return lat.intersects(rect.lat) && lng.intersects(rect.lng)
  }

  /// Reports whether this Rect contains the other Rect.
  public func contains(_ rect: S2Rect) -> Bool {
    return lat.contains(rect.lat) && lng.contains(rect.lng)
  }

  /// Reports whether the given Cell is contained by this Rect.
  public func contains(_ cell: Cell) -> Bool {
    // A latitude-longitude rectangle contains a cell if and only if it contains
    // the cell's bounding rectangle. This test is exact from a mathematical
    // point of view, assuming that the bounds returned by Cell.RectBound()
    // are tight. However, note that there can be a loss of precision when
    // converting between representations -- for example, if an s2.Cell is
    // converted to a polygon, the polygon's bounding rectangle may not contain
    // the cell's bounding rectangle. This has some slightly unexpected side
    // effects; for instance, if one creates an s2.Polygon from an s2.Cell, the
    // polygon will contain the cell, but the polygon's bounding box will not.
    return contains(cell.rectBound())
  }

  /// Reports whether the given LatLng is within the Rect.
  public func contains(_ latLng: LatLng) -> Bool {
    if !latLng.isValid {
      return false
    }
    return lat.contains(latLng.lat) && lng.contains(latLng.lng)
  }

  /// Reports whether the given Point is within the Rect.
  func contains(_ point: S2Point) -> Bool {
    return contains(LatLng(point: point))
  }

  /// Reports if the edge AB intersects the given edge of constant
  /// latitude. Requires the points to have unit length.
  func intersectsLatEdge(a: S2Point, b: S2Point, lat: Double, lng: S1Interval) -> Bool {
    // Unfortunately, lines of constant latitude are curves on
    // the sphere. They can intersect a straight edge in 0, 1, or 2 points.
    // First, compute the normal to the plane AB that points vaguely north.
    var z = a.pointCross(b)
    if z.z < 0 {
      let zv = z.v.mul(-1)
      z = S2Point(raw: zv)
    }
    // Extend this to an orthonormal frame (x,y,z) where x is the direction
    // where the great circle through AB achieves its maximium latitude.
    let y = z.pointCross(S2Point(x: 0, y: 0, z: 1))
    let x = S2Point(raw: y.v.cross(z.v))
    // Compute the angle "theta" from the x-axis (in the x-y plane defined
    // above) where the great circle intersects the given line of latitude.
    let sinLat = sin(lat)
    if abs(sinLat) >= x.z {
      // The great circle does not reach the given latitude.
      return false
    }
    let cosTheta = sinLat / x.z
    let sinTheta = sqrt(1.0 - cosTheta * cosTheta)
    let theta = atan2(sinTheta, cosTheta)
    // The candidate intersection points are located +/- theta in the x-y
    // plane. For an intersection to be valid, we need to check that the
    // intersection point is contained in the interior of the edge AB and
    // also that it is contained within the given longitude interval "lng".
    // Compute the range of theta values spanned by the edge AB.
    let abTheta = S1Interval(p1: atan2(a.v.dot(y.v), a.v.dot(x.v)), p2: atan2(b.v.dot(y.v), b.v.dot(x.v)))
    if abTheta.contains(theta) {
      // Check if the intersection point is also in the given lng interval.
      let isect = x.v.mul(cosTheta).add(y.v.mul(sinTheta))
      if lng.contains(atan2(isect.y, isect.x)) {
        return true
      }
    }
    if abTheta.contains(-theta) {
      // Check if the other intersection point is also in the given lng interval.
      let isect = x.v.mul(cosTheta).sub(y.v.mul(sinTheta))
      if lng.contains(atan2(isect.y, isect.x)) {
        return true
      }
    }
    return false
  }

  /// Reports if the edge AB intersects the given edge of constant
  /// longitude. Requires the points to have unit length.
  func intersectsLngEdge(a: S2Point, b: S2Point, lat: R1Interval, lng: Double) -> Bool {
    // The nice thing about edges of constant longitude is that
    // they are straight lines on the sphere (geodesics).
    let c = LatLng(lat: lat.lo, lng: lng).toPoint()
    let d = LatLng(lat: lat.hi, lng: lng).toPoint()
    return simpleCrossing(a, b: b, c: c, d: d)
  }

  /// Reports whether this rectangle intersects the given cell. This is an
  /// exact test and may be fairly expensive.
  public func intersects(_ cell: Cell) -> Bool {
    if isEmpty {
      return false
    }
    // First we eliminate the cases where one region completely contains the
    // rect. Once these are disposed of, then the regions will intersect
    // if and only if their boundaries intersect.
    if contains(cell.id.point()) {
      return true
    }
    if cell.contains(center.toPoint()) {
      return true
    }
    // Quick rejection test. This is optimistic speed up, and not required for correctness.
    if !intersects(cell.rectBound()) {
      return false
    }
    // Precompute the cell vertices as points and latitude-longitudes. We also
    // check whether the Cell contains any corner of the rectangle, or
    // vice-versa, since the edge-crossing tests only check the edge interiors.
    let vertices = (0..<4).map { cell.vertex($0) }
    var latlngs = [LatLng]() // 4
    for i in 0..<vertices.count {
      latlngs[i] = LatLng(point: vertices[i])
      if contains(latlngs[i]) {
        return true
      }
      if cell.contains(vertex(i).toPoint()) {
        return true
      }
    }
    // Now check whether the boundaries intersect. Unfortunately, a
    // latitude-longitude rectangle does not have straight edges: two edges
    // are curved, and at least one of them is concave.
    for i in 0..<vertices.count {
      let edgeLng = S1Interval(p1: latlngs[i].lng, p2: latlngs[(i+1)&3].lng)
      if !lng.intersects(edgeLng) {
        continue
      }
      // a->b defines the edge of the cell
      let a = vertices[i]
      let b = vertices[(i+1)&3]
      if edgeLng.contains(lng.lo) && intersectsLngEdge(a: a, b: b, lat: lat, lng: lng.lo) {
        return true
      }
      if edgeLng.contains(lng.hi) && intersectsLngEdge(a: a, b: b, lat: lat, lng: lng.hi) {
        return true
      }
      if intersectsLatEdge(a: a, b: b, lat: lat.lo, lng: lng) {
        return true
      }
      if intersectsLatEdge(a: a, b: b, lat: lat.hi, lng: lng) {
        return true
      }
    }
    return false
  }

}

extension S2Rect: CustomStringConvertible {

  public var description: String {
    return "[Lo \(lo()), Hi \(hi())]"
  }

}
