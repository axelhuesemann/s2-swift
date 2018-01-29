//
//  S2Cap.swift
//  s2-swift
//

import Foundation

/// Represents a disc-shaped region defined by a center and radius.
/// Technically this shape is called a "spherical cap" (rather than disc)
/// because it is not planar; the cap represents a portion of the sphere that
/// has been cut off by a plane. The boundary of the cap is the circle defined
/// by the intersection of the sphere and the plane. For containment purposes,
/// the cap is a closed set, i.e. it contains its boundary.
///
/// For the most part, you can use a spherical cap wherever you would use a
/// disc in planar geometry. The radius of the cap is measured along the
/// surface of the sphere (rather than the straight-line distance through the
/// interior). Thus a cap of radius π/2 is a hemisphere, and a cap of radius
/// π covers the entire sphere.
///
/// The center is a point on the surface of the unit sphere. (Hence the need for
/// it to be of unit length.)
///
/// Internally, the cap is represented by its center and "height". The height
/// is the distance from the center point to the cutoff plane. This
/// representation is much more efficient for containment tests than the
/// (center, radius) representation. There is also support for "empty" and
/// "full" caps, which contain no points and all points respectively.
///
/// The zero value of Cap is an invalid cap. Use EmptyCap to get a valid empty cap.
public struct S2Cap: S2Region {
  
  // MARK: constants
  
  // with radius 1.0, 0.0 means 0.0 size, and 2.0 means the entire globe
  static let zeroHeight = 0.0
  static let fullHeight = 2.0
  // negative is invalid and used as a marker for empty, which is different from zero
  static let emptyHeight = -1.0
  
  // TODO check if this is a good replacement value
  // static let roundUp = 1.0 + 1.0 / Double(1 << 52)
  static let roundUp = nextafter(1.0, 2.0)
  
  // centerPoint is the default center for S2Caps
  static let centerPoint = S2Point(x: 1.0, y: 0, z: 0)

  // MARK: state 
  
  /// Center of the cap represented as unit vector
  public let center: S2Point
  /// Distance from the plane that creates the cap
  let height: Double

  // MARK: inits / factory 
  
  /// Constructs a cap with the given center and height. A
  /// negative height yields an empty cap; a height of 2 or more yields a full cap.
  init(centerNormalized: S2Point, height: Double) {
    self.center = centerNormalized
    self.height = height
  }
  
  /// Constructs a cap with the given center and height. A
  /// negative height yields an empty cap; a height of 2 or more yields a full cap.
  init(center: S2Point, height: Double) {
    self.center = center
    self.height = height
  }
  
  /// Constructs a cap containing a single point.
  init(point: S2Point) {
    self.init(center: point, height: S2Cap.zeroHeight)
  }

  /// Constructs a cap with the given center and angle.
  init(centerNormalized: S2Point, angle: S1Angle) {
    let height = S2Cap.radiusToHeight(angle)
    self.init(centerNormalized: centerNormalized, height: height)
  }
  
  /// Constructs a cap with the given center and angle.
  init(center: S2Point, angle: S1Angle) {
    let height = S2Cap.radiusToHeight(angle)
    self.init(center: center, height: height)
  }
  
  /// Constructs a cap with the given center and surface area.
  /// Note that the area can also be interpreted as the solid angle subtended by the
  /// cap (because the sphere has unit radius). A negative area yields an empty cap;
  /// an area of 4*π or more yields a full cap.
  init(center: S2Point, area: Double) {
    let height = area / (.pi * 2.0)
    self.init(center: center, height: height)
  }

  // EmptyCap returns a cap that contains no points.
  public static let empty = S2Cap(center: centerPoint, height: emptyHeight)

  // FullCap returns a cap that contains all points.
  public static let full = S2Cap(center: centerPoint, height: fullHeight)

  // MARK: tests
  
  // IsValid reports whether the Cap is considered valid.
  // Heights are normalized so that they do not exceed 2.
  var isValid: Bool {
    return center.isUnit && height <= S2Cap.fullHeight
  }

  // IsEmpty reports whether the cap is empty, i.e. it contains no points.
  public var isEmpty: Bool {
    return height < S2Cap.zeroHeight
  }

  // IsFull reports whether the cap is full, i.e. it contains all points.
  public var isFull: Bool {
    return height == S2Cap.fullHeight
  }
  
  // MARK: contain / intersect 
  
  // Contains reports whether this cap contains the other.
  public func contains(_ cap: S2Cap) -> Bool {
    // In a set containment sense, every cap contains the empty cap.
    if isFull || cap.isEmpty {
      return true
    }
    return radius() >= center.distance(cap.center) + cap.radius()
  }
  
  // Intersects reports whether this cap intersects the other cap.
  // i.e. whether they have any points in common.
  public func intersects(_ cap: S2Cap) -> Bool {
    if isEmpty || cap.isEmpty {
      return false
    }
    return radius() + cap.radius() >= center.distance(cap.center)
  }
  
  // InteriorIntersects reports whether this caps interior intersects the other cap.
  public func interiorIntersects(_ cap: S2Cap) -> Bool {
    // Make sure this cap has an interior and the other cap is non-empty.
    if height <= S2Cap.zeroHeight || cap.isEmpty {
      return false
    }
    return radius() + cap.radius() > center.distance(cap.center)
  }
  
  /// Reports whether this cap contains the point.
  public func contains(_ point: S2Point) -> Bool {
    return center.v.sub(point.v).norm2 <= 2 * height
  }
  
  /// Reports whether the point is within the interior of this cap.
  public func interiorContains(_ point: S2Point) -> Bool {
    return isFull || center.v.sub(point.v).norm2 < 2 * height
  }

  // MARK: contains / intersects Cell
  
  /// Reports whether the cap contains the given cell.
  public func contains(_ cell: Cell) -> Bool {
    // If the cap does not contain all cell vertices, return false.
    var vertices = [S2Point]()
    for k in 0..<4 {
      let vertex = cell.vertex(k)
      vertices.append(vertex)
      if !contains(vertex) {
        return false
      }
    }
    // Otherwise, return true if the complement of the cap does not intersect the cell.
    return !complement().intersects(cell, vertices: vertices)
  }
  
  /// Reports whether the cap intersects the cell.
  public func intersects(_ cell: Cell) -> Bool {
    // If the cap contains any cell vertex, return true.
    var vertices = [S2Point]()
    for k in 0..<4 {
      let vertex = cell.vertex(k)
      vertices.append(vertex)
      if contains(vertex) {
        return true
      }
    }
    return intersects(cell, vertices: vertices)
  }
  
  /// Reports whether the cap intersects any point of the cell excluding
  /// its vertices (which are assumed to already have been checked).
  func intersects(_ cell: Cell, vertices: [S2Point]) -> Bool {
    // If the cap is a hemisphere or larger, the cell and the complement of the cap
    // are both convex. Therefore since no vertex of the cell is contained, no other
    // interior point of the cell is contained either.
    if height >= 1 {
      return false
    }
    // We need to check for empty caps due to the center check just below.
    if isEmpty {
      return false
    }
    // Optimization: return true if the cell contains the cap center. This allows half
    // of the edge checks below to be skipped.
    if cell.contains(center) {
      return true
    }
    // At this point we know that the cell does not contain the cap center, and the cap
    // does not contain any cell vertex. The only way that they can intersect is if the
    // cap intersects the interior of some edge.
    let sin2Angle = height * (2 - height)
    for k in 0..<4 {
      let edge = cell.edge(k)
      let dot = center.v.dot(edge.v)
      if dot > 0.0 {
        // The center is in the interior half-space defined by the edge. We do not need
        // to consider these edges, since if the cap intersects this edge then it also
        // intersects the edge on the opposite side of the cell, because the center is
        // not contained with the cell.
        continue
      }
      // The norm2 factor is necessary because "edge" is not normalized.
      if dot * dot > sin2Angle * edge.v.norm2 {
        return false
      }
      // Otherwise, the great circle containing this edge intersects the interior of the cap. We just
      // need to check whether the point of closest approach occurs between the two edge endpoints.
      let dir = edge.v.cross(center.v)
      if dir.dot(vertices[k].v) < 0 && dir.dot(vertices[(k+1) & 3].v) > 0 {
        return true
      }
    }
    return false
  }
  
  // MARK: computed members
  
  /// Returns the cap's radius.
  func radius() -> S1Angle {
    if isEmpty {
      return S2Cap.emptyHeight
    }
    // This could also be computed as acos(1 - height_), but the following
    // formula is much more accurate when the cap height is small. It
    // follows from the relationship h = 1 - cos(r) = 2 sin^2(r/2).
    return 2.0 * asin(sqrt(0.5 * height))
  }

  /// Returns the surface area of the Cap on the unit sphere.
  func area() -> Double {
    return 2.0 * .pi * max(S2Cap.zeroHeight, height)
  }

  /// Returns the complement of the interior of the cap. A cap and its
  /// complement have the same boundary but do not share any interior points.
  /// The complement operator is not a bijection because the complement of a
  /// singleton cap (containing a single point) is the same as the complement
  /// of an empty cap.
  func complement() -> S2Cap {
    let height: Double
    if isFull {
      height = S2Cap.emptyHeight
    } else {
      height = S2Cap.fullHeight - max(self.height, S2Cap.zeroHeight)
    }
    let antiCenter = center.inverse()
    return S2Cap(center: antiCenter, height: height)
  }

  /// Returns a bounding spherical cap. This is not guaranteed to be exact.
  public func capBound() -> S2Cap {
    return self
  }

  /// Returns a bounding latitude-longitude rectangle.
  /// The bounds are not guaranteed to be tight.
  public func rectBound() -> S2Rect {
    if isEmpty {
      return S2Rect.empty
    }
    //
    let capAngle = radius()
    let midAngle = center.lat
    var allLongitudes = false
    // Check whether cap includes the south pole.
    var latLo = midAngle - capAngle
    if latLo <= -.pi/2 {
      latLo = -.pi / 2
      allLongitudes = true
    }
    // Check whether cap includes the north pole.
    var latHi = midAngle + capAngle
    if latHi >= .pi/2 {
      latHi = .pi / 2
      allLongitudes = true
    }
    let lat = R1Interval(lo: latLo, hi: latHi)
    //
    var lng = S1Interval.full
    if !allLongitudes {
      // Compute the range of longitudes covered by the cap. We use the law
      // of sines for spherical triangles. Consider the triangle ABC where
      // A is the north pole, B is the center of the cap, and C is the point
      // of tangency between the cap boundary and a line of longitude. Then
      // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
      // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
      // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
      // minus the latitude). This formula also works for negative latitudes.
      //
      // The formula for sin(a) follows from the relationship h = 1 - cos(a).
      let sinA = sqrt(height * (2 - height))
      let sinC = cos(center.lat)
      if sinA <= sinC {
        let angleA = asin(sinA / sinC)
        let lngLo = (center.lng - angleA).truncatingRemainder(dividingBy: .pi*2)
        let lngHi = (center.lng + angleA).truncatingRemainder(dividingBy: .pi*2)
        lng = S1Interval(lo: lngLo, hi: lngHi)
      }
    }
    return S2Rect(lat: lat, lng: lng)
  }
  
  /// Computes a covering of the Cap. In general the covering
  /// consists of at most 4 cells except for very large caps, which may need
  /// up to 6 cells. The output is not sorted.
  func cellUnionBound() -> CellUnion {
    // TODO(roberts): The covering could be made quite a bit tighter by mapping
    // the cap to a rectangle in (i,j)-space and finding a covering for that.
    // Find the maximum level such that the cap contains at most one cell vertex
    // and such that CellID.AppendVertexNeighbors() can be called.
    let level = S2CellMetric.minWidth.maxLevel(radius()) - 1
    // If level < 0, more than three face cells are required.
    if level < 0 {
      let cellIds = (0..<6).map { CellId(face: $0) }
      return CellUnion(cellIds: cellIds)
    }
    // The covering consists of the 4 cells at the given level that share the
    // cell vertex that is closest to the cap center.
    let cellIds = CellId(point: center).vertexNeighbors(level)
    return CellUnion(cellIds: cellIds)
  }
  
  /// Increases the cap if necessary to include the given point. If this cap is empty,
  /// then the center is set to the point with a zero height. p must be unit-length.
  public func add(_ point: S2Point) -> S2Cap {
    if isEmpty {
      return S2Cap(center: point, height: S2Cap.zeroHeight)
    }
    // To make sure that the resulting cap actually includes this point,
    // we need to round up the distance calculation.  That is, after
    // calling cap.AddPoint(p), cap.Contains(p) should be true.
    let dist2 = center.v.sub(point.v).norm2
    let height = max(self.height, S2Cap.roundUp * 0.5 * dist2)
    return S2Cap(centerNormalized: center, height: height)
  }

  // MARK: arithmetic 
  
  /// Increases the cap height if necessary to include the other cap. If this cap is empty,
  /// it is set to the other cap.
  public func add(_ cap: S2Cap) -> S2Cap {
    if isEmpty {
      return cap
    }
    if cap.isEmpty {
      return self
    }
    //
    let radius = center.angle(cap.center) + cap.radius()
    let height = max(self.height, S2Cap.roundUp * S2Cap.radiusToHeight(radius))
    return S2Cap(centerNormalized: center, height: height)
  }

  /// Returns a new cap expanded by the given angle. If the cap is empty,
  /// it returns an empty cap.
  public func expanded(_ distance: S1Angle) -> S2Cap {
    if isEmpty {
      return S2Cap.empty
    }
    return S2Cap(centerNormalized: center, angle: radius() + distance)
  }

  /// Converts an S1Angle into the height of the cap.
  static func radiusToHeight(_ radius: S1Angle) -> Double {
    if radius < 0 {
      return S2Cap.emptyHeight
    }
    if radius >= .pi {
      return S2Cap.fullHeight
    }
    // The height of the cap can be computed as 1 - cos(r), but this isn't very
    // accurate for angles close to zero (where cos(r) is almost 1). The
    // formula below has much better precision.
    let d = sin(0.5 * radius)
    return 2 * d * d
  }

}

extension S2Cap: CustomStringConvertible, Approximatable {
  
  static let epsilon = 1e-14
  
  public var description: String {
    return "[Center=\(center.latLng), Radius=\(radius() * toDegrees)]"
  }
  
  /// Reports if this caps' center and height are within
  /// a reasonable epsilon from the other cap.
  public func approxEquals(_ cap: S2Cap) -> Bool {
    return center.approxEquals(cap.center) &&
      abs(height-cap.height) <= S2Cap.epsilon ||
      isEmpty && cap.height <= S2Cap.epsilon ||
      cap.isEmpty && height <= S2Cap.epsilon ||
      isFull && cap.height >= 2 - S2Cap.epsilon ||
      cap.isFull && height >= 2 - S2Cap.epsilon
  }
  
}
