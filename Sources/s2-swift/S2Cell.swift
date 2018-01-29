//
//  S2Cell.swift
//  s2-swift
//

import Foundation


// Cell is an S2 region object that represents a cell. Unlike CellIDs,
// it supports efficient containment and intersection tests. However, it is
// also a more expensive representation.
public struct Cell: S2Region {
  
  // TODO: move these package private variables to a more appropriate location.
  static let dblEpsilon = nextafter(1.0, 2.0) - 1.0
  static let poleMinLat = asin(sqrt(1.0 / 3.0)) - 0.5 * dblEpsilon

  //
  let face: UInt8
  let level: UInt8
  let orientation: UInt8
  let id: CellId
  let uv: R2Rect

  // MARK: inits / factory
  
  private init(face: UInt8, level: UInt8, orientation: UInt8, id: CellId, uv: R2Rect) {
    self.face = face
    self.level = level
    self.orientation = orientation
    self.id = id
    self.uv = uv
  }
  
  public init(id: CellId) {
    self.id = id
    level = UInt8(id.level())
    let (f, i, j, o) = id.faceIJOrientation()
    face = UInt8(f)
    orientation = UInt8(o)
    uv = CellId.ijLevelToBoundUV(i: i, j: j, level: Int(level))
  }

  public init(point: S2Point) {
    let cellId = CellId(point: point)
    self.init(id: cellId)
  }

  public init(latLng: LatLng) {
    let cellId = CellId(latLng: latLng)
    self.init(id: cellId)
  }

  // MARK: tests
  
  /// Returns whether this Cell is a leaf or not.
  func isLeaf() -> Bool {
    return level == UInt8(CellId.maxLevel)
  }

  // MARK: computed members
  
  /// Returns the CellID value for the cells level.
  func sizeIJ() -> Int {
    return CellId.sizeIJ(Int(level))
  }

  /// Returns the edge length of this cell in (s,t)-space.
  func sizeST() -> Double {
    return id.sizeST(level: Int(level))
  }
  
  /// Returns the bounds of this cell in (u,v)-space.
  func boundUV() -> R2Rect {
    return uv
  }

  /// Center returns the direction vector corresponding to the center in
  /// (s,t)-space of the given cell. This is the point at which the cell is
  /// divided into four subcells; it is not necessarily the centroid of the
  /// cell in (u,v)-space or (x,y,z)-space
  func center() -> S2Point {
    return S2Point(raw: id.rawPoint())
  }
  
  /// Returns the k-th vertex of the cell (k = [0,3]) in CCW order
  /// (lower left, lower right, upper right, upper left in the UV plane).
  func vertex(_ k: Int) -> S2Point {
    let u = uv.vertices[k].x
    let v = uv.vertices[k].y
    return S2Point(raw: S2Cube(face: Int(face), u: u, v: v).vector())
  }

  /// Returns the inward-facing normal of the great circle passing through
  /// the CCW ordered edge from vertex k to vertex k+1 (mod 4).
  func edge(_ k: Int) -> S2Point {
    switch k {
    case 0:
      return S2Cube.vNorm(face: Int(face), v: uv.y.lo, invert: false) // Bottom
    case 1:
      return S2Cube.uNorm(face: Int(face), u: uv.x.hi, invert: false) // Right
    case 2:
      return S2Cube.vNorm(face: Int(face), v: uv.y.hi, invert: true) // Top
    default:
      return S2Cube.uNorm(face: Int(face), u: uv.x.lo, invert: true) // Left
    }
  }

  /// Returns the area of this cell as accurately as possible.
  func exactArea() -> Double {
    let (v0, v1, v2, v3) = (vertex(0), vertex(1), vertex(2), vertex(3))
    return S2Point.pointArea(v0, v1, v2) + S2Point.pointArea(v0, v2, v3)
  }

  /// Returns the approximate area of this cell. This method is accurate
  /// to within 3% percent for all cell sizes and accurate to within 0.1% for cells
  /// at level 5 or higher (i.e. squares 350km to a side or smaller on the Earth's
  /// surface). It is moderately cheap to compute.
  func approxArea() -> Double {
    // All cells at the first two levels have the same area.
    if level < 2 {
      return averageArea()
    }
    // First, compute the approximate area of the cell when projected
    // perpendicular to its normal. The cross product of its diagonals gives
    // the normal, and the length of the normal is twice the projected area.
    let d0 = vertex(2).v - vertex(0).v
    let d1 = vertex(3).v - vertex(1).v
    let flatArea = 0.5 * d0.cross(d1).norm
    // Now, compensate for the curvature of the cell surface by pretending
    // that the cell is shaped like a spherical cap. The ratio of the
    // area of a spherical cap to the area of its projected disc turns out
    // to be 2 / (1 + sqrt(1 - r*r)) where r is the radius of the disc.
    // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
    // Here we set Pi*r*r == flatArea to find the equivalent disc.
    return flatArea * 2 / (1.0 + sqrt(1 - min(1.0 / .pi * flatArea, 1)))
  }

  /// Returns the average area of cells at the level of this cell.
  /// This is accurate to within a factor of 1.7.
  func averageArea() -> Double {
    return S2CellMetric.avgArea.value(Int(level))
  }

  // MARK: derive lat/lng from uv
  
  /// Returns the latitude of the cell vertex given by (i,j), where "i" and "j" are either 0 or 1.
  func latitude(i: Int, j: Int) -> Double {
//    let u: Double
//    let v: Double
//    switch (i, j) {
//    case (0, 0):
//      u = uv.x.lo
//      v = uv.y.lo
//    case (0, 1):
//      u = uv.x.lo
//      v = uv.y.hi
//    case (1, 0):
//      u = uv.x.hi
//      v = uv.y.lo
//    case (1, 1):
//      u = uv.x.hi
//      v = uv.y.hi
//    default:
//      fatalError("i and/or j is out of bound")
//    }
    let u = i==0 ? uv.x.lo : uv.x.hi
    let v = j==0 ? uv.y.lo : uv.y.hi
    let p = S2Point(raw: S2Cube(face: Int(face), u: u, v: v).vector())
    return p.lat
  }

  /// Returns the longitude of the cell vertex given by (i,j), where "i" and "j" are either 0 or 1.
  func longitude(i: Int, j: Int) -> Double {
//    let u: Double
//    let v: Double
//    switch (i, j) {
//    case (0, 0):
//      u = uv.x.lo
//      v = uv.y.lo
//    case (0, 1):
//      u = uv.x.lo
//      v = uv.y.hi
//    case (1, 0):
//      u = uv.x.hi
//      v = uv.y.lo
//    case (1, 1):
//      u = uv.x.hi
//      v = uv.y.hi
//    default:
//      fatalError("i and/or j is out of bound")
//    }
    let u = i==0 ? uv.x.lo : uv.x.hi
    let v = j==0 ? uv.y.lo : uv.y.hi
    let p = S2Point(raw: S2Cube(face: Int(face), u: u, v: v).vector())
    return p.lng
  }

  /// Returns the bounding rectangle of this cell.
  public func rectBound() -> S2Rect {
    if level > 0 {
      // Except for cells at level 0, the latitude and longitude extremes are
      // attained at the vertices.  Furthermore, the latitude range is
      // determined by one pair of diagonally opposite vertices and the
      // longitude range is determined by the other pair.
      //
      // We first determine which corner (i,j) of the cell has the largest
      // absolute latitude.  To maximize latitude, we want to find the point in
      // the cell that has the largest absolute z-coordinate and the smallest
      // absolute x- and y-coordinates.  To do this we look at each coordinate
      // (u and v), and determine whether we want to minimize or maximize that
      // coordinate based on the axis direction and the cell's (u,v) quadrant.
      let u = uv.x.lo + uv.x.hi
      let v = uv.y.lo + uv.y.hi
      var i = 0
      var j = 0
      if S2Cube.uAxis(face: Int(face)).z == 0 {
        if u < 0 {
          i = 1
        }
      } else if u > 0 {
        i = 1
      }
      if S2Cube.vAxis(face: Int(face)).z == 0 {
        if v < 0 {
          j = 1
        }
      } else if v > 0 {
        j = 1
      }
      let lat = R1Interval(point: latitude(i: i, j: j)).add(latitude(i: 1-i, j: 1-j))
      let lng = S1Interval.empty.add(longitude(i: i, j: 1-j)).add(longitude(i: 1-i, j: j))
      // We grow the bounds slightly to make sure that the bounding rectangle
      // contains LatLngFromPoint(P) for any point P inside the loop L defined by the
      // four *normalized* vertices.  Note that normalization of a vector can
      // change its direction by up to 0.5 * dblEpsilon radians, and it is not
      // enough just to add Normalize calls to the code above because the
      // latitude/longitude ranges are not necessarily determined by diagonally
      // opposite vertex pairs after normalization.
      //
      // We would like to bound the amount by which the latitude/longitude of a
      // contained point P can exceed the bounds computed above.  In the case of
      // longitude, the normalization error can change the direction of rounding
      // leading to a maximum difference in longitude of 2 * dblEpsilon.  In
      // the case of latitude, the normalization error can shift the latitude by
      // up to 0.5 * dblEpsilon and the other sources of error can cause the
      // two latitudes to differ by up to another 1.5 * dblEpsilon, which also
      // leads to a maximum difference of 2 * dblEpsilon.
      return S2Rect(lat: lat, lng: lng).expanded(LatLng(lat: 2 * Cell.dblEpsilon, lng: 2 * Cell.dblEpsilon)).polarClosure()
    }
    // The 4 cells around the equator extend to +/-45 degrees latitude at the
    // midpoints of their top and bottom edges.  The two cells covering the
    // poles extend down to +/-35.26 degrees at their vertices.  The maximum
    // error in this calculation is 0.5 * dblEpsilon.
    let bound: S2Rect
    switch face {
    case 0:
      bound = S2Rect(lat: R1Interval(lo: -.pi / 4, hi: .pi / 4), lng: S1Interval(lo: -.pi / 4, hi: .pi / 4))
    case 1:
      bound = S2Rect(lat: R1Interval(lo: -.pi / 4, hi: .pi / 4), lng: S1Interval(lo: .pi / 4, hi: 3 * .pi / 4))
    case 2:
      bound = S2Rect(lat: R1Interval(lo: Cell.poleMinLat, hi: .pi / 2), lng: S1Interval.full)
    case 3:
      bound = S2Rect(lat: R1Interval(lo: -.pi / 4, hi: .pi / 4), lng: S1Interval(lo: 3 * .pi / 4, hi: -3 * .pi / 4))
    case 4:
      bound = S2Rect(lat: R1Interval(lo: -.pi / 4, hi: .pi / 4), lng: S1Interval(lo: -3 * .pi / 4, hi: -.pi / 4))
    default:
      bound = S2Rect(lat: R1Interval(lo: -.pi / 2, hi: -Cell.poleMinLat), lng: S1Interval.full)
    }
    // Finally, we expand the bound to account for the error when a point P is
    // converted to a LatLng to test for containment. (The bound should be
    // large enough so that it contains the computed LatLng of any contained
    // point, not just the infinite-precision version.) We don't need to expand
    // longitude because longitude is calculated via a single call to math.Atan2,
    // which is guaranteed to be semi-monotonic.
    return bound.expanded(LatLng(lat: Cell.dblEpsilon, lng: 0.0))
  }

  /// Returns the bounding cap of this cell.
  public func capBound() -> S2Cap {
    // We use the cell center in (u,v)-space as the cap axis.  This vector is very close
    // to GetCenter() and faster to compute.  Neither one of these vectors yields the
    // bounding cap with minimal surface area, but they are both pretty close.
    let p = S2Point(raw: S2Cube(face: Int(face), u: uv.center.x, v: uv.center.y).vector())
    var cap = S2Cap(point: p)
    for k in 0..<4 {
      cap = cap.add(vertex(k))
    }
    return cap
  }

  // MARK: contains / intersects
  
  /// Reports whether the intersection of this cell and the other cell is not nil.
  public func intersects(_ cell: Cell) -> Bool {
    return id.intersects(cell.id)
  }
  
  // ContainsCell reports whether this cell contains the other cell.
  public func contains(_ cell: Cell) -> Bool {
    return id.contains(cell.id)
  }
  
  /// Reports whether this cell contains the given point. Note that
  /// unlike Loop/Polygon, a Cell is considered to be a closed set. This means
  /// that a point on a Cell's edge or vertex belong to the Cell and the relevant
  /// adjacent Cells too.
  /// If you want every point to be contained by exactly one Cell,
  /// you will need to convert the Cell to a Loop.
  public func contains(_ point: S2Point) -> Bool {
    guard let cube = S2Cube(point: point, face: Int(face)) else {
      return false
    }
    let uv2 = R2Point(x: cube.u, y: cube.v)
    // Expand the (u,v) bound to ensure that CellFromPoint(p).ContainsPoint(p)
    // is always true. To do this, we need to account for the error when
    // converting from (u,v) coordinates to (s,t) coordinates. In the
    // normal case the total error is at most dblEpsilon.
    return uv.expanded(Cell.dblEpsilon).contains(uv2)
  }

  // MARK: iterate
  
  /// Children returns the four direct children of this cell in traversal order
  /// and returns true. If this is a leaf cell, or the children could not be created,
  /// false is returned.
  /// The C++ method is called Subdivide.
  func children() -> [Cell]? {
    if id.isLeaf() {
      return nil
    }
    // Compute the cell midpoint in uv-space.
    let uvMid = id.centerUV()
    // Create four children with the appropriate bounds.
    var cid = id.childBegin()
    let children = (0..<4).map { (pos: Int) -> Cell in
      // We want to split the cell in half in u and v. To decide which
      // side to set equal to the midpoint value, we look at cell's (i,j)
      // position within its parent. The index for i is in bit 1 of ij.
      let ij = CellId.posToIJ[Int(orientation)][pos]
      let i = ij >> 1
      let j = ij & 1
      let xLo = (i == 1) ? uvMid.x : uv.x.lo
      let xHi = (i == 1) ? uv.x.hi : uvMid.x
      let x = R1Interval(lo: xLo, hi: xHi)
      let yLo = (j == 1) ? uvMid.y : uv.y.lo
      let yHi = (j == 1) ? uv.y.hi : uvMid.y
      let y = R1Interval(lo: yLo, hi: yHi)
      let cuv = R2Rect(x: x, y: y)
      let o = orientation ^ UInt8(CellId.posToOrientation[pos])
      let child = Cell(face: face, level: level + 1, orientation: o, id: cid, uv: cuv)
      cid = cid.next()
      return child
    }
    return children
  }

  // MARK: edge methods
  
  /// Returns the squared chord distance from point P to the
  /// given corner vertex specified by the Hi or Lo values of each.
  func vertexChordDist2(p: S2Point, xHi: Bool, yHi: Bool) -> Double {
    let x = xHi ? uv.x.hi : uv.x.lo
    let y = yHi ? uv.y.hi : uv.y.lo
    return p.sub(S2Point(x: x, y: y, z: 1)).norm2
  }
  
  // uEdgeIsClosest reports whether a point P is closer to the interior of the specified
  // Cell edge (either the lower or upper edge of the Cell) or to the endpoints.
  func uEdgeIsClosest(p: S2Point, vHi: Bool) -> Bool {
    let u0 = uv.x.lo
    let u1 = uv.x.hi
    let v = vHi ? uv.y.hi : uv.y.lo
    // These are the normals to the planes that are perpendicular to the edge
    // and pass through one of its two endpoints.
    let dir0 = R3Vector(x: v*v + 1, y: -u0 * v, z: -u0)
    let dir1 = R3Vector(x: v*v + 1, y: -u1 * v, z: -u1)
    return p.v.dot(dir0) > 0 && p.v.dot(dir1) < 0
  }
  
  // vEdgeIsClosest reports whether a point P is closer to the interior of the specified
  // Cell edge (either the right or left edge of the Cell) or to the endpoints.
  func vEdgeIsClosest(p: S2Point, uHi: Bool) -> Bool {
    let v0 = uv.y.lo
    let v1 = uv.y.hi
    let u = uHi ? uv.x.hi : uv.x.lo
    let dir0 = R3Vector(x: -u * v0, y: u*u + 1, z: -v0)
    let dir1 = R3Vector(x: -u * v1, y: u*u + 1, z: -v1)
    return p.v.dot(dir0) > 0 && p.v.dot(dir1) < 0
  }
  
  // edgeDistance reports the distance from a Point P to a given Cell edge. The point
  // P is given by its dot product, and the uv edge by its normal in the
  // given coordinate value.
  func edgeDistance(ij: Double, uv: Double) -> S1ChordAngle {
    // Let P by the target point and let R be the closest point on the given
    // edge AB.  The desired distance PR can be expressed as PR^2 = PQ^2 + QR^2
    // where Q is the point P projected onto the plane through the great circle
    // through AB.  We can compute the distance PQ^2 perpendicular to the plane
    // from "dirIJ" (the dot product of the target point P with the edge
    // normal) and the squared length the edge normal (1 + uv**2).
    let pq2  = (ij * ij) / (1 + uv * uv)
    // We can compute the distance QR as (1 - OQ) where O is the sphere origin,
    // and we can compute OQ^2 = 1 - PQ^2 using the Pythagorean theorem.
    // (This calculation loses accuracy as angle POQ approaches Pi/2.)
    let qr = 1 - sqrt(1 - pq2)
    return S1ChordAngle(squaredLength: pq2 + qr * qr)
  }
  
  // distanceInternal reports the distance from the given point to the interior of
  // the cell if toInterior is true or to the boundary of the cell otherwise.
  func distanceInternal(targetXYZ: S2Point, toInterior: Bool) -> S1ChordAngle {
    // All calculations are done in the (u,v,w) coordinates of this cell's face.
    let target = S2Cube.faceXYZtoUVW(face: Int(face), point: targetXYZ)
    // Compute dot products with all four upward or rightward-facing edge
    // normals. dirIJ is the dot product for the edge corresponding to axis
    // I, endpoint J. For example, dir01 is the right edge of the Cell
    // (corresponding to the upper endpoint of the u-axis).
    let dir00 = target.x - target.z * uv.x.lo
    let dir01 = target.x - target.z * uv.x.hi
    let dir10 = target.y - target.z * uv.y.lo
    let dir11 = target.y - target.z * uv.y.hi
    var inside = true
    if dir00 < 0 {
      inside = false // Target is to the left of the cell
      if vEdgeIsClosest(p: target, uHi: false) {
        return edgeDistance(ij: -dir00, uv: uv.x.lo)
      }
    }
    if dir01 > 0 {
      inside = false // Target is to the right of the cell
      if vEdgeIsClosest(p: target, uHi: true) {
        return edgeDistance(ij: dir01, uv: uv.x.hi)
      }
    }
    if dir10 < 0 {
      inside = false // Target is below the cell
      if uEdgeIsClosest(p: target, vHi: false) {
        return edgeDistance(ij: -dir10, uv: uv.y.lo)
      }
    }
    if dir11 > 0 {
      inside = false // Target is above the cell
      if uEdgeIsClosest(p: target, vHi: true) {
        return edgeDistance(ij: dir11, uv: uv.y.hi)
      }
    }
    if inside {
      if toInterior {
        return S1ChordAngle.zero
      }
      // Although you might think of Cells as rectangles, they are actually
      // arbitrary quadrilaterals after they are projected onto the sphere.
      // Therefore the simplest approach is just to find the minimum distance to
      // any of the four edges.
      return min(edgeDistance(ij: -dir00, uv: uv.x.lo),
                 edgeDistance(ij: dir01, uv: uv.x.hi),
                 edgeDistance(ij: -dir10, uv: uv.y.lo),
                 edgeDistance(ij: dir11, uv: uv.y.hi))
    }
    // Otherwise, the closest point is one of the four cell vertices. Note that
    // it is *not* trivial to narrow down the candidates based on the edge sign
    // tests above, because (1) the edges don't meet at right angles and (2)
    // there are points on the far side of the sphere that are both above *and*
    // below the cell, etc.
    let chordDist2 = min(vertexChordDist2(p: target, xHi: false, yHi: false),
                         vertexChordDist2(p: target, xHi: true, yHi: false),
                         vertexChordDist2(p: target, xHi: false, yHi: true),
                         vertexChordDist2(p: target, xHi: true, yHi: true))
    return S1ChordAngle(squaredLength: chordDist2)
  }
  
  // Distance reports the distance from the cell to the given point. Returns zero if
  // the point is inside the cell.
  func distance(target: S2Point) -> S1ChordAngle {
    return distanceInternal(targetXYZ: target, toInterior: true)
  }
  
  // BoundaryDistance reports the distance from the cell boundary to the given point.
  func boundaryDistance(target: S2Point) -> S1ChordAngle {
    return distanceInternal(targetXYZ: target, toInterior: false)
  }
  
}

extension Cell {
  
  // DistanceToEdge returns the minimum distance from the cell to the given edge AB. Returns
  // zero if the edge intersects the cell interior.
  func distanceToEdge(a: S2Point, b: S2Point) -> S1ChordAngle {
    // Possible optimizations:
    //  - Currently the (cell vertex, edge endpoint) distances are computed
    //    twice each, and the length of AB is computed 4 times.
    //  - To fix this, refactor GetDistance(target) so that it skips calculating
    //    the distance to each cell vertex. Instead, compute the cell vertices
    //    and distances in this function, and add a low-level UpdateMinDistance
    //    that allows the XA, XB, and AB distances to be passed in.
    //  - It might also be more efficient to do all calculations in UVW-space,
    //    since this would involve transforming 2 points rather than 4.
    // First, check the minimum distance to the edge endpoints A and B.
    // (This also detects whether either endpoint is inside the cell.)
    var minDist = min(distance(target: a), distance(target: b))
    if minDist == .zero {
      return minDist
    }
    // Otherwise, check whether the edge crosses the cell boundary.
    var crosser = EdgeCrosser(a: a, b: b, c: vertex(3))
    for i in 0..<4 {
      if crosser.chainCrossingSign(d: vertex(i)) != .doNotCross {
        return .zero
      }
    }
    // Finally, check whether the minimum distance occurs between a cell vertex
    // and the interior of the edge AB. (Some of this work is redundant, since
    // it also checks the distance to the endpoints A and B again.)
    //
    // Note that we don't need to check the distance from the interior of AB to
    // the interior of a cell edge, because the only way that this distance can
    // be minimal is if the two edges cross (already checked above).
    for i in 0..<4 {
      (minDist, _) = updateMinDistance(x: vertex(i), a: a, b: b, minDist: minDist)
    }
    return minDist
  }
  
}

extension Cell: Equatable, Comparable {
  
  public static func ==(lhs: Cell, rhs: Cell) -> Bool {
    return lhs.id == rhs.id
  }
  
  public static func <(lhs: Cell, rhs: Cell) -> Bool {
    return lhs.id < rhs.id
  }
  
}

