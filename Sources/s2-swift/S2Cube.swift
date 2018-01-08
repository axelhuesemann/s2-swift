//
//  S2Cube.swift
//  s2-swift
//

import Foundation


/// Helper class for Call and CellId. Needs not be published.
/// Represents a projection of S2 coordinates onto a unit cube with 6 faces (of course)
/// and an edge length of 1.
struct S2Cube {
  
  // projection of a point on the sphere onto a bounding cube.
  // faces 0,1,3,4 are around the equator, 2,5 are at the poles
  let face: Int
  // u and v are the coordinates on the respective cube face, and should be in [-1, 1]x[-1, 1]
  let u: Double
  let v: Double

  // MARK: inits/ factory
  
  init(face: Int, u: Double, v: Double) {
    self.face = face
    self.u = u
    self.v = v
  }
  
  init(point: S2Point) {
    let face = S2Cube.face(point: point)
    let (u, v) = S2Cube.validFaceXYZToUV(face: face, point: point)
    self.init(face: face, u: u, v: v)
  }
  
  init?(point: S2Point, face: Int) {
    let toCheck: Double
    switch face {
    case 0: toCheck = point.x
    case 1: toCheck = point.y
    case 2: toCheck = point.z
    case 3: toCheck = -point.x
    case 4: toCheck = -point.y
    case 5: toCheck = -point.z
    default: return nil
    }
    if toCheck <= 0 {
      return nil
    }
    let (u, v) = S2Cube.validFaceXYZToUV(face: face, point: point)
    self.init(face: face, u: u, v: v)
  }

  // MARK: non-linear projection
  
  // stToUV converts an s or t value to the corresponding u or v value.
  // This is a non-linear transformation from [-1,1] to [-1,1] that
  // attempts to make the cell sizes more uniform.
  // This uses what the C++ version calls 'the quadratic transform'.
  static func stToUV(_ s: Double) -> Double {
    if s >= 0.5 {
      return (1.0 / 3.0) * (4 * s * s - 1)
    }
    return (1.0 / 3.0) * (1 - 4 * (1 - s) * (1 - s))
  }
  
  // uvToST is the inverse of the stToUV transformation. Note that it
  // is not always true that uvToST(stToUV(x)) == x due to numerical
  // errors.
  static func uvToST(_ u: Double) -> Double {
    if u >= 0 {
      return 0.5 * sqrt(1 + 3 * u)
    }
    return 1 - 0.5 * sqrt(1 - 3 * u)
  }
  
  // MARK: face
  
  // face returns face ID from 0 to 5 containing the r. For points on the
  // boundary between faces, the result is arbitrary but deterministic.
  static func face(point r: S2Point) -> Int {
    // find largest component
    var id = 0
    var value = r.x
    if abs(r.y) > abs(value) {
      id = 1
      value = r.y
    }
    if abs(r.z) > abs(value) {
      id = 2
      value = r.z
    }
    // negative means the opposite face
    if value < 0.0 {
      id += 3
    }
    //
    return id
  }
  
  // validFaceXYZToUV given a valid face for the given point r (meaning that
  // dot product of r with the face normal is positive), returns
  // the corresponding u and v values, which may lie outside the range [-1,1].
  static func validFaceXYZToUV(face: Int, point r: S2Point) -> (Double, Double) {
    switch face {
    case 0: return (r.y / r.x, r.z / r.x)
    case 1: return (-r.x / r.y, r.z / r.y)
    case 2: return (-r.x / r.z, -r.y / r.z)
    case 3: return (r.z / r.x, r.y / r.x)
    case 4: return (r.z / r.y, -r.x / r.y)
    default: return (-r.y / r.z, -r.x / r.z)
    }
  }
  
  // xyzToFaceUV converts a direction S2Point (not necessarily unit length) to
  // (face, u, v) coordinates.
  static func xyzToFaceUV(r: S2Point) -> (Int, Double, Double) {
    let f = face(point: r)
    let (u, v) = validFaceXYZToUV(face: f, point: r)
    return (f, u, v)
  }
  
  // faceXYZToUV returns the u and v values (which may lie outside the range
  // [-1, 1]) if the dot product of the point p with the given face normal is positive.
  static func faceXYZToUV(face: Int, p: S2Point) -> (Double, Double)? {
    let toCheck: Double
    switch face {
    case 0: toCheck = p.x
    case 1: toCheck = p.y
    case 2: toCheck = p.z
    case 3: toCheck = -p.x
    case 4: toCheck = -p.y
    case 5: toCheck = -p.z
    default: return nil
    }
    if toCheck <= 0 {
      return nil
    }
    return validFaceXYZToUV(face: face, point: p)
  }

  // faceUVToXYZ turns face and UV coordinates into an unnormalized 3 R3Vector.
  func vector() -> R3Vector {
    switch face {
    case 0: return R3Vector(x: 1, y: u, z: v)
    case 1: return R3Vector(x: -u, y: 1, z: v )
    case 2: return R3Vector(x: -u, y: -v, z: 1)
    case 3: return R3Vector(x: -1, y: -v, z: -u)
    case 4: return R3Vector(x: v, y: -1, z: -u)
    default: return R3Vector(x: v, y: u, z: -1)
    }
  }
  
  // faceUVToXYZ turns face and UV coordinates into an unnormalized 3 R3Vector.
  static func faceUVToXYZ(face: Int, u: Double, v: Double) -> R3Vector {
    switch face {
    case 0: return R3Vector(x: 1, y: u, z: v)
    case 1: return R3Vector(x: -u, y: 1, z: v )
    case 2: return R3Vector(x: -u, y: -v, z: 1)
    case 3: return R3Vector(x: -1, y: -v, z: -u)
    case 4: return R3Vector(x: v, y: -1, z: -u)
    default: return R3Vector(x: v, y: u, z: -1)
    }
  }
  
  // faceXYZtoUVW transforms the given point P to the (u,v,w) coordinate frame of the given
  // face where the w-axis represents the face normal.
  static func faceXYZtoUVW(face: Int, point p: S2Point) -> S2Point {
    // The result coordinates are simply the dot products of P with the (u,v,w)
    // axes for the given face (see faceUVWAxes).
    switch face {
    case 0: return S2Point(x: p.y, y: p.z, z: p.x)
    case 1: return S2Point(x: -p.x, y: p.z, z: p.y)
    case 2: return S2Point(x: -p.x, y: -p.y, z: p.z)
    case 3: return S2Point(x: -p.z, y: -p.y, z: -p.x)
    case 4: return S2Point(x: -p.z, y: p.x, z: -p.y)
    default: return S2Point(x: p.y, y: p.x, z: -p.z)
    }
  }
  
  // MARK: u and v normals
  
  // uNorm returns the right-handed normal (not necessarily unit length) for an
  // edge in the direction of the positive v-axis at the given u-value on
  // the given face.  (This S2Point is perpendicular to the plane through
  // the sphere origin that contains the given edge.)
  static func uNorm(face: Int, u: Double, invert: Bool) -> S2Point {
    let one = invert ? -1.0 : 1.0
    let u = invert ? -u : u
    switch face {
    case 0: return S2Point(x: u, y: -one, z: 0)
    case 1: return S2Point(x: one, y: u, z: 0)
    case 2: return S2Point(x: one, y: 0, z: u)
    case 3: return S2Point(x: -u, y: 0, z: one)
    case 4: return S2Point(x: 0, y: -u, z: one)
    default: return S2Point(x: 0, y: -one, z: -u)
    }
  }
  
  // vNorm returns the right-handed normal (not necessarily unit length) for an
  // edge in the direction of the positive u-axis at the given v-value on
  // the given face.
  static func vNorm(face: Int, v: Double, invert: Bool) -> S2Point {
    let one = invert ? -1.0 : 1.0
    let v = invert ? -v : v
    switch face {
    case 0: return S2Point(x: -v, y: 0, z: one)
    case 1: return S2Point(x: 0, y: -v, z: one)
    case 2: return S2Point(x: 0, y: -one, z: -v)
    case 3: return S2Point(x: v, y: -one, z: 0)
    case 4: return S2Point(x: one, y: v, z: 0)
    default: return S2Point(x: one, y: 0, z: v)
    }
  }
  
  // MARK: face axes
  
  // faceUVWAxes are the U, V, and W axes for each face.
  static var faceUVWAxes = [
    [S2Point(x: 0, y: +1, z: 0), S2Point(x: 0, y: 0, z: +1), S2Point(x: +1, y: 0, z: 0)],
    [S2Point(x: -1, y: 0, z: 0), S2Point(x: 0, y: 0, z: +1), S2Point(x: 0, y: +1, z: 0)],
    [S2Point(x: -1, y: 0, z: 0), S2Point(x: 0, y: -1, z: 0), S2Point(x: 0, y: 0, z: +1)],
    [S2Point(x: 0, y: 0, z: -1), S2Point(x: 0, y: -1, z: 0), S2Point(x: -1, y: 0, z: 0)],
    [S2Point(x: 0, y: 0, z: -1), S2Point(x: +1, y: 0, z: 0), S2Point(x: 0, y: -1, z: 0)],
    [S2Point(x: 0, y: +1, z: 0), S2Point(x: +1, y: 0, z: 0), S2Point(x: 0, y: 0, z: -1)]]
  
  // uvwAxis returns the given axis of the given face.
  static func uvwAxis(face: Int, axis: Int) -> S2Point {
    return faceUVWAxes[face][axis]
  }
  
  // uAxis returns the u-axis for the given face.
  static func uAxis(face: Int) -> S2Point {
    return uvwAxis(face: face, axis: 0)
  }
  
  // vAxis returns the v-axis for the given face.
  static func vAxis(face: Int) -> S2Point {
    return uvwAxis(face: face, axis: 1)
  }
  
  // Return the unit-length normal for the given face.
  static func unitNorm(face: Int) -> S2Point {
    return uvwAxis(face: face, axis: 2)
  }
  
}
