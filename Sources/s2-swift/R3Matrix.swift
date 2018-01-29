//
//  R3Matrix.swift
//  s2-swift
//

import Foundation


/// Represents a traditional 3x3 matrix of floating point values.
/// This is not a full fledged matrix. It only contains the pieces needed
/// to satisfy the computations done within the s2 package.
public struct R3Matrix {
  
  let rows = 3
  let columns = 3
  var m: [Double]
  
  public init() {
    m = [Double](repeating: 0.0, count: rows * columns)
  }

  public init(r1: [Double], r2: [Double], r3: [Double]) {
    assert(r1.count == 3)
    assert(r2.count == 3)
    assert(r3.count == 3)
    m = r1 + r2 + r3
  }
  
  // MARK: subscript
  
  func indexIsValidForRow(_ row: Int, column: Int) -> Bool {
    return row >= 0 && row < rows && column >= 0 && column < columns
  }
  
  subscript(row: Int, column: Int) -> Double {
    get {
      assert(indexIsValidForRow(row, column: column), "Index out of range")
      return m[row * columns + column]
    }
    set {
      assert(indexIsValidForRow(row, column: column), "Index out of range")
      m[row * columns + column] = newValue
    }
  }

  // MARK: accessors
  
  /// Returns the given column as a Point.
  func col(_ col: Int) -> S2Point {
    return S2Point(x: self[0, col], y: self[1, col], z: self[2, col])
  }

  /// Returns the given row as a Point.
  func row(_ row: Int) -> S2Point {
    return S2Point(x: self[row, 0], y: self[row, 1], z: self[row, 2])
  }

  /// Multiplies the matrix by the given value.
  func scale(scalar f: Double) -> R3Matrix {
    return R3Matrix(
      r1: [f * self[0, 0], f * self[0, 1], f * self[0, 2]],
      r2: [f * self[1, 0], f * self[1, 1], f * self[1, 2]],
      r3: [f * self[2, 0], f * self[2, 1], f * self[2, 2]])
  }

  /// Returns the multiplication of m by the Point p and converts the
  /// resulting 1x3 matrix into a Point.
  func mul(point p: S2Point) -> S2Point {
    let x = self[0, 0]*p.x + self[0, 1]*p.y + self[0, 2]*p.z
    let y = self[1, 0]*p.x + self[1, 1]*p.y + self[1, 2]*p.z
    let z = self[2, 0]*p.x + self[2, 1]*p.y + self[2, 2]*p.z
    return S2Point(x: x, y: y, z: z)
  }

  /// Returns the determinant of this matrix.
  func det() -> Double {
    //      | a  b  c |
    //  det | d  e  f | = aei + bfg + cdh - ceg - bdi - afh
    //      | g  h  i |
    let aei = self[0, 0] * self[1, 1] * self[2, 2]
    let bfg = self[0, 1] * self[1, 2] * self[2, 0]
    let cdh = self[0, 2] * self[1, 0] * self[2, 1]
    let ceg = self[0, 2] * self[1, 1] * self[2, 0]
    let bdi = self[0, 1] * self[1, 0] * self[2, 2]
    let afh = self[0, 0] * self[1, 2] * self[2, 1]
    return aei + bfg + cdh - ceg - bdi - afh
  }

  /// Reflects the matrix along its diagonal and returns the result.
  func transpose() -> R3Matrix {
    return R3Matrix(
      r1: [self[0, 0], self[1, 0], self[2, 0]],
      r2: [self[0, 1], self[1, 1], self[2, 1]],
      r3: [self[0, 2], self[1, 2], self[2, 2]])
  }

}

extension R3Matrix: CustomStringConvertible {

  public var description: String {
    return String(format: "[ %0.4f %0.4f %0.4f ] [ %0.4f %0.4f %0.4f ] [ %0.4f %0.4f %0.4f ]",
      self[0, 0], self[0, 1], self[0, 2],
      self[1, 0], self[1, 1], self[1, 2],
      self[2, 0], self[2, 1], self[2, 2])
  }

}
