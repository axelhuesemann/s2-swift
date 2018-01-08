//
//  R3VectorTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class R3R3VectorTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func v(_ x: Double, _ y: Double, _ z: Double) -> R3Vector {
    return R3Vector(x: x, y: y, z: z)
  }
  
  func testNorm() {
    XCTAssertEqual(v(0, 0, 0).norm, 0.0, accuracy: 1e-14)
    XCTAssertEqual(v(0, 1, 0).norm, 1.0, accuracy: 1e-14)
    XCTAssertEqual(v(3, -4, 12).norm, 13.0, accuracy: 1e-14)
    XCTAssertEqual(v(1, 1e-16, 1e-32).norm, 1.0, accuracy: 1e-14)
  }

  func testNorm2() {
    XCTAssertEqual(v(0, 0, 0).norm2, 0.0, accuracy: 1e-14)
    XCTAssertEqual(v(0, 1, 0).norm2, 1.0, accuracy: 1e-14)
    XCTAssertEqual(v(1, 1, 1).norm2, 3.0, accuracy: 1e-14)
    XCTAssertEqual(v(1, 2, 3).norm2, 14.0, accuracy: 1e-14)
    XCTAssertEqual(v(3, -4, 12).norm2, 169.0, accuracy: 1e-14)
    XCTAssertEqual(v(1, 1e-16, 1e-32).norm2, 1.0, accuracy: 1e-14)
  }

  func testnormalized() {
    let vectors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 1e-16, 1e-32), (12.34, 56.78, 91.01)]
    for v1 in vectors {
      let v = R3Vector(x: v1.0, y: v1.1, z: v1.2)
      let nv = v.normalized()
      XCTAssertEqual(v.x*nv.y, v.y*nv.x, accuracy: 1e-14)
      XCTAssertEqual(v.x*nv.z, v.z*nv.x, accuracy: 1e-14)
      XCTAssertEqual(nv.norm, 1.0, accuracy: 1e-14)
    }
  }
  
  func testIsUnit() {
    let tests = [
      (v(0, 0, 0), false),
      (v(0, 1, 0), true),
      (v(1 + 2.0 * 1e-15, 0, 0), true),
      (v(1 * (1.0 + 1e-15), 0, 0), true),
      (v(1, 1, 1), false),
      (v(1, 1e-16, 1e-32), true)]
    for (v, want) in tests {
      XCTAssertEqual(v.isUnit, want)
    }
  }

  func testDot() {
    let tests = [
      (v(1, 0, 0), v(1, 0, 0), 1),
      (v(1, 0, 0), v(0, 1, 0), 0),
      (v(1, 0, 0), v(0, 1, 1), 0),
      (v(1, 1, 1), v(-1, -1, -1), -3),
      (v(1, 2, 2), v(-0.3, 0.4, -1.2), -1.9)]
    for (v1, v2, want) in tests {
      let v1 = v(v1.x, v1.y, v1.z)
      let v2 = v(v2.x, v2.y, v2.z)
      XCTAssertEqual(v1.dot(v2), want, accuracy: 1e-15)
      XCTAssertEqual(v2.dot(v1), want, accuracy: 1e-15)
    }
  }
      
  func testCross() {
    let tests = [
      (v(1, 0, 0), v(1, 0, 0), v(0, 0, 0)),
      (v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)),
      (v(0, 1, 0), v(1, 0, 0), v(0, 0, -1)),
      (v(1, 2, 3), v(-4, 5, -6), v(-27, -6, 13))]
    for (v1, v2, want) in tests {
      XCTAssert(v1.cross(v2).approxEquals(want))
    }
  }
  
  func testAdd() {
    let tests = [
      (v(0, 0, 0), v(0, 0, 0), v(0, 0, 0)),
      (v(1, 0, 0), v(0, 0, 0), v(1, 0, 0)),
      (v(1, 2, 3), v(4, 5, 7), v(5, 7, 10)),
      (v(1, -3, 5), v(1, -6, -6), v(2, -9, -1))]
    for (v1, v2, want) in tests {
      XCTAssert(v1.add(v2).approxEquals(want))
    }
  }
  
  func testSub() {
    let tests = [
      (v(0, 0, 0), v(0, 0, 0), v(0, 0, 0)),
      (v(1, 0, 0), v(0, 0, 0), v(1, 0, 0)),
      (v(1, 2, 3), v(4, 5, 7), v(-3, -3, -4)),
      (v(1, -3, 5), v(1, -6, -6), v(0, 3, 11))]
    for (v1, v2, want) in tests {
      XCTAssert(v1.sub(v2).approxEquals(want))
    }
  }
  
  func testDistance() {
    let tests = [
      (v(1, 0, 0), v(1, 0, 0), 0),
      (v(1, 0, 0), v(0, 1, 0), 1.41421356237310),
      (v(1, 0, 0), v(0, 1, 1), 1.73205080756888),
      (v(1, 1, 1), v(-1, -1, -1), 3.46410161513775),
      (v(1, 2, 2), v(-0.3, 0.4, -1.2), 3.80657326213486)]
    for (v1, v2, want) in tests {
      let v1 = v(v1.x, v1.y, v1.z)
      let v2 = v(v2.x, v2.y, v2.z)
      XCTAssertEqual(v1.distance(v2), want, accuracy: 1e-13)
      XCTAssertEqual(v1.distance(v2), want, accuracy: 1e-13)
    }
  }
  
  func testMul() {
    let tests = [
      (v(0, 0, 0), 3.0, v(0, 0, 0)),
      (v(1, 0, 0), 1, v(1, 0, 0)),
      (v(1, 0, 0), 0, v(0, 0, 0)),
      (v(1, 0, 0), 3, v(3, 0, 0)),
      (v(1, -3, 5), -1, v(-1, 3, -5)),
      (v(1, -3, 5), 2, v(2, -6, 10))]
    for (v, m, want) in tests {
      XCTAssert(v.mul(m).approxEquals(want))
    }
  }
  
  func testAngle() {
    let tests = [
      (v(1, 0, 0), v(1, 0, 0), 0),
      (v(1, 0, 0), v(0, 1, 0), .pi / 2),
      (v(1, 0, 0), v(0, 1, 1), .pi / 2),
      (v(1, 0, 0), v(-1, 0, 0), .pi),
      (v(1, 2, 3), v(2, 3, -1), 1.2055891055045298)]
    for (v1, v2, want) in tests {
      XCTAssertEqual(v1.angle(v2), want, accuracy: 1e-15)
      XCTAssertEqual(v2.angle(v1), want, accuracy: 1e-15)
    }
  }
  
  func testOrtho() {
    let vectors = [
      v(1, 0, 0),
      v(1, 1, 0),
      v(1, 2, 3),
      v(1, -2, -5),
      v(0.012, 0.0053, 0.00457),
      v(-0.012, -1, -0.00457)]
    for v in vectors {
      XCTAssertEqual(v.dot(v.ortho()), 0, accuracy: 1e-15)
      XCTAssertEqual(v.ortho().norm, 1, accuracy: 1e-15)
    }
  }
  
  func testIdentities() {
    let tests = [
      (v(0, 0, 0), v(0, 0, 0)),
      (v(0, 0, 0), v(0, 1, 2)),
      (v(1, 0, 0), v(0, 1, 0)),
      (v(1, 0, 0), v(0, 1, 1)),
      (v(1, 1, 1), v(-1, -1, -1)),
      (v(1, 2, 2), v(-0.3, 0.4, -1.2))]
    for (v1, v2) in tests {
      let a1 = v1.angle(v2)
      let a2 = v2.angle(v1)
      let c1 = v1.cross(v2)
      let c2 = v2.cross(v1)
      let d1 = v1.dot(v2)
      let d2 = v2.dot(v1)
      // Angle commutes
      XCTAssertEqual(a1, a2, accuracy: 1e-15)
      // Dot commutes
      XCTAssertEqual(d1, d2, accuracy: 1e-15)
      // Cross anti-commutes
      XCTAssert(c1.approxEquals(c2.mul(-1.0)))
      // Cross is orthogonal to original vectors
      XCTAssertEqual(v1.dot(c1), 0.0, accuracy: 1e-15)
      XCTAssertEqual(v2.dot(c1), 0.0, accuracy: 1e-15)
    }
  }

}
