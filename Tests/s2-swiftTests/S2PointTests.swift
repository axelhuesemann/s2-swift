//
//  S2PointTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2PointTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func p(_ x: Double, _ y: Double, _ z: Double) -> S2Point {
    return S2Point(x: x, y: y, z: z)
  }
  
  func p(_ lat: Double, _ lng: Double) -> S2Point {
    return llDegrees(lat, lng).toPoint()
  }
  
  func testOriginPoint() {
    XCTAssertEqual(S2Point.origin.v.norm, 1.0, accuracy: 1e-16)
  }
  
  func testPointCross() {
    let tests = [
      (1.0, 0.0, 0.0, 1.0, 0.0, 0.0),
      (1, 0, 0, 0, 1, 0),
      (0, 1, 0, 1, 0, 0),
      (1, 2, 3, -4, 5, -6)]
    for (p1x, p1y, p1z, p2x, p2y, p2z) in tests {
      let p1 = S2Point(x: p1x, y: p1y, z: p1z)
      let p2 = S2Point(x: p2x, y: p2y, z: p2z)
      let result = p1.pointCross(p2)
      XCTAssertEqual(result.v.norm, 1, accuracy: 1e-15)
      XCTAssertEqual(result.v.dot(p1.v), 0, accuracy: 1e-15)
      XCTAssertEqual(result.v.dot(p2.v), 0, accuracy: 1e-15)
    }
  }

  func testSign() {
    let tests: [(Double, Double, Double, Double, Double, Double, Double, Double, Double, Bool)] = [
      (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, true),
      (0, 1, 0, 0, 0, 1, 1, 0, 0, true),
      (0, 0, 1, 1, 0, 0, 0, 1, 0, true),
      (1, 1, 0, 0, 1, 1, 1, 0, 1, true),
      (-3, -1, 4, 2, -1, -3, 1, -2, 0, true),
      (-3, -1, 0, -2, 1, 0, 1, -2, 0, false),
      (-6, 3, 3, -4, 2, -1, -2, 1, 4, false),
      (0, -1, -1, 0, 1, -2, 0, 2, 1, false),
      (-1, 2, 7, 2, 1, -4, 4, 2, -8, false),
      (-4, -2, 7, 2, 1, -4, 4, 2, -8, false),
      (0, -5, 7, 0, -4, 8, 0, -2, 4, false),
      (-5, -2, 7, 0, 0, -2, 0, 0, -1, false),
      (0, -2, 7, 0, 0, 1, 0, 0, 2, false)]
    for (p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, want) in tests {
      let p1 = S2Point(x: p1x, y: p1y, z: p1z)
      let p2 = S2Point(x: p2x, y: p2y, z: p2z)
      let p3 = S2Point(x: p3x, y: p3y, z: p3z)
      let result = S2Point.sign(p1, b: p2, c: p3)
      XCTAssertEqual(result, want)
      // For these cases we can test the reversibility condition
      if want {
        XCTAssertNotEqual(S2Point.sign(p3, b: p2, c: p1), want)
      }
    }
  }
  
  // Points used in the various RobustSign tests.
  let x = S2Point(x: 1, y: 0, z: 0)
  let y = S2Point(x: 0, y: 1, z: 0)
  let z = S2Point(x: 0, y: 0, z: 1)
  
  // The following points happen to be *exactly collinear* along a line that it
  // approximate tangent to the surface of the unit sphere. In fact, C is the
  // exact midpoint of the line segment AB. All of these points are close
  // enough to unit length to satisfy r3.v.IsUnit.
  let poA = S2Point(x: 0.72571927877036835, y: 0.46058825605889098, z: 0.51106749730504852, normalize: false)
  let poB = S2Point(x: 0.7257192746638208, y: 0.46058826573818168, z: 0.51106749441312738, normalize: false)
  let poC = S2Point(x: 0.72571927671709457, y: 0.46058826089853633, z: 0.51106749585908795, normalize: false)
  
  // The points "x1" and "x2" are exactly proportional, i.e. they both lie
  // on a common line through the origin. Both points are considered to be
  // normalized, and in fact they both satisfy (x == x.normalized()).
  // Therefore the triangle (x1, x2, -x1) consists of three distinct points
  // that all lie on a common line through the origin.
  let x1 = S2Point(x: 0.99999999999999989, y: 1.4901161193847655e-08, z: 0)
  let x2 = S2Point(x: 1, y: 1.4901161193847656e-08, z: 0)
  
  // Here are two more points that are distinct, exactly proportional, and
  // that satisfy (x == x.normalized()).
  let x3 = S2Point(x: 1, y: 1, z: 1)
  let x4 = S2Point(x: 1, y: 1, z: 1).v.mul(0.99999999999999989).s2
  
  // The following three points demonstrate that normalized() is not idempotent, i.e.
  // y0.normalized() != y0.normalized().normalized(). Both points are exactly proportional.
  let y0 = S2Point(x: 1, y: 1, z: 0)
  let y1 = S2Point(x: 1, y: 1, z: 0)
  let y2 = S2Point(x: 1, y: 1, z: 0).v.normalized().s2
  
  // TODO: This test is missing the actual RobustSign() parts of the checks from C++
  // test method RobustSign::CollinearPoints.
  func testRobustSignEqualities() {
    XCTAssertEqual(poC.v.sub(poA.v), poB.v.sub(poC.v))
    XCTAssertEqual(poC.v.sub(poA.v).s2, poB.v.sub(poC.v).s2)
    XCTAssertEqual(x1, x1.v.normalized().s2)
    XCTAssertEqual(x2, x2.v.normalized().s2)
    XCTAssertEqual(x3, x3.v.normalized().s2)
    XCTAssertEqual(x4, x4.v.normalized().s2)
    XCTAssertNotEqual(x3, x4)
    XCTAssertNotEqual(y1, y2)
    XCTAssertEqual(y2, y2.v.normalized().s2)
  }

  func testRobustSign() {
    let tests = [
      // Simple collinear points test cases.
      // a == b != c
      (x, x, z, Direction.indeterminate),
      // a != b == c
      (x, y, y, .indeterminate),
      // c == a != b
      (z, x, z, .indeterminate),
      // CCW
      (x, y, z, .counterClockwise),
      // CW
      (z, y, x, .clockwise),
      // Edge cases:
      // The following points happen to be *exactly collinear* along a line that it
      // approximate tangent to the surface of the unit sphere. In fact, C is the
      // exact midpoint of the line segment AB. All of these points are close
      // enough to unit length to satisfy S2::IsUnitlength.
      // Until we get ExactSign, this will only return Indeterminate.
      
      // It should be Clockwise.
      (poA, poB, poC, .indeterminate),
      // The points "x1" and "x2" are exactly proportional, i.e. they both lie
      // on a common line through the origin. Both points are considered to be
      // normalized, and in fact they both satisfy (x == x.normalized()).
      // Therefore the triangle (x1, x2, -x1) consists of three distinct points
      // that all lie on a common line through the origin.
      // Until we get ExactSign, this will only return Indeterminate.
      // It should be CounterClockwise.
      (x1, x2, x1.inverse(), .indeterminate),
      // Here are two more points that are distinct, exactly proportional, and
      // that satisfy (x == x.normalized()).
      // Until we get ExactSign, this will only return Indeterminate.
      // It should be Clockwise.
      (x3, x4, x3.inverse(), .indeterminate),
      // The following points demonstrate that normalized() is not idempotent,
      // i.e. y0.normalized() != y0.normalized().normalized(). Both points satisfy
      // S2::IsNormalized(), though, and the two points are exactly proportional.
      // Until we get ExactSign, this will only return Indeterminate.
      // It should be CounterClockwise.
      (y1, y2, y1.inverse(), .indeterminate)]
    for (p1, p2, p3, want) in tests {
      let result = S2Point.robustSign(p1, p2, p3)
      XCTAssertEqual(result, want)
      // Test RobustSign(b,c,a) == RobustSign(a,b,c) for all a,b,c
      let rotated = S2Point.robustSign(p2, p3, p1)
      XCTAssertEqual(rotated, result)
      // Test RobustSign(c,b,a) == -RobustSign(a,b,c) for all a,b,c
      var want = Direction.clockwise
      if result == .clockwise {
        want = .counterClockwise
      } else if result == .indeterminate {
        want = .indeterminate
      }
      let reversed = S2Point.robustSign(p3, p2, p1)
      XCTAssertEqual(reversed, want)
    }
  }
  
  func testPointDistance() {
    let tests = [
      (1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
      (1, 0, 0, 0, 1, 0, .pi / 2),
      (1, 0, 0, 0, 1, 1, .pi / 2),
      (1, 0, 0, -1, 0, 0, .pi),
      (1, 2, 3, 2, 3, -1, 1.2055891055045298)]
    for (x1, y1, z1, x2, y2, z2, want) in tests {
      let p1 = S2Point(x: x1, y: y1, z: z1)
      let p2 = S2Point(x: x2, y: y2, z: z2)
      XCTAssertEqual(p1.distance(p2), want, accuracy: 1e-15)
      XCTAssertEqual(p2.distance(p1), want, accuracy: 1e-15)
    }
  }

  func testApproxEqual() {
    let epsilon = 1e-14
    let tests: [(Double, Double, Double, Double, Double, Double, Bool)] = [
      (1.0, 0.0, 0.0, 1.0, 0.0, 0.0, true),
      (1, 0, 0, 0, 1, 0, false),
      (1, 0, 0, 0, 1, 1, false),
      (1, 0, 0, -1, 0, 0, false),
      (1, 2, 3, 2, 3, -1, false),
      (1, 0, 0, 1 * (1 + epsilon), 0, 0, true),
      (1, 0, 0, 1 * (1 - epsilon), 0, 0, true),
      (1, 0, 0, 1 + epsilon, 0, 0, true),
      (1, 0, 0, 1 - epsilon, 0, 0, true),
      (1, 0, 0, 1, epsilon, 0, true),
      (1, 0, 0, 1, epsilon, epsilon, false),
      (1, epsilon, 0, 1, -epsilon, epsilon, false)]
    for (x1, y1, z1, x2, y2, z2, want) in tests {
      let p1 = S2Point(x: x1, y: y1, z: z1)
      let p2 = S2Point(x: x2, y: y2, z: z2)
      XCTAssertEqual(p1.approxEquals(p2), want)
    }
  }

  static func pp(_ x: Double, _ y: Double, _ z: Double) -> S2Point {
    return S2Point(x: x, y: y, z: z)
  }
  
  let pz = S2PointTests.pp(0, 0, 1)
  let p000 = S2PointTests.pp(1, 0, 0)
  let p045 = S2PointTests.pp(1, 1, 0)
  let p090 = S2PointTests.pp(0, 1, 0)
  let p180 = S2PointTests.pp(-1, 0, 0)
  // Degenerate triangles.
  let pr = S2PointTests.pp(0.257, -0.5723, 0.112)
  let pq = S2PointTests.pp(-0.747, 0.401, 0.2235)
  // For testing the Girard area fall through case.
  let g1 = S2PointTests.pp(1, 1, 1)
  let g2 = S2PointTests.pp(1, 1, 1).v.add(S2PointTests.pp(0.257, -0.5723, 0.112).v.mul(1e-15)).s2
  let g3 = S2PointTests.pp(1, 1, 1).v.add(S2PointTests.pp(-0.747, 0.401, 0.2235).v.mul(1e-15)).s2

  func testPointArea() {
    let epsilon = 1e-10
    let tests = [
      (p000, p090, pz, .pi / 2.0, 0),
      // This test case should give 0 as the epsilon, but either Go or C++'s value for Pi,
      // or the accuracy of the multiplications along the way, cause a difference ~15 decimal
      // places into the result, so it is not quite a difference of 0.
      (p045, pz, p180, 3.0 * .pi / 4.0, 1e-14),
      // Make sure that Area has good *relative* accuracy even for very small areas.
      (p(epsilon, 0, 1), p(0, epsilon, 1), pz, 0.5 * epsilon * epsilon, 1e-14),
      // Make sure that it can handle degenerate triangles.
      (pr, pr, pr, 0.0, 0),
      (pr, pq, pr, 0.0, 1e-15),
      (p000, p045, p090, 0.0, 0),
      // Try a very long and skinny triangle.
      (p000, p(1, 1, epsilon), p090, 5.8578643762690495119753e-11, 1e-9),
      // TODO:
      // C++ includes a 10,000 loop of perterbations to test out the Girard area
      // computation is less than some noise threshold.
      // Do we need that many? Will one or two suffice?
      (g1, g2, g3, 0.0, 1e-15)]
    for (a, b, c, want, nearness) in tests {
      XCTAssertEqual(S2Point.pointArea(a, b, c), want, accuracy: nearness)
    }
  }

  func testPointAreaQuarterHemisphere() {
    let tests = [
      // Triangles with near-180 degree edges that sum to a quarter-sphere.
      (p(1, 0.1*epsilon, epsilon), p000, p045, p180, pz, Double.pi),
      // Four other triangles that sum to a quarter-sphere.
      (p(1, 1, epsilon), p000, p045, p180, pz, Double.pi)]
      // TODO:
      // C++ Includes a loop of 100 perturbations on a hemisphere for more tests.
    for (a, b, c, d, e, want) in tests {
      let area =
        S2Point.pointArea(a, b, c) +
        S2Point.pointArea(a, c, d) +
        S2Point.pointArea(a, d, e) +
        S2Point.pointArea(a, e, b)
      XCTAssertEqual(area, want, accuracy: 1e-15)
    }
  }

  func testPlanarCentroid() {
    let tests = [
      (name: "xyz axis", p0: p(0, 0, 1), p1: p(0, 1, 0), p2: p(1, 0, 0), want: p(1.0/3, 1.0/3, 1.0/3)),
      (name: "Same point", p0: p(1, 0, 0), p1: p(1, 0, 0), p2: p(1, 0, 0), want: p(1, 0, 0))]
    for (_, p0, p1, p2, want) in tests {
      XCTAssertTrue(S2Point.planarCentroid(p0, p1, p2).approxEquals(want))
    }
  }

  func testTrueCentroid() {
    // Test TrueCentroid with very small triangles. This test assumes that
    // the triangle is small enough so that it is nearly planar.
    // The centroid of a planar triangle is at the intersection of its
    // medians, which is two-thirds of the way along each median.
    for _ in 0..<100 {
      let f = randomFrame()
      let p = f.col(0)
      let x = f.col(1)
      let y = f.col(2)
      let d = 1e-4 * pow(1e-4, randomFloat64())
      // Make a triangle with two equal sides.
      let p0 = p.v.sub(x.v.mul(d)).s2
      let p1 = p.v.add(x.v.mul(d)).s2
      let p2 = p.v.add(y.v.mul(d * 3)).s2
      let want = p.v.add(y.v.mul(d)).s2
      //
      XCTAssert(S2Point.trueCentroid(p0, p1, p2).s2.distance(want) < 2e-8)
      // Make a triangle with a right angle.
      let p0_ = p
      let p1_ = p.v.add(x.v.mul(d * 3)).s2
      let p2_ = p.v.add(y.v.mul(d * 6)).s2
      let want_ = p.v.add(x.v.add(y.v.mul(2)).mul(d)).s2
      XCTAssert(S2Point.trueCentroid(p0_, p1_, p2_).s2.distance(want_) < 2e-8)
    }
  }

  let N = 100
  
  func benchmarkPointArea() {
    for _ in 0..<N {
      let _ = S2Point.pointArea(p000, p090, pz)
    }
  }
  
  func benchmarkPointAreaGirardCase() {
    for _ in 0..<N {
      let _ = S2Point.pointArea(g1, g2, g3)
    }
  }
  
  func benchmarkSign() {
    let p1 = p(-3, -1, 4)
    let p2 = p(2, -1, -3)
    let p3 = p(1, -2, 0)
    for _ in 0..<N {
      let _ = S2Point.sign(p1, b: p2, c: p3)
    }
  }
  
  // BenchmarkRobustSignSimple runs the benchmark for points that satisfy the first
  // checks in RobustSign to compare the performance to that of Sign().
  func benchmarkRobustSignSimple() {
    let p1 = p(-3, -1, 4)
    let p2 = p(2, -1, -3)
    let p3 = p(1, -2, 0)
    for _ in 0..<N {
      let _ = S2Point.robustSign(p1, p2, p3)
    }
  }
  
  // BenchmarkRobustSignNearCollinear runs the benchmark for points that are almost but not
  // quite collinear, so the tests have to use most of the calculations of RobustSign
  // before getting to an answer.
  func benchmarkRobustSignNearCollinear() {
    for _ in 0..<N {
      let _ = S2Point.robustSign(poA, poB, poC)
    }
  }

  // MARK: ST UV tests
  
  func testSTUV() {
    XCTAssertEqual(S2Cube.stToUV(S2Cube.uvToST(0.125)), 0.125)
    XCTAssertEqual(S2Cube.uvToST(S2Cube.stToUV(0.125)), 0.125)
  }
  
  func testUVNorms() {
    for face in 0..<6 {
      for x in stride(from: -1.0, through: 1.0, by: 1.0 / 1024.0) {
        let x1 = S2Cube(face: face, u: x, v: -1).vector()
        let x2 = S2Cube(face: face, u: x, v: 1).vector()
        let uNorm1 = S2Cube.uNorm(face: face, u: x, invert: false)
        XCTAssertEqual(x1.cross(x2).angle(uNorm1.v), 0.0, accuracy: 1e-15)
        let y1 = S2Cube(face: face, u: -1, v: x).vector()
        let y2 = S2Cube(face: face, u: 1, v: x).vector()
        let vNorm1 = S2Cube.vNorm(face: face, v: x, invert: false)
        XCTAssertEqual(y1.cross(y2).angle(vNorm1.v), 0.0, accuracy: 1e-15)
      }
    }
  }
  
  func testFaceXYZToUV() {
    let point = p(1.1, 1.2, 1.3)
    let pointNeg = p(-1.1, -1.2, -1.3)
    let tests: [(Int, S2Point, Double, Double, Bool)] = [
      (0, point, 1.0 + (1.0 / 11), 1 + (2.0 / 11), true),
      (0, pointNeg, 0, 0, false),
      (1, point, -11.0 / 12, 1 + (1.0 / 12), true),
      (1, pointNeg, 0, 0, false),
      (2, point, -11.0 / 13, -12.0 / 13, true),
      (2, pointNeg, 0, 0, false),
      (3, point, 0, 0, false),
      (3, pointNeg, 1 + (2.0 / 11), 1 + (1.0 / 11), true),
      (4, point, 0, 0, false),
      (4, pointNeg, 1 + (1.0 / 12), -(11.0 / 12), true),
      (5, point, 0, 0, false),
      (5, pointNeg, -12.0 / 13, -11.0 / 13, true)]
    for (face, point, u, v, ok) in tests {
      guard let cube = S2Cube(point: point, face: face)
        else {
        if ok {
          XCTFail()
        }
        continue
      }
      XCTAssertEqual(cube.u, u, accuracy: 1e-15)
      XCTAssertEqual(cube.v, v, accuracy: 1e-15)
    }
  }
  
  func testFaceXYZtoUVW() {
//    let origin = p(0, 0, 0)
    let posX = p(1, 0, 0)
    let negX = p(-1, 0, 0)
    let posY = p(0, 1, 0)
    let negY = p(0, -1, 0)
    let posZ = p(0, 0, 1)
    let negZ = p(0, 0, -1)
    for face in 0..<6 {
//      XCTAssertEqual(S2Cube.faceXYZtoUVW(face, p: origin), origin)
      XCTAssertEqual(S2Cube.faceXYZtoUVW(face: face, point: S2Cube.uAxis(face: face)), posX)
      XCTAssertEqual(S2Cube.faceXYZtoUVW(face: face, point: S2Cube.uAxis(face: face).inverse()), negX)
      XCTAssertEqual(S2Cube.faceXYZtoUVW(face: face, point: S2Cube.vAxis(face: face)), posY)
      XCTAssertEqual(S2Cube.faceXYZtoUVW(face: face, point: S2Cube.vAxis(face: face).inverse()), negY)
      XCTAssertEqual(S2Cube.faceXYZtoUVW(face: face, point: S2Cube.unitNorm(face: face)), posZ)
      XCTAssertEqual(S2Cube.faceXYZtoUVW(face: face, point: S2Cube.unitNorm(face: face).inverse()), negZ)
    }
  }
  
  func testUVWAxis() {
    for face in 0..<6 {
      // Check that the axes are consistent with faceUVtoXYZ.
      XCTAssertEqual(S2Cube(face: face, u: 1, v: 0).vector().sub(S2Cube(face: face, u: 0, v: 0).vector()), S2Cube.uAxis(face: face).v)
      XCTAssertEqual(S2Cube(face: face, u: 0, v: 1).vector().sub(S2Cube(face: face, u: 0, v: 0).vector()), S2Cube.vAxis(face: face).v)
      XCTAssertEqual(S2Cube(face: face, u: 0, v: 0).vector(), S2Cube.unitNorm(face: face).v)
      // Check that every face coordinate frame is right-handed.
      XCTAssertEqual(S2Cube.uAxis(face: face).v.cross(S2Cube.vAxis(face: face).v).dot(S2Cube.unitNorm(face: face).v), 1)
      // Check that GetUVWAxis is consistent with GetUAxis, GetVAxis, GetNorm.
      XCTAssertEqual(S2Cube.uAxis(face: face), S2Cube.uvwAxis(face: face, axis: 0))
      XCTAssertEqual(S2Cube.vAxis(face: face), S2Cube.uvwAxis(face: face, axis: 1))
      XCTAssertEqual(S2Cube.unitNorm(face: face), S2Cube.uvwAxis(face: face, axis: 2))
    }
  }
  
}
