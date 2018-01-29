//
//  S2CapTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2CapTests: XCTestCase {
  
  let tinyRad = 1e-10
  
  let empty = S2Cap.empty
  let full = S2Cap.full
  let defaultCap = S2Cap.empty
  
  let xAxisPt = S2Point(x: 1, y: 0, z: 0)
  let yAxisPt = S2Point(x: 0, y: 1, z: 0)
  
  let xAxis = S2Cap(point: S2Point(x: 1, y: 0, z: 0))
  let yAxis = S2Cap(point: S2Point(x: 0, y: 1, z: 0))
  let xComp = S2Cap(point: S2Point(x: 1, y: 0, z: 0)).complement()
  
  let hemi = S2Cap(center: S2Point(x: 1, y: 0, z: 1), height: 1.0)
  let concave = S2Cap(center: LatLng(lat: 80.0 * toRadians, lng: 10.0 * toRadians).toPoint(), angle: 150.0 * toRadians)
  let tiny = S2Cap(center: S2Point(x: 1, y: 2, z: 3), angle: 1e-10)

  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func testCapBasicEmptyFullValid() {
    let tests = [
//      (S2Cap(), false, false, false),
      (empty, true, false, true),
      (empty.complement(), false, true, true),
      (full, false, true, true),
      (full.complement(), true, false, true),
      (defaultCap, true, false, true),
      (xComp, false, true, true),
      (xComp.complement(), true, false, true),
      (tiny, false, false, true),
      (concave, false, false, true),
      (hemi, false, false, true),
      (tiny, false, false, true)]
    for (got, empty, full, valid) in tests {
      XCTAssertEqual(got.isEmpty, empty, "\(got), \(empty)")
      XCTAssertEqual(got.isFull, full)
      XCTAssertEqual(got.isValid, valid, "\(got), \(valid)")
    }
  }
  
  func testCapCenterHeightRadius() {
    XCTAssertTrue(xAxis.approxEquals(xAxis.complement().complement()))
    XCTAssertEqual(full.height, S2Cap.fullHeight)
    XCTAssertEqual(full.radius(), 180.0 * toRadians)
    XCTAssertEqual(empty.center, defaultCap.center)
    XCTAssertEqual(empty.height, defaultCap.height)
    XCTAssertEqual(yAxis.height, S2Cap.zeroHeight)
    XCTAssertEqual(xAxis.height, S2Cap.zeroHeight)
    XCTAssertEqual(xAxis.radius(), S2Cap.zeroHeight)
    XCTAssertEqual(hemi.center.inverse(), hemi.complement().center)
    XCTAssertEqual(hemi.height, 1.0)
  }
  
  func testCapContains() {
    XCTAssertTrue(empty.contains(empty))
    XCTAssertTrue(full.contains(empty))
    XCTAssertTrue(full.contains(full))
    XCTAssertFalse(empty.contains(xAxis))
    XCTAssertTrue(full.contains(xAxis))
    XCTAssertFalse(xAxis.contains(full))
    XCTAssertTrue(xAxis.contains(xAxis))
    XCTAssertTrue(xAxis.contains(empty))
    XCTAssertTrue(hemi.contains(tiny))
    XCTAssertTrue(hemi.contains(S2Cap(center: xAxisPt, angle: .pi / 4 - S2Cap.epsilon)))
    XCTAssertFalse(hemi.contains(S2Cap(center: xAxisPt, angle: .pi / 4 + S2Cap.epsilon)))
    XCTAssertTrue(concave.contains(hemi))
    XCTAssertFalse(concave.contains(S2Cap(center: concave.center.inverse(), height: 0.1)))
  }
  
  func testCapContainsPoint() {
    // We don't use the standard epsilon in this test due different compiler
    // math optimizations that are permissible (FMA vs no FMA) that yield
    // slightly different floating point results between gccgo and gc.
    let tangent = tiny.center.v.cross(S2Point(x: 3, y: 2, z: 1).v)
    XCTAssertTrue(xAxis.contains(xAxisPt))
    XCTAssertFalse(xAxis.contains(S2Point(x: 1, y: 1e-20, z: 0)))
    XCTAssertFalse(yAxis.contains(xAxis.center))
    XCTAssertTrue(xComp.contains(xAxis.center))
    XCTAssertFalse(xComp.complement().contains(xAxis.center))
    XCTAssertTrue(tiny.contains(tiny.center.v.add(tangent.mul(tinyRad * 0.99)).s2))
    XCTAssertFalse(tiny.contains(tiny.center.v.add(tangent.mul(tinyRad * 1.81)).s2))
    XCTAssertTrue(hemi.contains(S2Point(x: 1, y: 0, z: -(1 - epsilon))))
    XCTAssertTrue(hemi.contains(xAxisPt))
    XCTAssertFalse(hemi.complement().contains(xAxisPt))
    XCTAssertTrue(concave.contains(llDegrees(-70 * (1 - S2Cap.epsilon), 10).toPoint()))
    XCTAssertFalse(concave.contains(llDegrees(-70 * (1 + S2Cap.epsilon), 10).toPoint()))
    // This test case is the one where the floating point values end up
    // different in the 15th place and beyond.
    XCTAssertTrue(concave.contains(llDegrees(-50 * (1 - S2Cap.epsilon), -170).toPoint()))
    XCTAssertFalse(concave.contains(llDegrees(-50 * (1 + S2Cap.epsilon), -170).toPoint()))
  }
  
  func testCapInteriorIntersects() {
    let tests = [
      (empty, empty, false),
      (empty, xAxis, false),
      (full, empty, false),
      (full, full, true),
      (full, xAxis, true),
      (xAxis, full, false),
      (xAxis, xAxis, false),
      (xAxis, empty, false),
      (concave, hemi.complement(), true)]
    for (c1, c2, want) in tests {
      XCTAssertEqual(c1.interiorIntersects(c2), want)
    }
  }
  
  func testCapInteriorContains() {
    XCTAssertFalse(hemi.interiorContains(S2Point(x: 1, y: 0, z: -(1 + S2Cap.epsilon).normalized())))
  }
  
  func testCapExpanded() {
    let cap50 = S2Cap(center: xAxisPt, angle: 50.0 * toRadians)
    let cap51 = S2Cap(center: xAxisPt, angle: 51.0 * toRadians)
    XCTAssertTrue(empty.expanded(S2Cap.fullHeight).isEmpty)
    XCTAssertTrue(full.expanded(S2Cap.fullHeight).isFull)
    XCTAssertTrue(cap50.expanded(0).approxEquals(cap50))
    XCTAssertTrue(cap50.expanded(1 * toRadians).approxEquals(cap51))
    XCTAssertFalse(cap50.expanded(129.99 * toRadians).isFull)
    XCTAssertTrue(cap50.expanded(130.01 * toRadians).isFull)
  }
  
  func testRadiusToHeight() {
    let tests = [
      // Above/below boundary checks.
      (-0.5, S2Cap.emptyHeight),
      (0, 0),
      (.pi, S2Cap.fullHeight),
      (2 * .pi, S2Cap.fullHeight),
      // Degree tests.
      (-7.0 * toRadians, S2Cap.emptyHeight),
      (-0.0 * toRadians, 0),
      (0.0 * toRadians, 0),
      (12.0 * toRadians, 0.0218523992661943),
      (30.0 * toRadians, 0.1339745962155613),
      (45.0 * toRadians, 0.2928932188134525),
      (90.0 * toRadians, 1.0),
      (179.99 * toRadians, 1.9999999847691292),
      (180.0 * toRadians, S2Cap.fullHeight),
      (270.0 * toRadians, S2Cap.fullHeight),
      // Radians tests.
      (-1.0, S2Cap.emptyHeight),
      (-0.0, 0),
      (0.0, 0),
      (1.0, 0.45969769413186),
      (.pi / 2.0, 1.0),
      (2.0, 1.4161468365471424),
      (3.0, 1.9899924966004454),
      (.pi, S2Cap.fullHeight),
      (4.0, S2Cap.fullHeight)]
    for (got, want) in tests {
      XCTAssertEqual(S2Cap.radiusToHeight(got), want, accuracy: 1e-14)
    }
  }
  
  func testCapGetRectBounds() {
    let epsilon = 1e-13
    let tests = [
      ("Cap that includes South Pole.", S2Cap(center: llDegrees(-45, 57).toPoint(), angle: 50 * toRadians), -90.0, 5.0, -180.0, 180.0, true),
      ("Cap that is tangent to the North Pole.", S2Cap(center: S2Point(x: 1, y: 0, z: 1), angle: .pi/4.0+1e-16), 0, 90, -180, 180, true),
      ("Cap that at 45 degree center that goes from equator to the pole.", S2Cap(center: S2Point(x: 1, y: 0, z: 1), angle: (45+5e-15) * toRadians), 0, 90, -180, 180, true),
      ("The eastern hemisphere.", S2Cap(center: S2Point(x: 0, y: 1, z: 0), angle: .pi/2+2e-16), -90, 90, -180, 180, true),
      ("A cap centered on the equator.", S2Cap(center: llDegrees(0, 50).toPoint(), angle: 20 * toRadians), -20, 20, 30, 70, false),
      ("A cap centered on the North Pole.", S2Cap(center: llDegrees(90, 123).toPoint(), angle: 10 * toRadians), 80, 90, -180, 180, true)]
  
    for (_, have, latLoDeg, latHiDeg, lngLoDeg, lngHiDeg, isFull) in tests {
      let r = have.rectBound()
      XCTAssertEqual(r.lat.lo, latLoDeg * toRadians, accuracy: epsilon)
      XCTAssertEqual(r.lat.hi, latHiDeg * toRadians, accuracy: epsilon)
      XCTAssertEqual(r.lng.lo, lngLoDeg * toRadians, accuracy: epsilon)
      XCTAssertEqual(r.lng.hi, lngHiDeg * toRadians, accuracy: epsilon)
      XCTAssertEqual(r.lng.isFull, isFull)
    }
    // Empty and full caps.
    XCTAssertTrue(S2Cap.empty.rectBound().isEmpty)
    
    XCTAssertTrue(S2Cap.full.rectBound().isFull)
  }
  
  func testCapAddPoint() {
    let tests = [
      // S2Cap plus its center equals itself.
      (xAxis, xAxisPt, xAxis),
      (yAxis, yAxisPt, yAxis),
      // Cap plus opposite point equals full.
      (xAxis, S2Point(x: -1, y: 0, z: 0), full),
      (yAxis, S2Point(x: 0, y: -1, z: 0), full),
      // Cap plus orthogonal axis equals half cap.
      (xAxis, S2Point(x: 0, y: 0, z: 1), S2Cap(center: xAxisPt, angle: .pi / 2.0)),
      (xAxis, S2Point(x: 0, y: 0, z: -1), S2Cap(center: xAxisPt, angle: .pi / 2.0)),
      // The 45 degree angled hemisphere plus some points.
      (hemi, S2Point(x: 0, y: 1, z: -1), S2Cap(center: S2Point(x: 1, y: 0, z: 1), angle: 120.0 * toRadians)),
      (hemi, S2Point(x: 0, y: -1, z: -1), S2Cap(center: S2Point(x: 1, y: 0, z: 1), angle: 120.0 * toRadians)),
      // This angle between this point and the center is acos(-sqrt(2/3))
      (hemi, S2Point(x: -1, y: -1, z: -1), S2Cap(center: S2Point(x: 1, y: 0, z: 1), angle: 2.5261129449194)),
      (hemi, S2Point(x: 0, y: 1, z: 1), hemi),
      (hemi, S2Point(x: 1, y: 0, z: 0), hemi)]
    for (have, p, want) in tests {
      let got = have.add(p)
      // NSLog("\n\(have) \(p)\n\(got)\n\(want)")
      XCTAssertTrue(got.approxEquals(want))
      XCTAssertTrue(got.contains(p))
    }
  }
  
  func testCapAddCap() {
    let tests = [
      // Identity cases.
      (empty, empty, empty),
      (full, full, full),
      // Anything plus empty equals itself.
      (full, empty, full),
      (empty, full, full),
      (xAxis, empty, xAxis),
      (empty, xAxis, xAxis),
      (yAxis, empty, yAxis),
      (empty, yAxis, yAxis),
      // Two halves make a whole.
      (xAxis, xComp, full),
      // Two zero-height orthogonal axis caps make a half-cap.
      (xAxis, yAxis, S2Cap(center: xAxisPt, angle: .pi/2.0))]
    for (have, other, want) in tests {
      let got = have.add(other)
      XCTAssertTrue(got.approxEquals(want))
    }
  }

  func testCapContainsCell() {
    // call this to initialize the lookup tables
    CellId.setup()
    //
    let eps = 1e-15
    let faceRadius = atan(sqrt(2.0))
    for face in 0..<6 {
      // The cell consisting of the entire face.
      let rootCell = Cell(id: CellId(face: face))
      // A leaf cell at the midpoint of the v=1 edge.
      let edgeCell = Cell(point: S2Cube(face: face, u: 0, v: 1-eps).vector().s2)
      // A leaf cell at the u=1, v=1 corner
      let cornerCell = Cell(point: S2Cube(face: face, u: 1-eps, v: 1-eps).vector().s2)
      // Quick check for full and empty caps.
      XCTAssertTrue(full.contains(rootCell), "face \(face)")
      // Check intersections with the bounding caps of the leaf cells that are adjacent to
      // cornerCell along the Hilbert curve.  Because this corner is at (u=1,v=1), the curve
      // stays locally within the same cube face.
      let first = cornerCell.id.advance(-3)
      let last = cornerCell.id.advance(4)
      var id = first
      while id.id < last.id {
        let c = Cell(id: id).capBound()
        XCTAssertEqual(c.contains(cornerCell), id == cornerCell.id)
        id = id.next()
      }
      //
      for capFace in 0..<6 {
        // A cap that barely contains all of capFace.
        let center = S2Cube.unitNorm(face: capFace)
        let covering = S2Cap(center: center, angle: faceRadius+eps)
        XCTAssertEqual(covering.contains(rootCell), capFace == face)
        XCTAssertEqual(covering.contains(edgeCell), center.v.dot(edgeCell.id.point().v) > 0.1)
        XCTAssertEqual(covering.contains(edgeCell), covering.intersects(edgeCell))
        XCTAssertEqual(covering.contains(cornerCell), capFace == face)
        // A cap that barely intersects the edges of capFace.
        let bulging = S2Cap(center: center, angle: .pi/4+eps)
        XCTAssertFalse(bulging.contains(rootCell))
        XCTAssertEqual(bulging.contains(edgeCell), capFace == face)
        XCTAssertFalse(bulging.contains(cornerCell))
      }
    }
  }
  
  func testCapIntersectsCell() {
    // call this to initialize the lookup tables
    CellId.setup()
    //
    let eps = 1e-15
    let faceRadius = atan(sqrt(2.0))
    for face in 0..<6 {
      // The cell consisting of the entire face.
      let rootCell = Cell(id: CellId(face: face))
      // A leaf cell at the midpoint of the v=1 edge.
      let edgeCell = Cell(point: S2Cube(face: face, u: 0, v: 1-eps).vector().s2)
      // A leaf cell at the u=1, v=1 corner
      let cornerCell = Cell(point: S2Cube(face: face, u: 1-eps, v: 1-eps).vector().s2)
      // Quick check for full and empty caps.
      XCTAssertFalse(empty.intersects(rootCell))
      // Check intersections with the bounding caps of the leaf cells that are adjacent to
      // cornerCell along the Hilbert curve.  Because this corner is at (u=1,v=1), the curve
      // stays locally within the same cube face.
      let first = cornerCell.id.advance(-3)
      let last = cornerCell.id.advance(4)
      var id = first
      while id.id < last.id {
        let c = Cell(id: id).capBound()
        XCTAssertEqual(c.intersects(cornerCell), id.immediateParent().contains(cornerCell.id))
        id = id.next()
      }
      let antiFace = (face + 3) % 6
      for capFace in 0..<6 {
        // A cap that barely contains all of capFace.
        let center = S2Cube.unitNorm(face: capFace)
        let covering = S2Cap(center: center, angle: faceRadius+eps)
        XCTAssertEqual(covering.intersects(rootCell), capFace != antiFace)
        XCTAssertEqual(covering.intersects(edgeCell), covering.contains(edgeCell))
        XCTAssertEqual(covering.intersects(cornerCell), center.v.dot(cornerCell.id.point().v) > 0)
        // A cap that barely intersects the edges of capFace.
        let bulging = S2Cap(center: center, angle: .pi/4+eps)
        XCTAssertEqual(bulging.intersects(rootCell), capFace != antiFace)
        XCTAssertEqual(bulging.intersects(edgeCell), center.v.dot(edgeCell.id.point().v) > 0.1)
        XCTAssertFalse(bulging.intersects(cornerCell))
        // A singleton cap.
        let singleton = S2Cap(center: center, angle: 0)
        XCTAssertEqual(singleton.intersects(rootCell), capFace == face)
        XCTAssertFalse(singleton.intersects(edgeCell))
        XCTAssertFalse(singleton.intersects(cornerCell))
      }
    }
  }

}

