//
//  S2PolylineTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2PolylineTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
  }
  
  override func tearDown() {
    super.tearDown()
  }
  
  func pl(_ coords: [(Double, Double)]) -> S2Polyline {
    let latLngs = coords.map { llDegrees($0, $1) }
    return S2Polyline(latLngs: latLngs)
  }
  
  func testCreate() {
    let coords = [(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)]
    let latLngs = coords.map { llDegrees($0, $1) }
    let poly = S2Polyline(latLngs: latLngs)
    XCTAssertNotNil(poly)
    let poly2 = pl([(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)])
    XCTAssertNotNil(poly2)
  }
  
  func testPolylineBasics() {
    let empty = S2Polyline(points: [])
    XCTAssert(empty.rectBound().isEmpty)
    XCTAssertEqual(empty.rectBound(), S2Rect.empty)
    XCTAssertEqual(empty.numVertices,  0, "empty Polyline should have no vertices")
    let emptyReversed = empty.reversed()
    XCTAssertEqual(emptyReversed.numVertices,  0, "empty Polyline should have no vertices")
    let latlngs = [llDegrees(0, 0), llDegrees(0, 90), llDegrees(0, 180)]
    let semiEquator = S2Polyline(latLngs: latlngs)
//    XCTAssert(semiEquator.interpolate(0.5).approxEquals(p(0, 1, 0)))
    let seReversed = semiEquator.reversed()
    XCTAssert(seReversed.vertex(2).approxEquals(p(1, 0, 0)))
  }

  func testPolylineShape() {
    let shape = makePolyline("0:0, 1:0, 1:1, 2:1")
    XCTAssertEqual(shape.numEdges(), 3)
    XCTAssertEqual(shape.numChains(), 1)
    XCTAssertEqual(shape.chain(0).start, 0)
    XCTAssertEqual(shape.chain(0).length, 3)
    let e = shape.edge(2)
    XCTAssert(e.v0.approxEquals(S2Point(latLng: llDegrees(1, 1))))
    XCTAssert(e.v1.approxEquals(S2Point(latLng: llDegrees(2, 1))))
    XCTAssert(!shape.hasInterior(), "polylines should not have an interior")
    XCTAssert(!shape.referencePoint().contained, "polylines should not contain their reference points")
    XCTAssertEqual(shape.dimension(), .polylineGeometry, "polylines should have PolylineGeometry")
    let empty = S2Polyline(points: [])
    XCTAssertEqual(empty.numEdges(), 0)
    XCTAssertEqual(empty.numChains(), 0)
  }
  
  func testPolylineLengthAndCentroid() {
    // Construct random great circles and divide them randomly into segments.
    // Then make sure that the length and centroid are correct.  Note that
    // because of the way the centroid is computed, it does not matter how
    // we split the great circle into segments.
    for _ in 0..<100 {
      // Choose a coordinate frame for the great circle.
      let f = randomFrame()
      var points: [S2Point] = []
      var theta = 0.0
      while theta < 2 * .pi {
        let p = S2Point(raw: f.row(0).mul(cos(theta)).add(f.row(1).mul(sin(theta))))
        if points.count == 0 || !p.approxEquals(points[points.count - 1]) {
          points.append(p)
        }
        theta += pow(randomFloat64(), 10)
      }
      // close the circle
      points.append(points[0])
      //
      let line = S2Polyline(points: points)
      let length = line.length
      XCTAssertEqual(length, 2 * .pi, accuracy: 2e-14)
      XCTAssert(line.centroid().norm < 2e-14)
    }
  }
  
  func testPolylineIntersectsCell() {
    let pline = S2Polyline(points: [
      S2Point(x: 1, y: -1.1, z: 0.8),
      S2Point(x: 1, y: -0.8, z: 1.1)])
    for face in 0..<6 {
      let cell = Cell(id: CellId(face: face))
      XCTAssertEqual(pline.intersects(cell), face & 1 == 0)
    }
  }
  
  func testPolylineSubsample() {
    let polyStr = "0:0, 0:1, -1:2, 0:3, 0:4, 1:4, 2:4.5, 3:4, 3.5:4, 4:4"
    let tests: [(String, Double, [Int])] = [
      ("", 1.0, []), // No vertices
      ("0:1", 1.0, [0]), // One vertex.
      ("10:10, 11:11", 5.0, [0, 1]), // Two vertices
      // Three points on a straight line. In theory, zero tolerance
      // should work, but in practice there are floating point errors.
      ("-1:0, 0:0, 1:0", 1e-15, [0, 2]),
      ("-1:0, 0:0, 1:1", 0.0, [0, 1, 2]), // Zero tolerance on a non-straight line.
      ("-1:0, 0:0, 1:1", -1.0, [0, 1, 2]), // Negative tolerance should return all vertices.
      ("0:1, 0:2, 0:3, 0:4, 0:5", 1.0, [0, 4]), // Non-zero tolerance with a straight line.
      // And finally, verify that we still do something
      // reasonable if the client passes in an invalid polyline
      // with two or more adjacent vertices.
      ("0:1, 0:1, 0:1, 0:2", 0.0, [0, 3]),
      // Simple examples
      (polyStr, 3.0, [0, 9]),
      (polyStr, 2.0, [0, 6, 9]),
      (polyStr, 0.9, [0, 2, 6, 9]),
      (polyStr, 0.4, [0, 1, 2, 3, 4, 6, 9]),
      (polyStr, 0, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
      ("10:10, 12:12, 10:10", 5.0, [0]), // Check that duplicate vertices are never generated.
      ("0:0, 1:1, 0:0, 0:120, 0:130", 5.0, [0, 3, 4]),
      // Check that points are not collapsed if they would create a line segment
      // longer than 90 degrees, and also that the code handles original polyline
      // segments longer than 90 degrees.
      ("90:0, 50:180, 20:180, -20:180, -50:180, -90:0, 30:0, 90:0", 5.0, [0, 2, 4, 5, 6, 7]),
      // Check that the output polyline is parametrically equivalent and not just
      // geometrically equivalent, i.e. that backtracking is preserved.  The
      // algorithm achieves this by requiring that the points must be encountered
      // in increasing order of distance along each output segment, except for
      // points that are within "tolerance" of the first vertex of each segment.
      ("10:10, 10:20, 10:30, 10:15, 10:40", 5.0, [0, 2, 3, 4]),
      ("10:10, 10:20, 10:30, 10:10, 10:30, 10:40", 5.0, [0, 2, 3, 5]),
      ("10:10, 12:12, 9:9, 10:20, 10:30", 5.0, [0, 4])]
    //
    for test in tests {
      print(test)
      let p = makePolyline(test.0)
      XCTAssertEqual(p.subsampleVertices(tolerance: S1Angle(degrees: test.1)), test.2)
    }
  }

}
