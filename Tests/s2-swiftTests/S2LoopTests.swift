//
//  S2LoopTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


func makeLoop(_ s: String) -> S2Loop {
  let points = parsePoints(s)
  return S2Loop(points: points)
//  return S2Loop.loopFromPoints(points)
}


class S2LoopTests: XCTestCase {

  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }

  // The northern hemisphere, defined using two pairs of antipodal points.
  let northHemi = makeLoop("0:-180, 0:-90, 0:0, 0:90")
  
  // The northern hemisphere, defined using three points 120 degrees apart.
  let northHemi3 = makeLoop("0:-180, 0:-60, 0:60")
  
  // The southern hemisphere, defined using two pairs of antipodal points.
  let southHemi = makeLoop("0:90, 0:0, 0:-90, 0:-180")
  
  // The western hemisphere, defined using two pairs of antipodal points.
  let westHemi = makeLoop("0:-180, -90:0, 0:0, 90:0")
  
  // The eastern hemisphere, defined using two pairs of antipodal points.
  let eastHemi = makeLoop("90:0, 0:0, -90:0, 0:-180")
  
  // The "near" hemisphere, defined using two pairs of antipodal points.
  let nearHemi = makeLoop("0:-90, -90:0, 0:90, 90:0")
  
  // The "far" hemisphere, defined using two pairs of antipodal points.
  let farHemi = makeLoop("90:0, 0:90, -90:0, 0:-90")
  
  // A spiral stripe that slightly over-wraps the equator.
  let candyCane = makeLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")
  
  // A small clockwise loop in the northern & eastern hemisperes.
  let smallNECW = makeLoop("35:20, 45:20, 40:25")
  
  // Loop around the north pole at 80 degrees.
  let arctic80 = makeLoop("80:-150, 80:-30, 80:90")
  
  // Loop around the south pole at 80 degrees.
  let antarctic80 = makeLoop("-80:120, -80:0, -80:-120")
  
  // A completely degenerate triangle along the equator that RobustCCW()
  // considers to be CCW.
  let lineTriangle = makeLoop("0:1, 0:2, 0:3")
  
  // A nearly-degenerate CCW chevron near the equator with very long sides
  // (about 80 degrees).  Its area is less than 1e-640, which is too small
  // to represent in double precision.
  let skinnyChevron = makeLoop("0:0, -1e-320:80, 0:1e-320, 1e-320:80")
  
  // A diamond-shaped loop around the point 0:180.
  let loopA = makeLoop("0:178, -1:180, 0:-179, 1:-180")
  
  // Like loopA, but the vertices are at leaf cell centers.
  let snappedLoopA = S2Loop(points: [
    CellId(latLng: parseLatLngs("0:178")[0]).point(),
    CellId(latLng: parseLatLngs("-1:180")[0]).point(),
    CellId(latLng: parseLatLngs("0:-179")[0]).point(),
    CellId(latLng: parseLatLngs("1:-180")[0]).point()])
  
  // A different diamond-shaped loop around the point 0:180.
  let loopB = makeLoop("0:179, -1:180, 0:-178, 1:-180")
  
  // The intersection of A and B.
  let aIntersectB = makeLoop("0:179, -1:180, 0:-179, 1:-180")
  
  // The union of A and B.
  let aUnionB = makeLoop("0:178, -1:180, 0:-178, 1:-180")
  
  // A minus B (concave).
  let aMinusB = makeLoop("0:178, -1:180, 0:179, 1:-180")
  
  // B minus A (concave).
  let bMinusA = makeLoop("0:-179, -1:180, 0:-178, 1:-180")
  
  // A shape gotten from A by adding a triangle to one edge, and
  // subtracting a triangle from the opposite edge.
  let loopC = makeLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180")
  
  // A shape gotten from A by adding a triangle to one edge, and
  // adding another triangle to the opposite edge.
  let loopD = makeLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180")
  
  //   3------------2
  //   |            |               ^
  //   |  7-8  b-c  |               |
  //   |  | |  | |  |      Latitude |
  //   0--6-9--a-d--1               |
  //   |  | |       |               |
  //   |  f-e       |               +----------->
  //   |            |                 Longitude
  //   4------------5
  //
  // Important: It is not okay to skip over collinear vertices when
  // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
  // uses symbolic perturbations to ensure that no three vertices are
  // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
  // collinear).  In other words, it is unpredictable (modulo knowing the
  // details of the symbolic perturbations) whether 0123 contains 06123
  // for example.
  
  // Loop E: 0,6,9,a,d,1,2,3
  // Loop F: 0,4,5,1,d,a,9,6
  // Loop G: 0,6,7,8,9,a,b,c,d,1,2,3
  // Loop H: 0,6,f,e,9,a,b,c,d,1,2,3
  // Loop I: 7,6,f,e,9,8
  let loopE = makeLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30")
  let loopF = makeLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34")
  let loopG = makeLoop("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
  let loopH = makeLoop("0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
  
  let loopI = makeLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36")

  func p(_ x: Double, _ y: Double, _ z: Double) -> S2Point {
      return S2Point(x: x, y: y, z: z)
  }
  
  func p(_ lat: Double, _ lng: Double) -> S2Point {
    return llDegrees(lat, lng).toPoint()
  }
  
  func testEmptyfulloops() {
    let emptyLoop = S2Loop.empty
    XCTAssertTrue(emptyLoop.isEmpty)
    XCTAssertFalse(emptyLoop.isFull)
    XCTAssertTrue(emptyLoop.isEmptyOrFull())
    let fulloop = S2Loop.full
    XCTAssertFalse(fulloop.isEmpty)
    XCTAssertTrue(fulloop.isFull)
    XCTAssertTrue(fulloop.isEmptyOrFull())
  }
  
  func r(_ x0: Double, _ x1: Double, _ y0: Double, _ y1: Double) -> S2Rect {
    let lat = R1Interval(lo: x0 * toRadians, hi: y0 * toRadians)
    let lng = S1Interval(lo: x1 * toRadians, hi: y1 * toRadians)
    return S2Rect(lat: lat, lng: lng)
  }
  
  func testLoopRectBound() {
    let _ = S2Loop(points: [S2Point(x: 0, y: 0, z: -1)])

    XCTAssertTrue(S2Loop.empty.rectBound().isEmpty)
    XCTAssertTrue(S2Loop.full.rectBound().isFull)
    XCTAssertTrue(candyCane.rectBound().lng.isFull)
    XCTAssert(candyCane.rectBound().lat.lo < -0.349066)
    XCTAssert(candyCane.rectBound().lat.hi > 0.174533)
    XCTAssertTrue(smallNECW.rectBound().isFull)
    XCTAssert(rectsApproxEqual(a: arctic80.rectBound(), b: r(80, -180, 90, 180), tolLat: rectErrorLat, tolLng: rectErrorLng))
    XCTAssert(rectsApproxEqual(a: antarctic80.rectBound(), b: r(-90, -180, -80, 180), tolLat: rectErrorLat, tolLng: rectErrorLng))
    XCTAssertTrue(southHemi.rectBound().lng.isFull)
    XCTAssert(southHemi.rectBound().lat.approxEquals(R1Interval(lo: -.pi/2, hi: 0)))
    
    // Create a loop that contains the complement of the arctic80 loop.
    let arctic80Inv = invert(l: arctic80)
    // The highest latitude of each edge is attained at its midpoint.
    let mid = arctic80Inv.vertices[0].v.add(arctic80Inv.vertices[1].v).mul(0.5)
    XCTAssertEqual(arctic80Inv.rectBound().lat.hi, Double(LatLng(point: S2Point(raw: mid)).lat), accuracy: 10 * Cell.dblEpsilon)
  }
  
  func testLoopCapBound() {
    XCTAssertTrue(S2Loop.empty.capBound().isEmpty)
    XCTAssertTrue(S2Loop.full.capBound().isFull)
    XCTAssertTrue(smallNECW.capBound().isFull)
    XCTAssert(arctic80.capBound().approxEquals(r(80, -180, 90, 180).capBound()))
    XCTAssert(antarctic80.capBound().approxEquals(r(-90, -180, -80, 180).capBound()))
  }
  
  func invert(l: S2Loop) -> S2Loop {
    return S2Loop(points: l.vertices.reversed())
  }
  
  func testOriginInside() {
    XCTAssertTrue(northHemi.originInside)
    XCTAssertTrue(northHemi3.originInside)
    XCTAssertFalse(southHemi.originInside)
    XCTAssertFalse(westHemi.originInside)
    XCTAssertTrue(eastHemi.originInside)
    XCTAssertTrue(nearHemi.originInside)
    XCTAssertFalse(farHemi.originInside)
    XCTAssertFalse(candyCane.originInside)
    XCTAssertTrue(smallNECW.originInside)
    XCTAssertFalse(arctic80.originInside)
    XCTAssertFalse(antarctic80.originInside)
    XCTAssertFalse(loopA.originInside)
  }
  
  func testLoopContainsPoint() {
    let north = S2Point(x: 0, y: 0, z: 1)
    let south = S2Point(x: 0, y: 0, z: -1)
    
    XCTAssertFalse(S2Loop.empty.contains(north))
    XCTAssertTrue(S2Loop.full.contains(south))
    
    let tests = [
      ("north hemisphere", northHemi, p(0, 0, 1), p(0, 0, -1)),
      ("south hemisphere", southHemi, p(0, 0, -1), p(0, 0, 1)),
      ("west hemisphere", westHemi, p(0, -1, 0), p(0, 1, 0)),
      ("east hemisphere", eastHemi, p(0, 1, 0), p(0, -1, 0)),
      ("candy cane", candyCane, p(5, 71), p(-8, 71))]
    for (_, l, inn, out) in tests {
      var l = l
      for _ in 0..<4 {
        XCTAssertTrue(l.contains(inn))
        XCTAssertFalse(l.contains(out))
        l = rotate(l: l)
      }
    }
  }

  func testVertex() {
    let tests = [
      (S2Loop.empty, 0, p(0, 0, 1)),
      (S2Loop.empty, 1, p(0, 0, 1)),
      (S2Loop.full, 0, p(0, 0, -1)),
      (S2Loop.full, 1, p(0, 0, -1)),
      (arctic80, 0, parsePoint("80:-150")),
      (arctic80, 1, parsePoint("80:-30")),
      (arctic80, 2, parsePoint("80:90")),
      (arctic80, 3, parsePoint("80:-150"))]
    for (loop, vertex, want) in tests {
      XCTAssert(pointsApproxEquals(a: loop.vertex(vertex), b: want, epsilon: epsilon))
    }
    // Check that wrapping is correct.
    XCTAssert(pointsApproxEquals(a: arctic80.vertex(2), b: arctic80.vertex(5), epsilon: epsilon))
    let loopAroundThrice = 2 + 3 * arctic80.vertices.count
    XCTAssert(pointsApproxEquals(a: arctic80.vertex(2), b: arctic80.vertex(loopAroundThrice), epsilon: epsilon))
  }

  func testNumEdges() {
    let tests = [
      (S2Loop.empty, 0),
      (S2Loop.full, 0),
      (farHemi, 4),
      (candyCane, 6),
      (smallNECW, 3),
      (arctic80, 3),
      (antarctic80, 3),
      (lineTriangle, 3),
      (skinnyChevron, 4)]
    for (loop, want) in tests {
      XCTAssertEqual(loop.numEdges(), want)
    }
  }

  func testEdge() {
    let tests = [
      (loop: farHemi, edge: 2, wantA: p(0, 0, -1), wantB: p(0, -1, 0)),
      (loop: candyCane, edge: 0, wantA: parsePoint("-20:150"), wantB: parsePoint("-20:-70")),
      (loop: candyCane, edge: 1, wantA: parsePoint("-20:-70"), wantB: parsePoint("0:70")),
      (loop: candyCane, edge: 2, wantA: parsePoint("0:70"), wantB: parsePoint("10:-150")),
      (loop: candyCane, edge: 3, wantA: parsePoint("10:-150"), wantB: parsePoint("10:70")),
      (loop: candyCane, edge: 4, wantA: parsePoint("10:70"), wantB: parsePoint("-10:-70")),
      (loop: candyCane, edge: 5, wantA: parsePoint("-10:-70"), wantB: parsePoint("-20:150")),
      (loop: skinnyChevron, edge: 2, wantA: parsePoint("0:1e-320"), wantB: parsePoint("1e-320:80")),
      (loop: skinnyChevron, edge: 3, wantA: parsePoint("1e-320:80"), wantB: parsePoint("0:0"))]
    for (loop, edge, wantA, wantB) in tests {
      let e = loop.edge(edge)
      XCTAssert(pointsApproxEquals(a: e.v0, b: wantA, epsilon: epsilon))
      XCTAssert(pointsApproxEquals(a: e.v1, b: wantB, epsilon: epsilon))
    }
  }

  func rotate(l: S2Loop) -> S2Loop {
    var vertices = [S2Point]()
    for i in 1..<l.vertices.count {
      vertices.append(l.vertices[i])
    }
    vertices.append(l.vertices[0])
    return S2Loop(points: vertices)
  }
  
  func testRegion() {
    let coords = [(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)]
    let latLngs = coords.map { llDegrees($0, $1) }
    let points = latLngs.map { S2Point(latLng: $0) }
    let loop = S2Loop(points: points)
    XCTAssertNotNil(loop)
    // print("coords \(latLngs)")
    // print("rectBound \(loop.rectBound())")
    // print("capBound \(loop.capBound())")
  }

}
