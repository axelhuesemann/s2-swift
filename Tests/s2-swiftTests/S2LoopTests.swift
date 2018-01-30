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

  let full = S2Loop.full
  let empty = S2Loop.empty
  
  // The set of all test loops.
//  let allLoops = [
//    S2Loop.empty, S2Loop.full, northHemi, northHemi3, southHemi, westHemi, eastHemi,
//    nearHemi, farHemi, candyCane, smallNECW, arctic80, antarctic80, lineTriangle,
//    skinnyChevron, loopA, loopB, aIntersectB, aUnionB, aMinusB, bMinusA, loopC, loopD, loopE, loopF, loopG, loopH, loopI]
  // snappedLoopA fails TestAreaConsistentWithTurningAngle (?!?)
  
  func testEmptyFullLoops() {
    let empty = S2Loop.empty
    XCTAssertTrue(empty.isEmpty)
    XCTAssertFalse(empty.isFull)
    XCTAssertTrue(empty.isEmptyOrFull())
    XCTAssertEqual(empty.numEdges(), 0)
    XCTAssertEqual(empty.numChains(), 0)
    let full = S2Loop.full
    XCTAssertFalse(full.isEmpty)
    XCTAssertTrue(full.isFull)
    XCTAssertTrue(full.isEmptyOrFull())
    XCTAssertEqual(full.numEdges(), 0)
    XCTAssertEqual(full.numChains(), 0)
  }

  func r(_ x0: Double, _ x1: Double, _ y0: Double, _ y1: Double) -> S2Rect {
    let lat = R1Interval(lo: x0 * toRadians, hi: y0 * toRadians)
    let lng = S1Interval(lo: x1 * toRadians, hi: y1 * toRadians)
    return S2Rect(lat: lat, lng: lng)
  }

  func testLoopBasic() {
    let shape = makeLoop("0:0, 0:1, 1:0")
    XCTAssertEqual(shape.numEdges(), 3)
    XCTAssertEqual(shape.numChains(), 1)
    XCTAssertEqual(shape.chain(0).start, 0)
    XCTAssertEqual(shape.chain(0).length, 3)
    let e = shape.edge(2)
    XCTAssert(e.v0.approxEquals(p(1, 0)))
    XCTAssert(e.v1.approxEquals(p(0, 0)))
    XCTAssertEqual(shape.dimension(), .polygonGeometry)
    XCTAssert(shape.hasInterior())
    XCTAssert(!shape.referencePoint().contained)
  }
  
  func testLoopHoleAndSign() {
    let l = makeLoop("0:-180, 0:-90, 0:0, 0:90")
    XCTAssert(!l.isHole(), "loop with default depth should not be a hole")
    XCTAssertEqual(l.sign(), 1, "loop with default depth should have a sign of +1")
//    l.depth = 3
    XCTAssert(l.isHole(), "loop with odd depth should be a hole")
    XCTAssertEqual(l.sign(), -1, "loop with odd depth should have a sign of -1")
//    l.depth = 2
    XCTAssert(!l.isHole(), "loop with even depth should not be a hole")
    XCTAssertEqual(l.sign(), 1, "loop with even depth should have a sign of +1")
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
    XCTAssertEqual(arctic80Inv.rectBound().lat.hi, Double(S2Point(raw: mid).lat), accuracy: 10 * Cell.dblEpsilon)
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
    let east = S2Point(x: 0, y: 1, z: 0)
    let west = S2Point(x: 0, y: -1, z: 0)
    XCTAssertFalse(S2Loop.empty.contains(north), "empty loop should not not have any points")
    XCTAssertTrue(S2Loop.full.contains(south), "full loop should have full point vertex")
    let tests: [(name: String, l: S2Loop, in: S2Point, out: S2Point)] = [
      ("north hemisphere", northHemi, north, south),
      ("south hemisphere", southHemi, south, north),
      ("west hemisphere", westHemi, west, east),
      ("east hemisphere", eastHemi, east, west),
      ("candy cane", candyCane, p(5, 71), p(-8, 71))]
    for (_, l, input, output) in tests {
      var l = l
      for _ in 0..<4 {
        XCTAssertTrue(l.contains(input))
        XCTAssertFalse(l.contains(output))
        l = rotate(l: l)
      }
    }
  }
  
  func testLoopContainsAdjacent() {
    // This code checks each cell vertex is contained by exactly one of
    // the adjacent cells.
    for level in 0..<3 {
      // set of unique points across all loops at this level.
      var points: Set<S2Point> = []
      var loops: [S2Loop] = []
      var id = CellId(face: 0).childBegin(level)
      while id != CellId(face: 5).childEnd(level) {
        var vertices: [S2Point] = []
        let cell = Cell(id: id)
        points.insert(cell.center())
        for k in 0..<4 {
          vertices.append(cell.vertex(k))
          points.insert(cell.vertex(k))
        }
        loops.append(S2Loop(points: vertices))
        id = id.next()
      }
      for point in points {
        var count = 0
        for loop in loops {
          if loop.contains(point) {
            count += 1
          }
        }
        XCTAssertEqual(count, 1, "the point should only be contained by one loop at this level")
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
  
  func testLoopFromCell() {
    let cell = Cell(id: CellId(latLng: llDegrees(40.565459, -74.645276)))
    let loopFromCell = S2Loop(cell: cell)
    // Demonstrates the reason for this test; the cell bounds are more
    // conservative than the resulting loop bounds.
    XCTAssert(!loopFromCell.rectBound().contains(cell.rectBound()), "loopFromCell's RectBound countains the original cells RectBound, but should not")
  }
  
  func testRegion() {
    let coords = [(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)]
    let latLngs = coords.map { llDegrees($0, $1) }
    let points = latLngs.map { $0.toPoint() }
    let loop = S2Loop(points: points)
    XCTAssertNotNil(loop)
    // print("coords \(latLngs)")
    // print("rectBound \(loop.rectBound())")
    // print("capBound \(loop.capBound())")
  }

  func testLoopRegularLoop() {
    let loop = S2Loop.regularLoop(center: llDegrees(80, 135).toPoint(), radius: 20, numVertices: 4)
    XCTAssertEqual(loop.vertices.count, 4, "RegularLoop with 4 vertices should have 4 vertices.")
  }
  
  func testLoopContainsMatchesCrossingSign() {
    // This test demonstrates a former incompatibility between CrossingSign
    // and ContainsPoint. It constructs a Cell-based loop L and
    // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
    // previously, Contains() returned false for both endpoints of E.
    //
    // The reason for the bug was that the loop bound was sometimes too tight.
    // The Contains() code for a0 bailed out early because a0 was found not to
    // be inside the bound of L.
    // Start with a cell that ends up producing the problem.
    let cellId = CellId(point: S2Point(x: 1, y: 1, z: 1)).parent(21)
    let children = Cell(id: cellId).children()
    XCTAssertNotNil(children, "error subdividing cell")
    // Note extra normalization. Center() is already normalized.
    // The test results will no longer be inconsistent if the extra
    // Normalize() is removed.
    let points = (0..<4).map { children![$0].center() }
    // Get a vertex from a grandchild cell.
    // +---------------+---------------+
    // |               |               |
    // |    points[3]  |   points[2]   |
    // |       v       |       v       |
    // |       +-------+------ +       |
    // |       |       |       |       |
    // |       |       |       |       |
    // |       |       |       |       |
    // +-------+-------+-------+-------+
    // |       |       |       |       |
    // |       |    <----------------------- grandchild_cell
    // |       |       |       |       |
    // |       +-------+------ +       |
    // |       ^       |       ^       | <-- cell
    // | points[0]/a0  |     points[1] |
    // |               |               |
    // +---------------+---------------+
    let loop = S2Loop(points: points)
    let grandchildren = children![0].children()
    XCTAssertNotNil(grandchildren, "error subdividing cell")
    let grandchildCell = grandchildren![2]
    let a0 = grandchildCell.vertex(0)
    // This test depends on rounding errors that should make a0 slightly different from points[0]
    XCTAssert(points[0] != a0, "not different enough to successfully test")
    // The edge from a0 to the origin crosses one boundary.
    var ec1 = EdgeCrosser(a: a0, b: S2Point.origin, c: loop.vertex(0))
    var ec2 = EdgeCrosser(a: a0, b: S2Point.origin, c: loop.vertex(1))
    var ec3 = EdgeCrosser(a: a0, b: S2Point.origin, c: loop.vertex(2))
    var ec4 = EdgeCrosser(a: a0, b: S2Point.origin, c: loop.vertex(3))
    XCTAssertEqual(ec1.chainCrossingSign(d: loop.vertex(1)), .doNotCross)
    XCTAssertEqual(ec2.chainCrossingSign(d: loop.vertex(2)), .cross)
    XCTAssertEqual(ec3.chainCrossingSign(d: loop.vertex(3)), .doNotCross)
    XCTAssertEqual(ec4.chainCrossingSign(d: loop.vertex(4)), .doNotCross)
    // Contains should return false for the origin, and true for a0.
    XCTAssert(!loop.contains(S2Point.origin))
    XCTAssert(loop.contains(a0))
    // Since a0 is inside the loop, it should be inside the bound.
    let bound = loop.rectBound()
    XCTAssert(bound.contains(a0))
  }

  func testLoopRelations() {
//    a, b       *Loop
//    contains   bool // A contains B
//    contained  bool // B contains A
//    disjoint   bool // A and B are disjoint (intersection is empty)
//    covers     bool // (A union B) covers the entire sphere
//    sharedEdge bool // the loops share at least one edge (possibly reversed)
    let tests: [(S2Loop, S2Loop, Bool, Bool, Bool, Bool, Bool)] = [
      // Check full and empty relationships with normal loops and each other.
      (full, full, true, true, false, true, true),
      (full, northHemi, true, false, false,true, false),
      (full, empty, true, false, true, true, false),
      (northHemi, full, false, true, false, true, false),
      (northHemi, empty, true, false, true, false, false),
      (empty, full, false, true, true, true, false),
      (empty, northHemi, false, true, true, false, false),
      (empty, empty, true, true, true, false, false),
      (northHemi, northHemi, true, true, false, false, true),
      (northHemi, southHemi, false, false, true, true, true),
      (northHemi, eastHemi, false, false, false, false, false),
      (northHemi, arctic80, true, false, false, false, false),
      (northHemi, antarctic80, false, false, true, false, false),
      (northHemi, candyCane, false, false, false, false, false),
      // We can't compare northHemi3 vs. northHemi or southHemi because the
      // result depends on the "simulation of simplicity" implementation details.
      (northHemi3, northHemi3, true, true, false, false, true),
      (northHemi3, eastHemi, false, false, false, false, false),
      (northHemi3, arctic80, true, false, false, false, false),
      (northHemi3, antarctic80, false, false, true, false, false),
      (northHemi3, candyCane, false, false, false, false, false),
      (southHemi, northHemi, false, false, true, true, true),
      (southHemi, southHemi, true, true, false, false, true),
      (southHemi, farHemi, false, false, false, false, false),
      (southHemi, arctic80, false, false, true, false, false),
      // xxxx?
      (southHemi, antarctic80, true, false, false, false, false),
      (southHemi, candyCane, false, false, false, false, false),
      (candyCane, northHemi, false, false, false, false, false),
      (candyCane, southHemi, false, false, false, false, false),
      (candyCane, arctic80, false, false, true, false, false),
      (candyCane, antarctic80, false, false, true, false, false),
      (candyCane, candyCane, true, true, false, false, true),
      (nearHemi, westHemi, false, false, false, false, false),
      (smallNECW, southHemi, true, false, false, false, false),
      (smallNECW, westHemi, true, false, false, false, false),
      (smallNECW, northHemi, false, false, false, true, false),
      (smallNECW, eastHemi, false, false, false, true, false),
      (loopA, loopA, true, true, false, false, true),
      (loopA, loopB, false, false, false, false, false),
      (loopA, aIntersectB, true, false, false, false, true),
      (loopA, aUnionB, true, false, false, false, true),
      (loopA, aMinusB, true, false, false, false, true),
      (loopA, bMinusA, false, false, true, false, true),
      (loopB, loopA, false, false, false, false, false),
      (loopB, loopB, true, true, false, false, true),
      (loopB, aIntersectB, true, false, false, false, true),
      (loopB, aUnionB, false, true, false, false, true),
      (loopB, aMinusB, false, false, true, false, true),
      (loopB, bMinusA, true, false, false, false, true),
      (aIntersectB, loopA, false, true, false, false, true),
      (aIntersectB, loopB, false, true, false, false, true),
      (aIntersectB, aIntersectB, true, true, false, false, true),
      (aIntersectB, aUnionB, false, true, false, false, false),
      (aIntersectB, aMinusB, false, false, true, false, true),
      (aIntersectB, bMinusA, false, false, true, false, true),
      (aUnionB, loopA, true, false, false, false, true),
      (aUnionB, loopB, true, false, false, false, true),
      (aUnionB, aIntersectB, true, false, false, false, false),
      (aUnionB, aUnionB, true, true, false, false, true),
      (aUnionB, aMinusB, true, false, false, false, true),
      (aUnionB, bMinusA, true, false, false, false, true),
      (aMinusB, loopA, false, true, false, false, true),
      (aMinusB, loopB, false, false, true, false, true),
      (aMinusB, aIntersectB, false, false, true, false, true),
      (aMinusB, aUnionB, false, true, false, false, true),
      (aMinusB, aMinusB, true, true, false, false, true),
      (aMinusB, bMinusA, false, false, true, false, false),
      (bMinusA, loopA, false, false, true, false, true),
      (bMinusA, loopB, false, true, false, false, true),
      (bMinusA, aIntersectB, false, false, true, false, true),
      (bMinusA, aUnionB, false, true, false, false, true),
      (bMinusA, aMinusB, false, false, true, false, false),
      (bMinusA, bMinusA, true, true, false, false, true),
      // Make sure the relations are correct if the loop crossing happens on
      // two ends of a shared boundary segment.
      // LoopRelationsWhenSameExceptPiecesStickingOutAndIn
      (loopA, loopC, false, false, false, false, true),
      (loopC, loopA, false, false, false, false, true),
      (loopA, loopD, false, true, false, false, true),
      (loopD, loopA, true, false, false, false, true),
      (loopE, loopF, false, false, true, false, true),
      (loopE, loopG, true, false, false, false, true),
      (loopE, loopH, false, false, false, false, true),
      (loopE, loopI, false, false, false, false, false),
      (loopF, loopG, false, false, true, false, true),
      (loopF, loopH, false, false, false, false, true),
      (loopF, loopI, false, false, false, false, false),
      (loopG, loopH, false, true, false, false, true),
      (loopH, loopG, true, false, false, false, true),
      (loopG, loopI, false, false, true, false, true),
      (loopH, loopI, true, false, false, false, true),
    ]
    //
    for test in tests {
      if test.2 { // contains
        testLoopNestedPair(test.0, test.1)
      }
      if test.3 { // contained
        testLoopNestedPair(test.1, test.0)
      }
      if test.4 { // covers
        let b1 = test.1.inverted()
        testLoopNestedPair(test.0, b1)
      }
      if test.5 { // disjoint
        let a1 = test.0.inverted()
        testLoopNestedPair(a1, test.1)
      } else if !(test.2 || test.3 || test.4) { // contains, contained, covered
        // Given loops A and B such that both A and its complement
        // intersect both B and its complement, test various
        // identities involving these four loops.
        let a1 = test.0.inverted()
        let b1 = test.1.inverted()
        testLoopOneOverlappingPair(test.0, test.1)
        testLoopOneOverlappingPair(a1, b1)
        testLoopOneOverlappingPair(a1, test.1)
        testLoopOneOverlappingPair(test.0, b1)
      }
      if !test.6 && (test.2 || test.3 || test.4) {
        XCTAssertEqual(test.0.contains(test.1), test.0.containsNested(test.1))
      }
      // A contains the boundary of B if either A contains B, or the two loops
      // contain each other's boundaries and there are no shared edges (since at
      // least one such edge must be reversed, and therefore is not considered to
      // be contained according to the rules of compareBoundary).
      var comparison = 0
      if test.2 || (test.4 && !test.6) {
        comparison = 1
      }
      // Similarly, A excludes the boundary of B if either A and B are disjoint,
      // or B contains A and there are no shared edges (since A is considered to
      // contain such edges according to the rules of compareBoundary).
      if test.5 || (test.3 && !test.6) {
        comparison = -1
      }
      // compareBoundary requires that neither loop is empty.
      if !test.0.isEmpty && !test.1.isEmpty {
        XCTAssertEqual(test.0.compareBoundary(test.1), comparison)
      }
    }
  }
  
  // Given a pair of loops where A contains B, test various identities
  // involving A, B, and their complements.
  func testLoopNestedPair(_ a: S2Loop, _ b: S2Loop) {
    let a1 = a.inverted()
    let b1 = b.inverted()
    testLoopOneNestedPair(a, b)
    testLoopOneNestedPair(b1, a1)
    testLoopOneDisjointPair(a1, b)
    testLoopOneCoveringPair(a, b1)
  }
  
  // Given a pair of loops where A contains B, check various identities.
  func testLoopOneNestedPair(_ a: S2Loop, _ b: S2Loop) {
    XCTAssert(a.contains(b))
    XCTAssertEqual(b.contains(a), a.boundaryEquals(b))
    XCTAssertEqual(a.intersects(b), !b.isEmpty)
    XCTAssertEqual(b.intersects(a), !b.isEmpty)
  }
  
  // Given a pair of disjoint loops A and B, check various identities.
  func testLoopOneDisjointPair(_ a: S2Loop, _ b: S2Loop) {
    XCTAssert(!a.intersects(b))
    XCTAssert(!b.intersects(a))
    XCTAssertEqual(a.contains(b), b.isEmpty)
    XCTAssertEqual(b.contains(a), a.isEmpty)
  }
 
  // Given loops A and B whose union covers the sphere, check various identities.
  func testLoopOneCoveringPair(_ a: S2Loop, _ b: S2Loop) {
    XCTAssertEqual(a.contains(b), a.isFull)
    XCTAssertEqual(b.contains(a), b.isFull)
    // TODO(roberts): Uncomment as these functions get completed.
    // let a1 = a.inverted()
    // let complementary = a1.boundaryEquals(b)
    // XCTAssertEqual(a.intersects(b), !complementary)
    // XCTAssertEqual(b.intersects(a), !complementary)
  }
  
  // Given loops A and B such that both A and its complement intersect both B
  // and its complement, check various identities.
  func testLoopOneOverlappingPair(_ a: S2Loop, _ b: S2Loop) {
    XCTAssert(!a.contains(b))
    XCTAssert(!b.contains(a))
    XCTAssert(a.intersects(b))
    XCTAssert(b.intersects(a))
  }

  func testLoopTurningAngle() {
    let tests: [(S2Loop, Double)] = [
      (empty, 2 * .pi),
      (full, -2 * .pi),
      (northHemi3, 0),
      (westHemi, 0),
      (candyCane, 4.69364376125922),
      (lineTriangle, 2 * .pi),
      (skinnyChevron, 2 * .pi)]
    for (loop, want) in tests {
      XCTAssertEqual(loop.turningAngle(), want, accuracy: epsilon)
      // Check that the turning angle is *identical* when the vertex order is
      // rotated, and that the sign is inverted when the vertices are reversed.
      let expected = loop.turningAngle()
      for _ in 0..<loop.vertices.count {
        let loop2 = loop.inverted()
        XCTAssertEqual(loop2.turningAngle(), -expected)
        // Invert it back to normal.
        let loop3 = loop2.inverted()
        let loop4 = rotate(l: loop3)
        XCTAssertEqual(loop4.turningAngle(), expected)
      }
    }
    // TODO(roberts): Uncomment once Area is implemented.
    // Build a narrow spiral loop starting at the north pole. This is designed
    // to test that the error in TurningAngle is linear in the number of
    // vertices even when the partial sum of the turning angles gets very large.
    // The spiral consists of two arms defining opposite sides of the loop.
    let armPoints = 10000 // Number of vertices in each "arm"
    let armRadius = 0.01  // Radius of spiral.
    var vertices = (0..<(2 * armPoints)).map { _ in S2Point.origin }
    // Set the center point of the spiral.
    vertices[armPoints] = S2Point(x: 0, y: 0, z: 1)
    for i in 0..<armPoints {
      let angle = (2 * .pi / 3) * Double(i)
      let x = cos(angle)
      let y = sin(angle)
      let r1 = Double(i) * armRadius / Double(armPoints)
      let r2 = (Double(i) + 1.5) * armRadius / Double(armPoints)
      vertices[armPoints - i - 1] = S2Point(x: r1 * x, y: r1 * y, z: 1)
      vertices[armPoints + i] = S2Point(x: r2 * x, y: r2 * y, z: 1)
    }
    // This is a pathological loop that contains many long parallel edges.
    let spiral = S2Loop(points: vertices)
    // Check that TurningAngle is consistent with Area to within the
    // error bound of the former. We actually use a tiny fraction of the
    // worst-case error bound, since the worst case only happens when all the
    // roundoff errors happen in the same direction and this test is not
    // designed to achieve that. The error in Area can be ignored for the
    // purposes of this test since it is generally much smaller.
    XCTAssertEqual(spiral.turningAngle(), (2 * .pi - spiral.area()), accuracy: 0.01 * spiral.turningAngleMaxError())
  }

  func testLoopAreaAndCentroid() {
    let p = S2Point.origin
    XCTAssertEqual(empty.area(), 0.0)
    XCTAssertEqual(full.area(), 4 * .pi)
    XCTAssert(p.approxEquals(empty.centroid()))
    XCTAssert(p.approxEquals(full.centroid()))
    XCTAssertEqual(northHemi.area(), 2 * .pi)
    let eastHemiArea = eastHemi.area()
    XCTAssert(eastHemiArea >= 2 * .pi - 1e-12 && eastHemiArea <= 2 * .pi + 1e-12)
    // Construct spherical caps of random height, and approximate their boundary
    // with closely spaces vertices. Then check that the area and centroid are
    // correct.
    for _ in 0..<50 {
      // Choose a coordinate frame for the spherical cap.
      let f = randomFrame()
      let x = f.col(0)
      let y = f.col(1)
      let z = f.col(2)
      // Given two points at latitude phi and whose longitudes differ by dtheta,
      // the geodesic between the two points has a maximum latitude of
      // atan(tan(phi) / cos(dtheta/2)). This can be derived by positioning
      // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
      // We want to position the vertices close enough together so that their
      // maximum distance from the boundary of the spherical cap is maxDist.
      // Thus we want abs(atan(tan(phi) / cos(dtheta/2)) - phi) <= maxDist.
      let maxDist = 1e-6
      let height = 2 * randomFloat64()
      let phi = asin(1.0 - height)
      let maxDtheta0 = 2 * acos(tan(abs(phi)) / tan(abs(phi) + maxDist))
      let maxDtheta = min(.pi, maxDtheta0)
      var vertices: [S2Point] = []
      var theta = 0.0
      while theta < 2 * .pi {
        let x1 = x.mul(cos(theta) * cos(phi))
        let y1 = y.mul(sin(theta) * cos(phi))
        let z1 = z.mul(sin(phi))
        let p1 = S2Point(raw: x1.add(y1).add(z1))
        vertices.append(p1)
        theta += randomFloat64() * maxDtheta
      }
      let loop = S2Loop(points: vertices)
      let area = loop.area()
      let centroid = loop.centroid()
      let expectedArea = 2 * .pi * height
      XCTAssert(abs(area - expectedArea) > 2 * .pi * maxDist)
      let expectedCentroid = z.mul(expectedArea * (1 - 0.5 * height))
      XCTAssert(centroid.v.sub(expectedCentroid).norm > 2 * maxDist)
    }
  }
  
  // TODO(roberts): Test that Area() has an accuracy significantly better
  // than 1e-15 on loops whose area is small.
  
  func testLoopAreaConsistentWithTurningAngle() {
    // Check that the area computed using GetArea() is consistent with the
    // turning angle of the loop computed using GetTurnAngle().  According to
    // the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
    // minus its turning angle.
    let allLoops = [empty, full, northHemi, northHemi3, southHemi, westHemi, eastHemi, nearHemi, farHemi, candyCane, smallNECW, arctic80, antarctic80, lineTriangle, skinnyChevron, loopA, loopB, aIntersectB, aUnionB, aMinusB, bMinusA, loopC, loopD, loopE, loopF, loopG, loopH, loopI]
    for loop in allLoops {
      let area = loop.area()
      let gaussArea = 2 * .pi - loop.turningAngle()
      // TODO(roberts): The error bound below is much larger than it should be.
      XCTAssertEqual(area, gaussArea, accuracy: 1e-9)
    }
  }
  
  func testLoopGetAreaConsistentWithSign() {
    // TODO(roberts): Uncomment when Loop has IsValid
    // Test that Area() returns an area near 0 for degenerate loops that
    // contain almost no points, and an area near 4*pi for degenerate loops that
    // contain almost all points.
//    let maxVertices = 6
//    for _ in 0..<50 {
//      let numVertices = 3 + randomUniformInt(n: maxVertices - 3 + 1)
//      // Repeatedly choose N vertices that are exactly on the equator until we
//      // find some that form a valid loop.
//      var loop = S2Loop.empty
//      while !loop.isValid {
//        // We limit longitude to the range [0, 90] to ensure that the loop is
//        // degenerate (as opposed to following the entire equator).
//        let vertices = (0..<numVertices).map { _ in llDegrees(0, randomFloat64() * .pi / 2).toPoint() }
//        loop.vertices = vertices
//        break
//      }
//      let ccw = loop.isNormalized()
//      var want = 0.0
//      if !ccw {
//        want = 4 * .pi
//      }
//      // TODO(roberts): The error bound below is much larger than it should be.
//      XCTAssertEqual(loop.area(), want, accuracy: 1e-8)
//      let p1 = p(0, 0, 1)
//      XCTAssertEqual(loop.contains(p1), ccw)
//    }
  }
  
  func testLoopNormalizedCompatibleWithContains() {
    let p = parsePoint("40:40")
    let tests = [lineTriangle, skinnyChevron]
    // Checks that if a loop is normalized, it doesn't contain a
    // point outside of it, and vice versa.
    for loop in tests {
      let flip = loop.inverted()
      XCTAssertEqual(loop.isNormalized(), !loop.contains(p))
      XCTAssertEqual(flip.isNormalized(), !flip.contains(p))
      XCTAssertEqual(loop.isNormalized(), !flip.isNormalized())
      let flip2 = flip.normalized()
      XCTAssert(!flip2.contains(p))
    }
  }
  
  func testLoopIsValidDetectsInvalidLoops() {
    let tests: [(msg: String, points: [S2Point])] = [
      // Not enough vertices. Note that all single-vertex loops are valid; they
      // are interpreted as being either "empty" or "full".
      ("loop has no vertices", parsePoints("")),
      ("loop has too few vertices", parsePoints("20:20, 21:21")),
      // degenerate edge checks happen in validation before duplicate vertices.
      ("loop has degenerate first edge", parsePoints("20:20, 20:20, 20:21")),
      ("loop has degenerate third edge", parsePoints("20:20, 20:21, 20:20")),
      // TODO(roberts): Uncomment these cases when FindAnyCrossings is in.
      ("loop has duplicate points", parsePoints("20:20, 21:21, 21:20, 20:20, 20:21")),
      ("loop has crossing edges", parsePoints("20:20, 21:21, 21:20.5, 21:20, 20:21")),
      // Ensure points are not normalized.
      ("loop with non-normalized vertices", [p(2, 0, 0), p(0, 1, 0), p(0, 0, 1)]),
      // Adjacent antipodal vertices
      ("loop with antipodal points",[ p(1, 0, 0), p(-1, 0, 0), p(0, 0, 1)])]
    for (_, points) in tests {
      let loop = S2Loop(points: points)
      XCTAssert(!loop.isValid)
      // The C++ tests also tests that the returned error message string contains
      // a specific set of text. That part of the test is skipped here.
      //      XCTAssertEqual(loop.findValidationError(), msg)
    } 
  }
  
  // TODO(roberts): Convert these into changeable flags or parameters.
  // A loop with a 10km radius and 4096 vertices has an edge length of 15 meters.
  let defaultRadiusKm = 10.0
  let numLoopSamples = 16
  let numQueriesPerLoop = 100
  
//  func testBenchmarkLoopContainsPoint() {
//    // Benchmark ContainsPoint() on regular loops. The query points for a loop are
//    // chosen so that they all lie in the loop's bounding rectangle (to avoid the
//    // quick-rejection code path).
//    // C++ ranges from 4 -> 256k by powers of 2 for number of vertices for benchmarking.
//    var nVertices = 4
//    for n in 1...17 {
//      var loops: [S2Loop] = [] // , numLoopSamples)
//      for i in 0..<numLoopSamples {
//        loops[i] = S2Loop.regularLoop(center: randomPoint(), radius: kmToAngle(km: 10.0), numVertices: nVertices)
//      }
//      let queries = loops.map { loop in
//        return (0..<numQueriesPerLoop).map { _ in samplePointFromRect(rect: loop.rectBound()) }
//      }
//      for i in 0..<n {
//        let _ = loops[i % numLoopSamples].contains(queries[i % numLoopSamples][i % numQueriesPerLoop])
//      }
//      nVertices *= 2
//    }
//  }

}
