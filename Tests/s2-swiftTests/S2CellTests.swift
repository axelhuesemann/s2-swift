//
//  S2CellTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2CellTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  
  func testCellObjectSize() {
    // maxCellSize is the upper bounds on the number of bytes we want the Cell object to ever be.
    let maxCellSize = 48
    XCTAssert(MemoryLayout.size(ofValue: Cell.self) <= maxCellSize)
  }
  
  func testCellFaces() {
    // call this to initialize the lookup tables
    CellId.setup()
    // hash counters
    var edgeCounts = [R3Vector: Int]()
    var vertexCounts = [R3Vector: Int]()
    //
    for face in 0..<6 {
      let id = CellId(face: face)
      let cell = Cell(id: id)
      XCTAssertEqual(cell.id, id)
      XCTAssertEqual(cell.face, UInt8(face))
      XCTAssertEqual(cell.level, 0)
      // Top-level faces have alternating orientations to get RHS coordinates.
//      XCTAssertEqual(cell.orientation, UInt8(face & CellId.swapMask))
      let edges = (0..<4).map { cell.edge($0).v }
      let vertices = (0..<4).map { cell.vertex($0).v }
      XCTAssertFalse(cell.isLeaf())
      for k in 0..<4 {
        edgeCounts[edges[k]] = (edgeCounts[edges[k]] ?? 0) + 1
        vertexCounts[vertices[k]] = (vertexCounts[vertices[k]] ?? 0) + 1
        XCTAssertEqual(vertices[k].dot(edges[k]), 0.0, accuracy: 1e-14)
        XCTAssertEqual(vertices[(k + 1) & 3].dot(edges[k]), 0.0, accuracy: 1e-14)
        XCTAssertEqual(vertices[k].cross(vertices[(k + 1) & 3]).normalized().dot(edges[k]), 1.0, accuracy: 1e-14)
      }
    }
    // Check that edges have multiplicity 2 and vertices have multiplicity 3.
    for (_, v) in edgeCounts {
      XCTAssertEqual(v, 2, "\(v)")
    }
    for (_, v) in vertexCounts {
      XCTAssertEqual(v, 3)
    }
  }
  
  func testAreas() {
    // relative error bounds for each type of area computation
    let exactError = log(1 + 1e-6)
    let approxError = log(1.03)
    let avgError = log(1 + 1e-15)
    // Test 1. Check the area of a top level cell.
    let level1Cell = CellId(id: 0x1000000000000000)
    let wantArea = 4.0 * .pi / 6.0
    XCTAssertEqual(Cell(id: level1Cell).exactArea(), wantArea, accuracy: 1e-14)
    // Test 2. Iterate inwards from this cell, checking at every level that
    // the sum of the areas of the children is equal to the area of the parent.
    var childIndex = 1
    var cell = CellId(id: 0x1000000000000000)
    while (cell.level() < 21) {
      var exactArea = 0.0
      var approxArea = 0.0
      var avgArea = 0.0
      for child in cell.children() {
        exactArea += Cell(id: child).exactArea()
        approxArea += Cell(id: child).approxArea()
        avgArea += Cell(id: child).averageArea()
      }
      XCTAssertEqual(Cell(id: cell).exactArea(), exactArea, accuracy: 1e-14)
      childIndex = (childIndex + 1) % 4
      // For ExactArea(), the best relative error we can expect is about 1e-6
      // because the precision of the unit vector coordinates is only about 1e-15
      // and the edge length of a leaf cell is about 1e-9.
      XCTAssert(abs(log(exactArea / Cell(id: cell).exactArea())) <= exactError)
      // For ApproxArea(), the areas are accurate to within a few percent.
      XCTAssert(abs(log(approxArea / Cell(id: cell).approxArea())) <= approxError)
      // For AverageArea(), the areas themselves are not very accurate, but
      // the average area of a parent is exactly 4 times the area of a child.
      XCTAssert(abs(log(avgArea / Cell(id: cell).averageArea())) <= avgError)
      //
      cell = cell.children()[childIndex]
    }
  }
  
  func testIntersectsCell() {
    let f0c2 = CellId(face: 0).childBegin(2)
    let tests = [
      (Cell(id: f0c2), Cell(id: f0c2), true),
      (Cell(id: f0c2), Cell(id: f0c2.childBegin(5)), true),
      (Cell(id: f0c2), Cell(id: f0c2.next()), false)]
    for (c, oc, want) in tests {
      XCTAssertEqual(c.intersects(oc), want)
    }
  }
  
  func testContainsCell() {
    let f0c2 = CellId(face: 0).childBegin(2)
    let tests = [
      (Cell(id: f0c2), Cell(id: f0c2), true),
      (Cell(id: f0c2), Cell(id: f0c2.childBegin(5)), true),
      (Cell(id: f0c2.childBegin(5)), Cell(id: f0c2), false),
      (Cell(id: f0c2.next()), Cell(id: f0c2), false),
      (Cell(id: f0c2), Cell(id: f0c2.next()), false)]
    for (c, oc, want) in tests {
      XCTAssertEqual(c.contains(oc), want)
    }
  }
  
  func testRectBound() {
    let tests = [
      (50.0, 50.0),
      (-50, 50),
      (50, -50),
      (-50, -50),
      (0, 0),
      (0, 180),
      (0, -179)]
    for (lat, lng) in tests {
      let c = Cell(latLng: llDegrees(lat, lng))
      let rect = c.rectBound()
      for i in 0..<4 {
        XCTAssertTrue(rect.contains(LatLng(vector: c.vertex(i).v)))
      }
    }
  }
  
  func testRectBoundAroundPoleMinLat() {
    // call this to initialize the lookup tables
    CellId.setup()
    //
    let f2 = Cell(id: CellId(face: 2, pos: 0, level: 0)).rectBound()
    XCTAssertFalse(f2.contains(llDegrees(3.0, 0)))
    XCTAssertTrue(f2.contains(llDegrees(50.0, 0)))
    let f5 = Cell(id: CellId(face: 5, pos: 0, level: 0)).rectBound()
    XCTAssertFalse(f5.contains(llDegrees(-3.0, 0)))
    XCTAssertTrue(f5.contains(llDegrees(-50.0, 0)))
  }
  
  func testCapBound() {
    let c = Cell(id: CellId(face: 0).childBegin(20))
    let s2Cap = c.capBound()
    for i in 0..<4 {
      XCTAssertTrue(s2Cap.contains(c.vertex(i)))
    }
  }
  
  func testContainsPoint() {
    let f0c2 = CellId(face: 0).childBegin(2)
    XCTAssertTrue(Cell(id: f0c2).contains(Cell(id: f0c2.childBegin(5)).vertex(1)))
    XCTAssertTrue(Cell(id: f0c2).contains(Cell(id: f0c2).vertex(1)))
    XCTAssertFalse(Cell(id: f0c2.childBegin(5)).contains(Cell(id: f0c2.next().childBegin(5)).vertex(1)))
  }
  
  func testContainsPointConsistentWithCellIdFromPoint() {
    // Construct many points that are nearly on a Cell edge, and verify that
    // Cell(id: cellIDFromPoint(p)).Contains(p) is always true.
    for _ in 0..<1000 {
      let cell = Cell(id: randomCellId())
      let i1 = randomUniformInt(n: 4)
      let i2 = (i1 + 1) & 3
      let v1 = cell.vertex(i1)
      let v2 = samplePointFromCap(c: S2Cap(center: cell.vertex(i2), angle: epsilon))
      let p = interpolate(t: randomFloat64(), a: v1, b: v2)
      XCTAssertTrue(Cell(id: CellId(point: p)).contains(p))
    }
  }
  
  func testContainsPointContainsAmbiguousPoint() {
    // This tests a case where S2CellId returns the "wrong" cell for a point
    // that is very close to the cell edge. (ConsistentWithS2CellIdFromPoint
    // generates more examples like this.)
    //
    // The Point below should have x = 0, but conversion from LatLng to
    // (x, y, z) gives x = ~6.1e-17. When xyz is converted to uv, this gives
    // u = -6.1e-17. However when converting to st, which has a range of [0, 1],
    // the low precision bits of u are lost and we wind up with s = 0.5.
    // cellIDFromPoint then chooses an arbitrary neighboring cell.
    //
    // This tests that Cell.ContainsPoint() expands the cell bounds sufficiently
    // so that the returned cell is still considered to contain p.
    let p = llDegrees(-2, 90).toPoint()
    let cell = Cell(id: CellId(point: p).parent(1))
    XCTAssertTrue(cell.contains(p))
  }

}
