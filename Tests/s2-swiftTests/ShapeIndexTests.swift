//
//  S2ShapeIndexTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


// testShape is a minimal implementation of the Shape interface for use in testing
// until such time as there are other s2 types that implement it.
struct TestShape: Shape {
  
  func referencePoint() -> ReferencePoint {
    return ReferencePoint(origin: true, contained: true)
  }
  
  func numChains() -> Int {
    return 1
  }
  
  func chain(_ chainId: Int) -> Chain {
    return Chain(start: 0, length: 1)
  }
  
  func chainEdge(chainId: Int, offset: Int) -> Edge {
    return Edge(v0: a, v1: b)
  }
  
  func chainPosition(_ edgeId: Int) -> ChainPosition {
    return ChainPosition(chainId: 0, offset: 0)
  }
  
  func dimension() -> ShapeDimension {
    return .pointGeometry
  }
  
  let a = S2Point(x: 0, y: 0, z: 0)
  let b = S2Point(x: 0, y: 0, z: 0)
  let edges = 1
  
  func numEdges() -> Int {
    return edges
  }
  
  func edge(_ i: Int) -> Edge {
    return Edge(v0: a, v1: b)
  }
  
  func hasInterior() -> Bool {
    return false
  }
  
  func containsOrigin() -> Bool {
    return false
  }
  
}

class S2ShapeIndexTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }

  func testShapeIndexBasics() {
    let si = ShapeIndex()
    XCTAssertEqual(si.count, 0)
    let s = TestShape()
    si.add(shape: s)
    // TODO: once an ID is available, use that rather than assuming the first one
    // is always 0.
//    XCTAssertEqual(si.at(0), s)
    si.reset()
    XCTAssertEqual(si.count, 0)
  }

  func testShapeIndexBasics2() {
    let index = ShapeIndex()
    XCTAssertEqual(index.count, 0, "initial index should be empty after creation")
    let s = EdgeVectorShape(edges: [])
    index.add(shape: s)
    XCTAssertNotEqual(index.count, 0, "index should not be empty after adding shape")
    index.reset()
    XCTAssertEqual(index.count, 0, "initial index should be empty after reset")
  }

  func e(_ x0: Double, _ y0: Double, _ z0: Double, _ x1: Double, _ y1: Double, _ z1: Double) -> Edge {
    return Edge(v0: S2Point(x: x0, y: y0, z: z0), v1: S2Point(x: x1, y: y1, z: z1))
  }
  
  func testShapeEdgeComparisons() {
    XCTAssert(e(-1, 0, 0, 0, 0, 0) < e(0, 0, 0, 0, 0, 0))
    XCTAssert(e(0, 2, 0, 0, 0, 5) == e(0, 2, 0, 0, 0, 5))
    XCTAssert(e(1, 0, 0, -6, 7, 8) > e(0, 0, 0, 1, 3, 5))
    XCTAssert(e(5, -2, -0.4, -1, 0, 0) < e(5, -2, -0.4, 0, -1, -1))
    XCTAssert(e(9, 8, 7, 12, 3, -4) == e(9, 8, 7, 12, 3, -4))
    XCTAssert(e(-11, 7.2, -4.6, 0, 1, 0) > e(-11, 7.2, -4.6, 0, 0, 0.9))
  }
  
  func testShapeIndexCellBasics() {
    var s = ShapeIndexCell(shapes: [])
    XCTAssertEqual(s.shapes.count, 0)
    // create some clipped shapes to add.
    let c1 = ClippedShape(shapeId: 1, containsCenter: true, edges: [])
    s.add(c1)
    let c2 = ClippedShape(shapeId: 7, containsCenter: true, edges: [])
    s.add(c2)
    let c3 = ClippedShape(shapeId: 3, containsCenter: true, edges: [])
    s.add(c3)
    // look up the element at a given index
    XCTAssertEqual(s.shapes[1], c2)
    // look for the clipped shape that is part of the given shape.
    XCTAssertEqual(s.find(shapeId: 7), c2)
  }
  
  // Determines whether or not the edge defined by points A and B should be
  // present in that CellID and verify that this matches hasEdge.
  func validateEdge(a: S2Point, b: S2Point, ci: CellId, hasEdge: Bool) {
    // Expand or shrink the padding slightly to account for errors in the
    // function we use to test for intersection (IntersectsRect).
    var padding = cellPadding
    var sign = 1.0
    if !hasEdge {
      sign = -1
    }
    padding += sign * intersectsRectErrorUVDist
    let bound = ci.boundUV().expanded(padding)
    let (aUV, bUV, ok) = clipToPaddedFace(a: a, b: b, f: ci.face(), padding: padding)
    XCTAssertEqual(ok && edgeIntersectsRect(a: aUV, b: bUV, r: bound), hasEdge)
  }
  
  // Tests if the given Shape contains the center of the given CellID,
  // and that this matches the expected value of indexContainsCenter.
  func validateInterior(shape: Shape?, ci: CellId, indexContainsCenter: Bool) {
    if let shape = shape {
      XCTAssertEqual(containsBruteForce(shape: shape, point: ci.point()), indexContainsCenter)
    } else if indexContainsCenter {
        XCTFail("was nil or does not have an interior, but should have")
    }
  }
  
  // Verifies that that every cell of the index contains the correct
  // edges, and that no cells are missing from the index.  The running time of this
  // function is quadratic in the number of edges.
  func quadraticValidate(index: ShapeIndex) {
    // Iterate through a sequence of nonoverlapping cell ids that cover the
    // sphere and include as a subset all the cell ids used in the index.  For
    // each cell id, verify that the expected set of edges is present.
    // "minCellID" is the first CellID that has not been validated yet.
    var minCellId = CellId(face: 0).childBegin(CellId.maxLevel)
    var it = index.iterator()
    while true {
      // Generate a list of CellIDs ("skipped cells") that cover the gap
      // between the last cell we validated and the next cell in the index.
      var skipped: CellUnion
      if !it.done() {
        let cellId = it.cellId()
        XCTAssert(cellId >= minCellId)
        skipped = CellUnion(begin: minCellId, end: cellId.rangeMin())
        minCellId = cellId.rangeMax().next()
      } else {
        // Validate the empty cells beyond the last cell in the index.
        let end = CellId(face: 5).childEnd(CellId.maxLevel)
        skipped = CellUnion(begin: minCellId, end: end)
      }
      // Iterate through all the shapes, simultaneously validating the current
      // index cell and all the skipped cells.
      var shortEdges = 0 // number of edges counted toward subdivision
      for (id, shape) in index.shapes {
        for j in 0..<skipped.count {
          validateInterior(shape: shape, ci: skipped[j], indexContainsCenter: false)
        }
        // First check that containsCenter() is set correctly.
        var clipped: ClippedShape? = nil
        if !it.done() {
          clipped = it.indexCell()?.find(shapeId: id)
          let containsCenter = clipped != nil && clipped!.containsCenter
          validateInterior(shape: shape, ci: it.cellId(), indexContainsCenter: containsCenter)
        }
//        // If this shape has been removed, it should not be present at all.
//        if shape == nil {
//          if clipped != nil {
//            XCTFail("clipped should be nil when shape is nil")
//          }
//          continue
//        }
        // Otherwise check that the appropriate edges are present.
        for e in 0..<shape.numEdges() {
          let edge = shape.edge(e)
          for j in 0..<skipped.count {
            validateEdge(a: edge.v0, b: edge.v1, ci: skipped[j], hasEdge: false)
          }
          if !it.done() {
            let hasEdge = clipped != nil && clipped!.containsEdge(id: e)
            validateEdge(a: edge.v0, b: edge.v1, ci: it.cellId(), hasEdge: hasEdge)
            if hasEdge && it.cellId().level() < ShapeIndex.maxLevel(edge: edge) {
              shortEdges += 1
            }
          }
        }
      }
      if shortEdges > index.maxEdgesPerCell {
        XCTFail("too many edges")
      }
      if it.done() {
        break
      }
      it.next()
    }
  }

  // copies the internal state of the given iterator to a new iterator.
  func copyIterator(_ i: ShapeIndexIterator) -> ShapeIndexIterator {
    return ShapeIndexIterator(index: i.index, position: i.position, id: i.id, cell: i.cell)
  }
  
  func iteratorMethods(index: ShapeIndex) {
    var it = index.iterator()
    if it.prev() {
      XCTFail("new iterator should not be able to go backwards")
    }
    it.end()
    if !it.done() {
      XCTFail("iterator positioned at end should report as done")
    }
    var ids: [CellId] = []
    // minCellID is the first CellID in a complete traversal.
    var minCellId = CellId(face: 0).childBegin(CellId.maxLevel)
    it.begin()
    while !it.done() {
      // Get the next cell in the iterator.
      let ci = it.cellId()
      let skipped = CellUnion(begin: minCellId, end: ci.rangeMin())
      var it2 = ShapeIndexIterator(index: index, pos: .end)
      for i in 0..<skipped.count {
        if it2.locate(point: skipped[i].point()) {
          XCTFail("iterator should not have been able to find the cell %v wihich was not in the index")
        }
        XCTAssertEqual(it2.locate(cellId: skipped[i]), .disjoint, "CellId location should be Disjoint for non-existent entry")
        it2.begin()
        it2.seek(target: skipped[i])
        XCTAssertEqual(ci, it2.cellId(), "seeking the current cell in the skipped list should match the current cellid")
      }
      if ids.count > 0 {
        let prevCell = ids[ids.count - 1]
        // C++ overloads operator= to clone the iterator. We can't
        // just assign directly since it2 will than change it when
        // it should not.
        var it2 = copyIterator(it)
        if !it2.prev() {
          XCTFail("should have been able to go back because there are cells")
        }
        XCTAssertEqual(prevCell, it2.cellId(), "ShapeIndexIterator should be positioned at the beginning and not equal to last entry")
        it2.next()
        XCTAssertEqual(ci, it2.cellId(), "advancing back one spot should give us the current cell")
        it2.seek(target: prevCell)
        XCTAssertEqual(prevCell, it2.cellId(), "seek from beginning for the first previous cell should not give us the current cell")
      }
      it2.begin()
      XCTAssertEqual(ci.point(), it.center(), "point at center of current position should equal center of the crrent CellID.")
      if !it2.locate(point: it.center()) {
        XCTFail("it.LocatePoint(it.Center()) should have been able to locate the point it is currently at")
      }
      XCTAssertEqual(ci, it2.cellId(), "CellID of the Point we just located should be equal.")
      it2.begin()
      XCTAssertEqual(it2.locate(cellId: ci), .indexed)
      XCTAssertEqual(ci, it2.cellId(), "CellID of the CellID we just located should match")
      if !ci.isFace() {
        it2.begin()
        XCTAssertEqual(it2.locate(cellId: ci.immediateParent()), .subdivided)
        if it2.cellId() > ci {
          XCTFail("CellID of the immediate parent should be above the current cell")
        }
        if it2.cellId() < ci.immediateParent().rangeMin() {
          XCTFail("CellID of the current position should fall below the RangeMin of the parent")
        }
      }
      if !ci.isLeaf() {
        for i in 0..<4 {
          it2.begin()
          XCTAssertEqual(it2.locate(cellId: ci.children()[i]), .indexed)
          XCTAssertEqual(ci, it2.cellId())
        }
      }
      // Add this cellID to the set of cells to examine.
      ids.append(ci)
      // Move the minimal CellID to the next CellID past our current position.
      minCellId = ci.rangeMax().next()
      //
      it.next()
    }
  }
  
  func testShapeIndexNoEdges() {
    let index = ShapeIndex()
    let iter = index.iterator()
    if !iter.done() {
      XCTFail("iterator for empty index should start at Done but did not")
    }
    iteratorMethods(index: index)
  }
  
  func testShapeIndexOneEdge() {
    let index = ShapeIndex()
    let a = S2Point(x: 1, y: 0, z: 0)
    let b = S2Point(x: 0, y: 1, z: 0)
    let e = EdgeVectorShape(a: a, b: b)
    index.add(shape: e)
    XCTAssertEqual(index.nextId, 0, "the first element added to the index should have id 0")
    quadraticValidate(index: index)
    iteratorMethods(index: index)
  }
  
  func testShapeIndexManyIdenticalEdges() {
    let numEdges = 100
    let a = p(0.99, 0.99, 1)
    let b = p(-0.99, -0.99, 1)
    let index = ShapeIndex()
    for i in 0..<numEdges {
      index.add(shape: EdgeVectorShape(a: a, b: b))
      XCTAssertEqual(index.nextId, Int32(i))
    }
    quadraticValidate(index: index)
    iteratorMethods(index: index)
    // Since all edges span the diagonal of a face, no subdivision should
    // have occurred (with the default index options).
    var it = index.iterator()
    while !it.done() {
      XCTAssertEqual(it.cellId().level(), 0)
      it.next()
    }
  }
  
  func testShapeIndexDegenerateEdge() {
    // This test verifies that degenerate edges are supported.  The following
    // point is a cube face vertex, and so it should be indexed in 3 cells.
    let a = p(1, 1, 1)
    let shape = EdgeVectorShape(a: a, b: a)
    let index = ShapeIndex()
    index.add(shape: shape)
    quadraticValidate(index: index)
    // Check that exactly 3 index cells contain the degenerate edge.
    var count = 0
    var it = index.iterator()
    while !it.done() {
      if !it.cellId().isLeaf() {
        XCTFail("the cell for this shape should be a leaf cell.")
      }
      XCTAssertEqual(it.indexCell()?.shapes.count, 1, "there should only be one shape stored in the index cell")
      XCTAssertEqual(it.indexCell()?.shapes[0].edges.count, 1, "point should only have one edge")
      count += 1
      it.next()
    }
    XCTAssertEqual(count, 3, "expected 3 index cells")
    it.next()
  }
  
  func testShapeIndexManyTinyEdges() {
    // This test adds many edges to a single leaf cell, to check that
    // subdivision stops when no further subdivision is possible.
    // Construct two points in the same leaf cell.
    let a = CellId(point: S2Point(x: 1, y: 0, z: 0)).point()
    let b = S2Point(raw: a.v.add(R3Vector(x: 0, y: 1e-12, z: 0)))
    var shape = EdgeVectorShape(a: a, b: b)
    for _ in 0..<100 {
      shape.add(a: a, b: b)
    }
    let index = ShapeIndex()
    index.add(shape: shape)
    quadraticValidate(index: index)
    // Check that there is exactly one index cell and that it is a leaf cell.
    var it = index.iterator()
    if it.done() {
      XCTFail("ShapeIndexIterator should not be positioned at the end")
      return
    }
    if !it.cellId().isLeaf() {
      XCTFail("there should be only one leaf cell in the index but it.CellID().IsLeaf() returned false")
    }
    it.next()
    if !it.done() {
      XCTFail("ShapeIndexIterator should be positioned at the end now since there should have been only one element")
    }
  }
  
  func testShapeIndexShrinkToFitOptimization() {
    // This used to trigger a bug in the ShrinkToFit optimization. The loop
    // below contains almost all of face 0 except for a small region in the
    // 0/00000 subcell. That subcell is the only one that contains any edges.
    // This caused the index to be built only in that subcell. However, all the
    // other cells on that face should also have index entries, in order to
    // indicate that they are contained by the loop.
    let loop = S2Loop.regularLoop(center: S2Point(x: 1, y: 0.5, z: 0.5), radius: S1Angle(degrees: 89), numVertices: 100)
    let index = ShapeIndex()
    index.add(shape: loop)
    quadraticValidate(index: index)
  }
  
  func testShapeIndexMixedGeometry() {
    // This test used to trigger a bug where the presence of a shape with an
    // interior could cause shapes that don't have an interior to suddenly
    // acquire one. This would cause extra ShapeIndex cells to be created
    // that are outside the bounds of the given geometry.
    let index = ShapeIndex()
    index.add(shape: makePolyline("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6"))
    index.add(shape: makePolyline("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6"))
    index.add(shape: makePolyline("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6"))
    let cell = Cell(id: CellId(face: 0).childBegin(CellId.maxLevel))
    let loop = S2Loop(cell: cell)
    index.add(shape: loop)
    var it = index.iterator()
    // No geometry intersects face 1, so there should be no index cells there.
    let c = CellId(face: 1)
    XCTAssertEqual(it.locate(cellId: c), .disjoint)
  }
  
  func testShapeIndexLoopSpanningThreeFaces() {
    let numEdges = 100
    // Construct two loops consisting of numEdges vertices each, centered
    // around the cube vertex at the start of the Hilbert curve.
    let polygon = concentricLoopsPolygon(center: S2Point(x: 1, y: -1, z: -1), numLoops: 2, verticesPerLoop: numEdges)
    let index = ShapeIndex()
    for l in polygon.loops {
      index.add(shape: l)
    }
    quadraticValidate(index: index)
    iteratorMethods(index: index)
  }
  
}
