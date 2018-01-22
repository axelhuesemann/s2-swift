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

}

