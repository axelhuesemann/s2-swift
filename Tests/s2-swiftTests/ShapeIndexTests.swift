//
//  S2ShapeIndexTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


// testShape is a minimal implementation of the Shape interface for use in testing
// until such time as there are other s2 types that implement it.
struct TestShape: Shape {
  
  let a = S2Point(x: 0, y: 0, z: 0)
  let b = S2Point(x: 0, y: 0, z: 0)
  let edges = 0
  
  func numEdges() -> Int {
    return edges
  }
  
  func edge(_ i: Int) -> (S2Point, S2Point) {
    return (a, b)
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
    si.add(s)
    // TODO: once an ID is available, use that rather than assuming the first one
    // is always 0.
//    XCTAssertEqual(si.at(0), s)
    si.reset()
    XCTAssertEqual(si.count, 0)
  }

}

