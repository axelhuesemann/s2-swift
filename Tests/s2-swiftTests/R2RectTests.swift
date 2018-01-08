//
//  R2RectTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class R2RectTests: XCTestCase {

  let sw = R2Point(x: 0, y: 0.25)
  let se = R2Point(x: 0.5, y: 0.25)
  let ne = R2Point(x: 0.5, y: 0.75)
  let nw = R2Point(x: 0, y: 0.75)
  
  let empty = R2Rect.empty
  let rect = R2Rect(p0: R2Point(x: 0, y: 0.25), p1: R2Point(x: 0.5, y: 0.75))
  let rectMid = R2Rect(p0: R2Point(x: 0.25, y: 0.5), p1: R2Point(x: 0.25, y: 0.5))
  let rectSW = R2Rect(p0: R2Point(x: 0, y: 0.25), p1: R2Point(x: 0, y: 0.25))
  let rectNE = R2Rect(p0: R2Point(x: 0.5, y: 0.75), p1: R2Point(x: 0.5, y: 0.75))

  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func p(_ x: Double, _ y: Double) -> R2Point {
    return R2Point(x: x, y: y)
  }
  
  func r(_ p0: R2Point, _ p1: R2Point) -> R2Rect {
    return R2Rect(p0: p0, p1: p1)
  }
  
  func testEmptyRect() {
    XCTAssertTrue(empty.isValid)
    XCTAssertTrue(empty.isEmpty)
  }
  
  func testFromVariousTypes() {
    let d1 = r(p(0.1, 0), p(0.25, 1))
    XCTAssert(R2Rect(center: p(0.3, 0.5), size: p(0.2, 0.4)).approxEquals(r(p(0.2, 0.3), p(0.4, 0.7))))
    XCTAssert(R2Rect(center: p(1, 0.1), size: p(0, 2)).approxEquals(r(p(1, -0.9), p(1, 1.1))))
    XCTAssert(d1.approxEquals(R2Rect(x: d1.x, y: d1.y)))
    // TODO hi before lo case in constructor
//    XCTAssert(r(p(0.15, 0.3), p(0.35, 0.9)).approxEquals(r(p(0.15, 0.9), p(0.35, 0.3))))
//    XCTAssert(r(p(0.12, 0), p(0.83, 0.5)).approxEquals(r(p(0.83, 0), p(0.12, 0.5))))
  }
    
  func testCenter() {
    XCTAssertEqual(empty.center, p(0.5, 0.5))
    XCTAssertEqual(rect.center, p(0.25, 0.5))
  }

  func testVertices() {
    let v = rect.vertices
    XCTAssertEqual(v[0], sw)
    XCTAssertEqual(v[1], se)
    XCTAssertEqual(v[2], ne)
    XCTAssertEqual(v[3], nw)
  }

  func testContainsPoint() {
    XCTAssertTrue(rect.contains(p(0.2, 0.4)))
    XCTAssertFalse(rect.contains(p(0.2, 0.8)))
    XCTAssertFalse(rect.contains(p(-0.1, 0.4)))
    XCTAssertFalse(rect.contains(p(0.6, 0.1)))
    XCTAssertTrue(rect.contains(p(rect.x.lo, rect.y.lo)))
    XCTAssertTrue(rect.contains(p(rect.x.hi, rect.y.hi)))
  }

  func testInteriorContainsPoint() {
    // Check corners are not contained.
    XCTAssertFalse(rect.interiorContains(sw))
    XCTAssertFalse(rect.interiorContains(ne))
    // Check a point on the border is not contained.
    XCTAssertFalse(rect.interiorContains(p(0, 0.5)))
    XCTAssertFalse(rect.interiorContains(p(0.25, 0.25)))
    XCTAssertFalse(rect.interiorContains(p(0.5, 0.5)))
    // Check points inside are contained.
    XCTAssertTrue(rect.interiorContains(p(0.125, 0.6)))
  }
              
  func testIntervalOps() {
    let tests = [
      (rect, rectMid, true, true, true, true, rect, rectMid),
      (rect, rectSW, true, false, true, false, rect, rectSW),
      (rect, rectNE, true, false, true, false, rect, rectNE),
      (rect, r(p(0.45, 0.1), p(0.75, 0.3)), false, false, true, true, r(p(0, 0.1), p(0.75, 0.75)), r(p(0.45, 0.25), p(0.5, 0.3))),
      (rect, r(p(0.5, 0.1), p(0.7, 0.3)), false, false, true, false, r(p(0, 0.1), p(0.7, 0.75)), r(p(0.5, 0.25), p(0.5, 0.3))),
      (rect, r(p(0.45, 0.1), p(0.7, 0.25)), false, false, true, false, r(p(0, 0.1), p(0.7, 0.75)), r(p(0.45, 0.25), p(0.5, 0.25))),
      (r(p(0.1, 0.2), p(0.1, 0.3)), r(p(0.15, 0.7), p(0.2, 0.8)), false, false, false, false, r(p(0.1, 0.2), p(0.2, 0.8)), R2Rect.empty),
      // Check that the intersection of two rectangles that overlap in x but not y is valid, and vice versa.
      (r(p(0.1, 0.2), p(0.4, 0.5)), r(p(0, 0), p(0.2, 0.1)), false, false, false, false, r(p(0, 0), p(0.4, 0.5)), R2Rect.empty),
      (r(p(0, 0), p(0.1, 0.3)), r(p(0.2, 0.1), p(0.3, 0.4)), false, false, false, false, r(p(0, 0), p(0.3, 0.4)), R2Rect.empty)]

    for (r1, r2, contains, intContains, intersects, intIntersects, wantUnion, wantIntersection) in tests {
      XCTAssertEqual(r1.contains(r2), contains)
      XCTAssertEqual(r1.interiorContains(r2), intContains)
      XCTAssertEqual(r1.intersects(r2), intersects)
      XCTAssertEqual(r1.interiorIntersects(r2), intIntersects)
      XCTAssertEqual(r1.union(r2).approxEquals(r1), r1.contains(r2))
      XCTAssertEqual(!r1.intersection(r2).isEmpty, r1.intersects(r2))
      XCTAssertEqual(r1.union(r2), wantUnion)
      XCTAssertEqual(r1.intersection(r2), wantIntersection)
      XCTAssertEqual(r1.add(r2), wantUnion)
    }
  }
                  
  func testAddPoint() {
    var r2 = R2Rect.empty
    r2 = r2.add(sw)
    r2 = r2.add(se)
    r2 = r2.add(nw)
    r2 = r2.add(p(0.1, 0.4))
    XCTAssertTrue(rect.approxEquals(r2))
  }

  func testClampPoint() {
    let r = R2Rect(x: R1Interval(lo: 0, hi: 0.5), y: R1Interval(lo: 0.25, hi: 0.75))
    XCTAssertEqual(r.clamp(p(-0.01, 0.24)), p(0, 0.25))
    XCTAssertEqual(r.clamp(p(-5.0, 0.48)), p(0, 0.48))
    XCTAssertEqual(r.clamp(p(-5.0, 2.48)), p(0, 0.75))
    XCTAssertEqual(r.clamp(p(0.19, 2.48)), p(0.19, 0.75))
    XCTAssertEqual(r.clamp(p(6.19, 2.48)), p(0.5, 0.75))
    XCTAssertEqual(r.clamp(p(6.19, 0.53)), p(0.5, 0.53))
    XCTAssertEqual(r.clamp(p(6.19, -2.53)), p(0.5, 0.25))
    XCTAssertEqual(r.clamp(p(0.33, -2.53)), p(0.33, 0.25))
    XCTAssertEqual(r.clamp(p(0.33, 0.37)), p(0.33, 0.37))
  }
                      
  func testExpandedEmpty() {
    XCTAssertTrue(R2Rect.empty.expanded(p(0.1, 0.3)).isEmpty)
    XCTAssertTrue(R2Rect.empty.expanded(p(-0.1, -0.3)).isEmpty)
    XCTAssertTrue(r(p(0.2, 0.4), p(0.3, 0.7)).expanded(p(-0.1, 0.3)).isEmpty)
    XCTAssertTrue(r(p(0.2, 0.4), p(0.3, 0.7)).expanded(p(0.1, -0.2)).isEmpty)
  }
          
  func testExpandedEquals() {
    XCTAssertTrue(r(p(0.2, 0.4), p(0.3, 0.7)).expanded(p(0.1, 0.3)).approxEquals(r(p(0.1, 0.1), p(0.4, 1.0))))
    XCTAssertTrue(r(p(0.2, 0.4), p(0.3, 0.7)).expanded(p(0.1, -0.1)).approxEquals(r(p(0.1, 0.5), p(0.4, 0.6))))
    XCTAssertTrue(r(p(0.2, 0.4), p(0.3, 0.7)).expanded(p(0.1, 0.1)).approxEquals(r(p(0.1, 0.3), p(0.4, 0.8))))
  }

}
