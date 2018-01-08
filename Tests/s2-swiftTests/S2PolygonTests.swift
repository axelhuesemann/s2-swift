//
//  S2PolygonTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2PolygonTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func makePolygon(_ coords: [(Double, Double)]) -> S2Polygon {
    let latLngs = coords.map { llDegrees($0, $1) }
    let points = latLngs.map { S2Point(latLng: $0) }
    let loop = S2Loop(points: points)
    return S2Polygon(loop: loop)
  }
  
  func makeLoop(_ coords: [(Double, Double)]) -> S2Loop {
    let latLngs = coords.map { llDegrees($0, $1) }
    let points = latLngs.map { S2Point(latLng: $0) }
    return S2Loop(points: points)
  }
  
  func testCreate() {
    let coords = [(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)]
    let latLngs = coords.map { llDegrees($0, $1) }
    let points = latLngs.map { S2Point(latLng: $0) }
    let loop = S2Loop(points: points)
    let poly = S2Polygon(loop: loop)
    XCTAssertNotNil(poly)
  }
  
}
