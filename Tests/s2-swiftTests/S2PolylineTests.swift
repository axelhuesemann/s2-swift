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
  
  func makePolyline(_ coords: [(Double, Double)]) -> S2Polyline {
    let latLngs = coords.map { llDegrees($0, $1) }
    return S2Polyline(latLngs: latLngs)
  }
  
  func testCreate() {
    let coords = [(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)]
    let latLngs = coords.map { llDegrees($0, $1) }
    let poly = S2Polyline(latLngs: latLngs)
    XCTAssertNotNil(poly)
    let poly2 = makePolyline([(37.5, 122.5), (37.5, 122), (37, 122), (37, 122.5)])
    XCTAssertNotNil(poly2)
  }
  
}
