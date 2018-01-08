//
//  S2MetricTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2MetricTests: XCTestCase {

    override func setUp() {
        super.setUp()
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }
    
    override func tearDown() {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
        super.tearDown()
    }

  func testMetric() {
    // This is not a thorough test.
    // TODO: Exercise this more.
    XCTAssertEqual(S2CellMetric.minWidth.maxLevel(0.001256), 9)
  }

}
