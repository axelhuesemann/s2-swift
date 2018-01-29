//
//  R1IntervalTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class R1IntervalTests: XCTestCase {

  // Some standard intervals for use throughout the tests.
  let unit = R1Interval(lo: 0, hi: 1)
  let negunit = R1Interval(lo: -1, hi: 0)
  let half = R1Interval(lo: 0.5, hi: 0.5)
  let empty = R1Interval.empty
  let zero = R1Interval(lo: 0, hi: 0)
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func i(_ lo: Double, _ hi: Double) -> R1Interval {
    return R1Interval(lo: lo, hi: hi)
  }
  
  func testisEmpty() {
    XCTAssertFalse(unit.isEmpty)
    XCTAssertFalse(half.isEmpty)
    XCTAssertTrue(empty.isEmpty)
    XCTAssertFalse(zero.isEmpty)
  }

  func testCenter() {
    XCTAssertEqual(unit.center, 0.5)
    XCTAssertEqual(negunit.center, -0.5)
    XCTAssertEqual(half.center, 0.5)
  }

  func testLength() {
    XCTAssertEqual(unit.length, 1.0)
    XCTAssertEqual(negunit.length, 1.0)
    XCTAssertEqual(half.length, 0.0)
    XCTAssert(empty.length < 0.0)
  }
  
  // TODO: Tests for Contains, InteriorContains, ContainsInterval, InteriorContainsInterval, Intersects, InteriorIntersects
  
  func testIntersection() {
    XCTAssertEqual(unit.intersection(half), half)
    XCTAssertEqual(unit.intersection(negunit), zero)
    XCTAssertEqual(negunit.intersection(half), empty)
    XCTAssertEqual(unit.intersection(empty), empty)
    XCTAssertEqual(empty.intersection(unit), empty)
  }

  func testUnion() {
    XCTAssertEqual(i(99, 100).union(empty), i(99, 100))
    XCTAssertEqual(empty.union(i(99, 100)), i(99, 100))
    XCTAssertEqual(i(5, 3).union(i(0, -2)), empty)
    XCTAssertEqual(i(0, -2).union(i(5, 3)), empty)
    XCTAssertEqual(unit.union(unit), unit)
    XCTAssertEqual(unit.union(negunit), i(-1, 1))
    XCTAssertEqual(negunit.union(unit), i(-1, 1))
    XCTAssertEqual(half.union(unit), unit)
  }

  func testAddPoint() {
    XCTAssertEqual(empty.add(5.0), i(5, 5))
    XCTAssertEqual(i(5, 5).add(-1.0), i(-1, 5))
    XCTAssertEqual(i(-1, 5).add(0.0), i(-1, 5))
    XCTAssertEqual(i(-1, 5).add(6.0), i(-1, 6))
  }
  
  func testClampPoint() {
    XCTAssertEqual(i(0.1, 0.4).clamp(0.3), 0.3)
    XCTAssertEqual(i(0.1, 0.4).clamp(-7.0), 0.1)
    XCTAssertEqual(i(0.1, 0.4).clamp(0.6), 0.4)
  }
  
  func testExpanded() {
    XCTAssertEqual(empty.expanded(0.45), empty)
    XCTAssertEqual(unit.expanded(0.5), i(-0.5, 1.5))
    XCTAssertEqual(unit.expanded(-0.5), i(0.5, 0.5))
    XCTAssertEqual(unit.expanded(-0.51), empty)
  }
  
  func testIntervalString() {
    XCTAssertEqual(i(2.0 * toRadians, 4.5 * toRadians).description, "[2.0000000, 4.5000000]")
  }
  
  func testApproxEqual() {
    // empty intervals
    XCTAssertTrue(empty.approxEquals(empty))
    XCTAssertTrue(zero.approxEquals(empty))
    XCTAssertTrue(empty.approxEquals(zero))
    XCTAssertTrue(i(1, 1).approxEquals(empty))
    XCTAssertTrue(empty.approxEquals(i(1, 1)))
    XCTAssertFalse(empty.approxEquals(i(0, 1)))
    XCTAssertTrue(empty.approxEquals(i(1, 1 + 2*epsilon)))
    // singleton intervals
    XCTAssertTrue(i(1, 1).approxEquals(i(1, 1)))
    XCTAssertTrue(i(1, 1).approxEquals(i(1 - epsilon, 1 - epsilon)))
    XCTAssertTrue(i(1, 1).approxEquals(i(1 + epsilon, 1 + epsilon)))
    XCTAssertFalse(i(1, 1).approxEquals(i(1 - 3*epsilon, 1)))
    XCTAssertFalse(i(1, 1).approxEquals(i(1, 1 + 3*epsilon)))
    XCTAssertTrue(i(1, 1).approxEquals(i(1 - epsilon, 1 + epsilon)))
    XCTAssertFalse(zero.approxEquals(i(1, 1)))
    // other intervals
    XCTAssertFalse(i(1 - R1Interval.epsilon, 2 + epsilon).approxEquals(i(1, 2)))
    XCTAssertTrue(i(1 + R1Interval.epsilon, 2 - epsilon).approxEquals(i(1, 2)))
    XCTAssertFalse(i(1 - 3 * R1Interval.epsilon, 2 + R1Interval.epsilon).approxEquals(i(1, 2)))
    XCTAssertFalse(i(1 + 3 * R1Interval.epsilon, 2 - R1Interval.epsilon).approxEquals(i(1, 2)))
    XCTAssertFalse(i(1 - R1Interval.epsilon, 2 + 3 * R1Interval.epsilon).approxEquals(i(1, 2)))
    XCTAssertFalse(i(1 + R1Interval.epsilon, 2 - 3 * R1Interval.epsilon).approxEquals(i(1, 2)))
  }

}
