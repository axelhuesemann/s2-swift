//
//  S1IntervalTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S1IntervalTests: XCTestCase {

  let empty = S1Interval.empty
  let full = S1Interval.full
  // Single-point intervals:
  let zero = S1Interval(p1: 0, p2: 0)
  let pi2 = S1Interval(p1: .pi/2, p2: .pi/2)
  let pi = S1Interval(p1: .pi, p2: .pi)
  let mipi = S1Interval(p1: -.pi, p2: -.pi) // same as pi after normalization
  let mipi2 = S1Interval(p1: -.pi/2, p2: -.pi/2)
  // Single quadrants:
  let quad1 = S1Interval(p1: 0, p2: .pi/2)
  let quad2 = S1Interval(p1: .pi/2, p2: -.pi) // equivalent to (pi/2, pi)
  let quad3 = S1Interval(p1: .pi, p2: -.pi/2)
  let quad4 = S1Interval(p1: -.pi/2, p2: 0)
  // Quadrant pairs:
  let quad12 = S1Interval(p1: 0, p2: -.pi)
  let quad23 = S1Interval(p1: .pi/2, p2: -.pi/2)
  let quad34 = S1Interval(p1: -.pi, p2: 0)
  let quad41 = S1Interval(p1: -.pi/2, p2: .pi/2)
  // Quadrant triples:
  let quad123 = S1Interval(p1: 0, p2: -.pi/2)
  let quad234 = S1Interval(p1: .pi/2, p2: 0)
  let quad341 = S1Interval(p1: .pi, p2: .pi/2)
  let quad412 = S1Interval(p1: -.pi/2, p2: -.pi)
  // Small intervals around the midpoints between quadrants,
  // such that the center of each interval is offset slightly CCW from the midpoint.
  let mid12 = S1Interval(p1: .pi/2-0.01, p2: .pi/2+0.02)
  let mid23 = S1Interval(p1: .pi-0.01, p2: -.pi+0.02)
  let mid34 = S1Interval(p1: -.pi/2-0.01, p2: -.pi/2+0.02)
  let mid41 = S1Interval(p1: -0.01, p2: 0.02)

  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func testConstructors() {
    // Check that [-π,-π] is normalized to [π,π].
    XCTAssertEqual(mipi.lo, .pi)
    XCTAssertEqual(mipi.hi, .pi)
//    XCTAssertTrue(S1Interval().isValid)
  }
  
  func testSimplePredicates() {
    XCTAssertFalse(!zero.isValid || zero.isEmpty || zero.isFull)
    XCTAssertFalse(!empty.isValid || !empty.isEmpty || empty.isFull)
    XCTAssertFalse(!empty.isInverted)
    XCTAssertFalse(!full.isValid || full.isEmpty || !full.isFull)
    XCTAssertFalse(!pi.isValid || pi.isEmpty || pi.isInverted)
    XCTAssertFalse(!mipi.isValid || mipi.isEmpty || mipi.isInverted)
  }
  
  func testCenter() {
    XCTAssertEqual(quad12.center, .pi / 2, accuracy: 1e-15)
    XCTAssertEqual(S1Interval(p1: 3.1, p2: 2.9).center, 3 - .pi, accuracy: 1e-15)
    XCTAssertEqual(S1Interval(p1: -2.9, p2: -3.1).center, .pi - 3, accuracy: 1e-15)
    XCTAssertEqual(S1Interval(p1: 2.1, p2: -2.1).center, .pi, accuracy: 1e-15)
    XCTAssertEqual(pi.center, .pi, accuracy: 1e-15)
    XCTAssertEqual(mipi.center, .pi, accuracy: 1e-15)
    XCTAssertEqual(quad23.center, .pi, accuracy: 1e-15)
    XCTAssertEqual(quad123.center, 0.75 * .pi, accuracy: 1e-15)
  }
  
  func testLength() {
    let tests = [
      (quad12, .pi),
      (pi, 0),
      (mipi, 0),
      // TODO: The C++ test for quad123 uses DOUBLE_EQ. Why?
      (quad123, 1.5 * .pi),
      // TODO: The C++ test for quad23 uses fabs. Why?
      (quad23, .pi),
      (full, 2 * .pi)]
    for (interval, want) in tests {
      XCTAssertEqual(interval.length, want)
    }
    XCTAssert(empty.length < 0)
  }
  
}

//func testContains() {
//  tests := []struct {
//    interval  Interval
//    in, out   []float64 // points that should be inside/outside the interval
//    iIn, iOut []float64 // points that should be inside/outside the interior
//  }{
//    (empty, nil, []float64{0, .pi, -.pi), nil, []float64{.pi, -.pi}),
//    (full, []float64{0, .pi, -.pi), nil, []float64{.pi, -.pi), nil),
//    (quad12, []float64{0, .pi, -.pi), nil,
//      []float64{.pi / 2), []float64{0, .pi, -.pi}),
//    (quad23, []float64{.pi / 2, -.pi / 2, .pi, -.pi), []float64{0),
//      []float64{.pi, -.pi), []float64{.pi / 2, -.pi / 2, 0}),
//    (pi, []float64{.pi, -.pi), []float64{0), nil, []float64{.pi, -.pi}),
//    (mipi, []float64{.pi, -.pi), []float64{0), nil, []float64{.pi, -.pi}),
//    (zero, []float64{0), nil, nil, []float64{0}),
//  }
//  for _, test := range tests {
//    for _, p := range test.in {
//      if !test.interval.Contains(p) {
//        t.Errorf("%v should contain %v", test.interval, p)
//      }
//    }
//    for _, p := range test.out {
//      if test.interval.Contains(p) {
//        t.Errorf("%v should not contain %v", test.interval, p)
//      }
//    }
//    for _, p := range test.iIn {
//      if !test.interval.InteriorContains(p) {
//        t.Errorf("interior of %v should contain %v", test.interval, p)
//      }
//    }
//    for _, p := range test.iOut {
//      if test.interval.InteriorContains(p) {
//        t.Errorf("interior %v should not contain %v", test.interval, p)
//      }
//    }
//  }
//}
//
//func testIntervalOperations() {
//  quad12eps := IntervalFromEndpoints(quad12.Lo, mid23.Hi)
//  quad2hi := IntervalFromEndpoints(mid23.Lo, quad12.Hi)
//  quad412eps := IntervalFromEndpoints(mid34.Lo, quad12.Hi)
//  quadeps12 := IntervalFromEndpoints(mid41.Lo, quad12.Hi)
//  quad1lo := IntervalFromEndpoints(quad12.Lo, mid41.Hi)
//  quad2lo := IntervalFromEndpoints(quad23.Lo, mid12.Hi)
//  quad3hi := IntervalFromEndpoints(mid34.Lo, quad23.Hi)
//  quadeps23 := IntervalFromEndpoints(mid12.Lo, quad23.Hi)
//  quad23eps := IntervalFromEndpoints(quad23.Lo, mid34.Hi)
//  quadeps123 := IntervalFromEndpoints(mid41.Lo, quad23.Hi)
//  
//  // This massive list of test cases is ported directly from the C++ test case.
//  tests := []struct {
//    x, y                               Interval
//    xContainsY, xInteriorContainsY     bool
//    xIntersectsY, xInteriorIntersectsY bool
//    wantUnion, wantIntersection        Interval
//  }{
//    // 0
//    (empty, empty, true, true, false, false, empty, empty),
//    (empty, full, false, false, false, false, full, empty),
//    (empty, zero, false, false, false, false, zero, empty),
//    (empty, pi, false, false, false, false, pi, empty),
//    (empty, mipi, false, false, false, false, mipi, empty),
//    
//    // 5
//    (full, empty, true, true, false, false, full, empty),
//    (full, full, true, true, true, true, full, full),
//    (full, zero, true, true, true, true, full, zero),
//    (full, pi, true, true, true, true, full, pi),
//    (full, mipi, true, true, true, true, full, mipi),
//    (full, quad12, true, true, true, true, full, quad12),
//    (full, quad23, true, true, true, true, full, quad23),
//    
//    // 12
//    (zero, empty, true, true, false, false, zero, empty),
//    (zero, full, false, false, true, false, full, zero),
//    (zero, zero, true, false, true, false, zero, zero),
//    (zero, pi, false, false, false, false, IntervalFromEndpoints(0, .pi), empty),
//    (zero, pi2, false, false, false, false, quad1, empty),
//    (zero, mipi, false, false, false, false, quad12, empty),
//    (zero, mipi2, false, false, false, false, quad4, empty),
//    (zero, quad12, false, false, true, false, quad12, zero),
//    (zero, quad23, false, false, false, false, quad123, empty),
//    
//    // 21
//    (pi2, empty, true, true, false, false, pi2, empty),
//    (pi2, full, false, false, true, false, full, pi2),
//    (pi2, zero, false, false, false, false, quad1, empty),
//    (pi2, pi, false, false, false, false, IntervalFromEndpoints(.pi/2, .pi), empty),
//    (pi2, pi2, true, false, true, false, pi2, pi2),
//    (pi2, mipi, false, false, false, false, quad2, empty),
//    (pi2, mipi2, false, false, false, false, quad23, empty),
//    (pi2, quad12, false, false, true, false, quad12, pi2),
//    (pi2, quad23, false, false, true, false, quad23, pi2),
//    
//    // 30
//    (pi, empty, true, true, false, false, pi, empty),
//    (pi, full, false, false, true, false, full, pi),
//    (pi, zero, false, false, false, false, IntervalFromEndpoints(.pi, 0), empty),
//    (pi, pi, true, false, true, false, pi, pi),
//    (pi, pi2, false, false, false, false, IntervalFromEndpoints(.pi/2, .pi), empty),
//    (pi, mipi, true, false, true, false, pi, pi),
//    (pi, mipi2, false, false, false, false, quad3, empty),
//    (pi, quad12, false, false, true, false, IntervalFromEndpoints(0, .pi), pi),
//    (pi, quad23, false, false, true, false, quad23, pi),
//    
//    // 39
//    (mipi, empty, true, true, false, false, mipi, empty),
//    (mipi, full, false, false, true, false, full, mipi),
//    (mipi, zero, false, false, false, false, quad34, empty),
//    (mipi, pi, true, false, true, false, mipi, mipi),
//    (mipi, pi2, false, false, false, false, quad2, empty),
//    (mipi, mipi, true, false, true, false, mipi, mipi),
//    (mipi, mipi2, false, false, false, false, IntervalFromEndpoints(-.pi, -.pi/2), empty),
//    (mipi, quad12, false, false, true, false, quad12, mipi),
//    (mipi, quad23, false, false, true, false, quad23, mipi),
//    
//    // 48
//    (quad12, empty, true, true, false, false, quad12, empty),
//    (quad12, full, false, false, true, true, full, quad12),
//    (quad12, zero, true, false, true, false, quad12, zero),
//    (quad12, pi, true, false, true, false, quad12, pi),
//    (quad12, mipi, true, false, true, false, quad12, mipi),
//    (quad12, quad12, true, false, true, true, quad12, quad12),
//    (quad12, quad23, false, false, true, true, quad123, quad2),
//    (quad12, quad34, false, false, true, false, full, quad12),
//    
//    // 56
//    (quad23, empty, true, true, false, false, quad23, empty),
//    (quad23, full, false, false, true, true, full, quad23),
//    (quad23, zero, false, false, false, false, quad234, empty),
//    (quad23, pi, true, true, true, true, quad23, pi),
//    (quad23, mipi, true, true, true, true, quad23, mipi),
//    (quad23, quad12, false, false, true, true, quad123, quad2),
//    (quad23, quad23, true, false, true, true, quad23, quad23),
//    (quad23, quad34, false, false, true, true, quad234, IntervalFromEndpoints(-.pi, -.pi/2)),
//    
//    // 64
//    (quad1, quad23, false, false, true, false, quad123, IntervalFromEndpoints(.pi/2, .pi/2)),
//    (quad2, quad3, false, false, true, false, quad23, mipi),
//    (quad3, quad2, false, false, true, false, quad23, pi),
//    (quad2, pi, true, false, true, false, quad2, pi),
//    (quad2, mipi, true, false, true, false, quad2, mipi),
//    (quad3, pi, true, false, true, false, quad3, pi),
//    (quad3, mipi, true, false, true, false, quad3, mipi),
//    
//    // 71
//    (quad12, mid12, true, true, true, true, quad12, mid12),
//    (mid12, quad12, false, false, true, true, quad12, mid12),
//    
//    // 73
//    (quad12, mid23, false, false, true, true, quad12eps, quad2hi),
//    (mid23, quad12, false, false, true, true, quad12eps, quad2hi),
//    
//    // This test checks that the union of two disjoint intervals is the smallest
//    // interval that contains both of them.  Note that the center of "mid34"
//    // slightly CCW of -Pi/2 so that there is no ambiguity about the result.
//    // 75
//    (quad12, mid34, false, false, false, false, quad412eps, empty),
//    (mid34, quad12, false, false, false, false, quad412eps, empty),
//    
//    // 77
//    (quad12, mid41, false, false, true, true, quadeps12, quad1lo),
//    (mid41, quad12, false, false, true, true, quadeps12, quad1lo),
//    
//    // 79
//    (quad23, mid12, false, false, true, true, quadeps23, quad2lo),
//    (mid12, quad23, false, false, true, true, quadeps23, quad2lo),
//    (quad23, mid23, true, true, true, true, quad23, mid23),
//    (mid23, quad23, false, false, true, true, quad23, mid23),
//    (quad23, mid34, false, false, true, true, quad23eps, quad3hi),
//    (mid34, quad23, false, false, true, true, quad23eps, quad3hi),
//    (quad23, mid41, false, false, false, false, quadeps123, empty),
//    (mid41, quad23, false, false, false, false, quadeps123, empty),
//  }
//  should := func(b bool) string {
//    if b {
//      return "should"
//    }
//    return "should not"
//  }
//  for _, test := range tests {
//    if test.x.ContainsInterval(test.y) != test.xContainsY {
//      t.Errorf("%v %s contain %v", test.x, should(test.xContainsY), test.y)
//    }
//    if test.x.InteriorContainsInterval(test.y) != test.xInteriorContainsY {
//      t.Errorf("interior of %v %s contain %v", test.x, should(test.xInteriorContainsY), test.y)
//    }
//    if test.x.Intersects(test.y) != test.xIntersectsY {
//      t.Errorf("%v %s intersect %v", test.x, should(test.xIntersectsY), test.y)
//    }
//    if test.x.InteriorIntersects(test.y) != test.xInteriorIntersectsY {
//      t.Errorf("interior of %v %s intersect %v", test.x, should(test.xInteriorIntersectsY), test.y)
//    }
//    if u := test.x.Union(test.y); u != test.wantUnion {
//      t.Errorf("%v ∪ %v was %v, want %v", test.x, test.y, u, test.wantUnion)
//    }
//    if u := test.x.Intersection(test.y); u != test.wantIntersection {
//      t.Errorf("%v ∩ %v was %v, want %v", test.x, test.y, u, test.wantIntersection)
//    }
//  }
//}
//
//func testAddPoint() {
//  tests := []struct {
//    interval Interval
//    points   []float64
//    want     Interval
//  }{
//    (empty, []float64{0), zero),
//    (empty, []float64{.pi), pi),
//    (empty, []float64{-.pi), mipi),
//    (empty, []float64{.pi, -.pi), pi),
//    (empty, []float64{-.pi, .pi), mipi),
//    (empty, []float64{mid12.Lo, mid12.Hi), mid12),
//    (empty, []float64{mid23.Lo, mid23.Hi), mid23),
//    
//    (quad1, []float64{-0.9 * .pi, -.pi / 2), quad123),
//    (full, []float64{0), full),
//    (full, []float64{.pi), full),
//    (full, []float64{-.pi), full),
//  }
//  for _, test := range tests {
//    got := test.interval
//    for _, point := range test.points {
//      got = got.AddPoint(point)
//    }
//    want := test.want
//    if math.Abs(got.Lo-want.Lo) > 1e-15 || math.Abs(got.Hi-want.Hi) > 1e-15 {
//      t.Errorf("%v.AddPoint(%v) = %v, want %v", test.interval, test.points, got, want)
//    }
//  }
//}
//
//func testExpanded() {
//  tests := []struct {
//    interval Interval
//    margin   float64
//    want     Interval
//  }{
//    (empty, 1, empty),
//    (full, 1, full),
//    (zero, 1, Interval{-1, 1}),
//    (mipi, 0.01, Interval{.pi - 0.01, -.pi + 0.01}),
//    (pi, 27, full),
//    (pi, .pi / 2, quad23),
//    (pi2, .pi / 2, quad12),
//    (mipi2, .pi / 2, quad34),
//    
//    (empty, -1, empty),
//    (full, -1, full),
//    (quad123, -27, empty),
//    (quad234, -27, empty),
//    (quad123, -.pi / 2, quad2),
//    (quad341, -.pi / 2, quad4),
//    (quad412, -.pi / 2, quad1),
//  }
//  for _, test := range tests {
//    if got, want := test.interval.Expanded(test.margin), test.want; math.Abs(got.Lo-want.Lo) > 1e-15 || math.Abs(got.Hi-want.Hi) > 1e-15 {
//      t.Errorf("%v.Expanded(%v) = %v, want %v", test.interval, test.margin, got, want)
//    }
//  }
//}
//
//func testIntervalString() {
//  if s, exp := pi.String(), "[3.1415927, 3.1415927]"; s != exp {
//    t.Errorf("pi.String() = %q, want %q", s, exp)
//  }
//}
