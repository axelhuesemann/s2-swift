//
//  S2LatLngTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2LatLngTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }

  func l(_ latDegrees: Double, _ lngDegrees: Double) -> LatLng {
    return llDegrees(latDegrees, lngDegrees)
  }
  
  func testLatLngNormalized() {
    let tests = [
      ("Valid lat/lng", l(21.8275043, 151.1979675), l(21.8275043, 151.1979675)),
      ("Valid lat/lng in the West", l(21.8275043, -151.1979675), l(21.8275043, -151.1979675)),
      ("Beyond the North pole", l(95, 151.1979675), l(90, 151.1979675)),
      ("Beyond the South pole", l(-95, 151.1979675), l(-90, 151.1979675)),
      ("At the date line (from East)", l(21.8275043, 180), l(21.8275043, 180)),
      ("At the date line (from West)", l(21.8275043, -180), l(21.8275043, -180)),
      ("Across the date line going East", l(21.8275043, 181.0012), l(21.8275043, -178.9988)),
      ("Across the date line going West", l(21.8275043, -181.0012), l(21.8275043, 178.9988)),
      ("All wrong", l(256, 256), l(90, -104))]
    for (desc, pos, want) in tests {
      let got = pos.normalized()
      // NSLog("\(got)")
      XCTAssertTrue(got.isValid, desc)
      XCTAssert(got.distance(want) < 1e-13 * toRadians, desc)
    }
  }
  
  func testLatLngString() {
    XCTAssertEqual(llDegrees(sqrt(2.0), -sqrt(5.0)).description, "[1.4142136, -2.2360680]")
  }
  
  func testLatLngPointConversion() {
    let eps = 1e-14
    // All test cases here have been verified against the C++ S2 implementation.
    let tests = [
      (0, 0, 1, 0, 0),
      (90, 0, 6.12323e-17, 0, 1),
      (-90, 0, 6.12323e-17, 0, -1),
      (0, 180, -1, 1.22465e-16, 0),
      (0, -180, -1, -1.22465e-16, 0),
      (90, 180, -6.12323e-17, 7.4988e-33, 1),
      (90, -180, -6.12323e-17, -7.4988e-33, 1),
      (-90, 180, -6.12323e-17, 7.4988e-33, -1),
      (-90, -180, -6.12323e-17, -7.4988e-33, -1),
      (-81.82750430354997, 151.19796752929685, -0.12456788151479525, 0.0684875268284729, -0.989844584550441)]
    for (lat, lng, x, y, z) in tests {
      let ll = llDegrees(lat, lng)
      let p = ll.toPoint()
      // TODO: Port Point.ApproxEquals, then use here.
      XCTAssertEqual(p.x, x, accuracy: eps)
      XCTAssertEqual(p.y, y, accuracy: eps)
      XCTAssertEqual(p.z, z, accuracy: eps)
      let ll2 = LatLng(point: p)
      // We need to be careful here, since if the latitude is +/- 90, any longitude
      // is now a valid conversion.
      let isPolar = (lat == 90 || lat == -90)
      XCTAssertEqual(ll2.lat, lat * toRadians, accuracy: eps)
      if !isPolar {
        XCTAssertEqual(ll2.lng, lng * toRadians, accuracy: eps)
      }
    }
  }
  
  func testLatLngDistance() {
    // Based on C++ S2LatLng::TestDistance.
    let tests = [
      (90.0, 0.0, 90.0, 0.0, 0.0, 0.0),
      (-37, 25, -66, -155, 77, 1e-13),
      (0, 165, 0, -80, 115, 1e-13),
      (47, -127, -47, 53, 180, 2e-6)]
    for (lat1, lng1, lat2, lng2, want, tolerance) in tests {
      let ll1 = llDegrees(lat1, lng1)
      let ll2 = llDegrees(lat2, lng2)
      let d = ll1.distance(ll2)
      XCTAssertEqual(d, want * toRadians, accuracy: tolerance)
    }
  }

}
