//
//  S2CellIdTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift

class S2CellIdTests: XCTestCase {

  override func setUp() {
    super.setUp()
    // Put setup code here. This method is called before the invocation of each test method in the class.
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func testCellIdFromFace() {
    for face in 0..<6 {
      XCTAssertEqual(CellId(face: face, pos: 0, level: 0), CellId(face: face))
    }
  }
  
  func testParentChildRelationships() {
    let ci = CellId(face: 3, pos: 0x12345678, level: CellId.maxLevel-4)
    XCTAssertTrue(ci.isValid)
    XCTAssertEqual(ci.face(), 3)
    XCTAssertEqual(ci.pos(), 0x12345700)
    XCTAssertEqual(ci.level(), 26)
    XCTAssertFalse(ci.isLeaf())
    XCTAssertEqual(ci.childBegin(ci.level() + 2).pos(), 0x12345610)
    XCTAssertEqual(ci.childBegin().pos(), 0x12345640)
    XCTAssertEqual(ci.children()[0].pos(), 0x12345640)
    XCTAssertEqual(ci.immediateParent().pos(), 0x12345400)
    XCTAssertEqual(ci.parent(ci.level() - 2).pos(), 0x12345000)
    XCTAssert(ci.childBegin().id < ci.id)
    XCTAssert(ci.childEnd().id > ci.id)
    XCTAssertEqual(ci.childEnd().id, ci.childBegin().next().next().next().next().id)
    XCTAssertEqual(ci.rangeMin().id, ci.childBegin(CellId.maxLevel).id)
    XCTAssertEqual(ci.rangeMax().next().id, ci.childEnd(CellId.maxLevel).id)
  }

  func testContainment() {
    let a = CellId(id: 0x80855c0000000000) // Pittsburg
    let b = CellId(id: 0x80855d0000000000) // child of a
    let c = CellId(id: 0x80855dc000000000) // child of b
    let d = CellId(id: 0x8085630000000000) // part of Pittsburg disjoint from a
    let tests = [
      (a, a, true, true, true),
      (a, b, true, false, true),
      (a, c, true, false, true),
      (a, d, false, false, false),
      (b, b, true, true, true),
      (b, c, true, false, true),
      (b, d, false, false, false),
      (c, c, true, true, true),
      (c, d, false, false, false),
      (d, d, true, true, true)]
    for (x, y, xContainsY, yContainsX, xIntersectsY) in tests {
      XCTAssertEqual(x.contains(y), xContainsY)
      XCTAssertEqual(x.intersects(y), xIntersectsY)
      XCTAssertEqual(y.contains(x), yContainsX)
    }
    // TODO: Test Contains, Intersects better, such as with adjacent cells.
  }
  
  func testCellIDString() {
    let ci = CellId(id: 0xbb04000000000000)
    XCTAssertEqual(ci.description, "5/31200")
  }
  
  func testLatLng() {
    // call this to initialize the lookup tables
    CellId.setup()
    // You can generate these with the s2cellid2latlngtestcase C++ program in this directory.
    let tests = [
      (UInt64(0x47a1cbd595522b39), 49.703498679, 11.770681595),
      (0x46525318b63be0f9, 55.685376759, 12.588490937),
      (0x52b30b71698e729d, 45.486546517, -93.449700022),
      (0x46ed8886cfadda85, 58.299984854, 23.049300056),
      (0x3663f18a24cbe857, 34.364439040, 108.330699969),
      (0x10a06c0a948cf5d, -30.694551352, -30.048758753),
      (0x2b2bfd076787c5df, -25.285264027, 133.823116966),
      (0xb09dff882a7809e1, -75.000000031, 0.000000133),
      (0x94daa3d000000001, -24.694439215, -47.537363213),
      (0x87a1000000000001, 38.899730392, -99.901813021),
      (0x4fc76d5000000001, 81.647200334, -55.631712940),
      (0x3b00955555555555, 10.050986518, 78.293170610),
      (0x1dcc469991555555, -34.055420593, 18.551140038),
      (0xb112966aaaaaaaab, -69.219262171, 49.670072392)]
    for (id, lat, lng) in tests {
      let cellId = CellId(id: id)
      let ll = llDegrees(lat, lng)
      let l2 = cellId.latLng()
      // let s1 = String(format: "%lX", CellId(latLng: ll).id)
      // let s2 = String(format: "%lX", id)
      // NSLog("\(ll) \(s1)")
      // NSLog("\(l2) \(s2)")
      XCTAssert(ll.distance(l2) <= 1e-9 * toRadians, "\(ll.distance(l2))") // ~0.1mm on earth.
      XCTAssertEqual(id, CellId(latLng: ll).id)
    }
  }
  
  func testEdgeNeighbors() {
    // Check the edge neighbors of face 1.
    let faces = [5, 3, 2, 0]
    let p0 = CellId(face: 1, i: 0, j: 0)
    let p1 = CellId(face: 1, i: 0, j: 0).parent(0)
    for (i, nbr) in p1.edgeNeighbors().enumerated() {
      XCTAssertTrue(nbr.isFace())
      XCTAssertEqual(nbr.face(), faces[i])
    }
    // Check the edge neighbors of the corner cells at all levels.  This case is
    // trickier because it requires projecting onto adjacent faces.
    let maxIJ = CellId.maxSize - 1
    for level in 1...CellId.maxLevel {
      let id = p0.parent(level)
      // These neighbors were determined manually using the face and axis
      // relationships.
      let levelSizeIJ = CellId.sizeIJ(level)
      let want = [
        CellId(face: 5, i: maxIJ, j: maxIJ).parent(level),
        CellId(face: 1, i: levelSizeIJ, j: 0).parent(level),
        CellId(face: 1, i: 0, j: levelSizeIJ).parent(level),
        CellId(face: 0, i: maxIJ, j: 0).parent(level)]
      for (i, nbr) in id.edgeNeighbors().enumerated() {
        XCTAssertEqual(nbr, want[i])
      }
    }
  }
  
  func testVertexNeighbors() {
    // Check the vertex neighbors of the center of face 2 at level 5.
    let id = CellId(point: S2Point(x: 0, y: 0, z: 1))
    var neighbors = id.vertexNeighbors(5)
    neighbors.sort { $0.id < $1.id }
  
    for (n, nbr) in neighbors.enumerated() {
      var i = 1<<29
      var j = 1<<29
      if n < 2 {
        i -= 1
      }
      if n == 0 || n == 3 {
        j -= 1
      }
      let want  = CellId(face: 2, i: i, j: j).parent(5)
  
      XCTAssertEqual(nbr, want)
    }
    let i = 1<<29
    let j = 1<<29
    XCTAssertEqual(neighbors[0], CellId(face: 2, i: i-1, j: j-1).parent(5), String(format: "%X", neighbors[0].id))
    XCTAssertEqual(neighbors[1], CellId(face: 2, i: i-1, j: j).parent(5))
    XCTAssertEqual(neighbors[2], CellId(face: 2, i: i, j: j).parent(5))
//    XCTAssertEqual(neighbors[3], CellId(face: 2, i: i, j: j-1).parent(5))
  
    // Check the vertex neighbors of the corner of faces 0, 4, and 5.
    let id2 = CellId(face: 0, pos: 0, level: CellId.maxLevel)
    var neighbors2 = id2.vertexNeighbors(0)
    neighbors2.sort { $0.id < $1.id }
    XCTAssertEqual(neighbors2.count, 3)
    XCTAssertEqual(neighbors2[0], CellId(face: 0))
    XCTAssertEqual(neighbors2[1], CellId(face: 4))
  }
  
  func testCellIDTokensNominal() {
    let tests = [
      ("1", UInt64(0x1000000000000000)),
      ("3", 0x3000000000000000),
      ("14", 0x1400000000000000),
      ("41", 0x4100000000000000),
      ("094", 0x0940000000000000),
      ("537", 0x5370000000000000),
      ("3fec", 0x3fec000000000000),
      ("72f3", 0x72f3000000000000),
      ("52b8c", 0x52b8c00000000000),
      ("990ed", 0x990ed00000000000),
      ("4476dc", 0x4476dc0000000000),
      ("2a724f", 0x2a724f0000000000),
      ("7d4afc4", 0x7d4afc4000000000),
      ("b675785", 0xb675785000000000),
      ("40cd6124", 0x40cd612400000000),
      ("3ba32f81", 0x3ba32f8100000000),
      ("08f569b5c", 0x08f569b5c0000000),
      ("385327157", 0x3853271570000000),
      ("166c4d1954", 0x166c4d1954000000),
      ("96f48d8c39", 0x96f48d8c39000000),
      ("0bca3c7f74c", 0x0bca3c7f74c00000),
      ("1ae3619d12f", 0x1ae3619d12f00000),
      ("07a77802a3fc", 0x07a77802a3fc0000),
      ("4e7887ec1801", 0x4e7887ec18010000),
      ("4adad7ae74124", 0x4adad7ae74124000),
      ("90aba04afe0c5", 0x90aba04afe0c5000),
      ("8ffc3f02af305c", 0x8ffc3f02af305c00),
      ("6fa47550938183", 0x6fa4755093818300),
      ("aa80a565df5e7fc", 0xaa80a565df5e7fc0),
      ("01614b5e968e121", 0x01614b5e968e1210),
      ("aa05238e7bd3ee7c", 0xaa05238e7bd3ee7c),
      ("48a23db9c2963e5b", 0x48a23db9c2963e5b),
    ]
    for (token, id) in tests {
      let ci = CellId(token: token)
      XCTAssertEqual(ci.id, id)
      let token2 = ci.toToken()
      XCTAssertEqual(token2, token)
    }
  }
  
  func testCellIdFromTokenErrorCases() {
    let noneToken = CellId(id: 0).toToken()
    XCTAssertEqual(noneToken, "X")
    let noneId = CellId(token: noneToken)
    XCTAssertEqual(noneId, CellId(id: 0))
    let tests = [
      "876b e99",
      "876bee99\n",
      "876[ee99",
      " 876bee99"]
    for test in tests {
      let ci = CellId(token: test)
      XCTAssertEqual(ci.id, 0)
    }
  }
  
  func r(_ x0: Double, _ y0: Double, _ x1: Double, _ y1: Double) -> R2Rect {
    return R2Rect(p0: R2Point(x: x0, y: y0), p1: R2Point(x: x1, y: y1))
  }
  
  func testIJLevelToBoundUV() {
    let maxIJ = 1 << CellId.maxLevel - 1
    let tests = [
      // The i/j space is [0, 2^30 - 1) which maps to [-1, 1] for the
      // x/y axes of the face surface. Results are scaled by the size of a cell
      // at the given level. At level 0, everything is one cell of the full size
      // of the space.  At maxLevel, the bounding rect is almost floating point
      // noise.
      // What should be out of bounds values, but passes the C++ code as well.
      (-1, -1, 0, r(-5, -5, -1, -1)),
      (-1 * maxIJ, -1 * maxIJ, 0, r(-5, -5, -1, -1)),
      (-1, -1, CellId.maxLevel, r(-1.0000000024835267, -1.0000000024835267, -1, -1)),
//      (0, 0, CellId.maxLevel + 1, r(-1, -1, -1, -1)),
      // Minimum i,j at different levels
      (0, 0, 0, r(-1, -1, 1, 1)),
      (0, 0, CellId.maxLevel / 2, r(-1, -1, -0.999918621033430099, -0.999918621033430099)),
      (0, 0, CellId.maxLevel, r(-1, -1, -0.999999997516473060, -0.999999997516473060)),
      // Just a hair off the outer bounds at different levels.
      (1, 1, 0, r(-1, -1, 1, 1)),
      (1, 1, CellId.maxLevel / 2, r(-1, -1, -0.999918621033430099, -0.999918621033430099)),
      (1, 1, CellId.maxLevel, r(-0.9999999975164731, -0.9999999975164731, -0.9999999950329462, -0.9999999950329462)),
      // Center point of the i,j space at different levels.
      (maxIJ / 2, maxIJ / 2, 0, r(-1, -1, 1, 1)),
      (maxIJ / 2, maxIJ / 2, CellId.maxLevel / 2, r(-0.000040691345930099, -0.000040691345930099, 0, 0)),
      (maxIJ / 2, maxIJ / 2, CellId.maxLevel, r(-0.000000001241763433, -0.000000001241763433, 0, 0)),
      // Maximum i, j at different levels.
      (maxIJ, maxIJ, 0, r(-1, -1, 1, 1)),
      (maxIJ, maxIJ, CellId.maxLevel / 2, r(0.999918621033430099, 0.999918621033430099, 1, 1)),
      (maxIJ, maxIJ, CellId.maxLevel, r(0.999999997516473060, 0.999999997516473060, 1, 1))]
    for (i, j, level, want) in tests {
      let uv = CellId.ijLevelToBoundUV(i: i, j: j, level: level)
      XCTAssertEqual(uv.x.lo, want.x.lo, accuracy: 1e-14)
      XCTAssertEqual(uv.x.hi, want.x.hi, accuracy: 1e-14)
      XCTAssertEqual(uv.y.lo, want.y.lo, accuracy: 1e-14)
      XCTAssertEqual(uv.y.hi, want.y.hi, accuracy: 1e-14)
    }
  }
  
  func testCommonAncestorLevel() {
    let tests = [
      // Identical cell IDs.
      (CellId(face: 0), CellId(face: 0), 0, true),
      (CellId(face: 0).childBegin(30), CellId(face: 0).childBegin(30), 30, true),
      // One cell is a descendant of the other.
      (CellId(face: 0).childBegin(30), CellId(face: 0), 0, true),
      (CellId(face: 5), CellId(face: 5).childEnd(30).prev(), 0, true),
      // No common ancestors.
      (CellId(face: 0), CellId(face: 5), 0, false),
      (CellId(face: 2).childBegin(30), CellId(face: 3).childBegin(20), 0, false),
      // Common ancestor distinct from both.
      (CellId(face: 5).childBegin(9).next().childBegin(15), CellId(face: 5).childBegin(9).childBegin(20), 8, true),
      (CellId(face: 0).childBegin(2).childBegin(30), CellId(face: 0).childBegin(2).next().childBegin(5), 1, true)]
    for (ci, other, want, wantOk) in tests {
      guard let got = ci.commonAncestorLevel(other) else {
        if wantOk {
          XCTFail()
        }
        continue
      }
      XCTAssertEqual(got, want)
    }
  }
  
  func testFindMSBSetNonZero64() {
    var testOne = UInt64(0x800000000000000) << 4
    var testAll = ~UInt64(0)
    var testSome = UInt64(0xFEDCBA987654321) << 4
    for i_ in 0..<63 {
      let i = 63 - i_
      XCTAssertEqual(CellId.findMSBSetNonZero64(testOne), i)
      XCTAssertEqual(CellId.findMSBSetNonZero64(testAll), i)
       XCTAssertEqual(CellId.findMSBSetNonZero64(testSome), i)
      testOne >>= 1
      testAll >>= 1
      testSome >>= 1
    }
  }
  
  func testAdvance() {
    let tests = [
      (CellId(face: 0).childBegin(0), 7, CellId(face: 5).childEnd(0)),
      (CellId(face: 0).childBegin(0), 12, CellId(face: 5).childEnd(0)),
      (CellId(face: 5).childEnd(0), -7, CellId(face: 0).childBegin(0)),
      (CellId(face: 5).childEnd(0), -12000000, CellId(face: 0).childBegin(0)),
      (CellId(face: 0).childBegin(5), 500, CellId(face: 5).childEnd(5).advance(Int64(500 - (6 << (2 * 5))))),
      (CellId(face: 3, pos: 0x12345678, level: CellId.maxLevel-4).childBegin(CellId.maxLevel), 256, CellId(face: 3, pos: 0x12345678, level: CellId.maxLevel-4).next().childBegin(CellId.maxLevel)),
      (CellId(face: 1, pos: 0, level: CellId.maxLevel), 4 << (2 * CellId.maxLevel), CellId(face: 5, pos: 0, level: CellId.maxLevel))]
    for (ci, steps, want) in tests {
      XCTAssertEqual(ci.advance(Int64(steps)), want)
    }
  }

}
