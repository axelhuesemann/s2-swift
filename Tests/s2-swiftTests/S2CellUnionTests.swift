//
//  CellUnionTests.swift
//  s2-swift
//

import XCTest
@testable import s2_swift


class S2CellUnionTests: XCTestCase {
  
  override func setUp() {
    super.setUp()
    // call this to initialize the lookup tables
    CellId.setup()
  }
  
  override func tearDown() {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    super.tearDown()
  }
  
  func testNormalization() {
    var cu = CellUnion(ids: [
      0x80855c0000000000, // A: a cell over Pittsburg CA
      0x80855d0000000000, // B, a child of A
      0x8085634000000000, // first child of X, disjoint from A
      0x808563c000000000, // second child of X
      0x80855dc000000000, // a child of B
      0x808562c000000000, // third child of X
      0x8085624000000000, // fourth child of X
      0x80855d0000000000]) // B again
    let exp = CellUnion(ids: [
      0x80855c0000000000, // A
      0x8085630000000000]) // X
    cu.normalize()
    XCTAssertEqual(cu, exp)

    // add a redundant cell
    /* TODO
    cu.Add(0x808562c000000000)
    if !reflect.DeepEqual(cu, exp) {
      t.Errorf("after redundant add, got %v, want %v", cu, exp)
    }
    */
  }

  func testCellUnionBasic() {
    var empty = CellUnion(ids: [])
    empty.normalize()
    XCTAssertEqual(empty.count, 0)
    let face1Id = CellId(face: 1)
    let face1Cell = Cell(id: face1Id)
    var face1Union = CellUnion(cellIds: [face1Id])
    face1Union.normalize()
    XCTAssertEqual(face1Union.count, 1)
    XCTAssertEqual(face1Id, face1Union[0])
    XCTAssertTrue(face1Union.contains(face1Cell))
    let face2Id = CellId(face: 2)
    let face2Cell = Cell(id: face2Id)
    var face2Union = CellUnion(cellIds: [face2Id])
    face2Union.normalize()
    XCTAssertEqual(face2Union.count, 1)
    XCTAssertEqual(face2Id, face2Union[0])
    XCTAssertFalse(face1Union.contains(face2Cell))
  }

  func testCellUnion() {
    let tests: [([UInt64], [UInt64], [UInt64], [UInt64])] = [
			// Single cell around NYC, and some simple nearby probes
      ([
        0x89c25c0000000000],
      [
        CellId(id: 0x89c25c0000000000).childBegin().id,
        CellId(id: 0x89c25c0000000000).childBegin(28).id],
			[
        CellId(id: 0x89c25c0000000000).immediateParent().id,
        CellId(face: CellId(id: 0x89c25c0000000000).face()).id], // the whole face
			[
        CellId(id: 0x89c25c0000000000).next().id,                       // Cell next to this one at same level
				CellId(id: 0x89c25c0000000000).next().childBegin(28).id, // Cell next to this one at deep level
				0x89c2700000000000,                                      // Big(er) neighbor cell
				0x89e9000000000000,                                      // Very big next door cell.
				0x89c1000000000000]),                                     // Very big cell, smaller value than probe
      // NYC and SFO:
			([
				0x89c25b0000000000, // NYC
				0x89c2590000000000, // NYC
				0x89c2f70000000000, // NYC
				0x89c2f50000000000, // NYC
				0x8085870000000000, // SFO
				0x8085810000000000, // SFO
				0x808f7d0000000000, // SFO
				0x808f7f0000000000], // SFO
			[
				0x808f7ef300000000, // SFO
				0x808f7e5cf0000000, // SFO
				0x808587f000000000, // SFO
				0x89c25ac000000000, // NYC
				0x89c259a400000000, // NYC
				0x89c258fa10000000, // NYC
				0x89c258f174007000], // NYC
			[
				0x808c000000000000, // Big SFO
				0x89c4000000000000], // Big NYC
			[
				0x89c15a4fcb1bb000, // outside NYC
				0x89c15a4e4aa95000, // outside NYC
				0x8094000000000000, // outside SFO (big)
				0x8096f10000000000, // outside SFO (smaller)

				0x87c0000000000000]), // Midwest very big
      ([
        // CellUnion with cells at many levels:
				0x8100000000000000, // starting around california
				0x8740000000000000, // adjacent cells at increasing
				0x8790000000000000, // levels, moving eastward.
				0x87f4000000000000,
				0x87f9000000000000, // going down across the midwest
				0x87ff400000000000,
				0x87ff900000000000,
				0x87fff40000000000,
				0x87fff90000000000,
				0x87ffff4000000000,
				0x87ffff9000000000,
				0x87fffff400000000,
				0x87fffff900000000,
				0x87ffffff40000000,
				0x87ffffff90000000,
				0x87fffffff4000000,
				0x87fffffff9000000,
				0x87ffffffff400000], // to a very small cell in Wisconsin
			[
				0x808f400000000000,
				0x80eb118b00000000,
				0x8136a7a11d000000,
				0x8136a7a11dac0000,
				0x876c7c0000000000,
				0x87f96d0000000000,
				0x87ffffffff400000],
			[
        CellId(id: 0x8100000000000000).immediateParent().id,
        CellId(id: 0x8740000000000000).immediateParent().id],
			[
				0x52aaaaaaab300000,
				0x52aaaaaaacd00000,
				0x87fffffffa100000,
				0x87ffffffed500000,
				0x87ffffffa0100000,
				0x87fffffed5540000,
				0x87fffffed6240000,
				0x52aaaacccb340000,
				0x87a0000400000000,
				0x87a000001f000000,
				0x87a0000029d00000,
				0x9500000000000000])]
    for (cells, contained, overlaps, disjoint) in tests {
      var union = CellUnion(ids: cells)
      union.normalize()
      // Ensure self-containment tests are correct.
      for id in cells {
        XCTAssertTrue(union.intersects(CellId(id: id)))
        XCTAssertTrue(union.contains(CellId(id: id)))
      }
      // Test for containment specified in test case.
      for id in contained {
        XCTAssertTrue(union.intersects(CellId(id: id)))
        XCTAssertTrue(union.contains(CellId(id: id)))
      }
      // Make sure the CellUnion intersect these cells but do not contain.
      for id in overlaps {
        XCTAssertTrue(union.intersects(CellId(id: id)))
        XCTAssertFalse(union.contains(CellId(id: id)))
      }
      // Negative cases make sure the CellUnion neither contain nor intersect these cells
      for id in disjoint {
        XCTAssertFalse(union.intersects(CellId(id: id)))
        XCTAssertFalse(union.contains(CellId(id: id)))
      }
    }
  }

  func addCells(id: CellId, selected: Bool, input: inout [CellId], expected: inout [CellId]) {
    // Decides whether to add "id" and/or some of its descendants to the test case.  If "selected"
    // is true, then the region covered by "id" *must* be added to the test case (either by adding
    // "id" itself, or some combination of its descendants, or both).  If cell ids are to the test
    // case "input", then the corresponding expected result after simplification is added to
    // "expected".
    if id.id == 0 {
      // Initial call: decide whether to add cell(s) from each face.
      for face in 0..<6 {
        addCells(id: CellId(face: face), selected: false, input: &input, expected: &expected)
      }
      return
    }
    var selected = selected
    if id.isLeaf() {
      // The oneIn() call below ensures that the parent of a leaf cell will always be selected (if
      // we make it that far down the hierarchy).
      XCTAssertTrue(selected)
      input.append(id)
      return
    }
    // The following code ensures that the probability of selecting a cell at each level is
    // approximately the same, i.e. we test normalization of cells at all levels.
    if !selected && oneIn(n: CellId.maxLevel-id.level()) {
      //  Once a cell has been selected, the expected output is predetermined.  We then make sure
      //  that cells are selected that will normalize to the desired output.
      expected.append(id)
      selected = true
    }
    // With the rnd.OneIn() constants below, this function adds an average
    // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
    // level at which the cell was first selected (level 15 on average).
    // Therefore the average number of input cells in a test case is about
    // (5/6 * 15 * 6) = 75.  The average number of output cells is about 6.
    // If a cell is selected, we add it to "input" with probability 5/6.
    var added = false
    if selected && !oneIn(n: 6) {
      input.append(id)
      added = true
    }
    var numChildren = 0
    var child = id.childBegin()
    while child != id.childEnd() {
      // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
      // This intentionally may result in a cell and some of its children
      // being included in the test case.
      //
      // If the cell is not selected, on average we recurse on one child.
      // We also make sure that we do not recurse on all 4 children, since
      // then we might include all 4 children in the input case by accident
      // (in which case the expected output would not be correct).
      let recurse = oneIn(n: selected ? 12 : 4)
      if recurse && numChildren < 3 {
        addCells(id: child, selected: selected, input: &input, expected: &expected)
        numChildren += 1
      }
      // If this cell was selected but the cell itself was not added, we
      // must ensure that all 4 children (or some combination of their
      // descendants) are added.
      if selected && !added {
        addCells(id: child, selected: selected, input: &input, expected: &expected)
      }
      child = child.next()
    }
  }

  func testNormalizePseudoRandom() {
    // Try a bunch of random test cases, and keep track of average statistics for normalization (to
    // see if they agree with the analysis above).
    var inSum = 0
    var outSum = 0
    let iters = 20
    //
    for _ in 0..<iters {
      var input = [CellId]()
      var expected = [CellId]()
      addCells(id: CellId(id: 0), selected: false, input: &input, expected: &expected)
      inSum += input.count
      outSum += expected.count
      var cellunion = CellUnion(cellIds: input)
      cellunion.normalize()
      //
      XCTAssertEqual(expected.count, cellunion.count)
      for j in input {
        XCTAssertTrue(cellunion.contains(j))
        XCTAssertTrue(cellunion.intersects(j))
        if !j.isFace() {
          XCTAssertTrue(cellunion.intersects(j.immediateParent()))
//            if j.Level() > 1 {
          
//              if cellunion.intersects(j.immediateParent().immediateParent()) == false)
//              if cellunion.intersects(j.Parent(0)) == false
//            }
        }
        //
        if !j.isLeaf() {
          XCTAssertTrue(cellunion.contains(j.childBegin()))
          XCTAssertTrue(cellunion.intersects(j.childBegin()))
          XCTAssertTrue(cellunion.contains(j.childEnd().prev()))
          XCTAssertTrue(cellunion.intersects(j.childEnd().prev()))
          XCTAssertTrue(cellunion.contains(j.childBegin(CellId.maxLevel)))
          XCTAssertTrue(cellunion.intersects(j.childBegin(CellId.maxLevel)))
        }
      }
      //
      for exp in expected {
        if !exp.isFace() {
          XCTAssertFalse(cellunion.contains(exp.parent(exp.level() - 1)))
          XCTAssertFalse(cellunion.contains(exp.parent(0)))
        }
      }
      //
      var test = [CellId]()
      var dummy = [CellId]()
      addCells(id: CellId(id: 0), selected: false, input: &test, expected: &dummy)
      for j in test {
        var intersects = false
        var contains = false
        for k in expected {
          if k.contains(j) {
            contains = true
          }
          if k.intersects(j) {
            intersects = true
          }
        }
        XCTAssertEqual(cellunion.contains(j), contains)
        XCTAssertEqual(cellunion.intersects(j), intersects)
      }
    }
  }

  func testDenormalized() {
    let tests = [
      ("not expanded, level mod == 1", 10, 1,
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(11),
          CellId(face: 2).childBegin(11),
          CellId(face: 3).childBegin(14),
          CellId(face: 0).childBegin(10)]),
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(11),
          CellId(face: 2).childBegin(11),
          CellId(face: 3).childBegin(14),
          CellId(face: 0).childBegin(10)])),
      ("not expanded, level mod > 1", 10, 2,
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(12),
          CellId(face: 2).childBegin(12),
          CellId(face: 3).childBegin(14),
          CellId(face: 0).childBegin(10)]),
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(12),
          CellId(face: 2).childBegin(12),
          CellId(face: 3).childBegin(14),
          CellId(face: 0).childBegin(10)])),
      ("expended, (level - min_level) is not multiple of level mod", 10, 3,
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(12),
          CellId(face: 5).childBegin(11)]),
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(12).children()[0],
          CellId(face: 2).childBegin(12).children()[1],
          CellId(face: 2).childBegin(12).children()[2],
          CellId(face: 2).childBegin(12).children()[3],
          CellId(face: 5).childBegin(11).children()[0].children()[0],
          CellId(face: 5).childBegin(11).children()[0].children()[1],
          CellId(face: 5).childBegin(11).children()[0].children()[2],
          CellId(face: 5).childBegin(11).children()[0].children()[3],
          CellId(face: 5).childBegin(11).children()[1].children()[0],
          CellId(face: 5).childBegin(11).children()[1].children()[1],
          CellId(face: 5).childBegin(11).children()[1].children()[2],
          CellId(face: 5).childBegin(11).children()[1].children()[3],
          CellId(face: 5).childBegin(11).children()[2].children()[0],
          CellId(face: 5).childBegin(11).children()[2].children()[1],
          CellId(face: 5).childBegin(11).children()[2].children()[2],
          CellId(face: 5).childBegin(11).children()[2].children()[3],
          CellId(face: 5).childBegin(11).children()[3].children()[0],
          CellId(face: 5).childBegin(11).children()[3].children()[1],
          CellId(face: 5).childBegin(11).children()[3].children()[2],
          CellId(face: 5).childBegin(11).children()[3].children()[3]])),
      ("expended, level < min_level", 10, 3,
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(9)]),
        CellUnion(cellIds: [
          CellId(face: 2).childBegin(9).children()[0],
          CellId(face: 2).childBegin(9).children()[1],
          CellId(face: 2).childBegin(9).children()[2],
          CellId(face: 2).childBegin(9).children()[3]]))]
      for (_, minL, lMod, cu, exp) in tests {
        var cu = cu
        cu.denormalize(minLevel: minL, levelMod: lMod)
        XCTAssertEqual(cu, exp)
      }
    }

  func testCellUnionRectBound() {
    let tests = [
      (CellUnion(ids: []), S2Rect.empty),
      (CellUnion(cellIds: [CellId(face: 1)]), S2Rect(
        lat: R1Interval(lo: -.pi / 4, hi: .pi / 4),
        lng: S1Interval(lo: .pi / 4, hi: 3 * .pi / 4))),
      (CellUnion(ids: [0x808c000000000000]), S2Rect( // Big SFO
        lat: R1Interval(lo: 34.644220547108482 * toRadians, hi: 38.011928357226651 * toRadians),
        lng: S1Interval(lo: -124.508522987668428 * toRadians, hi: -121.628309835221216 * toRadians))),
      (CellUnion(ids: [0x89c4000000000000]),S2Rect( // Big NYC
          lat: R1Interval(lo: 38.794595155857657 * toRadians, hi: 41.747046884651063 * toRadians),
          lng: S1Interval(lo: -76.456308667788633 * toRadians, hi: -73.465162142654819 * toRadians))),
      (CellUnion(ids: [0x89c4000000000000, 0x808c000000000000]), S2Rect( // Big NYC, Big SFO
          lat: R1Interval(lo: 34.644220547108482 * toRadians, hi: 41.747046884651063 * toRadians),
          lng: S1Interval(lo: -124.508522987668428 * toRadians, hi: -73.465162142654819 * toRadians)))]
      for (cu, want) in tests {
        XCTAssert(rectsApproxEqual(a: cu.rectBound(), b: want, tolLat: epsilon, tolLng: epsilon))
      }
    }

}
