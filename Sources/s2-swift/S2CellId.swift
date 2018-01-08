//
//  S2CellId.swift
//  s2-swift
//

import Foundation


prefix func -(id: UInt64) -> UInt64 {
  return UInt64(bitPattern: -Int64(bitPattern: id))
}

prefix func -(id: UInt32) -> UInt32 {
  return UInt32(bitPattern: -Int32(bitPattern: id))
}

/// CellId uniquely identifies a cell in the S2 cell decomposition.
/// The most significant 3 bits encode the face number (0-5). The
/// remaining 61 bits encode the position of the center of this cell
/// along the Hilbert curve on that face. The zero value and the value
/// (1<<64)-1 are invalid cell IDs. The first compares less than any
/// valid cell ID, the second as greater than any valid cell ID.
/// The major differences from the C++ version is that barely anything is implemented.
public struct CellId {

  // the world is divided into 6 faces and each face into a quadtree-like structure with 30 levels
  static let numFaces = 6
  public static let maxLevel = 30
  // we use the upper 3 bits for the face, and 61 bits to encode the quadtree position
  static let faceBits = 3
  static let posBits = 2 * maxLevel + 1
  // since we use hilbert, we can wrap
  static let wrapOffset = UInt64(numFaces) << UInt64(posBits)
  static let maxSize = 1 << maxLevel

  // Constants related to the bit mangling in the Cell ID.
  static let lookupBits = 4
  static let swapMask = 0x01
  static let invertMask = 0x02
  
  //
  let id: UInt64

  // MARK: inits / factory
  
  public init(id: UInt64) {
    self.id = id
  }
  
  /// Returns a cell given its face in the range
  /// [0,5], the 61-bit Hilbert curve position pos within that face, and
  /// the level in the range [0,maxLevel]. The position in the cell ID
  /// will be truncated to correspond to the Hilbert curve position at
  /// the center of the returned cell.
  public init(face: Int, pos: UInt64, level: Int) {
    let id = UInt64(face) << UInt64(CellId.posBits) + pos | 1
    self.init(id: id)
    self = self.parent(level)
  }

  /// Returns the cell corresponding to a given S2 cube face.
  public init(face: Int) {
    let id = UInt64(face) << UInt64(CellId.posBits) + CellId.lsb(0)
    self.init(id: id)
  }

  /// Returns the leaf cell containing ll.
  public init(latLng: LatLng) {
    let point = latLng.toPoint()
    self.init(point: point)
  }

  public init(token: String) {
    // add trailing zeros
    let token2 = token.padding(toLength: 16, withPad: "0", startingAt: 0)
    // deal with failure
    guard let id = UInt64(token2, radix: 16) else {
      self.init(id: 0)
      return
    }
    self.init(id: id)
  }
  
  // MARK: string token
  
  func toToken() -> String {
    // pad leading zeros
    var pad = ""
    for i in 0..<16 {
      if id & (UInt64(0xf) << UInt64(60 - i * 4)) == 0 {
        pad += "0"
      } else {
        break
      }
    }
    //
    if id == 0 {
      return "X"
    }
    // strip trailing zeros
    var id2 = id
    while id2 & 0xf == 0 {
      id2 >>= 4
    }
    // hex encoded string
    let s = String(format: "%llx", id2)
    // return hex encoded string
    return pad + s
  }
  
  // MARK: tests
  
  /// Reports whether ci represents a valid cell.
  var isValid: Bool {
    return face() < CellId.numFaces && (lsb() & 0x1555555555555555 != 0)
  }

  /// Returns whether this cell ID is at the deepest level;
  /// that is, the level at which the cells are smallest.
  func isLeaf() -> Bool {
    return id & 1 != 0
  }
  
  /// Returns whether this is a top-level (face) cell.
  func isFace() -> Bool {
    return id & (CellId.lsb(0) - 1) == 0
  }
  
  // MARK: computed members
  
  /// Returns the cube face for this cell ID, in the range [0,5].
  func face() -> Int {
    return Int(id >> UInt64(CellId.posBits))
  }

  /// Returns the position along the Hilbert curve of this cell ID, in the range [0,2^posBits-1].
  func pos() -> UInt64 {
    return id & (UInt64.max >> UInt64(CellId.faceBits))
  }

  /// Returns the subdivision level of this cell ID, in the range [0, maxLevel].
  func level() -> Int {
    // Fast path for leaf cells.
    if isLeaf() {
      return CellId.maxLevel
    }
    var x = UInt32(id & 0xffffffff)
    var level = -1
    if x != 0 {
      level += 16
    } else {
      x = UInt32(id >> 32)
    }
    // Only need to look at even-numbered bits for valid cell IDs.
    x &= -x // remove all but the LSB.
    if x&0x00005555 != 0 {
      level += 8
    }
    if x&0x00550055 != 0 {
      level += 4
    }
    if x&0x05050505 != 0 {
      level += 2
    }
    if x&0x11111111 != 0 {
      level += 1
    }
    return level
  }

  /// Returns the least significant bit that is set.
  func lsb() -> UInt64 {
    // this is the tricky one from the job interviews :-)
    return id & -id
  }
  
  /// Returns the lowest-numbered bit that is on for cells at the given level.
  static func lsb(_ level: Int) -> UInt64 {
    return 1 << UInt64(2 * (CellId.maxLevel - level))
  }
  
  // MARK: edge and vertex neighbors
  
  static func sizeIJ(_ level: Int) -> Int {
    return 1 << Int(CellId.maxLevel - level)
  }

  /// Returns the four cells that are adjacent across the cell's four edges.
  /// Edges 0, 1, 2, 3 are in the down, right, up, left directions in the face space.
  /// All neighbors are guaranteed to be distinct.
  func edgeNeighbors() -> [CellId] {
    let level = self.level()
    let size = CellId.sizeIJ(level)
    let (f, i, j, _) = faceIJOrientation()
    return [
      CellId(face: f, i: i, j: j-size, wrapped: true).parent(level),
      CellId(face: f, i: i+size, j: j, wrapped: true).parent(level),
      CellId(face: f, i: i, j: j+size, wrapped: true).parent(level),
      CellId(face: f, i: i-size, j: j, wrapped: true).parent(level)]
  }

  /// Returns the neighboring cellIds with vertex closest to this cell at the given level.
  /// Normally there are four neighbors, but the closest vertex may only have three neighbors
  /// if it is one of the 8 cube vertices.
  func vertexNeighbors(_ level: Int) -> [CellId] {
    let halfSize = CellId.sizeIJ(level + 1)
    let size = halfSize << 1
    let (f, i, j, _) = faceIJOrientation()
    //
    var isame: Bool
    var jsame: Bool
    var ioffset: Int
    var joffset: Int
    if i & halfSize != 0 {
      ioffset = size
      isame = (i + size) < CellId.maxSize
    } else {
      ioffset = -size
      isame = (i - size) >= 0
    }
    if j & halfSize != 0 {
      joffset = size
      jsame = (j + size) < CellId.maxSize
    } else {
      joffset = -size
      jsame = (j - size) >= 0
    }
    //
    var results = [
      parent(level),
      CellId(face: f, i: i+ioffset, j: j, sameFace: isame).parent(level),
      CellId(face: f, i: i, j: j+joffset, sameFace: jsame).parent(level)]
    //
    if isame || jsame {
      results.append(CellId(face: f, i: i+ioffset, j: j+joffset, sameFace: isame && jsame).parent(level))
    }
    //
    return results
  }

  // MARK: contain / intersect
  
  /// Returns the minimum CellId that is contained within this cell.
  public func rangeMin() -> CellId {
    return CellId(id: id - (lsb() - 1))
  }

  /// Returns the maximum CellId that is contained within this cell.
  public func rangeMax() -> CellId {
    return CellId(id: id + (lsb() - 1))
  }

  /// Returns true iff the CellId contains oci.
  public func contains(_ cellId: CellId) -> Bool {
    return rangeMin().id <= cellId.id && cellId.id <= rangeMax().id
  }

  /// Returns true iff the CellId intersects oci.
  public func intersects(_ cellId: CellId) -> Bool {
    return cellId.rangeMin().id <= rangeMax().id && cellId.rangeMax().id >= rangeMin().id
  }

  // MARK: conversions
  
  /// Returns the center of the s2 cell on the sphere as a Point.
  public func point() -> S2Point {
    return S2Point(raw: rawPoint())
  }

  /// Returns the center of the s2 cell on the sphere as a LatLng.
  public func latLng() -> LatLng {
    let point = self.point()
    return LatLng(point: point)
  }

  // MARK: cell hierarchy
  
  // Returns the cell at the given level, which must be no greater than the current level.
  public func parent(_ level: Int) -> CellId {
    let lsb = CellId.lsb(level)
    return CellId(id: (id & -lsb) | lsb)
  }
  
  /// Is cheaper than Parent, but assumes !isFace().
  public func immediateParent() -> CellId {
    assert(!isFace())
    let nlsb = lsb() << 2
    return CellId(id: (id & -nlsb) | nlsb)
  }
  
  /// Returns the child position (0..3) of this cell's
  /// ancestor at the given level, relative to its parent.  The argument
  /// should be in the range 1..kMaxLevel.  For example,
  /// ChildPosition(1) returns the position of this cell's level-1
  /// ancestor within its top-level face cell.
  public func childPosition(_ level: Int) -> Int {
    return Int(id >> UInt64(2 * (CellId.maxLevel-level) + 1)) & 3
  }
  
  // MARK: enumerate
  
  /// Returns the four immediate children of this cell.
  /// If ci is a leaf cell, it returns four identical cells that are not the children.
  public func children() -> [CellId] {
  //    var ch = [CellId]()
  //    var lsb = self.lsb()
  //    ch.append(CellId(id: id - lsb + lsb>>2))
  //    lsb >>= 1
  //    ch.append(CellId(id: ch[0].id + lsb))
  //    ch.append(CellId(id: ch[1].id + lsb))
  //    ch.append(CellId(id: ch[2].id + lsb))
  //    return ch
    let lsb = self.lsb()
    let lsb1 = lsb >> 1
    // move lsb marker
    // ....1000...
    // ....0010...
    // and increment
    // ....0110...
    // ....1010...
    // ....1110...
    let id0 = id &- lsb &+ lsb >> 2
    let id1 = id0 + lsb1
    let id2 = id1 + lsb1
    let id3 = id2 + lsb1
    return [CellId(id: id0), CellId(id: id1), CellId(id: id2), CellId(id: id3)]
  }
  
  // iterating the whole world needs to deal with the 6 faces
  static func iterate(level: Int) -> AnyIterator<CellId> {
    let levelLsb = CellId.lsb(level)
    let startCell = CellId(face: 0)
    let endCell = CellId(face: 5)
    var cellId = CellId(id: startCell.id - startCell.lsb() + levelLsb)
    let end = CellId(id: endCell.id + endCell.lsb() + levelLsb)
    return AnyIterator<CellId> {
      if cellId == end { return nil }
      let value = cellId
      let id = cellId.id
      cellId = CellId(id: id + (id & -id) << 1)
      return value
    }
  }
  
  func iterate(level: Int) -> AnyIterator<CellId> {
    let levelLsb = CellId.lsb(level)
    var cellId = CellId(id: id - lsb() + levelLsb)
    let end = CellId(id: id + lsb() + levelLsb)
    return AnyIterator<CellId> {
      if cellId == end { return nil }
      let value = cellId
      let id = cellId.id
      cellId = CellId(id: id + (id & -id) << 1)
      return value
    }
  }
  
  func iterate() -> AnyIterator<CellId> {
    let ol = lsb()
    var cellId = CellId(id: id - ol + ol>>2)
    let end = CellId(id: id + ol + ol>>2)
    return AnyIterator<CellId> {
      if cellId == end { return nil }
      let n = cellId
      let id = cellId.id
      cellId = CellId(id: id + (id & -id)<<1)
      return n
    }
  }
  
  // ChildBegin returns the first child in a traversal of the children of this cell, in Hilbert curve order.
  //
  //    for ci := c.ChildBegin(); ci != c.ChildEnd(); ci = Next() {
  //        ...
  //    }
  func childBegin() -> CellId {
    let ol = lsb()
    return CellId(id: id - ol + ol>>2)
  }

  // ChildEnd returns the first cell after a traversal of the children of this cell in Hilbert curve order.
  // The returned cell may be invalid.
  func childEnd() -> CellId {
    let ol = lsb()
    return CellId(id: id + ol + ol>>2)
  }

  // ChildBeginAtLevel returns the first cell in a traversal of children a given level deeper than this cell, in
  // Hilbert curve order. The given level must be no smaller than the cell's level.
  func childBegin(_ level: Int) -> CellId {
    return CellId(id: id - lsb() + CellId.lsb(level))
  }
  
  // ChildEndAtLevel returns the first cell after the last child in a traversal of children a given level deeper
  // than this cell, in Hilbert curve order.
  // The given level must be no smaller than the cell's level.
  // The returned cell may be invalid.
  func childEnd(_ level: Int) -> CellId {
    return CellId(id: id + lsb() + CellId.lsb(level))
  }

  /// Returns the next cell along the Hilbert curve.
  /// This is expected to be used with ChildStart and ChildEnd.
  func next() -> CellId {
    return CellId(id: id + lsb()<<1)
  }

  /// Returns the previous cell along the Hilbert curve.
  func prev() -> CellId {
    return CellId(id: id - lsb() << 1)
  }
  
  func children(level: Int?) -> AnyIterator<CellId> {
    let lsbParent = lsb()
    let lsbChild = (level != nil) ? CellId.lsb(level!) : lsbParent >> 2
    var current = CellId(id: id - lsbParent + lsbChild)
    let end = CellId(id: id + lsbParent + lsbChild)
    return AnyIterator {
      if current == end { return nil }
      let id = current.id
      let n = current
      current = CellId(id: id + (id & -id) << 1)
      return n
    }
  }

//  func children(level: Int) -> CellIdSequence {
//    return CellIdSequence(parent: self, level: level)
//  }
  
  /// Returns an unnormalized r3 vector from the origin through the center
  /// of the s2 cell on the sphere.
  func rawPoint() -> R3Vector {
    let (face, si, ti) = faceSiTi()
    let u = S2Cube.stToUV((0.5 / Double(CellId.maxSize)) * Double(si))
    let v = S2Cube.stToUV((0.5 / Double(CellId.maxSize)) * Double(ti))
    return S2Cube(face: face, u: u, v: v).vector()
  }

  // faceSiTi returns the Face/Si/Ti coordinates of the center of the cell.
  func faceSiTi() -> (Int, Int, Int) {
    let (face, i, j, _) = faceIJOrientation()
    // patching with delta values
    var delta = 0
    if isLeaf() {
      delta = 1
    } else {
      if (i^(Int(id>>2)))&1 != 0 {
        delta = 2
      }
    }
    //
    return (face, 2*i + delta, 2*j + delta)
  }

  /// Uses the global lookupIJ table to unfiddle the bits of ci.
  func faceIJOrientation() -> (Int, Int, Int, Int) {
    let f = face()
    let (i, j, orientation) = CellId.ijOrientation(face: f, id: id, lsb: lsb())
    return (f, i, j, orientation)
  }
  
  /// Uses the global lookupIJ table to unfiddle the bits of ci.
  static func ijOrientation(face: Int, id: UInt64, lsb: UInt64) -> (Int, Int, Int) {
    var i = 0
    var j = 0
    var bits = face & CellId.swapMask
    var nbits = CellId.maxLevel - 7 * CellId.lookupBits // first iteration
    for k_ in 0...7 {
      let k = 7 - k_
      let add = Int(id >> UInt64(k * 2 * CellId.lookupBits + 1)) & ((1 << (2 * nbits)) - 1)
      bits += add << 2
      bits = CellId.lookupIJ[bits]
      i += (bits >> (CellId.lookupBits + 2)) << (k * CellId.lookupBits)
      j += ((bits >> 2) & ((1 << CellId.lookupBits) - 1)) << (k * CellId.lookupBits)
      bits &= (CellId.swapMask | CellId.invertMask)
      nbits = CellId.lookupBits // following iterations
    }
    if lsb & 0x1111111111111110 != 0 {
      bits ^= CellId.swapMask
    }
    return (i, j, bits)
  }
  
  /// Contructs a leaf cell given its cube face (range 0..5) and IJ coordinates.
  init(face: Int, i: Int, j: Int) {
    // Note that this value gets shifted one bit to the left at the end
    // of the function.
    var n = face << (CellId.posBits - 1)
    // Alternating faces have opposite Hilbert curve orientations; this
    // is necessary in order for all faces to have a right-handed
    // coordinate system.
    var bits = face & CellId.swapMask
    // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
    // curve position.  The lookup table transforms a 10-bit key of the form
    // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
    // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
    // Hilbert curve orientation respectively.
    for k_ in 0...7 {
      let k = 7 - k_
      let mask = (1 << CellId.lookupBits) - 1
      bits += ((i >> (k * CellId.lookupBits)) & mask) << (CellId.lookupBits + 2)
      bits += ((j >> (k * CellId.lookupBits)) & mask) << 2
      bits = CellId.lookupPos[bits]
      n |= (bits >> 2) << (k * 2 * CellId.lookupBits)
      bits &= CellId.swapMask | CellId.invertMask
    }
    self.init(id: UInt64(n) * 2 + 1)
  }
  
  /// Returns a leaf cell given its cube face (range 0..5) and IJ coordinates.
  static func idFrom(_ face: Int, i: Int, j: Int) -> UInt64 {
    // Note that this value gets shifted one bit to the left at the end
    // of the function.
    var n = face << (CellId.posBits - 1)
    // Alternating faces have opposite Hilbert curve orientations; this
    // is necessary in order for all faces to have a right-handed
    // coordinate system.
    var bits = face & CellId.swapMask
    // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
    // curve position.  The lookup table transforms a 10-bit key of the form
    // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
    // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
    // Hilbert curve orientation respectively.
    for k_ in 0...7 {
      let k = 7 - k_
      let mask = (1 << CellId.lookupBits) - 1
      bits += ((i >> (k * CellId.lookupBits)) & mask) << (CellId.lookupBits + 2)
      bits += ((j >> (k * CellId.lookupBits)) & mask) << 2
      bits = CellId.lookupPos[bits]
      n |= (bits >> 2) << (k * 2 * CellId.lookupBits)
      bits &= CellId.swapMask | CellId.invertMask
    }
    return UInt64(n) * 2 + 1
  }
  
  init(face: Int, i: Int, j: Int, wrapped: Bool) {
    // Convert i and j to the coordinates of a leaf cell just beyond the
    // boundary of this face.  This prevents 32-bit overflow in the case
    // of finding the neighbors of a face cell.
    let i_ = clamp(i, min: -1, max: CellId.maxSize)
    let j_ = clamp(j, min: -1, max: CellId.maxSize)
    // We want to wrap these coordinates onto the appropriate adjacent face.
    // The easiest way to do this is to convert the (i,j) coordinates to (x,y,z)
    // (which yields a point outside the normal face boundary), and then call
    // xyzToFaceUV to project back onto the correct face.
    //
    // The code below converts (i,j) to (si,ti), and then (si,ti) to (u,v) using
    // the linear projection (u=2*s-1 and v=2*t-1).  (The code further below
    // converts back using the inverse projection, s=0.5*(u+1) and t=0.5*(v+1).
    // Any projection would work here, so we use the simplest.)  We also clamp
    // the (u,v) coordinates so that the point is barely outside the
    // [-1,1]x[-1,1] face rectangle, since otherwise the reprojection step
    // (which divides by the new z coordinate) might change the other
    // coordinates enough so that we end up in the wrong leaf cell.
    let scale = 1.0 / Double(CellId.maxSize)
    let limit = nextafter(1.0, 2.0)
    let u = max(-limit, min(limit, scale * Double((i_ << 1) + 1 - CellId.maxSize)))
    let v = max(-limit, min(limit, scale * Double((j_ << 1) + 1 - CellId.maxSize)))
    // Find the leaf cell coordinates on the adjacent face, and convert
    // them to a cell id at the appropriate level.
    let raw = S2Cube(face: face, u: u, v: v).vector()
    let cube = S2Cube(point: S2Point(raw: raw))
    let i__ = CellId.stToIJ(0.5 * (cube.u+1))
    let j__ = CellId.stToIJ(0.5 * (cube.v+1))
    //
    self.init(face: cube.face, i: i__, j: j__)
  }
  
  init(face: Int, i: Int, j: Int, sameFace: Bool) {
    if sameFace {
      self.init(face: face, i: i, j: j)
    } else {
      self.init(face: face, i: i, j: j, wrapped: true)
    }
  }

  /// Converts the i- or j-index of a leaf cell to the minimum corresponding
  // s- or t-value contained by that cell. The argument must be in the range
  // [0..2**30], i.e. up to one position beyond the normal range of valid leaf
  // cell indices.
  static func ijToSTMin(_ i: Int) -> Double {
    return Double(i) / Double(CellId.maxSize)
  }

  /// Converts value in ST coordinates to a value in IJ coordinates.
  static func stToIJ(_ s: Double) -> Int {
    let ij = Int(floor(Double(CellId.maxSize) * s))
    return clamp(ij, min: 0, max: CellId.maxSize - 1)
  }

  /// Returns a leaf cell containing point p. Usually there is
  /// exactly one such cell, but for points along the edge of a cell, any
  /// adjacent cell may be (deterministically) chosen. This is because
  /// s2.CellIds are considered to be closed sets. The returned cell will
  /// always contain the given point, i.e. CellFromPoint(p).ContainsPoint(p)
  /// is always true.
  init(point: S2Point) {
    let cube = S2Cube(point: point)
    let i = CellId.stToIJ(S2Cube.uvToST(cube.u))
    let j = CellId.stToIJ(S2Cube.uvToST(cube.v))
    self.init(face: cube.face, i: i, j: j)
  }

  /// Returns the bounds in (u,v)-space for the cell at the given
  /// level containing the leaf cell with the given (i,j)-coordinates.
  static func ijLevelToBoundUV(i: Int, j: Int, level: Int) -> R2Rect {
    let cellSize = sizeIJ(level)
    let xLo = i & -cellSize
    let yLo = j & -cellSize
    let x = R1Interval(lo: S2Cube.stToUV(CellId.ijToSTMin(xLo)), hi: S2Cube.stToUV(CellId.ijToSTMin(xLo + cellSize)))
    let y = R1Interval(lo: S2Cube.stToUV(CellId.ijToSTMin(yLo)), hi: S2Cube.stToUV(CellId.ijToSTMin(yLo + cellSize)))
    return R2Rect(x: x, y: y)
  }

  // MARK: Hilbert curve 
  
  static let posToIJ = [
    [0, 1, 3, 2], // canonical order:    (0,0), (0,1), (1,1), (1,0)
    [0, 2, 3, 1], // axes swapped:       (0,0), (1,0), (1,1), (0,1)
    [3, 2, 0, 1], // bits inverted:      (1,1), (1,0), (0,0), (0,1)
    [3, 1, 0, 2]] // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
  static let posToOrientation = [CellId.swapMask, 0, 0, invertMask | swapMask]
  static var lookupIJ = [Int](repeating: 0, count: 1 << (2*CellId.lookupBits + 2))
  static var lookupPos = [Int](repeating: 0, count: 1 << (2*CellId.lookupBits + 2))

  // TODO call this once
  static func setup() {
    if lookupIJ[0] != 0 {
      return
    }
    initLookupCell(level: 0, i: 0, j: 0, origOrientation: 0, pos: 0, orientation: 0)
    initLookupCell(level: 0, i: 0, j: 0, origOrientation: swapMask, pos: 0, orientation: swapMask)
    initLookupCell(level: 0, i: 0, j: 0, origOrientation: invertMask, pos: 0, orientation: invertMask)
    initLookupCell(level: 0, i: 0, j: 0, origOrientation: swapMask|invertMask, pos: 0, orientation: swapMask|invertMask)
  }

  /// Initializes the lookupIJ table at init time.
  static func initLookupCell(
    level: Int, i: Int, j: Int, origOrientation: Int, pos: Int, orientation: Int) {
    if level == lookupBits {
      let ij = (i << lookupBits) + j
      lookupPos[(ij<<2)+origOrientation] = (pos << 2) + orientation
      lookupIJ[(pos<<2)+origOrientation] = (ij << 2) + orientation
      return
    }
    let level = level + 1
    let i = i << 1
    let j = j << 1
    let pos = pos << 2
    let r = posToIJ[orientation]
    initLookupCell(level: level, i: i+(r[0]>>1), j: j+(r[0]&1), origOrientation: origOrientation, pos: pos, orientation: orientation^posToOrientation[0])
    initLookupCell(level: level, i: i+(r[1]>>1), j: j+(r[1]&1), origOrientation: origOrientation, pos: pos+1, orientation: orientation^posToOrientation[1])
    initLookupCell(level: level, i: i+(r[2]>>1), j: j+(r[2]&1), origOrientation: origOrientation, pos: pos+2, orientation: orientation^posToOrientation[2])
    initLookupCell(level: level, i: i+(r[3]>>1), j: j+(r[3]&1), origOrientation: origOrientation, pos: pos+3, orientation: orientation^posToOrientation[3])
  }

  /// Returns the level of the common ancestor of the two S2 CellIds.
  func commonAncestorLevel(_ cellId: CellId) -> Int? {
    var bits = UInt64(id ^ cellId.id)
    if bits < lsb() {
      bits = lsb()
    }
    if bits < cellId.lsb() {
      bits = cellId.lsb()
    }
    let msbPos = CellId.findMSBSetNonZero64(bits)
    if msbPos > 60 {
      return nil
    }
    return (60 - msbPos) >> 1
  }

  /// Returns the index (between 0 and 63) of the most
  /// significant set bit.
  static func findMSBSetNonZero64(_ bits: UInt64) -> Int {
    let val: [UInt64] = [0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000, 0xFFFFFFFF00000000]
    let shift: [UInt64] = [1, 2, 4, 8, 16, 32]
    var msbPos = UInt64(0)
    var bits_ = bits
    for i_ in 0...5 {
      let i = 5 - i_
      if bits_ & val[i] != 0 {
        bits_ >>= shift[i]
        msbPos |= shift[i]
      }
    }
    return Int(msbPos)
  }

  /// Advances or retreats the indicated number of steps along the
  /// Hilbert curve at the current level, and returns the new position. The
  /// position is never advanced past End() or before Begin().
  func advance(_ steps: Int64) -> CellId {
    if steps == 0 {
      return self
    }
    var steps_ = steps
    // We clamp the number of steps if necessary to ensure that we do not
    // advance past the End() or before the Begin() of this level. Note that
    // minSteps and maxSteps always fit in a signed 64-bit integer.
    let stepShift = UInt64(2 * (CellId.maxLevel - level()) + 1)
    if steps_ < 0 {
      let minSteps = -Int64(id >> stepShift)
      if steps_ < minSteps {
        steps_ = minSteps
      }
    } else {
      let maxSteps = Int64((CellId.wrapOffset + lsb() - id) >> stepShift)
      if steps_ > maxSteps {
        steps_ = maxSteps
      }
    }
    //
    let s = UInt64(bitPattern: steps_ << Int64(stepShift))
    return CellId(id: id &+ s)
  }

}

extension CellId: Equatable, CustomStringConvertible, Hashable {
  
  public static func ==(lhs:CellId, rhs: CellId) -> Bool {
    return lhs.id == rhs.id
  }
  
  public var description: String {
    guard isValid else {
      return "Invalid: " + String(id, radix: 16)
    }
    let pos = (1...level()).map { "\(childPosition($0))" }
    // looks like "1/3210"
    return "\(face())/\(pos.joined())"
  }
  
  public var hashValue: Int {
    return id.hashValue
  }
  
}
