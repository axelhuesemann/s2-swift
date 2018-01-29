//
//  PaddedCell.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/11/18.
//

import Foundation


// PaddedCell represents a Cell whose (u,v)-range has been expanded on
// all sides by a given amount of "padding". Unlike Cell, its methods and
// representation are optimized for clipping edges against Cell boundaries
// to determine which cells are intersected by a given set of edges.
struct PaddedCell {
  let id: CellId
  let padding: Double
  let bound: R2Rect
  let middle: R2Rect // A rect in (u, v)-space that belongs to all four children.
  let iLo: Int
  let jLo: Int     // Minimum (i,j)-coordinates of this cell before padding
  let orientation: Int     // Hilbert curve orientation of this cell.
  let level: Int
}

extension PaddedCell {

  /// Constructs a padded cell with the given padding.
  init(id: CellId, padding: Double) {
    self.id = id
    self.padding = padding
    // Fast path for constructing a top-level face (the most common case).
    if id.isFace() {
      let limit = padding + 1
      bound = R2Rect(x: R1Interval(lo: -limit, hi: limit), y: R1Interval(lo: -limit, hi: limit))
      middle = R2Rect(x: R1Interval(lo: -padding, hi: padding), y: R1Interval(lo: -padding, hi: padding))
      iLo = 0
      jLo = 0
      orientation = id.face() & 1
      level = 0 // face
    } else {
      let (_, iLo, jLo, orientation) = id.faceIJOrientation()
      level = id.level()
      let bound = CellId.ijLevelToBoundUV(i: iLo, j: jLo, level: level)
      self.bound = bound // .expanded(margin: padding)
      let ijSize = CellId.sizeIJ(level)
      self.iLo = iLo & -ijSize
      self.jLo = jLo & -ijSize
      let u = S2Cube.stToUV(CellId.siTiToST(si: UInt32(2 * self.iLo + ijSize)))
      let v = S2Cube.stToUV(CellId.siTiToST(si: UInt32(2 * self.jLo + ijSize)))
      middle = R2Rect(x: R1Interval(lo: u - padding, hi: u + padding), y: R1Interval(lo: v - padding, hi: v + padding))
      self.orientation = orientation
    }
  }
  
  /// Constructs the child of parent with the given (i,j) index.
  /// The four child cells have indices of (0,0), (0,1), (1,0), (1,1), where the i and j
  /// indices correspond to increasing u- and v-values respectively.
  init(parent: PaddedCell, i: Int, j: Int) {
    // Compute the position and orientation of the child incrementally from the
    // orientation of the parent.
    let pos = CellId.ijToPos[parent.orientation][2 * i + j]
    id = parent.id.children()[pos]
    padding = parent.padding
    let pbound = parent.bound
    orientation = parent.orientation ^ CellId.posToOrientation[pos]
    level = parent.level + 1
    let ijSize = CellId.sizeIJ(level)
    iLo = parent.iLo + i * ijSize
    jLo = parent.jLo + j * ijSize
    // We could compute middle lazily because it is not needed the majority of the
    // time (i.e., for cells where the recursion terminates).
    let u = S2Cube.stToUV(CellId.siTiToST(si: UInt32(2 * iLo + ijSize)))
    let v = S2Cube.stToUV(CellId.siTiToST(si: UInt32(2 * jLo + ijSize)))
    middle = R2Rect(x: R1Interval(lo: u - padding, hi: u + padding), y: R1Interval(lo: v - padding, hi: v + padding))
    // For each child, one corner of the bound is taken directly from the parent
    // while the diagonally opposite corner is taken from middle().
    let pmiddle = parent.middle
    let x = R1Interval(lo: (i == 1) ? pmiddle.x.lo : pbound.x.lo, hi: (i == 1) ? pbound.x.hi : pmiddle.x.hi)
    let y = R1Interval(lo: (j == 1) ? pmiddle.y.lo : pbound.x.lo, hi: (j == 1) ? pbound.y.hi : pmiddle.y.hi)
    bound = R2Rect(x: x, y: y)
  }
  
  // MARK:
  
  /// Returns the center of this cell.
  func center() -> S2Point {
    let ijSize = CellId.sizeIJ(level)
    let si = uint32(2 * iLo + ijSize)
    let ti = uint32(2 * jLo + ijSize)
    return CellId.faceSiTiToXYZ(face: id.face(), si: si, ti: ti)
  }
  
  /// Returns the (i,j) coordinates for the child cell at the given traversal
  /// position. The traversal position corresponds to the order in which child
  /// cells are visited by the Hilbert curve.
  func childIJ(pos: Int) -> (i: Int, j: Int) {
    let ij = CellId.posToIJ[orientation][pos]
    return (ij >> 1, ij & 1)
  }
  
  /// Return the vertex where the space-filling curve enters this cell.
  func entryVertex() -> S2Point {
    // The curve enters at the (0,0) vertex unless the axis directions are
    // reversed, in which case it enters at the (1,1) vertex.
    var i = iLo
    var j = jLo
    if orientation & CellId.invertMask != 0 {
      let ijSize = CellId.sizeIJ(level)
      i += ijSize
      j += ijSize
    }
    return CellId.faceSiTiToXYZ(face: id.face(), si: UInt32(2 * i), ti: UInt32(2 * j))
  }
  
  /// Returns the vertex where the space-filling curve exits this cell.
  func exitVertex() -> S2Point {
    // The curve exits at the (1,0) vertex unless the axes are swapped or
    // inverted but not both, in which case it exits at the (0,1) vertex.
    var i = iLo
    var j = jLo
    let ijSize = CellId.sizeIJ(level)
    if orientation == 0 || orientation == CellId.swapMask + CellId.invertMask {
      i += ijSize
    } else {
      j += ijSize
    }
    return CellId.faceSiTiToXYZ(face: id.face(), si: UInt32(2 * i), ti: UInt32(2 * j))
  }
  
  /// Returns the smallest CellID that contains all descendants of this
  /// padded cell whose bounds intersect the given rect. For algorithms that use
  /// recursive subdivision to find the cells that intersect a particular object, this
  /// method can be used to skip all of the initial subdivision steps where only
  /// one child needs to be expanded.
  ///
  /// Note that this method is not the same as returning the smallest cell that contains
  /// the intersection of this cell with rect. Because of the padding, even if one child
  /// completely contains rect it is still possible that a neighboring child may also
  /// intersect the given rect.
  ///
  /// The provided Rect must intersect the bounds of this cell.
  func shrinkToFit(rect: R2Rect) -> CellId {
    // Quick rejection test: if rect contains the center of this cell along
    // either axis, then no further shrinking is possible.
    if level == 0 {
      // Fast path (most calls to this function start with a face cell).
      if rect.x.contains(0) || rect.y.contains(0) {
        return id
      }
    }
    let ijSize = CellId.sizeIJ(level)
    let stX = CellId.siTiToST(si: UInt32(2 * iLo + ijSize))
    if rect.x.contains(S2Cube.stToUV(stX)) {
      return id
    }
    let stY = CellId.siTiToST(si: UInt32(2 * jLo + ijSize))
    if rect.y.contains(S2Cube.stToUV(stY)) {
      return id
    }
    // Otherwise we expand rect by the given padding on all sides and find
    // the range of coordinates that it spans along the i- and j-axes. We then
    // compute the highest bit position at which the min and max coordinates
    // differ. This corresponds to the first cell level at which at least two
    // children intersect rect.
    // Increase the padding to compensate for the error in uvToST.
    // (The constant below is a provable upper bound on the additional error.)
    let padded = rect.expanded(padding + 1.5 * S1Interval.dblEpsilon)
    var iMin = iLo
    var jMin = jLo // Min i- or j- coordinate spanned by padded
    var iXor: Int
    var jXor: Int         // XOR of the min and max i- or j-coordinates
    if iMin < CellId.stToIJ(S2Cube.uvToST(padded.x.lo)) {
      iMin = CellId.stToIJ(S2Cube.uvToST(padded.x.lo))
    }
    let a1 = iLo + ijSize - 1
    let b1 = CellId.stToIJ(S2Cube.uvToST(padded.x.hi))
    if a1 <= b1 {
      iXor = iMin ^ a1
    } else {
      iXor = iMin ^ b1
    }
    if jMin < CellId.stToIJ(S2Cube.uvToST(padded.y.lo)) {
      jMin = CellId.stToIJ(S2Cube.uvToST(padded.y.lo))
    }
    let a2 = jLo + ijSize - 1
    let b2 = CellId.stToIJ(S2Cube.uvToST(padded.y.hi))
    if a2 <= b2 {
      jXor = jMin ^ a2
    } else {
      jXor = jMin ^ b2
    }
    // Compute the highest bit position where the two i- or j-endpoints differ,
    // and then choose the cell level that includes both of these endpoints. So
    // if both pairs of endpoints are equal we choose maxLevel; if they differ
    // only at bit 0, we choose (maxLevel - 1), and so on.
    let levelMsb = UInt64(((iXor | jXor) << 1) + 1)
    let pLevel = CellId.maxLevel - Int(CellId.findMSBSetNonZero64(levelMsb))
    if pLevel <= level {
      return id
    }
    return CellId(face: id.face(), i: iMin, j: jMin).parent(pLevel)
  }
  
}
