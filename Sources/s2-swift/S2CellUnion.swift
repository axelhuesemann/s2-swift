//
//  S2CellUnion.swift
//  s2-swift
//

import Foundation


/// A collection of CellIDs.
/// It is normalized if it is sorted and does not contain redundancy.
/// Specifically, it may not contain the same CellId twice, nor a CellId that is contained by another,
/// nor the four sibling CellIds that are children of a single higher level CellId.
/// CellUnions are not required to be normalized, but certain operations will
/// return different results if they are not (e.g. Contains).
public struct CellUnion: S2RegionType {

  var cellIds: [CellId]
  
  // MARK: inits 
  
  public init(cellIds: [CellId]) {
    self.cellIds = cellIds
  }
  
  public init(ids: [UInt64]) {
    cellIds = ids.map { CellId(id: $0) }
  }
  
  /// init as union of given CellUnions
  init(union: [CellUnion]) {
    cellIds = []
    for cu in union {
      cellIds.append(contentsOf: cu.cellIds)
    }
    normalize()
  }
  
  /// init as intersection
  // CellUnionFromIntersection creates a CellUnion from the intersection of the given CellUnions.
  init(intersecting x: CellUnion, with y: CellUnion) {
    cellIds = []
    // This is a fairly efficient calculation that uses binary search to skip
    // over sections of both input vectors. It takes constant time if all the
    // cells of x come before or after all the cells of y in CellID order.
    var i = 0
    var j = 0
    while i < x.count && j < y.count {
      let iMin = x[i].rangeMin().id
      let jMin = y[j].rangeMin().id
      if iMin > jMin {
      // either j.Contains(i) or the two cells are disjoint.
        if x[i].id <= y[j].rangeMax().id {
          cellIds.append(x[i])
          i += 1
        } else {
          // Advance j to the first cell possibly contained by x[i].
          j = y.lowerBound(begin: j+1, end: y.count, id: iMin)
          // The previous cell y[j-1] may now contain x[i].
          if x[i].id <= y[j-1].rangeMax().id {
            j -= 1
          }
        }
      } else if jMin > iMin {
        // Identical to the code above with i and j reversed.
        if y[j].id <= x[i].rangeMax().id {
          cellIds.append(y[j])
          j += 1
        } else {
          i = x.lowerBound(begin: i+1, end: x.count, id: jMin)
          if y[j].id <= x[i-1].rangeMax().id {
            i -= 1
          }
        }
      } else {
        // i and j have the same RangeMin(), so one contains the other.
        if x[i].id < y[j].id {
          cellIds.append(x[i])
          i += 1
        } else {
          cellIds.append(y[j])
          j += 1
        }
      }
    }
    self.normalize()
  }
  
  /// Creates a CellUnion from the intersection of a CellUnion with the given CellID.
  /// This can be useful for splitting a CellUnion into chunks.
  init (intersecting x: CellUnion, with id: CellId) {
    cellIds = []
    if x.contains(id) {
      cellIds = [id]
      self.normalize()
      return
    }
    let idMax = id.rangeMax().id
    let begin = x.lowerBound(begin: 0, end: x.count, id: id.rangeMin().id)
    for i in begin..<x.count {
      if x[i].id > idMax { break }
      cellIds.append(x[i])
    }
    normalize()
  }
  
  // MARK: protocols
  
  public subscript(i: Int) -> CellId {
    get {
      return cellIds[i]
    }
    set(newValue) {
      cellIds[i] = newValue
    }
  }
  
  // MARK: arithmetic
  
  public mutating func add(_ cellId: CellId) {
    cellIds.append(cellId)
  }
  
  public func truncate(_ length: Int) -> CellUnion {
    let s = Array(cellIds[0..<length])
    return CellUnion(cellIds: s)
  }
  
  // CellUnionFromDifference creates a CellUnion from the difference (x - y)
  // of the given CellUnions.
  func difference(_ y: CellUnion) -> CellUnion {
    // TODO(roberts): This is approximately O(N*log(N)), but could probably
    // use similar techniques as CellUnionFromIntersectionWithCellID to be more efficient.
    var cu2 = CellUnion(cellIds: cellIds)
    for xid in cellIds {
      cu2.differenceInternal(id: xid, other: y)
    }
    // The output is generated in sorted order, and there should not be any
    // cells that can be merged (provided that both inputs were normalized).
    return cu2
  }

  // cellUnionDifferenceInternal adds the difference between the CellID and the union to
  // the result CellUnion. If they intersect but the difference is non-empty, it divides
  // and conquers.
  mutating func differenceInternal(id: CellId, other: CellUnion) {
    if !other.intersects(id) {
      cellIds.append(id)
      return
    }
    if !other.contains(id) {
      for child in id.children() {
        differenceInternal(id: child, other: other)
      }
    }
  }
  
  /// Creates a CellUnion that covers the half-open range
  /// of leaf cells [begin, end). If begin == end the resulting union is empty.
  /// This requires that begin and end are both leaves, and begin <= end.
  /// To create a closed-ended range, pass in end.Next().
  init(begin: CellId, end: CellId) {
    cellIds = []
    // We repeatedly add the largest cell we can.
    var id = begin.maxTile(limit: end)
    while id != end {
      cellIds.append(id)
      id = id.next().maxTile(limit: end)
    }
    // The output is normalized because the cells are added in order by the iteration.
  }

  /// Normalizes the CellUnion.
  public mutating func normalize() {
    cellIds.sort { $0.id < $1.id }
    var output = [CellId]() // the list of accepted cells
    // Loop invariant: output is a sorted list of cells with no redundancy.
    for ci in cellIds {
      // The first two passes here either ignore this new candidate,
      // or remove previously accepted cells that are covered by this candidate.
      // Ignore this cell if it is contained by the previous one.
      // We only need to check the last accepted cell. The ordering of the
      // cells implies containment (but not the converse), and output has no redundancy,
      // so if this candidate is not contained by the last accepted cell
      // then it cannot be contained by any previously accepted cell.
      if output.count > 0 && output[output.count - 1].contains(ci) {
        continue
      }
      // Discard any previously accepted cells contained by this one.
      // This could be any contiguous trailing subsequence, but it can't be
      // a discontiguous subsequence because of the containment property of
      // sorted S2 cells mentioned above.
      var j = output.count - 1 // last index to keep
      while j >= 0 {
        if !ci.contains(output[j]) {
          break
        }
        j -= 1
      }
      output = Array(output.prefix(j+1))
      // See if the last three cells plus this one can be collapsed.
      // We loop because collapsing three accepted cells and adding a higher level cell
      // could cascade into previously accepted cells.
      var ci_ = ci
      while output.count >= 3 {
        // check if those are siblings
        // if !areSiblings(output[len(output)-3], output[len(output)-2], output[len(output)-1], ci) { break }
        // last three
        let fin = Array(output.suffix(3))
        // fast XOR test; a necessary but not sufficient condition
        if fin[0].id ^ fin[1].id ^ fin[2].id ^ ci_.id != 0 {
          break
        }
        // more expensive test; exact.
        // Compute the two bit mask for the encoded child position,
        // then see if they all agree.
        var mask = ci_.lsb() << 1
        mask = ~(mask + mask << 1)
        let should = ci_.id & mask
        if fin[0].id & mask != should || fin[1].id & mask != should || fin[2].id & mask != should || ci_.isFace() {
          break
        }
        // replace four children by their parent cell.
        output = Array(output.prefix(output.count-3))
        ci_ = ci_.immediateParent() // checked !ci.isFace above
      }
      output.append(ci_)
    }
    cellIds = output
  }

  /// Replaces this CellUnion with an expanded version of the
  /// CellUnion where any cell whose level is less than minLevel or where
  /// (level - minLevel) is not a multiple of levelMod is replaced by its
  /// children, until either both of these conditions are satisfied or the
  /// maximum level is reached.
  public mutating func denormalize(minLevel: Int, levelMod: Int) {
    var denorm = [CellId]()
    for id in cellIds {
      let level = id.level()
      var newLevel = level
      if newLevel < minLevel {
        newLevel = minLevel
      }
      if levelMod > 1 {
        newLevel += (CellId.maxLevel - (newLevel - minLevel)) % levelMod
        if newLevel > CellId.maxLevel {
          newLevel = CellId.maxLevel
        }
      }
      if newLevel == level {
        denorm.append(id)
      } else {
        let end = id.childEnd(newLevel)
        var ci = id.childBegin(newLevel)
        while ci.id != end.id {
          denorm.append(ci)
          ci = ci.next()
        }
      }
    }
    cellIds = denorm
  }

  // MARK: computed members 
  
  var count: Int {
    return cellIds.count
  }
  
  /// Returns an S2Rect that bounds this entity.
  public func rectBound() -> S2Rect {
    var bound = S2Rect.empty
    for cellId in cellIds {
      let c = Cell(id: cellId)
      bound = bound.union(c.rectBound())
    }
    return bound
  }

  public func capBound() -> S2Cap {
    if cellIds.count == 0 {
      return S2Cap.empty
    }
    // Compute the approximate centroid of the region. This won't produce the
    // bounding cap of minimal area, but it should be close enough.
    let zero = R3Vector(x: 0.0, y: 0.0, z: 0.0)
    var centroid3 = zero
    for ci in cellIds {
      let area = S2CellMetric.avgArea.value(ci.level())
      centroid3 = centroid3.add(ci.point().mul(area))
    }
    let centroid: S2Point
    if centroid3 == zero {
      centroid = S2Point(x: 1, y: 0, z: 0)
    } else {
      centroid = S2Point(raw: centroid3.normalized())
    }
    // Use the centroid as the cap axis, and expand the cap angle so that it
    // contains the bounding caps of all the individual cells.  Note that it is
    // *not* sufficient to just bound all the cell vertices because the bounding
    // cap may be concave (i.e. cover more than one hemisphere).
    var c = S2Cap(point: centroid)
    for ci in cellIds {
      let bound = Cell(id: ci).capBound()
      c = c.add(bound)
    }
    return c
  }

  /// Computes a covering of the CellUnion.
//  func cellUnionBound() -> CellUnion {
//    return capBound().cellUnionBound()
//  }
  
  // MARK: tests
  
  /// Reports whether the cell union is valid, meaning that the CellIDs are
  /// valid, non-overlapping, and sorted in increasing order.
  var isValid: Bool {
    for i in 0..<cellIds.count {
      let id = cellIds[i]
      if !id.isValid {
        return false
      }
      if i > 0 && cellIds[i-1].rangeMax().id >= id.rangeMin().id {
        return false
      }
    }
    return true
  }
  
  /// Reports whether the cell union is normalized, meaning that it is
  /// satisfies IsValid and that no four cells have a common parent.
  /// Certain operations such as Contains will return a different
  /// result if the cell union is not normalized.
  var isNormalized: Bool {
    for i in 0..<cellIds.count {
      let cid = cellIds[i]
      if !cid.isValid {
        return false
      }
      if i == 0 {
        continue
      }
      if cellIds[i-1].rangeMax().id >= cid.rangeMin().id {
        return false
      }
      if i < 3 {
        continue
      }
      if areSiblings(a: cellIds[i-3], b: cellIds[i-2], c: cellIds[i-1], d: cid) {
        return false
      }
    }
    return true
  }
  
  // MARK: intersects / contains
  
  /// Reports whether this cell union intersects the given cell ID.
  /// This method assumes that the CellUnion has been normalized.
  func intersects(_ cellId: CellId) -> Bool {
    // Find index of array item that occurs directly after our probe cell:
    var i = cellIds.count
    for (index, ci) in cellIds.enumerated() {
      if cellId.id < ci.id {
        i = index
        break
      }
    }
    if i != cellIds.count && cellIds[i].rangeMin().id <= cellId.rangeMax().id {
      return true
    }
    return i != 0 && cellIds[i-1].rangeMax().id >= cellId.rangeMin().id
  }

  /// Reports whether the cell union contains the given cell ID.
  /// Containment is defined with respect to regions, e.g. a cell contains its 4 children.
  /// This method assumes that the CellUnion has been normalized.
  func contains(_ cellId: CellId) -> Bool {
    // Find index of array item that occurs directly after our probe cell:
    var i = cellIds.count
    for (index, ci) in cellIds.enumerated() {
      if cellId.id < ci.id {
        i = index
        break
      }
    }
    if i != cellIds.count && cellIds[i].rangeMin().id <= cellId.id {
      return true
    }
    return i != 0 && cellIds[i-1].rangeMax().id >= cellId.id
  }
  
  /// Reports whether this cell union contains the given cell.
  public func contains(_ cell: Cell) -> Bool {
    return contains(cell.id)
  }
  
  /// Reports whether this cell union intersects the given cell.
  public func intersects(_ cell: Cell) -> Bool {
    return intersects(cell.id)
  }

  // ContainsPoint reports whether this cell union contains the given point.
  func contains(_ p: S2Point) -> Bool {
    let cell = Cell(point: p)
    return contains(cell)
  }
  
  /// Reports whether this CellUnion contains all of the CellIDs of the given CellUnion.
  func contains(_ o: CellUnion) -> Bool {
    // TODO(roberts): Investigate alternatives such as divide-and-conquer
    // or alternating-skip-search that may be significantly faster in both
    // the average and worst case. This applies to Intersects as well.
    for id in o.cellIds {
      if !contains(id) {
        return false
      }
    }
    return true
  }
  
  /// Reports whether this CellUnion intersects any of the CellIDs of the given CellUnion.
  func intersects(_ o: CellUnion) -> Bool {
    for c in cellIds {
      if o.contains(c) {
        return true
      }
    }
    return false
  }
  
  // MARK: not sure where these go
  
  /// Returns the index in this CellUnion to the first element whose value
  /// is not considered to go before the given cell id. (i.e., either it is equivalent
  /// or comes after the given id.) If there is no match, then end is returned.
  func lowerBound(begin: Int, end: Int, id: UInt64) -> Int {
    for i in begin..<end {
      if self[i].id >= id {
        return i
      }
    }
    return end
  }
  
  /// Returns true if the given four cells have a common parent.
  /// This requires that the four CellIDs are distinct.
  func areSiblings(a: CellId, b: CellId, c: CellId, d: CellId) -> Bool {
    // A necessary (but not sufficient) condition is that the XOR of the
    // four cell IDs must be zero. This is also very fast to test.
    if (a.id ^ b.id ^ c.id) != d.id {
      return false
    }
    // Now we do a slightly more expensive but exact test. First, compute a
    // mask that blocks out the two bits that encode the child position of
    // "id" with respect to its parent, then check that the other three
    // children all agree with "mask".
    let mask0 = UInt64(d.lsb() << 1)
    let mask = ~(mask0 + (mask0 << 1))
    let idMasked = d.id & mask
    return (
      (a.id & mask) == idMasked &&
      (b.id & mask) == idMasked &&
      (c.id & mask) == idMasked &&
      !d.isFace())
  }

  /// Reports the number of leaf cells covered by this cell union.
  /// This will be no more than 6*2^60 for the whole sphere.
  func leafCellsCovered() -> UInt64 {
    var numLeaves: UInt64 = 0
    for c in cellIds {
      numLeaves += 1 << UInt64((CellId.maxLevel - c.level()) << 1)
    }
    return numLeaves
  }
  
  // MARK: areas
  
  /// Returns the average area of this CellUnion.
  /// This is accurate to within a factor of 1.7.
  func averageArea() -> Double {
    return S2CellMetric.avgArea.value(CellId.maxLevel) * Double(leafCellsCovered())
  }
  
  /// Returns the approximate area of this CellUnion. This method is accurate
  /// to within 3% percent for all cell sizes and accurate to within 0.1% for cells
  /// at level 5 or higher within the union.
  func approxArea() -> Double {
    var area = 0.0
    for id in cellIds {
      area += Cell(id: id).approxArea()
    }
    return area
  }
  
  /// Returns the area of this CellUnion as accurately as possible.
  func exactArea() -> Double {
    var area = 0.0
    for id in cellIds {
      area += Cell(id: id).exactArea()
    }
    return area
  }
 
  // MARK: expand
  
  /// Expands this CellUnion by adding a rim of cells at expandLevel
  /// around the unions boundary.
  /// For each cell c in the union, we add all cells at level
  /// expandLevel that abut c. There are typically eight of those
  /// (four edge-abutting and four sharing a vertex). However, if c is
  /// finer than expandLevel, we add all cells abutting
  /// c.Parent(expandLevel) as well as c.Parent(expandLevel) itself,
  /// as an expandLevel cell rarely abuts a smaller cell.
  ///
  /// Note that the size of the output is exponential in
  /// expandLevel. For example, if expandLevel == 20 and the input
  /// has a cell at level 10, there will be on the order of 4000
  /// adjacent cells in the output. For most applications the
  /// ExpandByRadius method below is easier to use.
  mutating func expand(level: Int) {
    var output = [CellId]()
    let levelLsb = CellId.lsb(level)
    for j in 0..<cellIds.count {
      var i = cellIds.count - j
      var id = cellIds[i]
      if id.lsb() < levelLsb {
        id = id.parent(level)
        // Optimization: skip over any cells contained by this one. This is
        // especially important when very small regions are being expanded.
        // TODO: this does not work in Swift loop
        while i > 0 && id.contains(cellIds[i - 1]) {
          i -= 1
        }
      }
      output.append(id)
      for n in id.allNeighbors(level: level) {
        output.append(n)
      }
    }
    cellIds = output
    normalize()
  }
  
  /// Expands this CellUnion such that it contains all points whose
  /// distance to the CellUnion is at most minRadius, but do not use cells that
  /// are more than maxLevelDiff levels higher than the largest cell in the input.
  /// The second parameter controls the tradeoff between accuracy and output size
  /// when a large region is being expanded by a small amount (e.g. expanding Canada
  /// by 1km). For example, if maxLevelDiff == 4 the region will always be expanded
  /// by approximately 1/16 the width of its largest cell. Note that in the worst case,
  /// the number of cells in the output can be up to 4 * (1 + 2 ** maxLevelDiff) times
  /// larger than the number of cells in the input.
  mutating func expandByRadius(minRadius: S1Angle, maxLevelDiff: Int) {
    var minLevel = CellId.maxLevel
    for cid in cellIds {
      minLevel = min(minLevel, cid.level())
    }
    // Find the maximum level such that all cells are at least "minRadius" wide.
    let radiusLevel = S2CellMetric.minWidth.maxLevel(minRadius)
    if radiusLevel == 0 && minRadius > S2CellMetric.minWidth.value(0) {
      // The requested expansion is greater than the width of a face cell.
      // The easiest way to handle this is to expand twice.
      expand(level: 0)
    }
    expand(level: min(minLevel + maxLevelDiff, radiusLevel))
  }
  
}

extension CellUnion: Equatable {

  public static func ==(lhs: CellUnion, rhs: CellUnion) -> Bool {
    return lhs.cellIds == rhs.cellIds
  }

}
