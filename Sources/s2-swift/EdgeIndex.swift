//
//  EdgeIndex3.swift
//  s2-swiftGo
//

import Foundation


class S2EdgeIndex {

  static let epsilon = 1e-14

  /**
   * Thicken the edge in all directions by roughly 1% of the edge length when
   * thickenEdge is true.
   */
  let Thickening = 0.01
  
  /**
   * Threshold for small angles, that help lenientCrossing to determine whether
   * two edges are likely to intersect.
   */
  let MaxDetError = 1e-14
  
  /**
   * The cell containing each edge, as given in the parallel array
   * <code>edges</code>.
   */
  var _cells: [UInt64]
  
  /**
   * The edge contained by each cell, as given in the parallel array
   * <code>cells</code>.
   */
  var _edges: [Int]
  
  /**
   * No cell strictly below this level appears in mapping. Initially leaf level,
   * that's the minimum level at which we will ever look for test edges.
   */
  
  /**
   * Has the index been computed already?
   */
  var _indexComputed: Bool
  var _minimumS2LevelUsed: Int
  
  /**
   * Number of queries so far
   */
  var _queryCount: Int
  
  // MARK: inits 
  
  init() {
    _cells = []
    _edges = []
    _indexComputed = false
    _minimumS2LevelUsed = 0
    _queryCount  = 0
  }
  
  /**
   * Empties the index in case it already contained something.
   */
  
  func Reset() {
    _minimumS2LevelUsed = CellId.maxLevel
    _indexComputed = false
    _queryCount = 0
    _cells = [] //nil
    _edges = [] //nil
  }
  
  /**
   * Compares [cell1, edge1] to [cell2, edge2], by cell first and edge second.
   *
   * @return -1 if [cell1, edge1] is less than [cell2, edge2], 1 if [cell1,
   *         edge1] is greater than [cell2, edge2], 0 otherwise.
   */
  
  static func Compare(cell1: UInt64, edge1: Int, cell2: UInt64, edge2: Int) -> Int {
    if (cell1 < cell2) {
      return -1
    } else if (cell1 > cell2) {
      return 1
    } else if (edge1 < edge2) {
      return -1
    } else if (edge1 > edge2) {
      return 1
    } else {
      return 0
    }
  }
  
  static func lessThan(cell1: UInt64, edge1: Int, cell2: UInt64, edge2: Int) -> Bool {
    if cell1 < cell2 {
      return true
    } else if cell1 == cell2 && edge1 < edge2 {
      return true
    } else {
      return false
    }
  }
  
  /** Computes the index (if it has not been previously done). */
  func ComputeIndex() {
    if (_indexComputed) {
      return
    }
    var cellList = [UInt64]()
    var edgeList = [Int]()
    for i in 0..<NumEdges() {
      let from = EdgeFrom(index: i)!
      let to = EdgeTo(index: i)!
      let (level, cover) = GetCovering(a: from, b: to, thickenEdge: true)
      _minimumS2LevelUsed = min(_minimumS2LevelUsed, level)
      for cellId in cover {
        cellList.append(cellId.id)
        edgeList.append(i)
      }
    }
    _cells = [UInt64](repeating: 0, count: _cells.count)
    _edges = [Int](repeating: 0, count: _cells.count)
    for i in 0..<_cells.count {
      _cells[i] = cellList[i]
      _edges[i] = edgeList[i]
    }
    SortIndex()
    _indexComputed = true
  }
  
  /** Sorts the parallel <code>cells</code> and <code>edges</code> arrays. */
  func SortIndex() {
    // create an array of indices and sort based on the values in the parallel
    // arrays at each index
    var indices = [Int](repeating: 0, count: _cells.count)
    for i in 0..<indices.count {
      indices[i] = i
    }
    indices.sort { _cells[$0] < _cells[$1] || (_cells[$0] == _cells[$1] && _edges[$0] < _edges[$1]) }
    // copy the cells and edges in the order given by the sorted list of indices
    var newCells = [UInt64](repeating: 0, count: _cells.count)
    var newEdges = [Int](repeating: 0, count: _cells.count)
    for i in 0..<indices.count {
      newCells[i] = _cells[indices[i]]
      newEdges[i] = _edges[indices[i]]
    }
    // replace the cells and edges with the sorted arrays
    _cells = newCells
    _edges = newEdges
  }
  
  var IsIndexComputed: Bool { return _indexComputed }
  
  /**
   * Tell the index that we just received a request for candidates. Useful
   * to compute when to switch to quad tree.
   */
  func IncrementQueryCount() {
    _queryCount += 1
  }
  
  /**
   * If the index hasn't been computed yet, looks at how much work has gone into
   * iterating using the brute force method, and how much more work is planned
   * as defined by 'cost'. If it were to have been cheaper to use a quad tree
   * from the beginning, then compute it now. This guarantees that we will never
   * use more than twice the time we would have used had we known in advance
   * exactly how many edges we would have wanted to test. It is the theoretical
   * best.
   *
   *  The value 'n' is the number of iterators we expect to request from this
   * edge index.
   *
   *  If we have m data edges and n query edges, then the brute force cost is m
   * * n * testCost where testCost is taken to be the cost of
   * EdgeCrosser.robustCrossing, measured to be about 30ns at the time of this
   * writing.
   *
   *  If we compute the index, the cost becomes: m * costInsert + n *
   * costFind(m)
   *
   *  - costInsert can be expected to be reasonably stable, and was measured at
   * 1200ns with the BM_QuadEdgeInsertionCost benchmark.
   *
   *  - costFind depends on the length of the edge . For m=1000 edges, we got
   * timings ranging from 1ms (edge the length of the polygon) to 40ms. The
   * latter is for very long query edges, and needs to be optimized. We will
   * assume for the rest of the discussion that costFind is roughly 3ms.
   *
   *  When doing one additional query, the differential cost is m * testCost -
   * costFind(m) With the numbers above, it is better to use the quad tree (if
   * we have it) if m >= 100.
   *
   *  If m = 100, 30 queries will give m*n*testCost = m * costInsert = 100ms,
   * while the marginal cost to find is 3ms. Thus, this is a reasonable thing to
   * do.
   */
  func PredictAdditionalCalls(n: Int) {
    if (_indexComputed) {
      return
    }
    if (NumEdges() > 100 && (_queryCount + n) > 30) {
      ComputeIndex()
    }
  }
  
  /**
   * Overwrite these functions to give access to the underlying data. The
   * function getNumEdges() returns the number of edges in the index, while
   * edgeFrom(index) and edgeTo(index) return the "from" and "to" endpoints of
   * the edge at the given index.
   */
  func NumEdges() -> Int { return 0 }
  func EdgeFrom(index: Int) ->  S2Point? { return nil }
  func EdgeTo(index: Int) ->  S2Point? { return nil }
  
  /**
   * Appends to "candidateCrossings" all edge references which may cross the
   * given edge. This is done by covering the edge and then finding all
   * references of edges whose coverings overlap this covering. Parent cells are
   * checked level by level. Child cells are checked all at once by taking
   * advantage of the natural ordering of CellIds.
   */
  func FindCandidateCrossings(a: S2Point, b: S2Point) -> [Int] {
    assert(_indexComputed)
    let (_, cover) = GetCovering(a: a, b: b, thickenEdge: false)
    // Edge references are inserted into the map once for each covering cell, so
    // absorb duplicates here
    var uniqueSet = GetEdgesInParentCells(cover: cover)
    // TODO(user): An important optimization for long query
    // edges (Contains queries): keep a bounding cap and clip the query
    // edge to the cap before starting the descent.
    GetEdgesInChildrenCells(a: a, b: b, cover: cover, candidateCrossings: &uniqueSet)
    let candidateCrossings = Array(uniqueSet)
    return candidateCrossings
  }
  
  /**
   * Returns the smallest cell containing all four points, or
   * {@link CellId#sentinel()} iffff they are not all on the same face. The
   * points don't need to be normalized.
   */
  static func ContainingCell(pa: S2Point, pb: S2Point, pc: S2Point, pd: S2Point) ->  CellId? {
    var a = CellId(point: pa)
    var b = CellId(point: pb)
    var c = CellId(point: pc)
    var d = CellId(point: pd)
    if (a.face() != b.face() || a.face() != c.face() || a.face() != d.face()) {
      return nil
    }
    while a != b || a != c || a != d {
      a = a.immediateParent()
      b = b.immediateParent()
      c = c.immediateParent()
      d = d.immediateParent()
    }
    return a
  }
  
  /**
   * Returns the smallest cell containing both points, or Sentinel if they are
   * not all on the same face. The points don't need to be normalized.
   */
  static func ContainingCell(pa: S2Point, pb: S2Point) -> CellId? {
    var a = CellId(point: pa)
    var b = CellId(point: pb)
    if a.face() != b.face() {
      return nil
    }
    while a != b {
      a = a.immediateParent()
      b = b.immediateParent()
    }
    return a
  }
  
  /**
   * Computes a cell covering of an edge. Clears edgeCovering and returns the
   * level of the s2 cells used in the covering (only one level is ever used for
   * each call).
   *
   *  If thickenEdge is true, the edge is thickened and extended by 1% of its
   * length.
   *
   *  It is guaranteed that no child of a covering cell will fully contain the
   * covered edge.
   */
  func GetCovering(a: S2Point, b: S2Point, thickenEdge: Bool) -> (level: Int, covering: [CellId]) {
    // Selects the ideal s2 level at which to cover the edge, this will be the
    // level whose S2 cells have a width roughly commensurate to the length of
    // the edge. We multiply the edge length by 2*THICKENING to guarantee the
    // thickening is honored (it's not a big deal if we honor it when we don't
    // request it) when doing the covering-by-cap trick.
    let edgeLength = a.angle(b)
    let idealLevel = S2CellMetric.minWidth.maxLevel(edgeLength * (1 + 2 * Thickening))
    var containingCellId: CellId?
    if !thickenEdge {
      containingCellId = S2EdgeIndex.ContainingCell(pa: a, pb: b)
    } else {
      if idealLevel == CellId.maxLevel {
        // If the edge is tiny, instabilities are more likely, so we
        // want to limit the number of operations.
        // We pretend we are in a cell much larger so as to trigger the
        // 'needs covering' case, so we won't try to thicken the edge.
        containingCellId = CellId(id: 0xFFF0).parent(3)
      } else {
        let pq = S2Point(raw: b.sub(a).mul(Thickening))
        let ortho = pq.cross(a).mul(edgeLength * Thickening).normalized()
        let p = a.sub(pq)
        let q = b.add(pq)
        // If p and q were antipodal, the edge wouldn't be lengthened,
        // and it could even flip! This is not a problem because
        // idealLevel != 0 here. The farther p and q can be is roughly
        // a quarter Earth away from each other, so we remain
        // Theta(THICKENING).
        containingCellId = S2EdgeIndex.ContainingCell(pa: p.sub(ortho).s2, pb: p.add(ortho).s2, pc: q.sub(ortho).s2, pd: q.add(ortho).s2)
      }
    }
    // Best case: edge is fully contained in a cell that's not too big.
    if let cellId = containingCellId, cellId.level() >= idealLevel - 2 {
      let edgeCovering = [cellId]
      return (cellId.level(), edgeCovering)
    }
    if idealLevel == 0 {
      // Edge is very long, maybe even longer than a face width, so the
      // trick below doesn't work. For now, we will add the whole S2 sphere.
      // TODO(user): Do something a tad smarter (and beware of the
      // antipodal case).
      var edgeCovering = [CellId]()
      for cellid in CellId.iterate(level: 0) {
        edgeCovering.append(cellid)
      }
      return (0, edgeCovering)
    }
    // TODO(user): Check trick below works even when vertex is at
    // interface
    // between three faces.
    // Use trick as in S2PolygonBuilder.PointIndex.findNearbyPoint:
    // Cover the edge by a cap centered at the edge midpoint, then cover
    // the cap by four big-enough cells around the cell vertex closest to the
    // cap center.
    let middle = S2Point(raw: a.add(b).mul(0.5))
    let actualLevel = min(idealLevel, CellId.maxLevel - 1)
    let edgeCovering = CellId(point: middle).vertexNeighbors(actualLevel)
    return (actualLevel, edgeCovering)
  }
  
  /**
   * Filters a list of entries down to the inclusive range defined by the given
   * cells, in <code>O(log N)</code> time.
   *
   * @param cell1 One side of the inclusive query range.
   * @param cell2 The other side of the inclusive query range.
   * @return An array of length 2, containing the start/end indices.
   */
  func GetEdges(cell1: UInt64, cell2: UInt64) -> [Int] {
    // ensure cell1 <= cell2
    let cell1 = min(cell1, cell2)
    let cell2 = max(cell1, cell2)
    // The binary search returns -N-1 to indicate an insertion point at index N,
    // if an exact match cannot be found. Since the edge indices queried for are
    // not valid edge indices, we will always get -N-1, so we immediately
    // convert to N.
    return [
      -1 - BinarySearch(cell: cell1, edge: Int.min),
      -1 - BinarySearch(cell: cell2, edge: Int.max)]
  }
  
  func BinarySearch(cell: UInt64, edge: Int) -> Int {
    var low = 0
    var high = _cells.count - 1
    while low <= high {
      let mid = (low + high) >> 1
      let cmp = S2EdgeIndex.Compare(cell1: _cells[mid], edge1: _edges[mid], cell2: cell, edge2: edge)
      if cmp < 0 {
        low = mid + 1
      } else if cmp > 0 {
        high = mid - 1
      } else {
        return mid
      }
    }
    return -(low + 1)
  }
  
  /**
   * Adds to candidateCrossings all the edges present in any ancestor of any
   * cell of cover, down to minimumS2LevelUsed. The cell->edge map is in the
   * variable mapping.
   */
  func GetEdgesInParentCells(cover: [CellId]) -> Set<Int> {
    // Find all parent cells of covering cells.
    var parentCells = Set<CellId>()
      for coverCell in cover {
      let n = coverCell.level() - 1 - _minimumS2LevelUsed
      for i in 0..<n {
        let parentLevel = coverCell.level() - 1 - i
        let (inserted, _) = parentCells.insert(coverCell.parent(parentLevel))
        if !inserted {
          break // cell is already in => parents are too.
        }
      }
    }
    // Put parent cell edge references into result.
    var candidateCrossings = Set<Int>()
    for parentCell in parentCells {
      var bounds = GetEdges(cell1: parentCell.id, cell2: parentCell.id)
      for i in bounds[0]..<bounds[1] {
        candidateCrossings.insert(_edges[i])
      }
    }
    return candidateCrossings
  }
  
  /**
   * Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
   */
  static func LenientCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
    // assert (S2.isUnitLength(a))
    // assert (S2.isUnitLength(b))
    // assert (S2.isUnitLength(c))
    let acb = a.cross(c).dot(b.v)
    let bda = b.cross(d).dot(a.v)
    if abs(acb) < S2EdgeIndex.epsilon || abs(bda) < S2EdgeIndex.epsilon {
      return true
    }
    if (acb * bda < 0) {
      return false
    }
    let cbd = c.cross(b).dot(d.v)
    let dac = c.cross(a).dot(c.v)
    if abs(cbd) < S2EdgeIndex.epsilon || abs(dac) < S2EdgeIndex.epsilon {
      return true
    }
    return (acb * cbd >= 0) && (acb * dac >= 0)
  }
  
  /**
   * Returns true if the edge and the cell (including boundary) intersect.
   */
  static func EdgeIntersectsCellBoundary(a: S2Point, b: S2Point, cell: Cell) -> Bool {
    let vertices = (0..<4).map { cell.vertex($0) }
    for i in 0..<4 {
      let fromPoint = vertices[i]
      let toPoint = vertices[(i + 1) % 4]
      if LenientCrossing(a: a, b: b, c: fromPoint, d: toPoint) {
        return true
      }
    }
    return false
  }
  
  /**
   * Appends to candidateCrossings the edges that are fully contained in an S2
   * covering of edge. The covering of edge used is initially cover, but is
   * refined to eliminate quickly subcells that contain many edges but do not
   * intersect with edge.
   */
  func GetEdgesInChildrenCells(a: S2Point, b: S2Point, cover: [CellId], candidateCrossings: inout Set<Int>) {
    // Put all edge references of (covering cells + descendant cells) into
    // result.
    // This relies on the natural ordering of CellIds.
    var cover = cover
    while cover.last != nil {
      let cell = cover.removeLast()
      var bounds = GetEdges(cell1: cell.rangeMin().id, cell2: cell.rangeMax().id)
      if (bounds[1] - bounds[0] <= 16) {
        for i in bounds[0]..<bounds[1] {
          candidateCrossings.insert(_edges[i])
        }
      } else {
        // Add cells at this level
        bounds = GetEdges(cell1: cell.id, cell2: cell.id)
        for i in bounds[0]..<bounds[1] {
          candidateCrossings.insert(_edges[i])
        }
        // Recurse on the children -- hopefully some will be empty.
        for child in cell.children() {
          // TODO(user): Do the check for the four cells at once,
          // as it is enough to check the four edges between the cells. At
          // this time, we are checking 16 edges, 4 times too many.
          //
          // Note that given the guarantee of AppendCovering, it is enough
          // to check that the edge intersect with the cell boundary as it
          // cannot be fully contained in a cell.
          let cell = Cell(id: child)
          if S2EdgeIndex.EdgeIntersectsCellBoundary(a: a, b: b, cell: cell) {
            cover.append(child)
          }
        }
      }
    }
  }
  
  /// <summary>
  /// An iterator on data edges that may cross a query edge (a,b). Create the
  /// iterator, call getCandidates(), then enumerating.
  ///
  /// The current edge in the iteration has index, goes between from()
  /// and to().
  /// </summary>
  /// <returns></returns>
  func GetIterator() -> DataEdgeIterator {
    return DataEdgeIterator(edgeIndex: self)
  }
  
}

/*
 * An iterator on data edges that may cross a query edge (a,b). Create the
 * iterator, call getCandidates(), then hasNext()/next() repeatedly.
 *
 * The current edge in the iteration has index index(), goes between from()
 * and to().
 */

class DataEdgeIterator {
  /**
   * The structure containing the data edges.
   */
  var _candidates: [Int]
  var _edgeIndex: S2EdgeIndex
  
  /**
   * Tells whether getCandidates() obtained the candidates through brute force
   * iteration or using the quad tree structure.
   */
  
  /**
   * Index of the current edge and of the edge before the last next() call.
   */
  var _currentIndex: Int = 0
  
  /**
   * Cache of edgeIndex.getNumEdges() so that hasNext() doesn't make an extra
   * call
   */
  
  /**
   * Index within array above. We have: currentIndex =
   * candidates.get(currentIndexInCandidates).
   */
  var _currentIndexInCandidates: Int = 0
  var _isBruteForce: Bool = false
  var _numEdges: Int = 0
  
  init(edgeIndex: S2EdgeIndex) {
    self._edgeIndex = edgeIndex
    _candidates = []
  }
  
  /**
   * Initializes the iterator to iterate over a set of candidates that may
   * cross the edge (a,b).
   */
  func GetCandidates(a: S2Point, b: S2Point) {
    _edgeIndex.PredictAdditionalCalls(n: 1)
    _isBruteForce = !_edgeIndex.IsIndexComputed
    if (_isBruteForce) {
      _edgeIndex.IncrementQueryCount()
      _currentIndex = 0
      _numEdges = _edgeIndex.NumEdges()
    } else {
      _candidates = _edgeIndex.FindCandidateCrossings(a: a, b: b)
      _currentIndexInCandidates = 0
      if let first = _candidates.first {
        _currentIndex = first
      }
    }
  }
  
  func GetEnumerator() -> AnyIterator<Int> {
    return AnyIterator {
      if self._isBruteForce && self._currentIndex < self._numEdges { return nil }
      if !self._isBruteForce && self._currentIndexInCandidates < self._candidates.count { return nil }
      let mark = self._currentIndex
      if self._isBruteForce {
        self._currentIndex += 1
      } else {
        self._currentIndexInCandidates += 1
        if self._currentIndexInCandidates < self._candidates.count {
          self._currentIndex = self._candidates[self._currentIndexInCandidates]
        }
      }
      return mark
    }
  }
  
}
