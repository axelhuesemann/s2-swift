//
//  S2RegionCoverer.swift
//  s2-swift
//

import Foundation

// BUG(akashagrawal): The differences from the C++ version FloodFill, SimpleCovering

// RegionCoverer allows arbitrary regions to be approximated as unions of cells (CellUnion).
// This is useful for implementing various sorts of search and precomputation operations.
//
// Typical usage:
//
//	rc := &s2.RegionCoverer{MaxLevel: 30, MaxCells: 5}
//	r := s2.Region(CapFromCenterArea(center, area))
//	covering := rc.Covering(r)
//
// This yields a CellUnion of at most 5 cells that is guaranteed to cover the
// given region (a disc-shaped region on the sphere).
//
// For covering, only cells where (level - MinLevel) is a multiple of LevelMod will be used.
// This effectively allows the branching factor of the S2 CellID hierarchy to be increased.
// Currently the only parameter values allowed are 0/1, 2, or 3, corresponding to
// branching factors of 4, 16, and 64 respectively.
//
// Note the following:
//
//  - MinLevel takes priority over MaxCells, i.e. cells below the given level will
//    never be used even if this causes a large number of cells to be returned.
//
//  - For any setting of MaxCells, up to 6 cells may be returned if that
//    is the minimum number of cells required (e.g. if the region intersects
//    all six face cells).  Up to 3 cells may be returned even for very tiny
//    convex regions if they happen to be located at the intersection of
//    three cube faces.
//
//  - For any setting of MaxCells, an arbitrary number of cells may be
//    returned if MinLevel is too high for the region being approximated.
//
//  - If MaxCells is less than 4, the area of the covering may be
//    arbitrarily large compared to the area of the original region even if
//    the region is convex (e.g. a Cap or Rect).
//
// The approximation algorithm is not optimal but does a pretty good job in
// practice. The output does not always use the maximum number of cells
// allowed, both because this would not always yield a better approximation,
// and because MaxCells is a limit on how much work is done exploring the
// possible covering as well as a limit on the final output size.
//
// Because it is an approximation algorithm, one should not rely on the
// stability of the output. In particular, the output of the covering algorithm
// may change across different versions of the library.
//
// One can also generate interior coverings, which are sets of cells which
// are entirely contained within a region. Interior coverings can be
// empty, even for non-empty regions, if there are no cells that satisfy
// the provided constraints and are contained by the region. Note that for
// performance reasons, it is wise to specify a MaxLevel when computing
// interior coverings - otherwise for regions with small or zero area, the
// algorithm may spend a lot of time subdividing cells all the way to leaf
// level to try to find contained cells.
public struct RegionCoverer {
  
  // MARK parameters
  
  let minLevel: Int // the minimum cell level to be used.
  let maxLevel: Int // the maximum cell level to be used.
  let levelMod: Int // the LevelMod to be used.
  let maxCells: Int // the maximum desired number of cells in the approximation.
  
  // MARK: coverer factory
  
  // newCoverer returns an instance of coverer.
  func newCoverer() -> Coverer {
    return Coverer(
      minLevel: max(0, min(CellId.maxLevel, minLevel)), // clamped
      maxLevel: max(0, min(CellId.maxLevel, maxLevel)), // clamped
      levelMod: max(1, min(3, levelMod)), // clamped
      maxCells: maxCells)
  }
  
  // MARK: region to cell union
  
  // Covering returns a CellUnion that covers the given region and satisfies the various restrictions.
  func covering(region: S2RegionType) -> CellUnion {
    let covering = cellUnion(region: region)
    covering.denormalize(minLevel: max(0, min(maxLevel, minLevel)), levelMod: max(1, min(3, levelMod)))
    return covering
  }
  
  // InteriorCovering returns a CellUnion that is contained within the given region and satisfies the various restrictions.
  func interiorCovering(region: S2RegionType) -> CellUnion {
    let covering = interiorCellUnion(region: region)
    covering.denormalize(minLevel: max(0, min(maxLevel, minLevel)), levelMod: max(1, min(3, levelMod)))
    return covering
  }
  
  // CellUnion returns a normalized CellUnion that covers the given region and
  // satisfies the restrictions except for minLevel and levelMod. These criteria
  // cannot be satisfied using a cell union because cell unions are
  // automatically normalized by replacing four child cells with their parent
  // whenever possible. (Note that the list of cell ids passed to the CellUnion
  // constructor does in fact satisfy all the given restrictions.)
  func cellUnion(region: S2RegionType) -> CellUnion {
    let coverer = newCoverer()
    coverer.coveringInternal(region: region)
    let union = coverer.result
    union.normalize()
    return union
  }
  
  // InteriorCellUnion returns a normalized CellUnion that is contained within the given region and
  // satisfies the restrictions except for minLevel and levelMod. These criteria
  // cannot be satisfied using a cell union because cell unions are
  // automatically normalized by replacing four child cells with their parent
  // whenever possible. (Note that the list of cell ids passed to the CellUnion
  // constructor does in fact satisfy all the given restrictions.)
  func interiorCellUnion(region: S2RegionType) -> CellUnion {
    let coverer = newCoverer()
    coverer.interiorCovering = true
    coverer.coveringInternal(region: region)
    let union = coverer.result
    union.normalize()
    return union
  }
  
  // FastCovering returns a CellUnion that covers the given region similar to Covering,
  // except that this method is much faster and the coverings are not as tight.
  // All of the usual parameters are respected (MaxCells, MinLevel, MaxLevel, and LevelMod),
  // except that the implementation makes no attempt to take advantage of large values of
  // MaxCells.  (A small number of cells will always be returned.)
  //
  // This function is useful as a starting point for algorithms that
  // recursively subdivide cells.
  func fastCovering(cap: S2Cap) -> CellUnion {
    let coverer = newCoverer()
    let union = coverer.rawFastCovering(cap: cap)
    return coverer.normalize(covering: union)
  }

}


// Nodes of a tree that dscribes the result
class Candidate: Comparable {
  
  let cell: Cell
  // Cell should not be expanded further
  var terminal: Bool
  // Number of children that intersect the region
  var numChildren = 0
  // Actual size may be 0, 4, 16, or 64 elements
  var children = [Candidate]()
  // Priority of the candiate
  var priority = 0
  
  // MARK: inits
  
  init(cell: Cell, terminal: Bool) {
    self.cell = cell
    self.terminal = terminal
  }
  
  // MARK: methods
  
  func addChild(_ child: Candidate) {
    children.append(child)
    numChildren += 1
  }
  
}

func ==(lhs: Candidate, rhs: Candidate) -> Bool {
  return lhs.priority == rhs.priority
}

func <(lhs: Candidate, rhs: Candidate) -> Bool {
  return lhs.priority < rhs.priority
}


class Coverer  {
  
  // MARK: parameters
  
  let minLevel: Int // the minimum cell level to be used.
  let maxLevel: Int // the maximum cell level to be used.
  let levelMod: Int // the LevelMod to be used.
  let maxCells: Int // the maximum desired number of cells in the approximation.
  
  // MARK: working variables
  
  var region: S2RegionType?
  var result: CellUnion
  var pq: PriorityQueue<Candidate>
  var interiorCovering: Bool
  
  // MARK: inits 
  
  init(minLevel: Int, maxLevel: Int, levelMod: Int, maxCells: Int) {
    self.minLevel = minLevel
    self.maxLevel = maxLevel
    self.levelMod = levelMod
    self.maxCells = maxCells
    //
    result = CellUnion(ids: [])
    pq = PriorityQueue<Candidate>()
    interiorCovering = false
  }
  
  // MARK: algorithms
  
  // newCandidate returns a new candidate with no children if the cell intersects the given region.
  // The candidate is marked as terminal if it should not be expanded further.
  func newCandidate(cell: Cell) -> Candidate? {
    guard let region = region else { return nil }
    if !region.intersects(cell) {
      return nil
    }
    var terminal = false
    let level = Int(cell.level)
    if level >= minLevel {
      if interiorCovering {
        if region.contains(cell) {
          terminal = true
        } else if level + levelMod > maxLevel {
          return nil
        }
      } else if level + levelMod > maxLevel || region.contains(cell) {
        terminal = true
      }
    }
    return Candidate(cell: cell, terminal: terminal)
  }
  
  // expandChildren populates the children of the candidate by expanding the given number of
  // levels from the given cell.  Returns the number of children that were marked "terminal".
  func expandChildren(candidate: Candidate, cell: Cell, numLevels: Int) -> Int {
    guard let region = region else { return 0 }
    let numLevels = numLevels - 1
    var numTerminals = 0
    let last = cell.id.childEnd()
    var ci = cell.id.childBegin()
    while ci.id != last.id {
      let childCell = Cell(id: ci)
      if numLevels > 0 {
        if region.intersects(childCell) {
          numTerminals += expandChildren(candidate: candidate, cell: childCell, numLevels: numLevels)
        }
        ci = ci.next()
        continue
      }
      let child = newCandidate(cell: childCell)
      if let child = child {
        candidate.addChild(child)
        if child.terminal {
          numTerminals += 1
        }
      }
      ci = ci.next()
    }
    return numTerminals
  }
  
  // addCandidate adds the given candidate to the result if it is marked as "terminal",
  // otherwise expands its children and inserts it into the priority queue.
  // Passing an argument of nil does nothing.
  func add(_ candidate: Candidate) {
    if candidate.terminal {
      result.add(candidate.cell.id)
      return
    }
    // Expand one level at a time until we hit minLevel to ensure that we don't skip over it.
    var numLevels = levelMod
    let level = Int(candidate.cell.level)
    if level < minLevel {
      numLevels = 1
    }
    let numTerminals = expandChildren(candidate: candidate, cell: candidate.cell, numLevels: numLevels)
    let maxChildrenShift = 2 * levelMod
    if candidate.numChildren == 0 {
      return
    } else if !interiorCovering && numTerminals == 1 << maxChildrenShift && level > minLevel {
      // Optimization: add the parent cell rather than all of its children.
      // We can't do this for interior coverings, since the children just
      // intersect the region, but may not be contained by it - we need to
      // subdivide them further.
      candidate.terminal = true
      add(candidate)
    } else {
      // We negate the priority so that smaller absolute priorities are returned
      // first. The heuristic is designed to refine the largest cells first,
      // since those are where we have the largest potential gain. Among cells
      // of the same size, we prefer the cells with the fewest children.
      // Finally, among cells with equal numbers of children we prefer those
      // with the smallest number of children that cannot be refined further.
      candidate.priority = -((level<<maxChildrenShift+candidate.numChildren)<<maxChildrenShift + numTerminals)
      pq.push(candidate)
    }
  }
  
  // adjustLevel returns the reduced "level" so that it satisfies levelMod. Levels smaller than minLevel
  // are not affected (since cells at these levels are eventually expanded).
  func adjustLevel(_ level: Int) -> Int {
    if levelMod > 1 && level > minLevel {
      return level - (level - minLevel) % levelMod
    }
    return level
  }
  
  // adjustCellLevels ensures that all cells with level > minLevel also satisfy levelMod,
  // by replacing them with an ancestor if necessary. Cell levels smaller
  // than minLevel are not modified (see AdjustLevel). The output is
  // then normalized to ensure that no redundant cells are present.
  func adjustCellLevels(_ cells: CellUnion) -> CellUnion {
    if levelMod == 1 {
      return cells
    }
    var out = 0
    for ci in cells.cellIds {
      let level = ci.level()
      let newLevel = adjustLevel(level)
      var ci = ci
      if newLevel != level {
        ci = ci.parent(newLevel)
      }
      if out > 0 && cells.cellIds[out-1].contains(ci) {
        continue
      }
      while out > 0 && ci.contains(cells.cellIds[out-1]) {
        out -= 1
      }
      cells.cellIds[out] = ci
      out += 1
    }
    return cells.truncate(out)
  }
  
  // initialCandidates computes a set of initial candidates that cover the given region.
  func initialCandidates() {
    // Optimization: start with a small (usually 4 cell) covering of the region's bounding cap.
    let temp = RegionCoverer(minLevel: 0, maxLevel: maxLevel, levelMod: 1, maxCells: min(4, maxCells))
    guard let region = region else { return }
    var cells = temp.fastCovering(cap: region.capBound())
    cells = adjustCellLevels(cells)
    for ci in cells.cellIds {
      if let candidate = newCandidate(cell: Cell(id: ci)) {
        add(candidate)
      }
    }
  }
  
  // coveringInternal generates a covering and stores it in result.
  // Strategy: Start with the 6 faces of the cube.  Discard any
  // that do not intersect the shape.  Then repeatedly choose the
  // largest cell that intersects the shape and subdivide it.
  //
  // result contains the cells that will be part of the output, while pq
  // contains cells that we may still subdivide further. Cells that are
  // entirely contained within the region are immediately added to the output,
  // while cells that do not intersect the region are immediately discarded.
  // Therefore pq only contains cells that partially intersect the region.
  // Candidates are prioritized first according to cell size (larger cells
  // first), then by the number of intersecting children they have (fewest
  // children first), and then by the number of fully contained children
  // (fewest children first).
  func coveringInternal(region: S2RegionType) {
    self.region = region
    initialCandidates()
    while pq.count > 0 && (!interiorCovering || result.count < maxCells) {
      guard let candidate = pq.pop() else { continue }
      // For interior covering we keep subdividing no matter how many children
      // candidate has. If we reach MaxCells before expanding all children,
      // we will just use some of them.
      // For exterior covering we cannot do this, because result has to cover the
      // whole region, so all children have to be used.
      // candidate.numChildren == 1 case takes care of the situation when we
      // already have more then MaxCells in result (minLevel is too high).
      // Subdividing of the candidate with one child does no harm in this case.
      if interiorCovering || Int(candidate.cell.level) < minLevel || candidate.numChildren == 1 || result.count + pq.count + candidate.numChildren <= maxCells {
        for child in candidate.children {
          if !interiorCovering || result.count < maxCells {
            add(child)
          }
        }
      } else {
        candidate.terminal = true
        add(candidate)
      }
    }
    pq = PriorityQueue<Candidate>()
    self.region = nil
  }
  
  // rawFastCovering computes a covering of the given cap. In general the covering consists of
  // at most 4 cells (except for very large caps, which may need up to 6 cells).
  // The output is not sorted.
  func rawFastCovering(cap: S2Cap) -> CellUnion {
    let covering = CellUnion(ids: [])
    // Find the maximum level such that the cap contains at most one cell vertex
    // and such that CellId.VertexNeighbors() can be called.
    let level = min(S2CellMetric.minWidth.maxLevel(2 * cap.radius()), CellId.maxLevel-1)
    if level == 0 {
      for face in 0..<6 {
        covering.add(CellId(face: face))
      }
    } else {
      for vn in CellId(point: cap.center).vertexNeighbors(level) {
        covering.add(vn)
      }
    }
    return covering
  }
  
  // normalizeCovering normalizes the "covering" so that it conforms to the current covering
  // parameters (MaxCells, minLevel, maxLevel, and levelMod).
  // This method makes no attempt to be optimal. In particular, if
  // minLevel > 0 or levelMod > 1 then it may return more than the
  // desired number of cells even when this isn't necessary.
  //
  // Note that when the covering parameters have their default values, almost
  // all of the code in this function is skipped.
  func normalize(covering: CellUnion) -> CellUnion {
    // If any cells are too small, or don't satisfy levelMod, then replace them with ancestors.
    if maxLevel < CellId.maxLevel || levelMod > 1 {
      for i in 0..<covering.count {
        let ci = covering[i]
        let level = ci.level()
        let newLevel = adjustLevel(min(level, maxLevel))
        if newLevel != level {
          covering[i] = ci.parent(newLevel)
        }
      }
    }
    // Sort the cells and simplify them.
    covering.normalize()
    // If there are still too many cells, then repeatedly replace two adjacent
    // cells in CellID order by their lowest common ancestor.
    while covering.count > maxCells {
      var bestIndex = -1
      var bestLevel = -1
      for i in 0..<covering.count-1 {
        guard var level = covering[i].commonAncestorLevel(covering[i+1]) else {
          continue
        }
        level = adjustLevel(level)
        if level > bestLevel {
          bestLevel = level
          bestIndex = i
        }
      }
      if bestLevel < minLevel {
        break
      }
      covering[bestIndex] = covering[bestIndex].parent(bestLevel)
      covering.normalize()
    }
    // Make sure that the covering satisfies minLevel and levelMod,
    // possibly at the expense of satisfying MaxCells.
    if minLevel > 0 || levelMod > 1 {
      covering.denormalize(minLevel: minLevel, levelMod: levelMod)
    }
    return covering
  }

}
