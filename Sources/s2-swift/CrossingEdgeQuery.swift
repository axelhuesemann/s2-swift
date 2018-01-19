//
//  CrossingEdgeQuery.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/11/18.
//

import Foundation


/// Stores a sorted set of edge ids for each shape.
//typealias EdgeMap = [Shape: Int]

/// CrossingEdgeQuery is used to find the Edge IDs of Shapes that are crossed by
/// a given edge(s).
///
/// Note that if you need to query many edges, it is more efficient to declare
/// a single CrossingEdgeQuery instance and reuse it.
///
/// If you want to find *all* the pairs of crossing edges, it is more efficient to
/// use the not yet implemented VisitCrossings in shapeutil.
struct CrossingEdgeQuery {
  var index: ShapeIndex
  // temporary values used while processing a query.
  var a: R2Point?
  var b: R2Point?
  var iter: ShapeIndexIterator
  // candidate cells generated when finding crossings.
  var cells: [ShapeIndexCell] = []
}

extension CrossingEdgeQuery {
 
  /// Creates a CrossingEdgeQuery for the given index.
  init(index: ShapeIndex) {
    self.index = index
    iter = index.iterator()
  }
  
  // Crossings returns the set of edge of the shape S that intersect the given edge AB.
  // If the CrossingType is Interior, then only intersections at a point interior to both
  // edges are reported, while if it is CrossingTypeAll then edges that share a vertex
  // are also reported.
  mutating func crossings(a: S2Point, b: S2Point, shape: Shape, crossType: CrossingType) -> [Int]? {
    var edges = candidates(a: a, b: b, shape: shape)
    if edges.count == 0 {
      return nil
    }
    var crosser = EdgeCrosser(a: a, b: b)
    var out = 0
    let n = edges.count
    for input in 0..<n {
      let b = shape.edge(edges[input])
      let sign = crosser.crossingSign(c: b.v0, d: b.v1)
      if crossType == .all && (sign == .maybeCross || sign == .cross) || crossType != .all && sign == .cross {
        edges[out] = edges[input]
        out += 1
      }
    }
    if out < n {
      return Array(edges[..<out])
    }
    return edges
  }

  // CrossingsEdgeMap returns the set of all edges in the index that intersect the given
  // edge AB. If crossType is CrossingTypeInterior, then only intersections at a
  // point interior to both edges are reported, while if it is CrossingTypeAll
  // then edges that share a vertex are also reported.
  //
  // The edges are returned as a mapping from shape to the edges of that shape
  // that intersect AB. Every returned shape has at least one crossing edge.
//  mutating func crossingsEdgeMap(a: S2Point, b: S2Point, crossType: CrossingType) -> EdgeMap? {
//    var edgeMap = candidatesEdgeMap(a: a, b: b)
//    if edgeMap.count == 0 {
//      return nil
//    }
//    var crosser = EdgeCrosser(a: a, b: b)
//    for (shape, edges) in edgeMap {
//      var out = 0
//      let n = edges.count
//      for input in 0..<n {
//        let edge = shape.edge(edges[input])
//        let sign = crosser.crossingSign(c: edge.v0, d: edge.v1)
//        if (crossType == .all && (sign == .maybeCross || sign == .cross)) || (crossType != .all && sign == .cross) {
//          edgeMap[shape][out] = edges[input]
//          out += 1
//        }
//      }
//      if out == 0 {
//        edgeMap.removeValue(forKey: shape)
//      } else {
//        if out < n {
//          edgeMap[shape] = Array(edgeMap[shape][..<out])
//        }
//      }
//    }
//    return edgeMap
//  }

  // candidates returns a superset of the edges of the given shape that intersect
  // the edge AB.
  mutating func candidates(a: S2Point, b: S2Point, shape: Shape) -> [Int] {
    // For small loops it is faster to use brute force. The threshold below was
    // determined using benchmarks.
    let maxBruteForceEdges = 27
    let maxEdges = shape.numEdges()
    if maxEdges <= maxBruteForceEdges {
      return (0..<maxEdges).map { $0 }
    }
    // Compute the set of index cells intersected by the query edge.
    getCellsForEdge(a: a, b: a)
    if cells.count == 0 {
      return []
    }
    // Gather all the edges that intersect those cells and sort them.
    // TODO(roberts): Shapes don't track their ID, so we need to range over
    // the index to find the ID manually.
    guard let shapeId = index.find(shape: shape) else { return [] }
    var edges: [Int] = []
    for cell in cells {
      guard let clipped = cell.find(shapeId: shapeId) else {
        continue
      }
      for j in clipped.edges {
        edges.append(j)
      }
    }
    return CrossingEdgeQuery.uniqueInts(input: edges)
  }

  // uniqueInts returns the sorted uniqued values from the given input.
  static func uniqueInts(input: [Int]) -> [Int] {
    if input.count <= 1 {
      return input
    }
    let edges = Set<Int>(input)
    return Array(edges)
//    var edges = [Int]()
//    var m = [Int: Bool]()
//    for i in input {
//      if m[i] != nil {
//        continue
//      }
//      m[i] = true
//      edges.append(i)
//    }
//    edges.sort()
//    return edges
  }

  // candidatesEdgeMap returns a map from shapes to the superse of edges for that
  // shape that intersect the edge AB.
  //
  // CAVEAT: This method may return shapes that have an empty set of candidate edges.
  // However the return value is non-empty only if at least one shape has a candidate edge.
//  mutating func candidatesEdgeMap(a: S2Point, b: S2Point) -> EdgeMap {
//    var edgeMap: EdgeMap = [:]
//    // If there are only a few edges then it's faster to use brute force. We
//    // only bother with this optimization when there is a single shape.
//    if let shape = index.shape(id: 0) {
//      // Typically this method is called many times, so it is worth checking
//      // whether the edge map is empty or already consists of a single entry for
//      // this shape, and skip clearing edge map in that case.
//      // Note that we leave the edge map non-empty even if there are no candidates
//      // (i.e., there is a single entry with an empty set of edges).
//      edgeMap[shape] = candidates(a, b, shape)
//      return edgeMap
//    }
//    // Compute the set of index cells intersected by the query edge.
//    getCellsForEdge(a: a, b: b)
//    if cells.count == 0 {
//      return edgeMap
//    }
//    // Gather all the edges that intersect those cells and sort them.
//    for cell in cells {
//      for clipped in cell.shapes {
//        if let s = index.shape(id: clipped.shapeId) {
//          for j in 0..<clipped.numEdges() {
//            edgeMap[s].append(clipped.edges[j])
//          }
//        }
//      }
//    }
//    if cells.count > 1 {
//      for (s, edges) in edgeMap {
//        edgeMap[s] = uniqueInts(edges)
//      }
//    }
//    return edgeMap
//  }
  
  // getCells returns the set of ShapeIndexCells that might contain edges intersecting
  // the edge AB in the given cell root. This method is used primarly by loop and shapeutil.
  mutating func getCells(a: S2Point, b: S2Point, root: PaddedCell) -> [ShapeIndexCell] {
    let (aUV, bUV, intersects) = clipToFace(a: a, b: b, face: root.id.face())
    if intersects {
      self.a = aUV
      self.b = bUV
      let edgeBound = R2Rect(p0: aUV, p1: bUV)
      if root.bound.intersects(edgeBound) {
        computeCellsIntersected(pcell: root, edgeBound: edgeBound)
      }
    }
    return cells
  }

  /// Populates the cells field to the set of index cells intersected by an edge AB.
  mutating func getCellsForEdge(a: S2Point, b: S2Point) {
    cells = []
    let segments = faceSegments(a: a, b: b)
    for segment in segments {
      self.a = segment.a
      self.b = segment.b
      // Optimization: rather than always starting the recursive subdivision at
      // the top level face cell, instead we start at the smallest S2CellId that
      // contains the edge (the edge root cell). This typically lets us skip
      // quite a few levels of recursion since most edges are short.
      let edgeBound = R2Rect(p0: segment.a, p1: segment.b)
      var pcell = PaddedCell(id: CellId(face: segment.face), padding: 0)
      let edgeRoot = pcell.shrinkToFit(rect: edgeBound)
      // Now we need to determine how the edge root cell is related to the cells
      // in the spatial index (cellMap). There are three cases:
      //
      //  1. edgeRoot is an index cell or is contained within an index cell.
      //     In this case we only need to look at the contents of that cell.
      //  2. edgeRoot is subdivided into one or more index cells. In this case
      //     we recursively subdivide to find the cells intersected by AB.
      //  3. edgeRoot does not intersect any index cells. In this case there
      //     is nothing to do.
      let relation = iter.locate(cellId: edgeRoot)
      if relation == .indexed {
        // edgeRoot is an index cell or is contained by an index cell (case 1).
        guard let indexCell = iter.indexCell() else { return }
        cells.append(indexCell)
      } else if relation == .subdivided {
        // edgeRoot is subdivided into one or more index cells (case 2). We
        // find the cells intersected by AB using recursive subdivision.
        if !edgeRoot.isFace() {
          pcell = PaddedCell(id: edgeRoot, padding: 0)
        }
        computeCellsIntersected(pcell: pcell, edgeBound: edgeBound)
      }
    }
  }

  // computeCellsIntersected computes the index cells intersected by the current
  // edge that are descendants of pcell and adds them to this queries set of cells.
  mutating func computeCellsIntersected(pcell: PaddedCell, edgeBound: R2Rect) {
    iter.seek(target: pcell.id.rangeMin())
    if iter.done() || iter.cellId() > pcell.id.rangeMax() {
      // The index does not contain pcell or any of its descendants.
      return
    }
    if iter.cellId() == pcell.id {
      // The index contains this cell exactly.
      cells.append(iter.indexCell()!)
      return
    }
    // Otherwise, split the edge among the four children of pcell.
    let center = pcell.middle.lo
    if edgeBound.x.hi < center.x {
      // Edge is entirely contained in the two left children.
      clipVAxis(edgeBound: edgeBound, center: center.y, i: 0, pcell: pcell)
      return
    } else if edgeBound.x.lo >= center.x {
      // Edge is entirely contained in the two right children.
      clipVAxis(edgeBound: edgeBound, center: center.y, i: 1, pcell: pcell)
      return
    }
    let childBounds = splitUBound(edgeBound: edgeBound, u: center.x)
    if edgeBound.y.hi < center.y {
      // Edge is entirely contained in the two lower children.
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: 0, j: 0), edgeBound: childBounds.0)
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: 1, j: 0), edgeBound: childBounds.1)
    } else if edgeBound.y.lo >= center.y {
      // Edge is entirely contained in the two upper children.
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: 0, j: 1), edgeBound: childBounds.0)
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: 1, j: 1), edgeBound: childBounds.1)
    } else {
      // The edge bound spans all four children. The edge itself intersects
      // at most three children (since no padding is being used).
      clipVAxis(edgeBound: childBounds.0, center: center.y, i: 0, pcell: pcell)
      clipVAxis(edgeBound: childBounds.1, center: center.y, i: 1, pcell: pcell)
    }
  }
  
  // clipVAxis computes the intersected cells recursively for a given padded cell.
  // Given either the left (i=0) or right (i=1) side of a padded cell pcell,
  // determine whether the current edge intersects the lower child, upper child,
  // or both children, and call c.computeCellsIntersected recursively on those children.
  // The center is the v-coordinate at the center of pcell.
  mutating func clipVAxis(edgeBound: R2Rect, center: Double, i: Int, pcell: PaddedCell) {
    if edgeBound.y.hi < center {
      // Edge is entirely contained in the lower child.
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: i, j: 0), edgeBound: edgeBound)
    } else if edgeBound.y.lo >= center {
      // Edge is entirely contained in the upper child.
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: i, j: 1), edgeBound: edgeBound)
    } else {
      // The edge intersects both children.
      let childBounds = splitVBound(edgeBound: edgeBound, v: center)
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: i, j: 0), edgeBound: childBounds.0)
      computeCellsIntersected(pcell: PaddedCell(parent: pcell, i: i, j: 1), edgeBound: childBounds.1)
    }
  }
  
  // splitUBound returns the bound for two children as a result of spliting the
  // current edge at the given value U.
  func splitUBound(edgeBound: R2Rect, u: Double) -> (R2Rect, R2Rect) {
    guard let a = a, let b = b else { return (R2Rect.empty, R2Rect.empty) }
    let value = interpolateFloat64(x: u, a: a.x, b: b.x, a1: a.y, b1: b.y)
    let v = edgeBound.y.clamp(value)
    // diag indicates which diagonal of the bounding box is spanned by AB:
    // it is 0 if AB has positive slope, and 1 if AB has negative slope.
    let diag = (a.x > b.x) != (a.y > b.y)
    return CrossingEdgeQuery.splitBound(edgeBound: edgeBound, uEnd: false, vEnd: diag, u: u, v: v)
  }
  
  // splitVBound returns the bound for two children as a result of spliting the
  // current edge into two child edges at the given value V.
  func splitVBound(edgeBound: R2Rect, v: Double) -> (R2Rect, R2Rect) {
    guard let a = a, let b = b else { return (R2Rect.empty, R2Rect.empty) }
    let value = interpolateFloat64(x: v, a: a.y, b: b.y, a1: a.x, b1: b.x)
    let u = edgeBound.x.clamp(value)
    let diag = (a.x > b.x) != (a.y > b.y)
    return CrossingEdgeQuery.splitBound(edgeBound: edgeBound, uEnd: diag, vEnd: false, u: u, v: v)
  }

  // splitBound returns the bounds for the two childrenn as a result of spliting
  // the current edge into two child edges at the given point (u,v). uEnd and vEnd
  // indicate which bound endpoints of the first child will be updated.
  static func splitBound(edgeBound: R2Rect, uEnd: Bool, vEnd: Bool, u: Double, v: Double) -> (R2Rect, R2Rect) {
    let x0 = R1Interval(lo: uEnd ? u : edgeBound.x.lo, hi: !uEnd ? u : edgeBound.x.hi)
    let x1 = R1Interval(lo: !uEnd ? u : edgeBound.x.lo, hi: uEnd ? u : edgeBound.x.hi)
    let y0 = R1Interval(lo: vEnd ? v : edgeBound.y.lo, hi: !vEnd ? v : edgeBound.y.hi)
    let y1 = R1Interval(lo: !vEnd ? v : edgeBound.y.lo, hi: vEnd ? v : edgeBound.y.hi)
    let bound0 = R2Rect(x: x0, y: y0)
    let bound1 = R2Rect(x: x1, y: y1)
    return (bound0, bound1)
  }

}
