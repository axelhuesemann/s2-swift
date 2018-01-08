//
//  S2Polygon.swift
//  s2-swift
//

import Foundation

// Polygon represents a sequence of zero or more loops; recall that the
// interior of a loop is defined to be its left-hand side (see Loop).
//
// When the polygon is initialized, the given loops are automatically converted
// into a canonical form consisting of "shells" and "holes".  Shells and holes
// are both oriented CCW, and are nested hierarchically.  The loops are
// reordered to correspond to a preorder traversal of the nesting hierarchy.
//
// Polygons may represent any region of the sphere with a polygonal boundary,
// including the entire sphere (known as the "full" polygon).  The full polygon
// consists of a single full loop (see Loop), whereas the empty polygon has no
// loops at all.
//
// Use FullPolygon() to construct a full polygon. The zero value of Polygon is
// treated as the empty polygon.
//
// Polygons have the following restrictions:
//
//  - Loops may not cross, i.e. the boundary of a loop may not intersect
//    both the interior and exterior of any other loop.
//
//  - Loops may not share edges, i.e. if a loop contains an edge AB, then
//    no other loop may contain AB or BA.
//
//  - Loops may share vertices, however no vertex may appear twice in a
//    single loop (see Loop).
//
//  - No loop may be empty.  The full loop may appear only in the full polygon.
public struct S2Polygon: S2RegionType {
  
  let loops: [S2Loop]
  
  // index is a spatial index of all the polygon loops.
  let index: ShapeIndex
  
  // hasHoles tracks if this polygon has at least one hole.
  let hasHoles: Bool
  
  // numVertices keeps the running total of all of the vertices of the contained loops.
  let numVertices: Int
  
  // bound is a conservative bound on all points contained by this loop.
  // If l.ContainsPoint(P), then l.bound.ContainsPoint(P).
  let bound: S2Rect
  
  // Since bound is not exact, it is possible that a loop A contains
  // another loop B whose bounds are slightly larger. subregionBound
  // has been expanded sufficiently to account for this error, i.e.
  // if A.Contains(B), then A.subregionBound.Contains(B.bound).
  let subregionBound: S2Rect
  
  // MARK: inits / factory
  
  private init(loops: [S2Loop], index: ShapeIndex, hasHoles: Bool, numVertices: Int, bound: S2Rect, subregionBound: S2Rect) {
    self.loops = loops
    self.index = index
    self.hasHoles = hasHoles
    self.numVertices = numVertices
    self.bound = bound
    self.subregionBound = subregionBound
  }
  
  // PolygonFromLoops constructs a polygon from the given hierarchically nested
  // loops. The polygon interior consists of the points contained by an odd
  // number of loops. (Recall that a loop contains the set of points on its
  // left-hand side.)
  //
  // This method figures out the loop nesting hierarchy and assigns every loop a
  // depth. Shells have even depths, and holes have odd depths.
  //
  // NOTE: this function is NOT YET IMPLEMENTED for more than one loop and will
  // panic if given a slice of length > 1.
//  init(loops: [S2Loop]) {
//    if loops.count > 1 {
//      fatalError("multiple loops are not yet implemented")
//    }
//    // TODO: Once multi-loop is supported, fix this.
//    self.init(loops: loops, index: ShapeIndex(), hasHoles: false, numVertices: loops[0].vertices.count, bound: S2Rect.empty, subregionBound: S2Rect.empty)
//  }

  init(loop: S2Loop) {
    self.init(loops: [loop], index: ShapeIndex(), hasHoles: false, numVertices: loop.vertices.count, bound: loop.bound, subregionBound: loop.subregionBound)
  }

  // FullPolygon returns a special "full" polygon.
  static let full = S2Polygon(loops: [S2Loop.full], index: ShapeIndex(), hasHoles: false, numVertices: S2Loop.full.vertices.count, bound:S2Rect.full, subregionBound: S2Rect.full)
  static let empty = S2Polygon(loops: [], index: ShapeIndex(), hasHoles: false, numVertices: 0, bound:S2Rect.empty, subregionBound: S2Rect.empty)
  
  // MARK: tests
  
  // IsEmpty reports whether this is the special "empty" polygon (consisting of no loops).
  var isEmpty: Bool {
    return loops.count == 0
  }

  // IsFull reports whether this is the special "full" polygon (consisting of a
  // single loop that encompasses the entire sphere).
  var isFull: Bool {
    return loops.count == 1 && loops[0].isFull
  }

  // CapBound returns a bounding spherical cap.
  public func capBound() -> S2Cap {
    return bound.capBound()
  }

  // RectBound returns a bounding latitude-longitude rectangle.
  public func rectBound() -> S2Rect {
    return bound
  }

  // ContainsCell reports whether the polygon contains the given cell.
  public func contains(_ cell: Cell) -> Bool {
    return loops.count == 1 && loops[0].contains(cell)
  }
  
  // IntersectsCell reports whether the polygon intersects the given cell.
  public func intersects(_ cell: Cell) -> Bool {
    return loops.count == 1 && loops[0].intersects(cell)
  }

}
