//
//  S2Loop.swift
//  s2-swift
//

import Foundation


/// Represents a simple spherical polygon. It consists of a sequence
/// of vertices where the first vertex is implicitly connected to the
/// last. All loops are defined to have a CCW orientation, i.e. the interior of
/// the loop is on the left side of the edges. This implies that a clockwise
/// loop enclosing a small area is interpreted to be a CCW loop enclosing a
/// very large area.
///
/// Loops are not allowed to have any duplicate vertices (whether adjacent or
/// not), and non-adjacent edges are not allowed to intersect. Loops must have
/// at least 3 vertices (except for the "empty" and "full" loops discussed
/// below).
///
/// There are two special loops: the "empty" loop contains no points and the
/// "full" loop contains all points. These loops do not have any edges, but to
/// preserve the invariant that every loop can be represented as a vertex
/// chain, they are defined as having exactly one vertex each (see EmptyLoop
/// and FullLoop).
public struct S2Loop: Shape, S2Region {
  
  var vertices: [S2Point]
  
  // originInside keeps a precomputed value whether this loop contains the origin
  // versus computing from the set of vertices every time.
  var originInside: Bool
  
  // bound is a conservative bound on all points contained by this loop.
  // If l.ContainsPoint(P), then l.bound.ContainsPoint(P).
  var bound: S2Rect
  
  // Since "bound" is not exact, it is possible that a loop A contains
  // another loop B whose bounds are slightly larger. subregionBound
  // has been expanded sufficiently to account for this error, i.e.
  // if A.Contains(B), then A.subregionBound.Contains(B.bound).
  var subregionBound: S2Rect
  
  // depth is the nesting depth of this Loop if it is contained by a Polygon
  // or other shape and is used to determine if this loop represents a hole
  // or a filled in portion.
  var depth: Int = 0
  
  // index is the spatial index for this Loop.
  var index: ShapeIndex

  // MARK: inits / factory
  
  init(vertices: [S2Point], originInside: Bool, bound: S2Rect, subregionBound: S2Rect) {
    self.vertices = vertices
    self.originInside = originInside
    self.bound = bound
    self.subregionBound = subregionBound
    // create a new index and add us to it.
    index = ShapeIndex()
    index.add(shape: self)
  }
  
  /// Constructs a loop from the given points
  /// 1 point is a special empty or full, depending on it being in the northern hemisphere or not
  /// 0 or 2 points are incomplete loops
  init(points: [S2Point]) {
    let originInside = S2Loop.computeOrigin(vertices: points)
    let (bound, subregionBound) = S2Loop.computeBounds(vertices: points, originInside: originInside)
    self.init(vertices: points, originInside: originInside, bound: bound, subregionBound: subregionBound)
  }
  
  /// Constructs a loop corresponding to the given cell.
  /// Note that the loop and cell *do not* contain exactly the same set of
  /// points, because Loop and Cell have slightly different definitions of
  /// point containment. For example, a Cell vertex is contained by all
  /// four neighboring Cells, but it is contained by exactly one of four
  /// Loops constructed from those cells. As another example, the cell
  /// coverings of cell and LoopFromCell(cell) will be different, because the
  /// loop contains points on its boundary that actually belong to other cells
  /// (i.e., the covering will include a layer of neighboring cells).
  init(cell: Cell) {
    let points = (0..<4).map { cell.vertex($0) }
    let originInside = S2Loop.computeOrigin(vertices: points)
    let (bound, subregionBound) = S2Loop.computeBounds(vertices: points, originInside: originInside)
    self.init(vertices: points, originInside: originInside, bound: bound, subregionBound: subregionBound)
  }

  static func computeOrigin(vertices: [S2Point]) -> Bool {
    switch vertices.count {
    case 1:
      // this is the special empty or full loop
      // the origin depends on if the vertex is in the southern hemisphere or not.
      return vertices[0].z < 0
    case 0, 2:
      // these are incomplete loops
      return false
    default:
      // Point containment testing is done by counting edge crossings starting
      // at a fixed point on the sphere (OriginPoint). We need to know whether
      // the reference point (OriginPoint) is inside or outside the loop before
      // we can construct the S2ShapeIndex. We do this by first guessing that
      // it is outside, and then seeing whether we get the correct containment
      // result for vertex 1. If the result is incorrect, the origin must be
      // inside the loop.
      // A loop with consecutive vertices A,B,C contains vertex B if and only if
      // the fixed vector R = B.Ortho is contained by the wedge ABC. The
      // wedge is closed at A and open at C, i.e. the point B is inside the loop
      // if A = R but not if C = R. This convention is required for compatibility
      // with VertexCrossing. (Note that we can't use OriginPoint
      // as the fixed vector because of the possibility that B == OriginPoint.)
      let ortho = S2Point(raw: vertices[1].v.ortho())
      let v1Inside = S2Point.orderedCCW(ortho, vertices[0], vertices[2], vertices[1])
      return v1Inside != S2Loop.contains(vertices[1], vertices: vertices, originInside: false)
    }
  }
  
  static func computeBounds(vertices: [S2Point], originInside: Bool) -> (bound: S2Rect, subregionBound: S2Rect) {
    //
    switch vertices.count {
    case 1:
      // this is the special empty or full loop
      let bound = originInside ? S2Rect.full : S2Rect.empty
      return (bound: bound, subregionBound: bound)
    case 0, 2:
      // these are incomplete loops
      return (bound: S2Rect.empty, subregionBound: S2Rect.empty)
    default: break
    }
    // We *must* call initBound before initIndex, because initBound calls
    // ContainsPoint(s2.Point), and ContainsPoint(s2.Point) does a bounds check whenever the
    // index is not fresh (i.e., the loop has been added to the index but the
    // index has not been updated yet).
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices. First, the maximal latitude may be
    // attained along the interior of an edge. Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe). Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    var bounder = RectBounder()
    for p in vertices {
      bounder.add(point: p)
    }
    bounder.add(point: vertices[0])
    var b = bounder.rectBound()
    if S2Loop.contains(S2Point(x: 0, y: 0, z: 1), vertices: vertices, originInside: originInside) {
      b = S2Rect(lat: R1Interval(lo: b.lat.lo, hi: .pi / 2.0), lng: S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.Lng.isFull due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.Lng.isFull.
    if b.lng.isFull && S2Loop.contains(S2Point(x: 0, y: 0, z: -1), vertices: vertices, originInside: originInside) {
      b = S2Rect(lat: R1Interval(lo: -.pi / 2.0, hi: b.lat.hi), lng: S1Interval.full)
    }
    return (bound: b, subregionBound: b.expandForSubregions())
  }
  
  static func xcomputeBounds(vertices: [S2Point]) -> (originInside: Bool, bound: S2Rect, subregionBound: S2Rect) {
    //
    switch vertices.count {
    case 1:
      // this is the special empty or full loop
      // the origin depends on if the vertex is in the southern hemisphere or not.
      let originInside = vertices[0].z < 0
      let bound = originInside ? S2Rect.full : S2Rect.empty
      return (originInside: originInside, bound: bound, subregionBound: bound)
    case 0, 2:
      // these are incomplete loops
      return (originInside: false, bound: S2Rect.empty, subregionBound: S2Rect.empty)
    default: break
    }
    // Point containment testing is done by counting edge crossings starting
    // at a fixed point on the sphere (OriginPoint). We need to know whether
    // the reference point (OriginPoint) is inside or outside the loop before
    // we can construct the S2ShapeIndex. We do this by first guessing that
    // it is outside, and then seeing whether we get the correct containment
    // result for vertex 1. If the result is incorrect, the origin must be
    // inside the loop.
    // A loop with consecutive vertices A,B,C contains vertex B if and only if
    // the fixed vector R = B.Ortho is contained by the wedge ABC. The
    // wedge is closed at A and open at C, i.e. the point B is inside the loop
    // if A = R but not if C = R. This convention is required for compatibility
    // with VertexCrossing. (Note that we can't use OriginPoint
    // as the fixed vector because of the possibility that B == OriginPoint.)
    let ortho = S2Point(raw: vertices[1].v.ortho())
    let v1Inside = S2Point.orderedCCW(ortho, vertices[0], vertices[2], vertices[1])
    let originInside = v1Inside != S2Loop.contains(vertices[1], vertices: vertices, originInside: false)
    // We *must* call initBound before initIndex, because initBound calls
    // ContainsPoint(s2.Point), and ContainsPoint(s2.Point) does a bounds check whenever the
    // index is not fresh (i.e., the loop has been added to the index but the
    // index has not been updated yet).
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices. First, the maximal latitude may be
    // attained along the interior of an edge. Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe). Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    var bounder = RectBounder()
    for p in vertices {
      bounder.add(point: p)
    }
    bounder.add(point: vertices[0])
    var b = bounder.rectBound()
    if S2Loop.contains(S2Point(x: 0, y: 0, z: 1), vertices: vertices, originInside: originInside) {
      b = S2Rect(lat: R1Interval(lo: b.lat.lo, hi: .pi / 2.0), lng: S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.Lng.isFull due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.Lng.isFull.
    if b.lng.isFull && S2Loop.contains(S2Point(x: 0, y: 0, z: -1), vertices: vertices, originInside: originInside) {
      b = S2Rect(lat: R1Interval(lo: -.pi / 2.0, hi: b.lat.hi), lng: S1Interval.full)
    }
    return (originInside: originInside, bound: b, subregionBound: b.expandForSubregions())
  }
  
  /// Returns true if the loop contains the point
  static func contains(_ point: S2Point, vertices: [S2Point], originInside: Bool) -> Bool {
    // TODO: Move to bruteForceContains and update with ShapeIndex when available.
    // Empty and full loops don't need a special case, but invalid loops with
    // zero vertices do, so we might as well handle them all at once.
    if vertices.count < 3 {
      return originInside
    }
    var originInside = originInside
    var crosser = EdgeCrosser(a: S2Point.origin, b: point, c: vertices[0])
    for i in 1..<vertices.count {
      originInside = originInside != crosser.isEdgeOrVertexChainCrossing(d: vertices[i])
    }
    // test the closing edge of the loop last
    originInside = originInside != crosser.isEdgeOrVertexChainCrossing(d: vertices[0])
    return originInside
  }
  
  // These two points are used for the special Empty and Full loops.
  let emptyLoopPoint = S2Point(x: 0, y: 0, z: 1)
  let fullLoopPoint = S2Point(x: 0, y: 0, z: -1)

  /// Returns a special "empty" loop.
  static let empty = S2Loop(points: [S2Point(x: 0, y: 0, z: 1)])

  /// Returns a special "full" loop.
  static let full = S2Loop(points: [S2Point(x: 0, y: 0, z: -1)])

  // MARK: tests
  
  /// Reports true if this loop contains s2.OriginPoint().
  public func containsOrigin() -> Bool {
    return originInside
  }

  /// Returns true because all loops have an interior.
  public func hasInterior() -> Bool {
    return true
  }

  /// Reports true if this is the special "empty" loop that contains no points.
  var isEmpty: Bool {
    return isEmptyOrFull() && !containsOrigin()
  }
  
  /// Reports true if this is the special "full" loop that contains all points.
  var isFull: Bool {
    return isEmptyOrFull() && containsOrigin()
  }
  
  /// Reports true if this loop is either the "empty" or "full" special loops.
  func isEmptyOrFull() -> Bool {
    return vertices.count == 1
  }
  
  /// Reports whether this is a valid loop or not.
  var isValid: Bool {
    // all vertices must be unit length.
    for v in vertices {
      if !v.isUnit {
        return false
      }
    }
    // Loops must have at least 3 vertices (except for empty and full).
    if vertices.count < 3 {
      if isEmptyOrFull() {
        // skip remaining tests
        return true
      }
      // non-empty, non-full loops must have at least 3 vertices
      return false
    }
    // Loops are not allowed to have any duplicate vertices or edge crossings.
    // We split this check into two parts. First we check that no edge is
    // degenerate (identical endpoints). Then we check that there are no
    // intersections between non-adjacent edges (including at vertices). The
    // second check needs the ShapeIndex, so it does not fall within the scope
    // of this method.
    for i in 0..<vertices.count {
      let v = vertex(i)
      let v1 = vertex(i + 1)
      if v == v1 {
        // edge i is degenerate (duplicate vertex)
        return false
      }
      // Antipodal vertices are not allowed.
      let antipode = S2Point(raw: v1.mul(-1))
      if v == antipode {
        // vertices v and v1 are antipodal
        return false
      }
    }
    return true
  }

  /// Reports whether this loop represents a hole in its containing polygon.
  func isHole() -> Bool {
    return depth & 1 != 0
  }
  
  /// Returns -1 if this Loop represents a hole in its containing polygon, and +1 otherwise.
  func sign() -> Int {
    return isHole() ? -1 : 1
  }
  
  /// Reports whether the two loops have the same boundary. This is
  /// true if and only if the loops have the same vertices in the same cyclic order
  /// (i.e., the vertices may be cyclically rotated). The empty and full loops are
  /// considered to have different boundaries.
  func boundaryEquals(_ o: S2Loop) -> Bool {
    if vertices.count != o.vertices.count {
      return false
    }
    // Special case to handle empty or full loops.  Since they have the same
    // number of vertices, if one loop is empty/full then so is the other.
    if isEmptyOrFull() {
      return isEmpty == o.isEmpty
    }
    // Loop through the vertices to find the first of ours that matches the
    // starting vertex of the other loop. Use that offset to then 'align' the
    // vertices for comparison.
    for (offset, v) in vertices.enumerated() {
      if v == o.vertex(0) {
        // There is at most one starting offset since loop vertices are unique.
        for i in 0..<vertices.count {
          if vertex(i + offset) != o.vertex(i) {
            return false
          }
        }
        return true
      }
    }
    return false
  }
  
  /// Returns +1 if this loop contains the boundary of the other loop,
  /// -1 if it excludes the boundary of the other, and 0 if the boundaries of the two
  /// loops cross. Shared edges are handled as follows:
  ///   If XY is a shared edge, define Reversed(XY) to be true if XY
  ///     appears in opposite directions in both loops.
  ///   Then this loop contains XY if and only if Reversed(XY) == the other loop is a hole.
  ///   (Intuitively, this checks whether this loop contains a vanishingly small region
  ///   extending from the boundary of the other toward the interior of the polygon to
  ///   which the other belongs.)
  /// This function is used for testing containment and intersection of
  /// multi-loop polygons. Note that this method is not symmetric, since the
  /// result depends on the direction of this loop but not on the direction of
  /// the other loop (in the absence of shared edges).
  /// This requires that neither loop is empty, and if other loop IsFull, then it must not
  /// be a hole.
  func compareBoundary(_ o: S2Loop) -> Int {
    // The bounds must intersect for containment or crossing.
    if !bound.intersects(o.bound) {
      return -1
    }
    // Full loops are handled as though the loop surrounded the entire sphere.
    if isFull {
      return 1
    }
    if o.isFull {
      return -1
    }
    // Check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    let relation = CompareBoundaryRelation(reverse: o.isHole())
    if hasCrossingRelation(a: self, b: o, relation: relation) {
      return 0
    }
    if relation.foundSharedVertex {
      if relation.containsEdge {
        return 1
      }
      return -1
    }
    // There are no edge intersections or shared vertices, so we can check
    // whether A contains an arbitrary vertex of B.
    if contains(o.vertex(0)) {
      return 1
    }
    return -1
  }
  
  /// Returns the reference point for this loop.
  public func referencePoint() -> ReferencePoint {
    return ReferencePoint(origin: true, contained: originInside)
  }
  
  // MARK: computed members
  
  /// Returns the number of edges in this shape.
  public func numEdges() -> Int {
    if isEmptyOrFull() {
      return 0
    }
    return vertices.count
  }

  /// Returns the endpoints for the given edge index.
  public func edge(_ i: Int) -> Edge {
    return Edge(v0: vertex(i), v1: vertex(i + 1))
  }
//  public func edge(_ i: Int) -> (S2Point, S2Point) {
//    return (vertex(i), vertex(i + 1))
//  }

  /// Reports the number of contiguous edge chains in the Loop.
  public func numChains() -> Int {
    if isEmptyOrFull() {
      return 0
    }
    return 1
  }
  
  /// Returns the i-th edge chain in the Shape.
  public func chain(_ chainId: Int) -> Chain {
    return Chain(start: 0, length: numEdges())
  }
  
  /// Returns the j-th edge of the i-th edge chain.
  public func chainEdge(chainId: Int, offset: Int) -> Edge {
    return Edge(v0: vertex(offset), v1: vertex(offset + 1))
  }
  
  // ChainPosition returns a ChainPosition pair (i, j) such that edgeID is the
  // j-th edge of the Loop.
  public func chainPosition(_ edgeId: Int) -> ChainPosition {
    return ChainPosition(chainId: 0, offset: edgeId)
  }
  
  // dimension returns the dimension of the geometry represented by this Loop.
  public func dimension() -> ShapeDimension {
    return .polygonGeometry
  }
  
  /// Returns a tight bounding rectangle. If the loop contains the point,
  /// the bound also contains it.
  public func rectBound() -> S2Rect {
    return bound
  }

  /// Returns a bounding cap that may have more padding than the corresponding
  /// RectBound. The bound is conservative such that if the loop contains a point P,
  /// the bound also contains it.
  public func capBound() -> S2Cap {
    return bound.capBound()
  }

  /// Returns the vertex for the given index. For convenience, the vertex indices
  /// wrap automatically for methods that do index math such as Edge.
  /// i.e., Vertex(NumEdges() + n) is the same as Vertex(n).
  func vertex(_ i: Int) -> S2Point {
    return vertices[i % vertices.count]
  }

  /// Returns the vertex in reverse order if the loop represents a polygon
  /// hole. For example, arguments 0, 1, 2 are mapped to vertices n-1, n-2, n-3, where
  /// n == len(vertices). This ensures that the interior of the polygon is always to
  /// the left of the vertex chain.
  ///This requires: 0 <= i < 2 * len(vertices)
  func orientedVertex(_ i: Int) -> S2Point {
    var j = i - vertices.count
    if j < 0 {
      j = i
    }
    if isHole() {
      j = vertices.count - 1 - j
    }
    return vertex(i)
  }
  
  // NumVertices returns the number of vertices in this loop.
  func numVertices() -> Int {
    return vertices.count
  }
  
  // MARK: intersects / contains
  
  /// Reports if the given point is contained by this loop.
  /// This method does not use the ShapeIndex, so it is only preferable below a certain
  /// size of loop.
  func bruteForceContains(_ p: S2Point) -> Bool {
    let origin = S2Point.origin
    var inside = originInside
    var crosser = EdgeCrosser(a: origin, b: p, c: vertex(0))
    for i in 1...vertices.count { // add vertex 0 twice
      inside = inside != crosser.isEdgeOrVertexChainCrossing(d: vertex(i))
    }
    return inside
  }
  
  /// Returns true if the loop contains the point
  func contains(_ point: S2Point) -> Bool {
    // empty and full loops don't need a special case, but invalid loops with
    // zero vertices do, so we might as well handle them all at once.
    if vertices.count < 3 {
      return originInside
    }
    // For small loops, and during initial construction, it is faster to just
    // check all the crossing.
    let maxBruteForceVertices = 32
    if vertices.count < maxBruteForceVertices { //} || index == nil {
      return bruteForceContains(point)
    }
    // Otherwise, look up the point in the index.
    var it = index.iterator()
    if !it.locate(point: point) {
      return false
    }
    return iteratorContains(iterator: it, point: point)
    // TODO: Move to bruteForceContains and update with ShapeIndex when available.
//    let origin = S2Point.origin
//    var inside = originInside
//    var crosser = EdgeCrosser(a: origin, b: point, c: vertices[0])
//    for i in 1..<vertices.count {
//      inside = inside != crosser.isEdgeOrVertexChainCrossing(d: vertices[i])
//    }
//    // test the closing edge of the loop last
//    inside = inside != crosser.isEdgeOrVertexChainCrossing(d: vertices[0])
//    return inside
  }

  // MARK: intersects / contains Cell
  
  /// Reports whether the loop contains the given cell
  public func contains(_ cell: Cell) -> Bool {
    var it = index.iterator()
    let relation = it.locate(cellId: cell.id)
    // If "target" is disjoint from all index cells, it is not contained.
    // Similarly, if "target" is subdivided into one or more index cells then it
    // is not contained, since index cells are subdivided only if they (nearly)
    // intersect a sufficient number of edges.  (But note that if "target" itself
    // is an index cell then it may be contained, since it could be a cell with
    // no edges in the loop interior.)
    if relation != .indexed {
      return false
    }
    // Otherwise check if any edges intersect "target".
    if boundaryApproxIntersects(iterator: it, cell: cell) {
      return false
    }
    // Otherwise check if the loop contains the center of "target".
    return iteratorContains(iterator: it, point: cell.center())
    // old implementation
//    let cellVertices = (0..<4).map { cell.vertex($0) }
//    // if the loop does not contain all cell vertices, return false
//    for k in 0..<4 {
//      if !contains(cellVertices[k]) {
//        return false
//      }
//    }
//    // if there are any edge crossing, it is not containing
//    for j in 0..<4 {
//      var crosser = EdgeCrosser(a: cellVertices[j], b: cellVertices[(j+1)&3], c: vertices[0])
//      for i in 1..<vertices.count {
//        if crosser.chainCrossingSign(d: vertices[i]) != .doNotCross {
//          // There is a proper crossing, or two vertices were the same.
//          return false
//        }
//      }
//    }
//    return true
  }

  /// Reports whether the loop intersects the cell
  public func intersects(_ cell: Cell) -> Bool {
    var it = index.iterator()
    let relation = it.locate(cellId: cell.id)
    // If target does not overlap any index cell, there is no intersection.
    if relation == .disjoint {
      return false
    }
    // If target is subdivided into one or more index cells, there is an
    // intersection to within the ShapeIndex error bound (see Contains).
    if relation == .subdivided {
      return true
    }
    // If target is an index cell, there is an intersection because index cells
    // are created only if they have at least one edge or they are entirely
    // contained by the loop.
    if it.cellId() == cell.id {
      return true
    }
    // Otherwise check if any edges intersect target.
    if boundaryApproxIntersects(iterator: it, cell: cell) {
      return true
    }
    // Otherwise check if the loop contains the center of target.
    return iteratorContains(iterator: it, point: cell.center())
    // old implementation
//    let cellVertices = (0..<4).map { cell.vertex($0) }
//    // intersects if the loop contains any cell vertex
//    for k in 0..<4 {
//      let vertex = cell.vertex(k)
//      if contains(vertex) {
//        return true
//      }
//    }
//    // intersects if any loop edge crosses with any cell edge
//    for j in 0..<4 {
//      var crosser = EdgeCrosser(a: cellVertices[j], b: cellVertices[(j+1)&3], c: vertices[0])
//      for i in 1..<vertices.count {
//        if crosser.chainCrossingSign(d: vertices[i]) != .doNotCross {
//          // There is a proper crossing, or two vertices were the same.
//          return true
//        }
//      }
//    }
//    // intersects if tightest rect is enclosed in cell
//    var bounder = RectBounder()
//    for p in vertices {
//      bounder.add(point: p)
//    }
//    bounder.add(point: vertices[0])
//    let rect = bounder.rectBound()
//    if cell.rectBound().contains(rect) {
//      return true
//    }
//    // otherwise fals
//    return false
  }
  
  /// Computes a covering of the Loop.
  func cellUnionBound() -> CellUnion {
    return capBound().cellUnionBound()
  }
  
  /// Reports if the loop's boundary intersects target.
  /// It may also return true when the loop boundary does not intersect target but
  /// some edge comes within the worst-case error tolerance.
  /// This requires that it.Locate(target) returned Indexed.
  func boundaryApproxIntersects(iterator it: ShapeIndexIterator, cell target: Cell) -> Bool {
    guard let aClipped = it.indexCell()?.find(shapeId: 0) else { return false }
    // If there are no edges, there is no intersection.
    if aClipped.edges.count == 0 {
      return false
    }
    // We can save some work if target is the index cell itself.
    if it.cellId() == target.id {
      return true
    }
    // Otherwise check whether any of the edges intersect target.
    let maxError = faceClipErrorUVCoord + intersectsRectErrorUVDist
    let bound = target.boundUV().expanded(maxError)
    for ai in aClipped.edges {
      let (v0, v1, ok) = clipToPaddedFace(a: vertex(ai), b: vertex(ai + 1), f: Int(target.face), padding: maxError)
      if ok && edgeIntersectsRect(a: v0, b: v1, r: bound) {
        return true
      }
    }
    return false
  }
  
  /// Reports if the iterator that is positioned at the ShapeIndexCell
  /// that may contain p, contains the point p.
  func iteratorContains(iterator it: ShapeIndexIterator, point p: S2Point) -> Bool {
    // Test containment by drawing a line segment from the cell center to the
    // given point and counting edge crossings.
    guard let aClipped = it.indexCell()?.find(shapeId: 0) else { return false }
    var inside = aClipped.containsCenter
    if aClipped.edges.count > 0 {
      let center = it.center()
      var crosser = EdgeCrosser(a: center, b: p)
      var aiPrev = -2
      for ai in aClipped.edges {
        if ai != aiPrev + 1 {
          crosser.restart(at: vertex(ai))
        }
        aiPrev = ai
        inside = inside != crosser.isEdgeOrVertexChainCrossing(d: vertex(ai + 1))
      }
    }
    return inside
  }
  
  // Contains reports whether the region contained by this loop is a superset of the
  // region contained by the given other loop.
  func contains(_ o: S2Loop) -> Bool {
    // For a loop A to contain the loop B, all of the following must
    // be true:
    //  (1) There are no edge crossings between A and B except at vertices.
    //  (2) At every vertex that is shared between A and B, the local edge
    //      ordering implies that A contains B.
    //  (3) If there are no shared vertices, then A must contain a vertex of B
    //      and B must not contain a vertex of A. (An arbitrary vertex may be
    //      chosen in each case.)
    // The second part of (3) is necessary to detect the case of two loops whose
    // union is the entire sphere, i.e. two loops that contains each other's
    // boundaries but not each other's interiors.
    if !subregionBound.contains(o.bound) {
      return false
    }
    // special cases to handle either loop being empty or full.
    if isEmptyOrFull() || o.isEmptyOrFull() {
      return isFull || o.isEmpty
    }
    // Check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    let relation = ContainsRelation()
    if hasCrossingRelation(a: self, b: o, relation: relation) {
      return false
    }
    // There are no crossings, and if there are any shared vertices then A
    // contains B locally at each shared vertex.
    if relation.foundSharedVertex {
      return true
    }
    // Since there are no edge intersections or shared vertices, we just need to
    // test condition (3) above. We can skip this test if we discovered that A
    // contains at least one point of B while checking for edge crossings.
    if !contains(o.vertex(0)) {
      return false
    }
    // We still need to check whether (A union B) is the entire sphere.
    // Normally this check is very cheap due to the bounding box precondition.
    if (o.subregionBound.contains(bound) || o.bound.union(bound).isFull) && o.contains(vertex(0)) {
      return false
    }
    return true
  }
  
  // Intersects reports whether the region contained by this loop intersects the region
  // contained by the other loop.
  func intersects(_ o: S2Loop) -> Bool {
    // Given two loops, A and B, A.Intersects(B) if and only if !A.Complement().Contains(B).
    // This code is similar to Contains, but is optimized for the case
    // where both loops enclose less than half of the sphere.
    if !bound.intersects(o.bound) {
      return false
    }
    // Check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    let relation = IntersectsRelation()
    if hasCrossingRelation(a: self, b: o, relation: relation) {
      return true
    }
    if relation.foundSharedVertex {
      return false
    }
    // Since there are no edge intersections or shared vertices, the loops
    // intersect only if A contains B, B contains A, or the two loops contain
    // each other's boundaries.  These checks are usually cheap because of the
    // bounding box preconditions.  Note that neither loop is empty (because of
    // the bounding box check above), so it is safe to access vertex(0).
    // Check whether A contains B, or A and B contain each other's boundaries.
    // (Note that A contains all the vertices of B in either case.)
    if subregionBound.contains(o.bound) || bound.union(o.bound).isFull {
      if contains(o.vertex(0)) {
        return true
      }
    }
    // Check whether B contains A.
    if o.subregionBound.contains(bound) {
      if o.contains(vertex(0)) {
        return true
      }
    }
    return false
  }
  
  /// Reports whether given two loops whose boundaries
  /// do not cross (see compareBoundary), if this loop contains the boundary of the
  /// other loop. If reverse is true, the boundary of the other loop is reversed
  /// first (which only affects the result when there are shared edges). This method
  /// is cheaper than compareBoundary because it does not test for edge intersections.
  ///
  /// This function requires that neither loop is empty, and that if the other is full,
  /// then reverse == false.
  func containsNonCrossingBoundary(_ other: S2Loop, reverse reverseOther: Bool) -> Bool {
    // The bounds must intersect for containment.
    if !bound.intersects(other.bound) {
      return false
    }
    // Full loops are handled as though the loop surrounded the entire sphere.
    if isFull {
      return true
    }
    if other.isFull {
      return false
    }
    guard let m = find(vertex: other.vertex(0)) else {
      // Since the other loops vertex 0 is not shared, we can check if this contains it.
      return contains(other.vertex(0))
    }
    // Otherwise check whether the edge (b0, b1) is contained by this loop.
    return wedgeContainsSemiwedge(a0: vertex(m-1), ab1: vertex(m), a2: vertex(m+1), b2: other.vertex(1), reverse: reverseOther)
  }

  /// Returns the index of the vertex at the given Point in the range
  /// 1..numVertices, and a boolean indicating if a vertex was found.
  func find(vertex p: S2Point) -> Int? {
    if vertices.count < 10 {
      // Exhaustive search for loops below a small threshold.
      for i in 1...vertices.count {
        if vertex(i) == p {
          return i
        }
      }
      return nil
    }
    var it = index.iterator()
    guard it.locate(point: p) else { return nil }
    guard let aClipped = it.indexCell()!.find(shapeId: 0) else { return nil }
    let nEdges = aClipped.numEdges()
    for j in 1...nEdges {
      let i = nEdges - j
      let ai = aClipped.edges[i]
      if vertex(ai) == p {
        if ai == 0 {
          return vertices.count
        }
      return ai
      }
      if vertex(ai + 1) == p {
        return ai + 1
      }
    }
    return nil
  }
 
  // MARK:
  
  /// Creates a loop with the given number of vertices, all
  /// located on a circle of the specified radius around the given center.
  static func regularLoop(center: S2Point, radius: S1Angle, numVertices: Int) -> S2Loop {
    return S2Loop.regularLoopForFrame(frame: S2Point.getFrame(center), radius: radius, numVertices: numVertices)
  }
  
  /// Creates a loop centered around the z-axis of the given
  /// coordinate frame, with the first vertex in the direction of the positive x-axis.
  static func regularLoopForFrame(frame: R3Matrix, radius: S1Angle, numVertices: Int) -> S2Loop {
    return S2Loop(points: S2Point.regularPointsForFrame(frame: frame, radius: radius, numVertices: numVertices))
  }

  // MARK:
  
  /// Returns a first index and a direction (either +1 or -1)
  /// such that the vertex sequence (first, first+dir, ..., first+(n-1)*dir) does
  /// not change when the loop vertex order is rotated or inverted. This allows the
  /// loop vertices to be traversed in a canonical order. The return values are
  /// chosen such that (first, ..., first+n*dir) are in the range [0, 2*n-1] as
  /// expected by the Vertex method.
  func canonicalFirstVertex() -> (firstIdx: Int, direction: Int) {
    var firstIdx = 0
    let n = vertices.count
    for i in 1..<n {
      if vertex(i) < vertex(firstIdx) {
        firstIdx = i
      }
    }
    // 0 <= firstIdx <= n-1, so (firstIdx+n*dir) <= 2*n-1.
    if vertex(firstIdx + 1) < vertex(firstIdx + n - 1) {
      return (firstIdx, 1)
    }
    // n <= firstIdx <= 2*n-1, so (firstIdx+n*dir) >= 0.
    firstIdx += n
    return (firstIdx, -1)
  }
  
  // TurningAngle returns the sum of the turning angles at each vertex. The return
  // value is positive if the loop is counter-clockwise, negative if the loop is
  // clockwise, and zero if the loop is a great circle. Degenerate and
  // nearly-degenerate loops are handled consistently with Sign. So for example,
  // if a loop has zero area (i.e., it is a very small CCW loop) then the turning
  // angle will always be negative.
  //
  // This quantity is also called the "geodesic curvature" of the loop.
  func turningAngle() -> Double {
    // For empty and full loops, we return the limit value as the loop area
    // approaches 0 or 4*Pi respectively.
    if isEmptyOrFull() {
      if containsOrigin() {
        return -2 * .pi
      }
      return 2 * .pi
    }
    // Don't crash even if the loop is not well-defined.
    if vertices.count < 3 {
      return 0
    }
    // To ensure that we get the same result when the vertex order is rotated,
    // and that the result is negated when the vertex order is reversed, we need
    // to add up the individual turn angles in a consistent order. (In general,
    // adding up a set of numbers in a different order can change the sum due to
    // rounding errors.)
    //
    // Furthermore, if we just accumulate an ordinary sum then the worst-case
    // error is quadratic in the number of vertices. (This can happen with
    // spiral shapes, where the partial sum of the turning angles can be linear
    // in the number of vertices.) To avoid this we use the Kahan summation
    // algorithm (http://en.wikipedia.org/wiki/Kahan_summation_algorithm).
    var n = vertices.count
    var (i, dir) = canonicalFirstVertex()
    var sum = S2Point.turnAngle(a: vertex((i + n - dir) % n), b: vertex(i), c: vertex((i + dir) % n))
    var compensation = S1Angle(0)
    while n - 1 > 0 {
      i += dir
      var angle = S2Point.turnAngle(a: vertex(i - dir), b: vertex(i), c: vertex(i + dir))
      let oldSum = sum
      angle += compensation
      sum += angle
      compensation = (oldSum - sum) + angle
      n -= 1
    }
    return Double(dir) * Double(sum + compensation)
  }
  
  // turningAngleMaxError return the maximum error in TurningAngle. The value is not
  // constant; it depends on the loop.
  func turningAngleMaxError() -> Double {
    // The maximum error can be bounded as follows:
    //   2.24 * dblEpsilon    for RobustCrossProd(b, a)
    //   2.24 * dblEpsilon    for RobustCrossProd(c, b)
    //   3.25 * dblEpsilon    for Angle()
    //   2.00 * dblEpsilon    for each addition in the Kahan summation
    //   ------------------
    //   9.73 * dblEpsilon
    let maxErrorPerVertex = 9.73 * S1Interval.dblEpsilon
    return maxErrorPerVertex * Double(vertices.count)
  }
  
  // IsNormalized reports whether the loop area is at most 2*pi. Degenerate loops are
  // handled consistently with Sign, i.e., if a loop can be
  // expressed as the union of degenerate or nearly-degenerate CCW triangles,
  // then it will always be considered normalized.
  func isNormalized() -> Bool {
    // Optimization: if the longitude span is less than 180 degrees, then the
    // loop covers less than half the sphere and is therefore normalized.
    if bound.lng.length < .pi {
      return true
    }
    // We allow some error so that hemispheres are always considered normalized.
    // TODO(roberts): This is no longer required by the Polygon implementation,
    // so alternatively we could create the invariant that a loop is normalized
    // if and only if its complement is not normalized.
    return turningAngle() >= -turningAngleMaxError()
  }
  
  // Normalize inverts the loop if necessary so that the area enclosed by the loop
  // is at most 2*pi.
  mutating func normalize() {
    if !isNormalized() {
      invert()
    }
  }
  
  // Invert reverses the order of the loop vertices, effectively complementing the
  // region represented by the loop. For example, the loop ABCD (with edges
  // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
  // Notice that the last edge is the same in both cases except that its
  // direction has been reversed.
  mutating func invert() {
    index.reset()
    if isEmptyOrFull() {
      if isFull {
        vertices[0] = emptyLoopPoint
      } else {
        vertices[0] = fullLoopPoint
      }
    } else {
      // For non-special loops, reverse the slice of vertices.
      vertices.reverse()
    }
    // originInside must be set correctly before building the ShapeIndex.
    originInside = !originInside
    if bound.lat.lo > -.pi / 2 && bound.lat.hi < .pi / 2 {
      // The complement of this loop contains both poles.
      bound = S2Rect.full
      subregionBound = bound
    } else {
      self.originInside = S2Loop.computeOrigin(vertices: vertices)
    }
    index.add(shape: self)
  }
  
  // ContainsNested reports whether the given loops is contained within this loop.
  // This function does not test for edge intersections. The two loops must meet
  // all of the Polygon requirements; for example this implies that their
  // boundaries may not cross or have any shared edges (although they may have
  // shared vertices).
  func containsNested(_ other: S2Loop) -> Bool {
    if !subregionBound.contains(other.bound) {
      return false
    }
    // Special cases to handle either loop being empty or full.  Also bail out
    // when B has no vertices to avoid heap overflow on the vertex(1) call
    // below.  (This method is called during polygon initialization before the
    // client has an opportunity to call IsValid().)
    if isEmptyOrFull() || other.numVertices() < 2 {
      return isFull || other.isEmpty
    }
    // We are given that A and B do not share any edges, and that either one
    // loop contains the other or they do not intersect.
    guard let m = find(vertex: other.vertex(1)) else {
      // Since other.vertex(1) is not shared, we can check whether A contains it.
      return contains(other.vertex(1))
    }
    // Check whether the edge order around other.Vertex(1) is compatible with
    // A containing B.
    return wedgeContains(a0: vertex(m - 1), ab1: vertex(m), a2: vertex(m + 1), b0: other.vertex(0), b2: other.vertex(2))
  }

  // MARK:
  
  /// Computes the oriented surface integral of some quantity f(x)
  /// over the loop interior, given a function f(A,B,C) that returns the
  /// corresponding integral over the spherical triangle ABC. Here "oriented
  /// surface integral" means:
  /// (1) f(A,B,C) must be the integral of f if ABC is counterclockwise,
  ///     and the integral of -f if ABC is clockwise.
  /// (2) The result of this function is *either* the integral of f over the
  ///     loop interior, or the integral of (-f) over the loop exterior.
  /// Note that there are at least two common situations where it easy to work
  /// around property (2) above:
  ///  - If the integral of f over the entire sphere is zero, then it doesn't
  ///    matter which case is returned because they are always equal.
  ///  - If f is non-negative, then it is easy to detect when the integral over
  ///    the loop exterior has been returned, and the integral over the loop
  ///    interior can be obtained by adding the integral of f over the entire
  ///    unit sphere (a constant) to the result.
  /// Any changes to this method may need corresponding changes to surfaceIntegralPoint as well.
  func surfaceIntegral(f: (S2Point, S2Point, S2Point) -> Double) -> Double {
    // We sum f over a collection T of oriented triangles, possibly
    // overlapping. Let the sign of a triangle be +1 if it is CCW and -1
    // otherwise, and let the sign of a point x be the sum of the signs of the
    // triangles containing x. Then the collection of triangles T is chosen
    // such that either:
    //  (1) Each point in the loop interior has sign +1, and sign 0 otherwise; or
    //  (2) Each point in the loop exterior has sign -1, and sign 0 otherwise.
    // The triangles basically consist of a fan from vertex 0 to every loop
    // edge that does not include vertex 0. These triangles will always satisfy
    // either (1) or (2). However, what makes this a bit tricky is that
    // spherical edges become numerically unstable as their length approaches
    // 180 degrees. Of course there is not much we can do if the loop itself
    // contains such edges, but we would like to make sure that all the triangle
    // edges under our control (i.e., the non-loop edges) are stable. For
    // example, consider a loop around the equator consisting of four equally
    // spaced points. This is a well-defined loop, but we cannot just split it
    // into two triangles by connecting vertex 0 to vertex 2.
    // We handle this type of situation by moving the origin of the triangle fan
    // whenever we are about to create an unstable edge. We choose a new
    // location for the origin such that all relevant edges are stable. We also
    // create extra triangles with the appropriate orientation so that the sum
    // of the triangle signs is still correct at every point.
    // The maximum length of an edge for it to be considered numerically stable.
    // The exact value is fairly arbitrary since it depends on the stability of
    // the function f. The value below is quite conservative but could be
    // reduced further if desired.
    let maxLength = .pi - 1e-5
    var sum = 0.0
    var origin = vertex(0)
    for i in 1..<vertices.count - 1 {
      // Let V_i be vertex(i), let O be the current origin, and let length(A,B)
      // be the length of edge (A,B). At the start of each loop iteration, the
      // "leading edge" of the triangle fan is (O,V_i), and we want to extend
      // the triangle fan so that the leading edge is (O,V_i+1).
      // Invariants:
      //  1. length(O,V_i) < maxLength for all (i > 1).
      //  2. Either O == V_0, or O is approximately perpendicular to V_0.
      //  3. "sum" is the oriented integral of f over the area defined by
      //     (O, V_0, V_1, ..., V_i).
      if vertex(i + 1).angle(origin) > maxLength {
        // We are about to create an unstable edge, so choose a new origin O'
        // for the triangle fan.
        let oldOrigin = origin
        if origin == vertex(0) {
          // The following point is well-separated from V_i and V_0 (and
          // therefore V_i+1 as well).
          origin = vertex(0).pointCross(vertex(i))
        } else if vertex(i).angle(vertex(0)) < maxLength {
          // All edges of the triangle (O, V_0, V_i) are stable, so we can
          // revert to using V_0 as the origin.
          origin = vertex(0)
        } else {
          // (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are
          // perpendicular. Therefore V_0.CrossProd(O) is approximately
          // perpendicular to all of {O, V_0, V_i, V_i+1}, and we can choose
          // this point O' as the new origin.
          origin = S2Point(raw: vertex(0).cross(oldOrigin))
          // Advance the edge (V_0,O) to (V_0,O').
          sum += f(vertex(0), oldOrigin, origin)
        }
        // Advance the edge (O,V_i) to (O',V_i).
        sum += f(oldOrigin, vertex(i), origin)
      }
      // Advance the edge (O,V_i) to (O,V_i+1).
      sum += f(origin, vertex(i), vertex(i + 1))
    }
    // If the origin is not V_0, we need to sum one more triangle.
    if origin != vertex(0) {
      // Advance the edge (O,V_n-1) to (O,V_0).
      sum += f(origin, vertex(vertices.count - 1), vertex(0))
    }
    return sum
  }
  
  /// Mirrors the surfaceIntegralFloat64 method but over Points;
  /// see that method for commentary. The C++ version uses a templated method.
  /// Any changes to this method may need corresponding changes to surfaceIntegralFloat64 as well.
  func surfaceIntegral(f: (S2Point, S2Point, S2Point) -> R3Vector) -> S2Point {
    let maxLength = .pi - 1e-5
    var sum = R3Vector(x: 0, y: 0, z: 0)
    var origin = vertex(0)
    for i in 1..<vertices.count - 1 {
      if vertex(i + 1).angle(origin) > maxLength {
        let oldOrigin = origin
        if origin == vertex(0) {
          origin = vertex(0).pointCross(vertex(i))
        } else if vertex(i).angle(vertex(0)) < maxLength {
          origin = vertex(0)
        } else {
          origin = S2Point(raw: vertex(0).cross(oldOrigin))
          sum = sum.add(f(vertex(0), oldOrigin, origin))
        }
        sum = sum.add(f(oldOrigin, vertex(i), origin))
      }
      sum = sum.add(f(origin, vertex(i), vertex(i + 1)))
    }
    if origin != vertex(0) {
      sum = sum.add(f(origin, vertex(vertices.count - 1), vertex(0)))
    }
    return S2Point(raw: sum)
  }
  
  /// Returns the area of the loop interior, i.e. the region on the left side of
  /// the loop. The return value is between 0 and 4*pi. (Note that the return
  /// value is not affected by whether this loop is a "hole" or a "shell".)
  func area() -> Double {
    // It is suprisingly difficult to compute the area of a loop robustly. The
    // main issues are (1) whether degenerate loops are considered to be CCW or
    // not (i.e., whether their area is close to 0 or 4*pi), and (2) computing
    // the areas of small loops with good relative accuracy.
    //
    // With respect to degeneracies, we would like Area to be consistent
    // with ContainsPoint in that loops that contain many points
    // should have large areas, and loops that contain few points should have
    // small areas. For example, if a degenerate triangle is considered CCW
    // according to s2predicates Sign, then it will contain very few points and
    // its area should be approximately zero. On the other hand if it is
    // considered clockwise, then it will contain virtually all points and so
    // its area should be approximately 4*pi.
    //
    // More precisely, let U be the set of Points for which IsUnitLength
    // is true, let P(U) be the projection of those points onto the mathematical
    // unit sphere, and let V(P(U)) be the Voronoi diagram of the projected
    // points. Then for every loop x, we would like Area to approximately
    // equal the sum of the areas of the Voronoi regions of the points p for
    // which x.ContainsPoint(p) is true.
    //
    // The second issue is that we want to compute the area of small loops
    // accurately. This requires having good relative precision rather than
    // good absolute precision. For example, if the area of a loop is 1e-12 and
    // the error is 1e-15, then the area only has 3 digits of accuracy. (For
    // reference, 1e-12 is about 40 square meters on the surface of the earth.)
    // We would like to have good relative accuracy even for small loops.
    //
    // To achieve these goals, we combine two different methods of computing the
    // area. This first method is based on the Gauss-Bonnet theorem, which says
    // that the area enclosed by the loop equals 2*pi minus the total geodesic
    // curvature of the loop (i.e., the sum of the "turning angles" at all the
    // loop vertices). The big advantage of this method is that as long as we
    // use Sign to compute the turning angle at each vertex, then
    // degeneracies are always handled correctly. In other words, if a
    // degenerate loop is CCW according to the symbolic perturbations used by
    // Sign, then its turning angle will be approximately 2*pi.
    //
    // The disadvantage of the Gauss-Bonnet method is that its absolute error is
    // about 2e-15 times the number of vertices (see turningAngleMaxError).
    // So, it cannot compute the area of small loops accurately.
    //
    // The second method is based on splitting the loop into triangles and
    // summing the area of each triangle. To avoid the difficulty and expense
    // of decomposing the loop into a union of non-overlapping triangles,
    // instead we compute a signed sum over triangles that may overlap (see the
    // comments for surfaceIntegral). The advantage of this method
    // is that the area of each triangle can be computed with much better
    // relative accuracy (using l'Huilier's theorem). The disadvantage is that
    // the result is a signed area: CCW loops may yield a small positive value,
    // while CW loops may yield a small negative value (which is converted to a
    // positive area by adding 4*pi). This means that small errors in computing
    // the signed area may translate into a very large error in the result (if
    // the sign of the sum is incorrect).
    //
    // So, our strategy is to combine these two methods as follows. First we
    // compute the area using the "signed sum over triangles" approach (since it
    // is generally more accurate). We also estimate the maximum error in this
    // result. If the signed area is too close to zero (i.e., zero is within
    // the error bounds), then we double-check the sign of the result using the
    // Gauss-Bonnet method. (In fact we just call IsNormalized, which is
    // based on this method.) If the two methods disagree, we return either 0
    // or 4*pi based on the result of IsNormalized. Otherwise we return the
    // area that we computed originally.
    if isEmptyOrFull() {
      if containsOrigin() {
        return 4 * .pi
      }
      return 0
    }
    var area = surfaceIntegral(f: S2Point.signedArea)
    // TODO(roberts): This error estimate is very approximate. There are two
    // issues: (1) SignedArea needs some improvements to ensure that its error
    // is actually never higher than GirardArea, and (2) although the number of
    // triangles in the sum is typically N-2, in theory it could be as high as
    // 2*N for pathological inputs. But in other respects this error bound is
    // very conservative since it assumes that the maximum error is achieved on
    // every triangle.
    let maxError = turningAngleMaxError()
    // The signed area should be between approximately -4*pi and 4*pi.
    if area < 0 {
      // We have computed the negative of the area of the loop exterior.
      area += 4 * .pi
    }
    if area > 4 * .pi {
      area = 4 * .pi
    }
    if area < 0 {
      area = 0
    }
    // If the area is close enough to zero or 4*pi so that the loop orientation
    // is ambiguous, then we compute the loop orientation explicitly.
    if area < maxError && !isNormalized() {
      return 4 * .pi
    } else if area > (4 * .pi - maxError) && isNormalized() {
      return 0
    }
    return area
  }

  
  /// Returns the true centroid of the loop multiplied by the area of the
  /// loop. The result is not unit length, so you may want to normalize it. Also
  /// note that in general, the centroid may not be contained by the loop.
  /// We prescale by the loop area for two reasons: (1) it is cheaper to
  /// compute this way, and (2) it makes it easier to compute the centroid of
  /// more complicated shapes (by splitting them into disjoint regions and
  /// adding their centroids).
  /// Note that the return value is not affected by whether this loop is a
  /// "hole" or a "shell".
  func centroid() -> S2Point {
    // surfaceIntegralPoint() returns either the integral of position over loop
    // interior, or the negative of the integral of position over the loop
    // exterior. But these two values are the same (!), because the integral of
    // position over the entire sphere is (0, 0, 0).
    return surfaceIntegral(f: S2Point.trueCentroid)
  }
  
}

// MARK: Loop Relations

/// Representing the possible crossing target cases for relations.
enum CrossingTarget: Int {
  case dontCare, dontCross, cross
}

/// Defines the interface for checking a type of relationship between two loops.
/// Some examples of relations are Contains, Intersects, or CompareBoundary.
protocol LoopRelation {
  // Optionally, aCrossingTarget and bCrossingTarget can specify an early-exit
  // condition for the loop relation. If any point P is found such that
  //   A.ContainsPoint(P) == aCrossingTarget() &&
  //   B.ContainsPoint(P) == bCrossingTarget()
  // then the loop relation is assumed to be the same as if a pair of crossing
  // edges were found. For example, the ContainsPoint relation has
  //   aCrossingTarget() == crossingTargetDontCross
  //   bCrossingTarget() == crossingTargetCross
  // because if A.ContainsPoint(P) == false and B.ContainsPoint(P) == true
  // for any point P, then it is equivalent to finding an edge crossing (i.e.,
  // since Contains returns false in both cases).
  // Loop relations that do not have an early-exit condition of this form
  // should return crossingTargetDontCare for both crossing targets.
  // aCrossingTarget reports whether loop A crosses the target point with
  // the given relation type.
  func aCrossingTarget() -> CrossingTarget
  // bCrossingTarget reports whether loop B crosses the target point with
  // the given relation type.
  func bCrossingTarget() -> CrossingTarget
  // wedgesCross reports if a shared vertex ab1 and the two associated wedges
  // (a0, ab1, b2) and (b0, ab1, b2) are equivalent to an edge crossing.
  // The loop relation is also allowed to maintain its own internal state, and
  // can return true if it observes any sequence of wedges that are equivalent
  // to an edge crossing.
  func wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool
}

/// Implements loopRelation for a contains operation. If
/// A.ContainsPoint(P) == false && B.ContainsPoint(P) == true, it is equivalent
/// to having an edge crossing (i.e., Contains returns false).
class ContainsRelation: LoopRelation {

  var foundSharedVertex = false

  func aCrossingTarget() -> CrossingTarget { return .dontCross }
  func bCrossingTarget() -> CrossingTarget { return .cross }
  func wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
    foundSharedVertex = true
    return !wedgeContains(a0: a0, ab1: ab1, a2: a2, b0: b0, b2: b2)
  }
  
}

/// Implements loopRelation for an intersects operation. Given
/// two loops, A and B, if A.ContainsPoint(P) == true && B.ContainsPoint(P) == true,
/// it is equivalent to having an edge crossing (i.e., Intersects returns true).
class IntersectsRelation: LoopRelation {
  
  var foundSharedVertex = false

  func aCrossingTarget() -> CrossingTarget { return .cross }
  func bCrossingTarget() -> CrossingTarget { return .cross }
  func wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
    foundSharedVertex = true
    return wedgeIntersects(a0: a0, ab1: ab1, a2: a2, b0: b0, b2: b2)
  }

}

// compareBoundaryRelation implements loopRelation for comparing boundaries.
// The compare boundary relation does not have a useful early-exit condition,
// so we return crossingTargetDontCare for both crossing targets.
// Aside: A possible early exit condition could be based on the following.
//   If A contains a point of both B and ~B, then A intersects Boundary(B).
//   If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
//   So if the intersections of {A, ~A} with {B, ~B} are all non-empty,
//   the return value is 0, i.e., Boundary(A) intersects Boundary(B).
// Unfortunately it isn't worth detecting this situation because by the
// time we have seen a point in all four intersection regions, we are also
// guaranteed to have seen at least one pair of crossing edges.
class CompareBoundaryRelation: LoopRelation {
  
  let reverse: Bool // True if the other loop should be reversed.
  // state
  var foundSharedVertex = false // True if any wedge was processed.
  var containsEdge = false // True if any edge of the other loop is contained by this loop.
  var excludesEdge = false // True if any edge of the other loop is excluded by this loop.

  init(reverse: Bool) {
    self.reverse = reverse
  }

  func aCrossingTarget() -> CrossingTarget { return .dontCare }
  func bCrossingTarget() -> CrossingTarget { return .dontCare }

  func wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
    // Because we don't care about the interior of the other, only its boundary,
    // it is sufficient to check whether this one contains the semiwedge (ab1, b2).
    foundSharedVertex = true
    if wedgeContainsSemiwedge(a0: a0, ab1: ab1, a2: a2, b2: b2, reverse: reverse) {
      containsEdge = true
    } else {
      excludesEdge = true
    }
    return containsEdge && excludesEdge
  }

}

// wedgeContainsSemiwedge reports whether the wedge (a0, ab1, a2) contains the
// "semiwedge" defined as any non-empty open set of rays immediately CCW from
// the edge (ab1, b2). If reverse is true, then substitute clockwise for CCW;
// this simulates what would happen if the direction of the other loop was reversed.
func wedgeContainsSemiwedge(a0: S2Point, ab1: S2Point, a2: S2Point, b2: S2Point, reverse: Bool) -> Bool {
  if b2 == a0 || b2 == a2 {
    // We have a shared or reversed edge.
    return (b2 == a0) == reverse
  }
  return S2Point.orderedCCW(a0, a2, b2, ab1)
}

// loopCrosser is a helper type for determining whether two loops cross.
// It is instantiated twice for each pair of loops to be tested, once for the
// pair (A,B) and once for the pair (B,A), in order to be able to process
// edges in either loop nesting order.
struct LoopCrosser {

  let a: S2Loop
  let b: S2Loop
  let relation: LoopRelation
  let swapped: Bool
  let aCrossingTarget: CrossingTarget
  let bCrossingTarget: CrossingTarget
  
  // state maintained by startEdge and edgeCrossesCell.
  var crosser: EdgeCrosser? = nil
  var aj = 0
  var bjPrev = -2
  
  // temporary data declared here to avoid repeated memory allocations.
  var bQuery: CrossingEdgeQuery
  var bCells: [ShapeIndexCell] = []

  /// Creates a loopCrosser from the given values. If swapped is true,
  /// the loops A and B have been swapped. This affects how arguments are passed to
  /// the given loop relation, since for example A.Contains(B) is not the same as
  /// B.Contains(A).
  init(a: S2Loop, b: S2Loop, relation: LoopRelation, swapped: Bool) {
    let act = relation.aCrossingTarget()
    let bct = relation.bCrossingTarget()
    self.a = a
    self.b = b
    self.relation = relation
    self.swapped = swapped
    aCrossingTarget = swapped ? act : bct
    bCrossingTarget = swapped ? bct : act
    bQuery = CrossingEdgeQuery(index: b.index)
  }

  /// Sets the crossers state for checking the given edge of loop A.
  mutating func startEdge(aj: Int) {
    crosser = EdgeCrosser(a: a.vertex(aj), b: a.vertex(aj + 1))
    self.aj = aj
    bjPrev = -2
  }

  /// Reports whether the current edge of loop A has any crossings with
  /// edges of the index cell of loop B.
  mutating func edgeCrossesCell(bClipped: ClippedShape) -> Bool {
    // Test the current edge of A against all edges of bClipped
    let bNumEdges = bClipped.numEdges()
    for j in 0..<bNumEdges {
      let bj = bClipped.edges[j]
      if bj != bjPrev + 1 {
        crosser!.restart(at: b.vertex(bj))
      }
      bjPrev = bj
      let crossing = crosser!.chainCrossingSign(d: b.vertex(bj + 1))
      if crossing == .doNotCross {
        continue
      } else if crossing == .cross {
        return true
      }
      // We only need to check each shared vertex once, so we only
      // consider the case where l.aVertex(l.aj+1) == l.b.Vertex(bj+1).
      if a.vertex(aj + 1) == b.vertex(bj + 1) {
        if swapped {
          if relation.wedgesCross(a0: b.vertex(bj), ab1: b.vertex(bj + 1), a2: b.vertex(bj + 2), b0: a.vertex(aj), b2: a.vertex(aj + 2)) {
            return true
          }
        } else {
          if relation.wedgesCross(a0: a.vertex(aj), ab1: a.vertex(aj + 1), a2: a.vertex(aj + 2), b0: b.vertex(bj), b2: b.vertex(bj + 2)) {
            return true
          }
        }
      }
    }
    return false
  }

  // cellCrossesCell reports whether there are any edge crossings or wedge crossings
  // within the two given cells.
  mutating func cellCrossesCell(aClipped: ClippedShape, bClipped: ClippedShape) -> Bool {
    // Test all edges of aClipped against all edges of bClipped.
    for edge in aClipped.edges {
      startEdge(aj: edge)
      if edgeCrossesCell(bClipped: bClipped) {
        return true
      }
    }
    return false
  }

  // cellCrossesAnySubcell reports whether given an index cell of A, if there are any
  // edge or wedge crossings with any index cell of B contained within bID.
  mutating func cellCrossesAnySubcell(aClipped: ClippedShape, bId: CellId) -> Bool {
    // Test all edges of aClipped against all edges of B. The relevant B
    // edges are guaranteed to be children of bID, which lets us find the
    // correct index cells more efficiently.
    let bRoot = PaddedCell(id: bId, padding: 0)
    for aj in aClipped.edges {
      // Use an CrossingEdgeQuery starting at bRoot to find the index cells
      // of B that might contain crossing edges.
      bCells = bQuery.getCells(a: a.vertex(aj), b: a.vertex(aj+1), root: bRoot)
      if bCells.count == 0 {
        continue
      }
      startEdge(aj: aj)
      for c in 0..<bCells.count {
        if edgeCrossesCell(bClipped: bCells[c].shapes[0]) {
          return true
        }
      }
    }
    return false
  }

  // hasCrossing reports whether given two iterators positioned such that
  // ai.cellID().ContainsCellID(bi.cellID()), there is an edge or wedge crossing
  // anywhere within ai.cellID(). This function advances bi only past ai.cellID().
  mutating func hasCrossing(ai: RangeIterator, bi: inout RangeIterator) -> Bool {
    // If ai.CellID() intersects many edges of B, then it is faster to use
    // CrossingEdgeQuery to narrow down the candidates. But if it intersects
    // only a few edges, it is faster to check all the crossings directly.
    // We handle this by advancing bi and keeping track of how many edges we
    // would need to test.
    let edgeQueryMinEdges = 20 // Tuned from benchmarks.
    var totalEdges = 0
    bCells = []
    while true {
      let n = bi.it.indexCell()!.shapes[0].numEdges()
      if n > 0 {
        totalEdges += n
        if totalEdges >= edgeQueryMinEdges {
          // There are too many edges to test them directly, so use CrossingEdgeQuery.
          if cellCrossesAnySubcell(aClipped: ai.it.indexCell()!.shapes[0], bId: ai.cellId()) {
            return true
          }
          bi.seek(beyond: ai)
          return false
        }
        if let indexCell = bi.indexCell() {
          bCells.append(indexCell)
        }
      }
      bi.next()
      if bi.cellId() > ai.rangeMax {
        break
      }
    }
    // Test all the edge crossings directly.
    for c in bCells {
      if cellCrossesCell(aClipped: ai.it.indexCell()!.shapes[0], bClipped: c.shapes[0]) {
        return true
      }
    }
    return false
  }

  /// Reports if the clippedShapes containsCenter boolean corresponds
  /// to the crossing target type given. (This is to work around C++ allowing false == 0,
  /// true == 1 type implicit conversions and comparisons)
  static func containsCenterMatches(a: ClippedShape, target: CrossingTarget) -> Bool {
    return (!a.containsCenter && target == .dontCross) || (a.containsCenter && target == .cross)
  }

  /// Reports whether given two iterators positioned such that
  /// ai.cellID().ContainsCellID(bi.cellID()), there is a crossing relationship
  /// anywhere within ai.cellID(). Specifically, this method returns true if there
  /// is an edge crossing, a wedge crossing, or a point P that matches both relations
  /// crossing targets. This function advances both iterators past ai.cellID.
  mutating func hasCrossingRelation(ai: inout RangeIterator, bi: inout RangeIterator) -> Bool {
    let aClipped = ai.it.indexCell()!.shapes[0]
    if aClipped.numEdges() != 0 {
      // The current cell of A has at least one edge, so check for crossings.
      if hasCrossing(ai: ai, bi: &bi) {
        return true
      }
      ai.next()
      return false
    }
    if LoopCrosser.containsCenterMatches(a: aClipped, target: aCrossingTarget) {
      // The crossing target for A is not satisfied, so we skip over these cells of B.
      bi.seek(beyond: ai)
      ai.next()
      return false
    }
    // All points within ai.cellID() satisfy the crossing target for A, so it's
    // worth iterating through the cells of B to see whether any cell
    // centers also satisfy the crossing target for B.
    while bi.cellId() <= ai.rangeMax {
      let bClipped = bi.it.indexCell()!.shapes[0]
      if LoopCrosser.containsCenterMatches(a: bClipped, target: bCrossingTarget) {
        return true
      }
      bi.next()
    }
    ai.next()
    return false
  }

}

// hasCrossingRelation checks all edges of loop A for intersection against all edges
// of loop B and reports if there are any that satisfy the given relation. If there
// is any shared vertex, the wedges centered at this vertex are sent to the given
// relation to be tested.
//
// If the two loop boundaries cross, this method is guaranteed to return
// true. It also returns true in certain cases if the loop relationship is
// equivalent to crossing. For example, if the relation is Contains and a
// point P is found such that B contains P but A does not contain P, this
// method will return true to indicate that the result is the same as though
// a pair of crossing edges were found (since Contains returns false in
// both cases).
//
// See Contains, Intersects and CompareBoundary for the three uses of this function.
func hasCrossingRelation(a: S2Loop, b: S2Loop, relation: LoopRelation) -> Bool {
  // We look for CellID ranges where the indexes of A and B overlap, and
  // then test those edges for crossings.
  var ai = RangeIterator(index: a.index)
  var bi = RangeIterator(index: b.index)
  var ab = LoopCrosser(a: a, b: b, relation: relation, swapped: false) // Tests edges of A against B
  var ba = LoopCrosser(a: b, b: a, relation: relation, swapped: true)  // Tests edges of B against A
  while !ai.done() || !bi.done() {
    if ai.rangeMax < bi.rangeMin {
      // The A and B cells don't overlap, and A precedes B.
      ai.seek(to: bi)
    } else if bi.rangeMax < ai.rangeMin {
      // The A and B cells don't overlap, and B precedes A.
      bi.seek(to: ai)
    } else {
      // One cell contains the other. Determine which cell is larger.
      let abRelation = Int64(ai.it.cellId().lsb() - bi.it.cellId().lsb())
      if abRelation > 0 {
        // A's index cell is larger.
        if ab.hasCrossingRelation(ai: &ai, bi: &bi) {
          return true
        }
      } else if abRelation < 0 {
        // B's index cell is larger.
        if ba.hasCrossingRelation(ai: &bi, bi: &ai) {
          return true
        }
      } else {
        // The A and B cells are the same. Since the two cells
        // have the same center point P, check whether P satisfies
        // the crossing targets.
        let aClipped = ai.it.indexCell()!.shapes[0]
        let bClipped = bi.it.indexCell()!.shapes[0]
        if LoopCrosser.containsCenterMatches(a: aClipped, target: ab.aCrossingTarget) && LoopCrosser.containsCenterMatches(a: bClipped, target: ab.bCrossingTarget) {
          return true
        }
        // Otherwise test all the edge crossings directly.
        if aClipped.numEdges() > 0 && bClipped.numEdges() > 0 && ab.cellCrossesCell(aClipped: aClipped, bClipped: bClipped) {
          return true
        }
        ai.next()
        bi.next()
      }
    }
  }
  return false
}
