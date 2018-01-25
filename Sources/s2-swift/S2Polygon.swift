//
//  S2Polygon.swift
//  s2-swift
//

import Foundation


/// Represents a sequence of zero or more loops; recall that the
/// interior of a loop is defined to be its left-hand side (see Loop).
/// When the polygon is initialized, the given loops are automatically converted
/// into a canonical form consisting of "shells" and "holes".  Shells and holes
/// are both oriented CCW, and are nested hierarchically.  The loops are
/// reordered to correspond to a preorder traversal of the nesting hierarchy.
/// Polygons may represent any region of the sphere with a polygonal boundary,
/// including the entire sphere (known as the "full" polygon).  The full polygon
/// consists of a single full loop (see Loop), whereas the empty polygon has no
/// loops at all.
/// Use FullPolygon() to construct a full polygon. The zero value of Polygon is
/// treated as the empty polygon.
/// Polygons have the following restrictions:
///  - Loops may not cross, i.e. the boundary of a loop may not intersect
///    both the interior and exterior of any other loop.
///  - Loops may not share edges, i.e. if a loop contains an edge AB, then
///    no other loop may contain AB or BA.
///  - Loops may share vertices, however no vertex may appear twice in a
///    single loop (see Loop).
///  - No loop may be empty.  The full loop may appear only in the full polygon.
public struct S2Polygon: S2RegionType, Shape {
  
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

  // numEdges tracks the total number of edges in all the loops in this polygon.
  let numEdges: Int
  
  // A slice where element i is the cumulative number of edges in the
  // preceding loops in the polygon. This field is used for polygons that
  // have a large number of loops, and may be empty for polygons with few loops.
  let cumulativeEdges: [Int]

  // MARK: inits / factory
  
  /// Performs the shape related initializations and adds the final
  /// polygon to the index.
  private static func computeEdges(loops: [S2Loop]) -> ([Int], Int) {
    // check if full
    if loops.count == 1 && loops[0].isFull {
      return ([], 0)
    }
    //    let maxLinearSearchLoops = 12 // Based on benchmarks.
    //    if loops.count > maxLinearSearchLoops {
    //      cumulativeEdges = []
    //    }
    var cumulativeEdges: [Int] = []
    var numEdges = 0
    for l in loops {
      cumulativeEdges.append(numEdges)
      numEdges += l.vertices.count
    }
    return (cumulativeEdges, numEdges)
  }
  
  private init(loops: [S2Loop], hasHoles: Bool, numVertices: Int, bound: S2Rect, subregionBound: S2Rect) {
    self.loops = loops
    self.hasHoles = hasHoles
    self.numVertices = numVertices
    self.bound = bound
    self.subregionBound = subregionBound
    //
    let (cumulativeEdges, numEdges) = S2Polygon.computeEdges(loops: loops)
    self.cumulativeEdges = cumulativeEdges
    self.numEdges = numEdges
    // has its own shapeindex for recursive loop thingies (?)
    self.index = ShapeIndex()
    self.index.add(shape: self)
  }
  
  /// Returns a special "full" polygon.
  static let full = S2Polygon(loops: [S2Loop.full], hasHoles: false, numVertices: S2Loop.full.vertices.count, bound:S2Rect.full, subregionBound: S2Rect.full)
  static let empty = S2Polygon(loops: [], hasHoles: false, numVertices: 0, bound:S2Rect.empty, subregionBound: S2Rect.empty)
  
  /// Set the properties for a polygon made of a single loop.
  init(loop: S2Loop) {
    guard loop.depth == 0 else { fatalError("nested loops not supported") }
    self.init(loops: [loop], hasHoles: false, numVertices: loop.vertices.count, bound: loop.bound, subregionBound: loop.subregionBound)
  }

  /// Constructs a polygon from the given set of loops. The polygon
  /// interior consists of the points contained by an odd number of loops. (Recall
  /// that a loop contains the set of points on its left-hand side.)
  ///
  /// This method determines the loop nesting hierarchy and assigns every loop a
  /// depth. Shells have even depths, and holes have odd depths.
  ///
  /// Note: The given set of loops are reordered by this method so that the hierarchy
  /// can be traversed using Parent, LastDescendant and the loops depths.
  init(loops: [S2Loop]) {
    if loops.count == 0 {
      self = S2Polygon.empty
      return
    }
    // simple
    if loops.count == 1 {
      self.init(loop: loops[0])
      return
    }
    // take the set of loops in this polygon and performs the nesting
    // computations to set the proper nesting and parent/child relationships.
    let lm = LoopMap(loops: loops)
    // Reorder the loops in depth-first traversal order.
    let loops2 = S2Polygon.computeLoops(lm)
    let (hasHoles, numVertices, bound, subregionBound) = S2Polygon.computeLoopProperties(loops: loops2)
    self.init(loops: loops2, hasHoles: hasHoles, numVertices: numVertices, bound: bound, subregionBound: subregionBound)
  }
  
  /// A map of a loop to its immediate children with respect to nesting.
  /// It is used to determine which loops are shells and which are holes.
  struct LoopMap {

    private let loopMap: [Int: [(Int, S2Loop)]]

    init(loops: [S2Loop]) {
      var loopMap: [Int: [(Int, S2Loop)]] = [:]
      var parent = (-1, S2Loop.full)
      for (index, loop) in loops.enumerated() {
        var children: [(Int, S2Loop)] = []
        // find parent by narrowing down
        while true {
          children = loopMap[parent.0] ?? []
          if let child = children.first(where: { $0.1.containsNested(other: loop) }) {
            parent = child
          } else {
            break
          }
        }
        // Now, we have found a parent for this loop, it may be that some of the
        // children of the parent of this loop may now be children of the new loop.
        var newChildren = children
        var parentsChildren = [(index, loop)]
        for child in children {
          if loop.containsNested(other: child.1) {
            newChildren.append(child)
          } else {
            parentsChildren.append(child)
          }
        }
        loopMap[index] = newChildren
        loopMap[parent.0] = parentsChildren
      }
      self.loopMap = loopMap
    }
    
    subscript(i: Int) -> [(Int, S2Loop)] {
      guard let children = loopMap[i] else { return [] }
      return children
    }
    
  }

  /// Like S2Polygon(loops), returns a Polygon from the
  /// given set of loops. It expects loops to be oriented such that the polygon
  /// interior is on the left-hand side of all loops. This implies that shells
  /// and holes should have opposite orientations in the input to this method.
  /// (During initialization, loops representing holes will automatically be
  /// inverted.)
  init(orientedLoops: [S2Loop]) {
    self.init(loops: loops)
    initOriented()
  }
  
  /// Takes the loops in this polygon and performs the nesting
  /// computations. It expects the loops to be oriented such that the polygon
  /// interior is on the left-hand side of all loops. This implies that shells
  /// and holes should have opposite orientations in the input to this method.
  /// (During initialization, loops representing holes will automatically be
  /// inverted.)
  func initOriented() {
    // Here is the algorithm:
    //
    // 1. Remember which of the given loops contain OriginPoint.
    //
    // 2. Invert loops as necessary to ensure that they are nestable (i.e., no
    //    loop contains the complement of any other loop). This may result in a
    //    set of loops corresponding to the complement of the given polygon, but
    //    we will fix that problem later.
    //
    //    We make the loops nestable by first normalizing all the loops (i.e.,
    //    inverting any loops whose turning angle is negative). This handles
    //    all loops except those whose turning angle is very close to zero
    //    (within the maximum error tolerance). Any such loops are inverted if
    //    and only if they contain OriginPoint(). (In theory this step is only
    //    necessary if there are at least two such loops.) The resulting set of
    //    loops is guaranteed to be nestable.
    //
    // 3. Build the polygon. This yields either the desired polygon or its
    //    complement.
    //
    // 4. If there is at least one loop, we find a loop L that is adjacent to
    //    OriginPoint() (where "adjacent" means that there exists a path
    //    connecting OriginPoint() to some vertex of L such that the path does
    //    not cross any loop). There may be a single such adjacent loop, or
    //    there may be several (in which case they should all have the same
    //    contains_origin() value). We choose L to be the loop containing the
    //    origin whose depth is greatest, or loop(0) (a top-level shell) if no
    //    such loop exists.
    //
    // 5. If (L originally contained origin) != (polygon contains origin), we
    //    invert the polygon. This is done by inverting a top-level shell whose
    //    turning angle is minimal and then fixing the nesting hierarchy. Note
    //    that because we normalized all the loops initially, this step is only
    //    necessary if the polygon requires at least one non-normalized loop to
    //    represent it.
    fatalError("initOriented not yet implemented")
  }
  
  // loopStack simplifies access to the loops while being initialized.
  struct LoopStack {
    
    var loops: [S2Loop] = []
  
    mutating func push(_ v: S2Loop) {
      loops.append(v)
    }
  
    mutating func pop() -> S2Loop? {
      return loops.popLast()
    }
  
    var isEmpty: Bool { return loops.count == 0 }
    
  }
  
  /// Walks the mapping of loops to all of their children, and adds them in
  /// order into to the polygons set of loops.
  func computeLoops(lm: LoopMap) -> [S2Loop] {
    var loops: [S2Loop] = []
    var stack = LoopStack()
    stack.push(S2Loop.full)
    var depth = -1
    while !stack.isEmpty {
      if let loop = stack.pop() {
        depth = loop.depth
        loops.append(loop)
      }
      var children = lm[loop]
      for j in 0..<children.count {
        let i = children.count - 1 - j
        let child = children[i]
        child.depth = depth + 1
        stack.push(child)
      }
    }
    return loops
  }
  
  /// Returns a Polygon from a single loop created from the given Cell.
  init(cell: Cell) {
    let loop = S2Loop(cell: cell)
    self.init(loops: [loop])
  }
  
  // MARK: tests
  
  /// Reports whether this is the special "empty" polygon (consisting of no loops).
  var isEmpty: Bool {
    return loops.count == 0
  }

  /// Reports whether this is the special "full" polygon (consisting of a
  /// single loop that encompasses the entire sphere).
  var isFull: Bool {
    return loops.count == 1 && loops[0].isFull
  }

  /// Returns a bounding spherical cap.
  public func capBound() -> S2Cap {
    return bound.capBound()
  }

  /// Returns a bounding latitude-longitude rectangle.
  public func rectBound() -> S2Rect {
    return bound
  }

  // Reports whether the polygon contains the given cell.
  public func contains(_ cell: Cell) -> Bool {
    return loops.count == 1 && loops[0].contains(cell)
  }
  
  /// Reports whether the polygon intersects the given cell.
  public func intersects(_ cell: Cell) -> Bool {
    return loops.count == 1 && loops[0].intersects(cell)
  }

  /// Sets the properties for polygons with multiple loops.
  static func initLoopProperties(loops: [S2Loop]) -> (Bool, Int, S2Rect, S2Rect) {
    // the loops depths are set by initNested/initOriented prior to this.
    var hasHoles = false
    var numVertices = 0
    var bound = S2Rect.empty
    for l in loops {
      if l.isHole() {
        hasHoles = true
      } else {
        bound = bound.union(l.rectBound())
      }
      numVertices += l.numVertices()
    }
    var subregionBound = expandForSubregions(bound)
    let (cumulativeEdges, numEdges) = S2Polygon.initEdgesAndIndex(loops: loops)
    return (hasHoles, numVertices, bound, subregionBound)
  }
  
  /// Returns a special "full" polygon.
  static func fullPolygon() -> S2Polygon {
    let full = S2Loop.full
    let (cumulativeEdges, numEdges) = computeEdges(loops: [full])
    var p = S2Polygon(loops: [full], numVertices: full.vertices.count, bound: S2Rect.full, subregionBound: S2Rect.full)
    return p
  }
  
  /// Returns the number of loops in this polygon.
  func numLoops() -> Int {
    return loops.count
  }
  
  /// Returns the loop at the given index. Note that during initialization,
  /// the given loops are reordered according to a pre-order traversal of the loop
  /// nesting hierarchy. This implies that every loop is immediately followed by
  /// its descendants. This hierarchy can be traversed using the methods Parent,
  /// LastDescendant, and Loop.depth.
  func loop(_ k: Int) -> S2Loop {
    return loops[k]
  }
  
  /// Returns the index of the parent of loop k.
  /// If the loop does not have a parent, ok=false is returned.
  func parent(k: Int) -> Int? {
    // See where we are on the depth hierarchy.
    let depth = loops[k].depth
    if depth == 0 {
      return nil
    }
    // There may be several loops at the same nesting level as us that share a
    // parent loop with us. (Imagine a slice of swiss cheese, of which we are one loop.
    // we don't know how many may be next to us before we get back to our parent loop.)
    // Move up one position from us, and then begin traversing back through the set of loops
    // until we find the one that is our parent or we get to the top of the polygon.
    var k = k - 1
    while k >= 0 && loops[k].depth <= depth {
      k -= 1
    }
    return k
  }
 
  /// Returns the index of the last loop that is contained within loop k.
  /// If k is negative, it returns the last loop in the polygon.
  /// Note that loops are indexed according to a pre-order traversal of the nesting
  /// hierarchy, so the immediate children of loop k can be found by iterating over
  /// the loops (k+1)..LastDescendant(k) and selecting those whose depth is equal
  /// to Loop(k).depth+1.
  func lastDescendant(k: Int) -> Int {
    if k < 0 {
      return loops.count - 1
    }
    let depth = loops[k].depth
    // Find the next loop immediately past us in the set of loops, and then start
    // moving down the list until we either get to the end or find the next loop
    // that is higher up the hierarchy than we are.
    var k = k + 1
    while k < loops.count && loops[k].depth > depth {
      k += 1
    }
    return k - 1
  }
  
  /// Reports whether the polygon contains the point.
  func contains(_ point: S2Point) -> Bool {
    // NOTE: A bounds check slows down this function by about 50%. It is
    // worthwhile only when it might allow us to delay building the index.
    if !index.isFresh() && !bound.contains(point) {
      return false
    }
    // For small polygons, and during initial construction, it is faster to just
    // check all the crossing.
    let maxBruteForceVertices = 32
    if numVertices < maxBruteForceVertices || index == nil {
      var inside = false
      for l in loops {
        // use loops bruteforce to avoid building the index on each loop.
        inside = inside != l.bruteForceContains(point)
      }
      return inside
    }
    // Otherwise, look up the ShapeIndex cell containing this point.
    var it = index.iterator()
    if !it.locate(point: point) {
      return false
    }
    return iteratorContains(iterator: it, point: point)
  }
  
  /// Reports whether the polygon contains the given cell.
  func contains(_ cell: Cell) -> Bool {
    var it = index.iterator()
    var relation = it.locate(cellId: cell.id)
    
    // If "cell" is disjoint from all index cells, it is not contained.
    // Similarly, if "cell" is subdivided into one or more index cells then it
    // is not contained, since index cells are subdivided only if they (nearly)
    // intersect a sufficient number of edges.  (But note that if "cell" itself
    // is an index cell then it may be contained, since it could be a cell with
    // no edges in the loop interior.)
    if relation != .indexed {
      return false
    }
    // Otherwise check if any edges intersect "cell".
    if boundaryApproxIntersects(iterator: it, cell: cell) {
      return false
    }
    // Otherwise check if the loop contains the center of "cell".
    return iteratorContains(iterator: it, point: cell.center())
  }
  
  /// Reports whether the polygon intersects the given cell.
  func intersects(_ cell: Cell) -> Bool {
    var it = index.iterator()
    var relation = it.locate(cellId: cell.id)
    // If cell does not overlap any index cell, there is no intersection.
    if relation == .disjoint {
      return false
    }
    // If cell is subdivided into one or more index cells, there is an
    // intersection to within the S2ShapeIndex error bound (see Contains).
    if relation == .subdivided {
      return true
    }
    // If cell is an index cell, there is an intersection because index cells
    // are created only if they have at least one edge or they are entirely
    // contained by the loop.
    if it.cellId() == cell.id {
      return true
    }
    // Otherwise check if any edges intersect cell.
    if boundaryApproxIntersects(iterator: it, cell: cell) {
      return true
    }
    // Otherwise check if the loop contains the center of cell.
    return iteratorContains(iterator: it, point: cell.center())
  }
  
  /// Computes a covering of the Polygon.
  func cellUnionBound() -> CellUnion {
    // TODO(roberts): Use ShapeIndexRegion when it's available.
    return capBound().cellUnionBound()
  }
  
  /// Reports whether the loop's boundary intersects cell.
  /// It may also return true when the loop boundary does not intersect cell but
  /// some edge comes within the worst-case error tolerance.
  /// This requires that it.Locate(cell) returned Indexed.
  func boundaryApproxIntersects(iterator it: ShapeIndexIterator, cell: Cell) -> Bool {
    guard let aClipped = it.indexCell()?.find(shapeId: 0) else { return false }
    // If there are no edges, there is no intersection.
    if aClipped.edges.count == 0 {
      return false
    }
    // We can save some work if cell is the index cell itself.
    if it.cellId() == cell.id {
      return true
    }
    // Otherwise check whether any of the edges intersect cell.
    let maxError = faceClipErrorUVCoord + intersectsRectErrorUVDist
    let bound = cell.boundUV().expanded(maxError)
    for e in aClipped.edges {
      if let edge = index.shape(id: 0)?.edge(e) {
        let (v0, v1, ok) = clipToPaddedFace(a: edge.v0, b: edge.v1, f: Int(cell.face), padding: maxError)
        if ok && edgeIntersectsRect(v0, v1, bound) {
          return true
        }
      }
    }
    return false
  }
  
  /// Reports whether the iterator that is positioned at the
  /// ShapeIndexCell that may contain p, contains the point p.
  func iteratorContains(iterator it: ShapeIndexIterator, point: S2Point) -> Bool {
    // Test containment by drawing a line segment from the cell center to the
    // given point and counting edge crossings.
    guard let aClipped = it.indexCell()?.find(shapeId: 0) else { return false }
    var inside = aClipped.containsCenter
    if aClipped.edges.count == 0 {
      return inside
    }
    // This block requires ShapeIndex.
    var crosser = EdgeCrosser(a: it.center(), b: point)
    if let shape = index.shape(id: 0) {
      for e in aClipped.edges {
        let edge = shape.edge(e)
        inside = inside != crosser.isEdgeOrVertexCrossing(c: edge.v0, d: edge.v1)
      }
    }
    return inside
  }
  
  // MARK: Shape Interface
  
  /// Returns endpoints for the given edge index.
  public func edge(_ e: Int) -> Edge {
    var e = e
    if cumulativeEdges.count > 0 {
      for i in cumulativeEdges {
        if i + 1 >= cumulativeEdges.count || e < cumulativeEdges[i + 1] {
          e -= cumulativeEdges[i]
          break
        }
      }
    } else {
      // When the number of loops is small, use linear search. Most often
      // there is exactly one loop and the code below executes zero times.
      var i = 0
      while e >= loop(i).vertices.count {
        e -= loop(i).vertices.count
        i += 1
      }
    }
    return edge(v0: loop(i).orientedVertex(e), v1: loop(i).orientedVertex(e + 1))
  }
  
  //// Reports whether this Polygon has an interior.
  public func hasInterior() -> Bool {
    return dimension() == .polygonGeometry
  }
  
  /// Returns the reference point for this polygon.
  public func referencePoint() -> ReferencePoint {
    var containsOrigin = false
    for l in loops {
      containsOrigin = containsOrigin != l.containsOrigin()
    }
    return ReferencePoint(origin: true, contained: containsOrigin)
  }
  
  /// Reports the number of contiguous edge chains in the Polygon.
  public func numChains() -> Int {
    if isFull {
      return 0
    }
    return numLoops()
  }
  
  /// Returns the i-th edge Chain (loop) in the Shape.
  public func chain(_ chainId: Int) -> Chain {
    if cumulativeEdges != nil {
      return Chain(start: cumulativeEdges[chainId], length: loop(chainId).vertices.count)
    }
    var e = 0
    for j in 0..<chainId {
      e += loop(j).vertices.count
    }
    return Chain(start: e, length: loop(chainId).vertices.count)
  }
  
  /// Returns the j-th edge of the i-th edge Chain (loop).
  public func chainEdge(chainId i: Int, offset j: Int) -> Edge {
    return Edge(v0: loop(i).orientedVertex(j), v1: loop(i).orientedVertex(j + 1))
  }
  
  /// Returns a pair (i, j) such that edgeID is the j-th edge
  /// of the i-th edge Chain.
  public func chainPosition(_ edgeId: Int) -> ChainPosition {
    var edgeId = edgeId
    var i = 0
    if cumulativeEdges.count > 0 {
      for i in cumulativeEdges {
        if i + 1 >= cumulativeEdges.count || edgeId < cumulativeEdges[i + 1] {
          edgeId -= cumulativeEdges[i]
          break
        }
      }
    } else {
      // When the number of loops is small, use linear search. Most often
      // there is exactly one loop and the code below executes zero times.
      i = 0
      while edgeId >= loop(i).vertices.count {
        edgeId -= loop(i).vertices.count
        i += 1
      }
    }
    // TODO(roberts): unify this and Edge since they are mostly identical.
    return ChainPosition(chainId: i, offset: edgeId)
  }
  
  /// Returns the dimension of the geometry represented by this Polygon.
  public func dimension() -> ShapeDimension {
    return .polygonGeometry
  }
  
  /// Reports whether this polygon contains the other polygon.
  /// Specifically, it reports whether all the points in the other polygon
  /// are also in this polygon.
  func contains(_ o: S2Polygon) -> Bool {
    // If both polygons have one loop, use the more efficient Loop method.
    // Note that Loop's Contains does its own bounding rectangle check.
    if loops.count == 1 && o.loops.count == 1 {
      return loops[0].contains(o.loops[0])
    }
    // Otherwise if neither polygon has holes, we can still use the more
    // efficient Loop's Contains method (rather than compareBoundary),
    // but it's worthwhile to do our own bounds check first.
    if !subregionBound.contains(o.bound) {
      // Even though Bound(A) does not contain Bound(B), it is still possible
      // that A contains B. This can only happen when union of the two bounds
      // spans all longitudes. For example, suppose that B consists of two
      // shells with a longitude gap between them, while A consists of one shell
      // that surrounds both shells of B but goes the other way around the
      // sphere (so that it does not intersect the longitude gap).
      if !bound.lng.union(o.bound.lng).isFull {
        return false
      }
    }
    if !hasHoles && !o.hasHoles {
      for l in o.loops {
        if !anyLoopContains(l) {
          return false
        }
      }
      return true
    }
    // Polygon A contains B iff B does not intersect the complement of A. From
    // the intersection algorithm below, this means that the complement of A
    // must exclude the entire boundary of B, and B must exclude all shell
    // boundaries of the complement of A. (It can be shown that B must then
    // exclude the entire boundary of the complement of A.) The first call
    // below returns false if the boundaries cross, therefore the second call
    // does not need to check for any crossing edges (which makes it cheaper).
    return containsBoundary(o) && o.excludesNonCrossingComplementShells(self)
  }
  
  /// Reports whether this polygon intersects the other polygon, i.e.
  /// if there is a point that is contained by both polygons.
  func intersects(_ o: S2Polygon) -> Bool {
    // If both polygons have one loop, use the more efficient Loop method.
    // Note that Loop Intersects does its own bounding rectangle check.
    if loops.count == 1 && o.loops.count == 1 {
      return loops[0].intersects(o.loops[0])
    }
    // Otherwise if neither polygon has holes, we can still use the more
    // efficient Loop.Intersects method. The polygons intersect if and
    // only if some pair of loop regions intersect.
    if !bound.intersects(o.bound) {
      return false
    }
    if !hasHoles && !o.hasHoles {
      for l in o.loops {
        if anyLoopIntersects(l) {
          return true
        }
      }
      return false
    }
    // Polygon A is disjoint from B if A excludes the entire boundary of B and B
    // excludes all shell boundaries of A. (It can be shown that B must then
    // exclude the entire boundary of A.) The first call below returns false if
    // the boundaries cross, therefore the second call does not need to check
    // for crossing edges.
    return !excludesBoundary(o) || !o.excludesNonCrossingShells(self)
  }
  
  /// Returns +1 if this polygon contains the boundary of B, -1 if A
  /// excludes the boundary of B, and 0 if the boundaries of A and B cross.
  func compareBoundary(_ o: S2Loop) -> Int {
    var result = -1
    var i = 0
    while i < loops.count && result != 0 {
      // If B crosses any loop of A, the result is 0. Otherwise the result
      // changes sign each time B is contained by a loop of A.
      result *= -loops[i].compareBoundary(o: o)
      i += 1
    }
    return result
  }
  
  /// Reports whether this polygon contains the entire boundary of B.
  func containsBoundary(_ o: S2Polygon) -> Bool {
    for l in o.loops {
      if compareBoundary(l) <= 0 {
        return false
      }
    }
    return true
  }
  
  /// Reports whether this polygon excludes the entire boundary of B.
  func excludesBoundary(_ o: S2Polygon) -> Bool {
    for l in o.loops {
      if compareBoundary(l) >= 0 {
        return false
      }
    }
    return true
  }
  
  /// Reports whether polygon A contains the boundary of
  /// loop B. Shared edges are handled according to the rule described in loops
  /// containsNonCrossingBoundary.
  func containsNonCrossingBoundary(_ o: S2Loop, reverse: Bool) -> Bool {
    var inside = false
    for l in loops {
      let x = l.containsNonCrossingBoundary(o, reverse: reverse)
      inside = (inside != x)
    }
    return inside
  }
  
  /// Reports wheterh given two polygons A and B such that the
  /// boundary of A does not cross any loop of B, if A excludes all shell boundaries of B.
  func excludesNonCrossingShells(_ o: S2Polygon) -> Bool {
    for l in o.loops {
      if l.isHole() {
        continue
      }
      if containsNonCrossingBoundary(l, reverse: false) {
        return false
      }
    }
    return true
  }
  
  /// Reports whether given two polygons A and B
  /// such that the boundary of A does not cross any loop of B, if A excludes all
  /// shell boundaries of the complement of B.
  func excludesNonCrossingComplementShells(_ o: S2Polygon) -> Bool {
    // Special case to handle the complement of the empty or full polygons.
    if o.isEmpty {
      return !isFull
    }
    if o.isFull {
      return true
    }
    // Otherwise the complement of B may be obtained by inverting loop(0) and
    // then swapping the shell/hole status of all other loops. This implies
    // that the shells of the complement consist of loop 0 plus all the holes of
    // the original polygon.
    for (j, l) in o.loops.enumerated() {
      if j > 0 && !l.isHole() {
        continue
      }
      // The interior of the complement is to the right of loop 0, and to the
      // left of the loops that were originally holes.
      if containsNonCrossingBoundary(l, reverse: j == 0) {
        return false
      }
    }
    return true
  }
  
  /// Reports whether any loop in this polygon contains the given loop.
  func anyLoopContains(_ o: S2Loop) -> Bool {
    for l in loops {
      if contains(o) {
        return true
      }
    }
    return false
  }
  
  /// Reports whether any loop in this polygon intersects the given loop.
  func anyLoopIntersects(_ o: S2Loop) -> Bool {
    for l in loops {
      if intersects(o) {
        return true
      }
    }
    return false
  }

}
