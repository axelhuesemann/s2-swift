//
//  S2ShapeIndex.swift
//  s2-swift
//

import Foundation


/// Describes the possible relationships between a target cell
/// and the cells of the ShapeIndex. If the target is an index cell or is
/// contained by an index cell, it is Indexed. If the target is subdivided
/// into one or more index cells, it is Subdivided. Otherwise it is Disjoint.
enum CellRelation: Int {
  case indexed = 0
  case subdivided = 1
  case disjoint = 2
}

// cellPadding defines the total error when clipping an edge which comes
// from two sources:
// (1) Clipping the original spherical edge to a cube face (the face edge).
//     The maximum error in this step is faceClipErrorUVCoord.
// (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
//     The maximum error in this step is edgeClipErrorUVCoord.
// Finally, since we encounter the same errors when clipping query edges, we
// double the total error so that we only need to pad edges during indexing
// and not at query time.
let cellPadding = 2.0 * (faceClipErrorUVCoord + edgeClipErrorUVCoord)

// cellSizeToLongEdgeRatio defines the cell size relative to the length of an
// edge at which it is first considered to be long. Long edges do not
// contribute toward the decision to subdivide a cell further. For example,
// a value of 2.0 means that the cell must be at least twice the size of the
// edge in order for that edge to be counted. There are two reasons for not
// counting long edges: (1) such edges typically need to be propagated to
// several children, which increases time and memory costs without much benefit,
// and (2) in pathological cases, many long edges close together could force
// subdivision to continue all the way to the leaf cell level.
let cellSizeToLongEdgeRatio = 1.0

/// Represents the part of a shape that intersects a Cell.
/// It consists of the set of edge IDs that intersect that cell and a boolean
/// indicating whether the center of the cell is inside the shape (for shapes
/// that have an interior).
///
/// Note that the edges themselves are not clipped; we always use the original
/// edges for intersection tests so that the results will be the same as the
/// original shape.
struct ClippedShape {
  /// shapeId is the index of the shape this clipped shape is a part of.
  let shapeId: Int32
  /// containsCenter indicates if the center of the CellID this shape has been
  /// clipped to falls inside this shape. This is false for shapes that do not
  /// have an interior.
  let containsCenter: Bool
  /// edges is the ordered set of ShapeIndex original edge IDs. Edges
  /// are stored in increasing order of edge ID.
  let edges: [Int]
}

extension ClippedShape {
  
  /// Returns the number of edges that intersect the CellID of the Cell this was clipped to.
  func numEdges() -> Int {
    return edges.count
  }
  
  // Reports if this clipped shape contains the given edge ID.
  func containsEdge(id: Int) -> Bool {
    // Linear search is fast because the number of edges per shape is typically
    // very small (less than 10).
    for e in edges {
      if e == id {
        return true
      }
    }
    return false
  }

}

extension ClippedShape: Equatable {
  
  public static func ==(lhs: ClippedShape, rhs: ClippedShape) -> Bool {
    return lhs.shapeId == rhs.shapeId && lhs.edges == rhs.edges
  }
  
}

/// Stores the index contents for a particular CellId.
struct ShapeIndexCell {
  var shapes: [ClippedShape]
}

extension ShapeIndexCell {
  
  // Creates a new cell that is sized to hold the given number of shapes.
  init(numShapes: Int) {
    self.init(shapes: [ClippedShape]())
  }

  // numEdges reports the total number of edges in all clipped shapes in this cell.
  func numEdges() -> Int {
    return shapes.reduce(0) { $0 + $1.numEdges() }
  }

  // add adds the given clipped shape to this index cell.
  mutating func add(_ c: ClippedShape) {
    shapes.append(c)
  }

  // findByShapeID returns the clipped shape that contains the given shapeID,
  // or nil if none of the clipped shapes contain it.
  func find(shapeId: Int32) -> ClippedShape? {
    // Linear search is fine because the number of shapes per cell is typically
    // very small (most often 1), and is large only for pathological inputs
    // (e.g. very deeply nested loops).
    return shapes.first(where: { $0.shapeId == shapeId })
  }

}

// faceEdge and clippedEdge store temporary edge data while the index is being
// updated.
//
// While it would be possible to combine all the edge information into one
// structure, there are two good reasons for separating it:
//
//  - Memory usage. Separating the two means that we only need to
//    store one copy of the per-face data no matter how many times an edge is
//    subdivided, and it also lets us delay computing bounding boxes until
//    they are needed for processing each face (when the dataset spans
//    multiple faces).
//
//  - Performance. UpdateEdges is significantly faster on large polygons when
//    the data is separated, because it often only needs to access the data in
//    clippedEdge and this data is cached more successfully.

/// Represents an edge that has been projected onto a given face,
struct FaceEdge {
  let shapeId: Int32 // The ID of shape that this edge belongs to
  let edgeId: Int // Edge ID within that shape
  let maxLevel: Int // Not desirable to subdivide this edge beyond this level
  let hasInterior: Bool // Belongs to a shape that has an interior
  let edge: Edge // The original edge.
  var a: R2Point // The edge endpoints, clipped to a given face
  var b: R2Point
}

//extension FaceEdge {
//  init(shapeId: Int32, hasInterior: Bool) {
//
//  }
//}

/// Represents the portion of that edge that has been clipped to a given Cell.
struct ClippedEdge {
  let faceEdge: FaceEdge // The original unclipped edge
  let bound: R2Rect   // Bounding box for the clipped portion
}

/// ShapeIndexIteratorPos defines the set of possible iterator starting positions. By
/// default iterators are unpositioned, since this avoids an extra seek in this
/// situation where one of the seek methods (such as Locate) is immediately called.
enum ShapeIndexIteratorPos {
  case begin
  case end
}

/// An iterator that provides low-level access to
/// the cells of the index. Cells are returned in increasing order of CellID.
///   for it := index.Iterator(); !it.Done(); it.Next() {
///     fmt.Print(it.CellID())
///   }
struct ShapeIndexIterator {
  let index: ShapeIndex
  var position: Int
  var id: CellId
  var cell: ShapeIndexCell?
}

let sentinelCellId = CellId(id: 0)

extension ShapeIndexIterator {
  
  /// Creates a new iterator for the given index. If a starting
  /// position is specified, the iterator is positioned at the given spot.
  init(index: ShapeIndex, pos: ShapeIndexIteratorPos) {
    self.index = index
    switch pos {
    case .begin: position = 0
    case .end: position = index.cells.count
    }
    if position < index.cells.count {
      id = index.cells[position]
      cell = index.cellMap[id]
    } else {
      id = sentinelCellId
      cell = nil
    }
  }

  /// Returns the CellID of the current index cell.
  /// If s.Done() is true, a value larger than any valid CellID is returned.
  func cellId() -> CellId {
    return id
  }

  /// Returns the current index cell.
  func indexCell() -> ShapeIndexCell? {
    // TODO(roberts): C++ has this call a virtual method to allow subclasses
    // of ShapeIndexIterator to do other work before returning the cell. Do
    // we need such a thing?
    return cell
  }

  /// Returns the Point at the center of the current position of the iterator.
  func center() -> S2Point {
    return id.point()
  }

  /// Positions the iterator at the beginning of the index.
  mutating func begin() {
    if !index.isFresh() {
      index.maybeApplyUpdates()
    }
    position = 0
    refresh()
  }

  /// Positions the iterator at the next index cell.
  mutating func next() {
    position += 1
    refresh()
  }

  /// Advances the iterator to the previous cell in the index and returns true to
  /// indicate it was not yet at the beginning of the index. If the iterator is at the
  /// first cell the call does nothing and returns false.
  mutating func prev() -> Bool {
    if position <= 0 {
      return false
    }
    position -= 1
    refresh()
    return true
  }

  /// Positions the iterator at the end of the index.
  mutating func end() {
    position = index.cells.count
    refresh()
  }

  /// Reports if the iterator is positioned at or after the last index cell.
  func done() -> Bool {
    return id == sentinelCellId
  }

  /// Updates the stored internal iterator values.
  mutating func refresh() {
    if position < index.cells.count {
      id = index.cells[position]
      cell = index.cellMap[id]
    } else {
      id = sentinelCellId
      cell = nil
    }
  }

  /// Positions the iterator at the first cell whose ID >= target, or at the
  /// end of the index if no such cell exists.
  mutating func seek(target: CellId) {
    position = 0
    // In C++, this relies on the lower_bound method of the underlying btree_map.
    // TODO(roberts): Convert this to a binary search since the list of cells is ordered.
    for (k, v) in index.cells.enumerated() {
      // We've passed the cell that is after us, so we are done.
      if v >= target {
        position = k
        break
      }
      // Otherwise, advance the position.
      position += 1
    }
    refresh()
  }

  /// Positions the iterator at the cell that contains the given Point.
  /// If no such cell exists, the iterator position is unspecified, and false is returned.
  /// The cell at the matched position is guaranteed to contain all edges that might
  /// intersect the line segment between target and the cell's center.
  mutating func locate(point: S2Point) -> Bool {
    // Let I = cellMap.LowerBound(T), where T is the leaf cell containing
    // point P. Then if T is contained by an index cell, then the
    // containing cell is either I or I'. We test for containment by comparing
    // the ranges of leaf cells spanned by T, I, and I'.
    let target = CellId(point: point)
    seek(target: target)
    if !done() && cellId().rangeMin() <= target {
      return true
    }
    if prev() && cellId().rangeMax() >= target {
      return true
    }
    return false
  }

  /// Attempts to position the iterator at the first matching index cell
  /// in the index that has some relation to the given CellID. Let T be the target CellID.
  /// If T is contained by (or equal to) some index cell I, then the iterator is positioned
  /// at I and returns Indexed. Otherwise if T contains one or more (smaller) index cells,
  /// then the iterator is positioned at the first such cell I and return Subdivided.
  /// Otherwise Disjoint is returned and the iterator position is undefined.
  mutating func locate(cellId target: CellId) -> CellRelation {
    // Let T be the target, let I = cellMap.LowerBound(T.RangeMin()), and
    // let I' be the predecessor of I. If T contains any index cells, then T
    // contains I. Similarly, if T is contained by an index cell, then the
    // containing cell is either I or I'. We test for containment by comparing
    // the ranges of leaf cells spanned by T, I, and I'.
    seek(target: target.rangeMin())
    if !done() {
      if cellId() >= target && cellId().rangeMin() <= target {
        return .indexed
      }
      if cellId() <= target.rangeMax() {
        return .subdivided
      }
    }
    if prev() && cellId().rangeMax() >= target {
      return .indexed
    }
    return .disjoint
  }

}

/// Keeps track of which shapes in a given set contain a particular point
/// (the focus). It provides an efficient way to move the focus from one point
/// to another and incrementally update the set of shapes which contain it. We use
/// this to compute which shapes contain the center of every CellID in the index,
/// by advancing the focus from one cell center to the next.
///
/// Initially the focus is at the start of the CellID space-filling curve. We then
/// visit all the cells that are being added to the ShapeIndex in increasing order
/// of CellID. For each cell, we draw two edges: one from the entry vertex to the
/// center, and another from the center to the exit vertex (where entry and exit
/// refer to the points where the space-filling curve enters and exits the cell).
/// By counting edge crossings we can incrementally compute which shapes contain
/// the cell center. Note that the same set of shapes will always contain the exit
/// point of one cell and the entry point of the next cell in the index, because
/// either (a) these two points are actually the same, or (b) the intervening
/// cells in CellID order are all empty, and therefore there are no edge crossings
/// if we follow this path from one cell to the other.
class Tracker {
  var isActive: Bool
  var a: S2Point
  var b: S2Point
  var nextCellId: CellId
  var crosser: EdgeCrosser
  var shapeIds: [Int32]
  // Shape ids saved by saveAndClearStateBefore. The state is never saved
  // recursively so we don't need to worry about maintaining a stack.
  var savedIds: [Int32]

  // newTracker returns a new tracker with the appropriate defaults.
  init() {
    // As shapes are added, we compute which ones contain the start of the
    // CellID space-filling curve by drawing an edge from OriginPoint to this
    // point and counting how many shape edges cross this edge.
    isActive = false
    a = Tracker.trackerOrigin()
    b = a
    nextCellId = CellId(face: 0).childBegin(CellId.maxLevel)
    crosser = EdgeCrosser(a: a, b: b, c: a)
    shapeIds = []
    savedIds = []
  }
  
}

extension Tracker {
  
  // trackerOrigin returns the initial focus point when the tracker is created
  // (corresponding to the start of the CellID space-filling curve).
  static func trackerOrigin() -> S2Point {
    // The start of the S2CellId space-filling curve.
    return S2Point(raw: S2Cube.faceUVToXYZ(face: 0, u: -1, v: -1).normalized())
  }

  // focus returns the current focus point of the tracker.
  func focus() -> S2Point {
    return b
  }

  /// Adds a shape whose interior should be tracked. containsOrigin indicates
  /// whether the current focus point is inside the shape. Alternatively, if
  /// the focus point is in the process of being moved (via moveTo/drawTo), you
  /// can also specify containsOrigin at the old focus point and call testEdge
  /// for every edge of the shape that might cross the current drawTo line.
  /// This updates the state to correspond to the new focus point.
  /// This requires shape.HasInterior
  func addShape(shapeId: Int32, containsFocus: Bool) {
    isActive = true
    if containsFocus {
      toggleShape(shapeId: shapeId)
    }
  }

  /// Moves the focus of the tracker to the given point. This method should
  /// only be used when it is known that there are no edge crossings between the old
  /// and new focus locations; otherwise use drawTo.
  func move(to b: S2Point) {
    self.b = b
  }

  /// Moves the focus of the tracker to the given point. After this method is
  /// called, testEdge should be called with all edges that may cross the line
  /// segment between the old and new focus locations.
  func draw(to b: S2Point) {
    a = self.b
    self.b = b
    // TODO: the edge crosser may need an in-place Init method if this gets expensive
    crosser = EdgeCrosser(a: a, b: b)
  }

  /// testEdge checks if the given edge crosses the current edge, and if so, then
  /// toggle the state of the given shapeID.
  /// This requires shape to have an interior.
  func testEdge(shapeId: Int32, edge: Edge) {
    if crosser.isEdgeOrVertexCrossing(c: edge.v0, d: edge.v1) {
      toggleShape(shapeId: shapeId)
    }
  }

  /// Used to indicate that the last argument to moveTo or drawTo
  /// was the entry vertex of the given CellID, i.e. the tracker is positioned at the
  /// start of this cell. By using this method together with isAt, the caller
  /// can avoid calling moveTo in cases where the exit vertex of the previous cell
  /// is the same as the entry vertex of the current cell.
  func setNextCellId(nextCellId: CellId) {
    self.nextCellId = nextCellId.rangeMin()
  }

  /// Reports if the focus is already at the entry vertex of the given
  /// CellID (provided that the caller calls setNextCellID as each cell is processed).
  func isAt(cellId: CellId) -> Bool {
    return cellId.rangeMin() == nextCellId
  }

  /// Adds or removes the given shapeID from the set of IDs it is tracking.
  func toggleShape(shapeId: Int32) {
    // Most shapeIDs slices are small, so special case the common steps.
    // If there is nothing here, add it.
    if shapeIds.count == 0 {
      shapeIds.append(shapeId)
      return
    }
    // If it's the first element, drop it from the slice.
    if shapeIds[0] == shapeId {
      shapeIds.removeFirst()
      return
    }
    for (i, s) in shapeIds.enumerated() {
      if s < shapeId {
        continue
      }
      // If it's in the set, cut it out.
      if s == shapeId {
        shapeIds.remove(at: i)
        return
      }
      // We've got to a point in the slice where we should be inserted.
      // (the given shapeID is now less than the current positions id.)
      shapeIds.insert(shapeId, at: i)
      return
    }
    // We got to the end and didn't find it, so add it to the list.
    shapeIds.append(shapeId)
  }

  /// Makes an internal copy of the state for shape ids below
  /// the given limit, and then clear the state for those shapes. This is used during
  /// incremental updates to track the state of added and removed shapes separately.
  func saveAndClearStateBefore(limitShapeId: Int32) {
    let limit = lowerBound(shapeId: limitShapeId)
    savedIds = Array(shapeIds[..<limit])
    shapeIds = Array(shapeIds[limit...])
  }

  // restoreStateBefore restores the state previously saved by saveAndClearStateBefore.
  // This only affects the state for shapeIDs below "limitShapeID".
  func restoreStateBefore(limitShapeId: Int32) {
    let limit = lowerBound(shapeId: limitShapeId)
    shapeIds = savedIds + Array(shapeIds[limit...])
    savedIds = []
  }

  // lowerBound returns the shapeID of the first entry x where x >= shapeID.
  func lowerBound(shapeId: Int32) -> Int {
    fatalError("not implemented")
  }

}

// There are three basic states the index can be in.
enum StateIndex {
  case stale // There are pending updates.
  case updating // Updates are currently being applied.
  case fresh // There are no pending updates.
}

/// Indexes a set of Shapes, where a Shape is some collection of
/// edges. A shape can be as simple as a single ed/ge, or as complex as a set of loops.
/// For Shapes that have interiors, the index make/s it very fast to determine which
/// Shape(s) that contain a given point or region./
/// The index can be updated incrementally by adding or removing shapes. It is
/// designed to handle up to hundreds of millions of edges. All data structures
/// are designed to be small, so the index is compact; generally it is smaller
/// than the underlying data being indexed. The index is also fast to construct.
///
/// Polygon, Loop, and Polyline implement Shape which allows these objects to
/// be indexed easily. You can find useful query methods in CrossingEdgeQuery
/// and ClosestEdgeQuery (Not yet implemented in Go).
///
/// Example showing how to build an index of Polylines:
///
///   index := NewShapeIndex()
///   for _, polyline := range polylines {
///       index.Add(polyline);
///   }
/// Now you can use a CrossingEdgeQuery or ClosestEdgeQuery here.
public class ShapeIndex {
  // shapes is a map of shape ID to shape.
  var shapes: [Int32: Shape]
  // The maximum number of edges per cell.
  // TODO(roberts): Update the comments when the usage of this is implemented.
  let maxEdgesPerCell: Int
  // nextID tracks the next ID to hand out. IDs are not reused when shapes
  // are removed from the index.
  var nextId: Int32
  // cellMap is a map from CellID to the set of clipped shapes that intersect that
  // cell. The cell IDs cover a set of non-overlapping regions on the sphere.
  // In C++, this is a BTree, so the cells are ordered naturally by the data structure.
  var cellMap: [CellId: ShapeIndexCell]
  // Track the ordered list of cell IDs.
  var cells: [CellId]
  // The current status of the index; accessed atomically.
  var status: StateIndex
  // Additions and removals are queued and processed on the first subsequent
  // query. There are several reasons to do this:
  //
  //  - It is significantly more efficient to process updates in batches if
  //    the amount of entities added grows.
  //  - Often the index will never be queried, in which case we can save both
  //    the time and memory required to build it. Examples:
  //     + Loops that are created simply to pass to an Polygon. (We don't
  //       need the Loop index, because Polygon builds its own index.)
  //     + Applications that load a database of geometry and then query only
  //       a small fraction of it.
  //
  // The main drawback is that we need to go to some extra work to ensure that
  // some methods are still thread-safe. Note that the goal is *not* to
  // make this thread-safe in general, but simply to hide the fact that
  // we defer some of the indexing work until query time.
  //
  // This mutex protects all of following fields in the index.
//  mu sync.RWMutex
  let lockQueue = DispatchQueue(label: "s2-swift.ShapeIndex.LockQueue")
  // pendingAdditionsPos is the index of the first entry that has not been processed
  // via applyUpdatesInternal.
  var pendingAdditionsPos: Int32
  // The set of shapes that have been queued for removal but not processed yet by
  
  init() {
    shapes = [:]
    maxEdgesPerCell = 10
    nextId = -1
    cellMap = [:]
    cells = []
    status = .fresh
    pendingAdditionsPos = -1
  }

}

extension ShapeIndex {
  
  // Returns an iterator for this index.
  func iterator() -> ShapeIndexIterator {
//    maybeApplyUpdates()
    return ShapeIndexIterator(index: self, pos: .begin)
  }
  
  /// Positions the iterator at the first cell in the index.
  func begin() -> ShapeIndexIterator {
    maybeApplyUpdates()
    return ShapeIndexIterator(index: self, pos: .begin)
  }
  
  // Positions the iterator at the last cell in the index.
  func end() -> ShapeIndexIterator {
    // TODO(roberts): It's possible that updates could happen to the index between
    // the time this is called and the time the iterators position is used and this
    // will be invalid or not the end. For now, things will be undefined if this
    // happens. See about referencing the IsFresh to guard for this in the future.
    maybeApplyUpdates()
    return ShapeIndexIterator(index: self, pos: .end)
  }
  
  /// Reports the number of Shapes in this index.
  var count: Int {
    return shapes.count
  }
  
  /// Resets the index to its original state.
  func reset() {
    shapes = [:]
    nextId = 0
    cellMap = [:]
    cells = []
    status = .fresh
  }
  
  /// Returns the number of edges in this index.
  func numEdges() -> Int {
    return shapes.values.reduce(0) { $0 + $1.numEdges() }
  }
  
  /// Returns the shape with the given ID, or nil if the shape has been removed from the index.
  func shape(id: Int32) -> Shape? {
    return shapes[id]
  }
  
  /// Adds the given shape to the index and returns the assigned ID..
  func add(shape: Shape) {
    shapes[nextId] = shape
    nextId += 1
    status = .stale
  }
  
  /// Reports if there are no pending updates that need to be applied.
  /// This can be useful to avoid building the index unnecessarily, or for
  /// choosing between two different algorithms depending on whether the index
  /// is available.
  ///
  /// The returned index status may be slightly out of date if the index was
  /// built in a different thread. This is fine for the intended use (as an
  /// efficiency hint), but it should not be used by internal methods.
  func isFresh() -> Bool {
    return status == .fresh
  }

  /// Reports if this is the first update to the index.
  func isFirstUpdate() -> Bool {
    // Note that it is not sufficient to check whether cellMap is empty, since
    // entries are added to it during the update process.
    return pendingAdditionsPos == 0
  }

  /// Reports if the shape with the given ID is currently slated for removal.
  func isShapeBeingRemoved(shapeId: Int32) -> Bool {
    // All shape ids being removed fall below the index position of shapes being added.
    return shapeId < pendingAdditionsPos
  }

  /// Checks if the index pieces have changed, and if so, applies pending updates.
  func maybeApplyUpdates() {
    // TODO(roberts): To avoid acquiring and releasing the mutex on every
    // query, we should use atomic operations when testing whether the status
    // is fresh and when updating the status to be fresh. This guarantees
    // that any thread that sees a status of fresh will also see the
    // corresponding index updates.
    if status != .fresh {
//      lockQueue.sync() {
        applyUpdatesInternal()
        status = .fresh
//      }
    }
  }

  /// Does the actual work of updating the index by applying all
  /// pending additions and removals. It does *not* update the indexes status.
  func applyUpdatesInternal() {
    // TODO(roberts): Building the index can use up to 20x as much memory per
    // edge as the final index memory size. If this causes issues, add in
    // batched updating to limit the amount of items per batch to a
    // configurable memory footprint overhead.
    let t = Tracker()
    // allEdges maps a Face to a collection of faceEdges.
    var allEdges = (0..<6).map { _ in [FaceEdge]() } // make([][]faceEdge, 6)
    for id in pendingAdditionsPos..<Int32(shapes.count) {
      addShapeInternal(shapeId: id, allEdges: &allEdges, t: t)
    }
    for face in 0..<6 {
      updateFaceEdges(face: face, faceEdges: allEdges[face], t: t)
    }
    pendingAdditionsPos = Int32(shapes.count)
    // It is the caller's responsibility to update the index status.
  }

  /// Clips all edges of the given shape to the six cube faces,
  /// adds the clipped edges to the set of allEdges, and starts tracking its
  /// interior if necessary.
  func addShapeInternal(shapeId: Int32, allEdges: inout [[FaceEdge]], t: Tracker) {
    guard let shape = shapes[shapeId] else {
      // This shape has already been removed.
      return
    }
    if shape.hasInterior() {
      let containsFocus = containsBruteForce(shape: shape, point: t.focus())
      t.addShape(shapeId: shapeId, containsFocus: containsFocus)
    }
    let numEdges = shape.numEdges()
    for e in 0..<numEdges {
      let edge = shape.edge(e)
      let maxLevel = ShapeIndex.maxLevel(edge: edge)
      addFaceEdge(shapeId: shapeId, edgeId: e, maxLevel: maxLevel, hasInterior: shape.hasInterior(), edge: edge, allEdges: &allEdges)
    }
  }

  /// Adds the given faceEdge into the collection of all edges.
  func addFaceEdge(shapeId: Int32, edgeId: Int, maxLevel: Int, hasInterior: Bool, edge: Edge, allEdges: inout [[FaceEdge]]) {
    let aFace = S2Cube.face(point: edge.v0)
    // See if both endpoints are on the same face, and are far enough from
    // the edge of the face that they don't intersect any (padded) adjacent face.
    if aFace == S2Cube.face(point: edge.v1) {
      let (xa, ya) = S2Cube.validFaceXYZToUV(face: aFace, point: edge.v0)
      let a = R2Point(x: xa, y: ya)
      let (xb, yb) = S2Cube.validFaceXYZToUV(face: aFace, point: edge.v1)
      let b = R2Point(x: xb, y: yb)
      let maxUV = 1 - cellPadding
      if fabs(a.x) <= maxUV && fabs(a.y) <= maxUV && fabs(b.x) <= maxUV && fabs(b.y) <= maxUV {
        let fe = FaceEdge(shapeId: shapeId, edgeId: edgeId, maxLevel: maxLevel, hasInterior: hasInterior, edge: edge, a: a, b: b)
        allEdges[aFace].append(fe)
        return
      }
    }
    // Otherwise, we simply clip the edge to all six faces.
    for face in 0..<6 {
      let (aClip, bClip, intersects) = clipToPaddedFace(a: edge.v0, b: edge.v1, f: face, padding: cellPadding)
      if intersects {
        let fe = FaceEdge(shapeId: shapeId, edgeId: edgeId, maxLevel: maxLevel, hasInterior: hasInterior, edge: edge, a: aClip, b: bClip)
        allEdges[face].append(fe)
      }
    }
    return
  }

  /// Adds or removes the various edges from the index.
  /// An edge is added if shapes[id] is not nil, and removed otherwise.
  func updateFaceEdges(face: Int, faceEdges: [FaceEdge], t: Tracker) {
    let numEdges = faceEdges.count
    if numEdges == 0 && t.shapeIds.count == 0 {
      return
    }
    // Create the initial clippedEdge for each faceEdge. Additional clipped
    // edges are created when edges are split between child cells. We create
    // two arrays, one containing the edge data and another containing pointers
    // to those edges, so that during the recursion we only need to copy
    // pointers in order to propagate an edge to the correct child.
    var clippedEdges = [ClippedEdge]() // numEdges
    var bound = R2Rect.empty
    for e in 0..<numEdges {
      let clippedBound = R2Rect(p0: faceEdges[e].a, p1: faceEdges[e].b)
      let clipped = ClippedEdge(faceEdge: faceEdges[e], bound: clippedBound)
      clippedEdges.append(clipped)
      bound = bound.add(clippedBound)
    }
    // Construct the initial face cell containing all the edges, and then update
    // all the edges in the index recursively.
    let faceId = CellId(face: face)
    var pcell = PaddedCell(id: faceId, padding: cellPadding)
    let disjointFromIndex = isFirstUpdate()
    if numEdges > 0 {
      let shrunkId = shrinkToFit(pcell: pcell, bound: bound)
      if shrunkId != pcell.id {
        // All the edges are contained by some descendant of the face cell. We
        // can save a lot of work by starting directly with that cell, but if we
        // are in the interior of at least one shape then we need to create
        // index entries for the cells we are skipping over.
        skipCellRange(begin: faceId.rangeMin(), end: shrunkId.rangeMin(), t: t, disjointFromIndex: disjointFromIndex)
        pcell = PaddedCell(id: shrunkId, padding: cellPadding)
        updateEdges(pcell: pcell, edges: &clippedEdges, t: t, disjointFromIndex: disjointFromIndex)
        skipCellRange(begin: shrunkId.rangeMax().next(), end: faceId.rangeMax().next(), t: t, disjointFromIndex: disjointFromIndex)
        return
      }
    }
    // Otherwise (no edges, or no shrinking is possible), subdivide normally.
    updateEdges(pcell: pcell, edges: &clippedEdges, t: t, disjointFromIndex: disjointFromIndex)
  }

  /// Shrinks the PaddedCell to fit within the given bounds.
  func shrinkToFit(pcell: PaddedCell, bound: R2Rect) -> CellId {
    var shrunkId = pcell.shrinkToFit(rect: bound)
    if !isFirstUpdate() && shrunkId != pcell.id {
      // Don't shrink any smaller than the existing index cells, since we need
      // to combine the new edges with those cells.
      var iter = iterator()
      if iter.locate(cellId: shrunkId) == .indexed {
        shrunkId = iter.cellId()
      }
    }
    return shrunkId
  }
  
  /// Skips over the cells in the given range, creating index cells if we are
  /// currently in the interior of at least one shape.
  func skipCellRange(begin: CellId, end: CellId, t: Tracker, disjointFromIndex: Bool) {
    // If we aren't in the interior of a shape, then skipping over cells is easy.
    if t.shapeIds.count == 0 {
      return
    }
    // Otherwise generate the list of cell ids that we need to visit, and create
    // an index entry for each one.
    let skipped = CellUnion(begin: begin, end: end)
    for cellId in skipped.cellIds {
      var clippedEdges: [ClippedEdge] = []
      let pcell = PaddedCell(id: cellId, padding: cellPadding)
      updateEdges(pcell: pcell, edges: &clippedEdges, t: t, disjointFromIndex: disjointFromIndex)
    }
  }
  
  /// Adds or removes the given edges whose bounding boxes intersect a
  /// given cell. disjointFromIndex is an optimization hint indicating that cellMap
  /// does not contain any entries that overlap the given cell.
  func updateEdges(pcell: PaddedCell, edges: inout [ClippedEdge], t: Tracker, disjointFromIndex: Bool) {
    // This function is recursive with a maximum recursion depth of 30 (maxLevel).
    // Incremental updates are handled as follows. All edges being added or
    // removed are combined together in edges, and all shapes with interiors
    // are tracked using tracker. We subdivide recursively as usual until we
    // encounter an existing index cell. At this point we absorb the index
    // cell as follows:
    //
    //   - Edges and shapes that are being removed are deleted from edges and
    //     tracker.
    //   - All remaining edges and shapes from the index cell are added to
    //     edges and tracker.
    //   - Continue subdividing recursively, creating new index cells as needed.
    //   - When the recursion gets back to the cell that was absorbed, we
    //     restore edges and tracker to their previous state.
    //
    // Note that the only reason that we include removed shapes in the recursive
    // subdivision process is so that we can find all of the index cells that
    // contain those shapes efficiently, without maintaining an explicit list of
    // index cells for each shape (which would be expensive in terms of memory).
    var indexCellAbsorbed = false
    var disjointFromIndex = disjointFromIndex
    if !disjointFromIndex {
      // There may be existing index cells contained inside pcell. If we
      // encounter such a cell, we need to combine the edges being updated with
      // the existing cell contents by absorbing the cell.
      var iter = iterator()
      let r = iter.locate(cellId: pcell.id)
      if r == .disjoint {
        disjointFromIndex = true
      } else if r == .indexed {
        // Absorb the index cell by transferring its contents to edges and
        // deleting it. We also start tracking the interior of any new shapes.
        absorbIndexCell(p: pcell, iter: iter, edges: &edges, t: t)
        indexCellAbsorbed = true
        disjointFromIndex = true
      } else {
        // DCHECK_EQ(SUBDIVIDED, r)
      }
    }
    // If there are existing index cells below us, then we need to keep
    // subdividing so that we can merge with those cells. Otherwise,
    // makeIndexCell checks if the number of edges is small enough, and creates
    // an index cell if possible (returning true when it does so).
    if !disjointFromIndex || !makeIndexCell(p: pcell, edges: edges, t: t) {
      // TODO(roberts): If it turns out to have memory problems when there
      // are 10M+ edges in the index, look into pre-allocating space so we
      // are not always appending.
      let x = [ClippedEdge]()
      var childEdges = [[x, x], [x, x]] // [i][j]
      // Compute the middle of the padded cell, defined as the rectangle in
      // (u,v)-space that belongs to all four (padded) children. By comparing
      // against the four boundaries of middle we can determine which children
      // each edge needs to be propagated to.
      let middle = pcell.middle
      // Build up a vector edges to be passed to each child cell. The (i,j)
      // directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
      // Note that the vast majority of edges are propagated to a single child.
      for edge in edges {
        if edge.bound.x.hi <= middle.x.lo {
          // Edge is entirely contained in the two left children.
          let (a, b) = clipVAxis(edge: edge, middle: middle.y)
          if let a = a {
            childEdges[0][0].append(a)
          }
          if let b = b {
            childEdges[0][1].append(b)
          }
        } else if edge.bound.x.lo >= middle.x.hi {
          // Edge is entirely contained in the two right children.
          let (a, b) = clipVAxis(edge: edge, middle: middle.y)
          if let a = a {
            childEdges[1][0].append(a)
          }
          if let b = b {
            childEdges[1][1].append(b)
          }
        } else if edge.bound.y.hi <= middle.y.lo {
          // Edge is entirely contained in the two lower children.
          let a = clipUBound(edge: edge, uEnd: 1, u: middle.x.hi)
          childEdges[0][0].append(a)
          let b = clipUBound(edge: edge, uEnd: 0, u: middle.x.lo)
          childEdges[1][0].append(b)
        } else if edge.bound.y.lo >= middle.y.hi {
          // Edge is entirely contained in the two upper children.
          let a = clipUBound(edge: edge, uEnd: 1, u: middle.x.hi)
          childEdges[0][1].append(a)
          let b = clipUBound(edge: edge, uEnd: 0, u: middle.x.lo)
          childEdges[1][1].append(b)
        } else {
          // The edge bound spans all four children. The edge
          // itself intersects either three or four padded children.
          let left = clipUBound(edge: edge, uEnd: 1, u: middle.x.hi)
          let (a1, b1) = clipVAxis(edge: left, middle: middle.y)
          if let a = a1 {
            childEdges[0][0].append(a)
          }
          if let b = b1 {
            childEdges[0][1].append(b)
          }
          let right = clipUBound(edge: edge, uEnd: 0, u: middle.x.lo)
          let (a2, b2) = clipVAxis(edge: right, middle: middle.y)
          if let a = a2 {
            childEdges[1][0].append(a)
          }
          if let b = b2 {
            childEdges[1][1].append(b)
          }
        }
      }
      // Now recursively update the edges in each child. We call the children in
      // increasing order of CellID smo that when the index is first constructed,
      // all insertions into cellMap are at the end (which is much faster).
      for pos in 0..<4 {
        let (i, j) = pcell.childIJ(pos: pos)
        if childEdges[i][j].count > 0 || t.shapeIds.count > 0 {
          updateEdges(pcell: PaddedCell(parent: pcell, i: i, j: j), edges: &childEdges[i][j], t: t, disjointFromIndex: disjointFromIndex)
        }
      }
    }
    if indexCellAbsorbed {
      // Restore the state for any edges being removed that we are tracking.
      t.restoreStateBefore(limitShapeId: pendingAdditionsPos)
    }
  }

  /// Builds an indexCell from the given padded cell and set of edges and adds
  /// it to the index. If the cell or edges are empty, no cell is added.
  func makeIndexCell(p: PaddedCell, edges: [ClippedEdge], t: Tracker) -> Bool {
    // If the cell is empty, no index cell is needed. (In most cases this
    // situation is detected before we get to this point, but this can happen
    // when all shapes in a cell are removed.)
    if edges.count == 0 && t.shapeIds.count == 0 {
      return true
    }
    // Count the number of edges that have not reached their maximum level yet.
    // Return false if there are too many such edges.
    var count = 0
    for ce in edges {
      if p.level < ce.faceEdge.maxLevel {
        count += 1
      }
      if count > maxEdgesPerCell {
        return false
      }
    }
    // Possible optimization: Continue subdividing as long as exactly one child
    // of the padded cell intersects the given edges. This can be done by finding
    // the bounding box of all the edges and calling ShrinkToFit:
    //
    // cellID = p.ShrinkToFit(RectBound(edges));
    //
    // Currently this is not beneficial; it slows down construction by 4-25%
    // (mainly computing the union of the bounding rectangles) and also slows
    // down queries (since more recursive clipping is required to get down to
    // the level of a spatial index cell). But it may be worth trying again
    // once containsCenter is computed and all algorithms are modified to
    // take advantage of it.
    // We update the InteriorTracker as follows. For every Cell in the index
    // we construct two edges: one edge from entry vertex of the cell to its
    // center, and one from the cell center to its exit vertex. Here entry
    // and exit refer the CellID ordering, i.e. the order in which points
    // are encountered along the 2 space-filling curve. The exit vertex then
    // becomes the entry vertex for the next cell in the index, unless there are
    // one or more empty intervening cells, in which case the InteriorTracker
    // state is unchanged because the intervening cells have no edges.
    // Shift the InteriorTracker focus point to the center of the current cell.
    if t.isActive && edges.count != 0 {
      if !t.isAt(cellId: p.id) {
        t.move(to: p.entryVertex())
      }
      t.draw(to: p.center())
      testAllEdges(edges: edges, t: t)
    }
    // Allocate and fill a new index cell. To get the total number of shapes we
    // need to merge the shapes associated with the intersecting edges together
    // with the shapes that happen to contain the cell center.
    var cshapeIds = t.shapeIds
    let numShapes = countShapes(edges: edges, shapeIds: cshapeIds)
    // To fill the index cell we merge the two sources of shapes: edge shapes
    // (those that have at least one edge that intersects this cell), and
    // containing shapes (those that contain the cell center). We keep track
    // of the index of the next intersecting edge and the next containing shape
    // as we go along. Both sets of shape ids are already sorted.
    var eNext = 0
    var cNextIdx = 0
    let cellShapes = (0..<numShapes).map { (i: Int) -> ClippedShape in
      // advance to next value base + i
      var eshapeId = Int32(shapes.count)
      var cshapeId = Int32(eshapeId) // Sentinels
      if eNext != edges.count {
        eshapeId = edges[eNext].faceEdge.shapeId
      }
      if cNextIdx != cshapeIds.count {
        cshapeId = cshapeIds[cNextIdx]
      }
      let eBegin = eNext
      if cshapeId < eshapeId {
        cNextIdx += 1
        // The entire cell is in the shape interior.
        return ClippedShape(shapeId: cshapeId, containsCenter: true, edges: [])
      }
      // Count the number of edges for this shape and allocate space for them.
      while eNext < edges.count && edges[eNext].faceEdge.shapeId == eshapeId {
        eNext += 1
      }
      let cedges = (eBegin..<eNext).map { edges[$0].faceEdge.edgeId }
      if cshapeId == eshapeId {
        cNextIdx += 1
      }
      let containsCenter = (cshapeId == eshapeId)
      return ClippedShape(shapeId: cshapeId, containsCenter: containsCenter, edges: cedges)
    }
    let cell = ShapeIndexCell(shapes: cellShapes)
    // Add this cell to the map.
    cellMap[p.id] = cell
    cells.append(p.id)
    // Shift the tracker focus point to the exit vertex of this cell.
    if t.isActive && edges.count != 0 {
      t.draw(to: p.exitVertex())
      testAllEdges(edges: edges, t: t)
      t.setNextCellId(nextCellId: p.id.next())
    }
    return true
  }

  /// Updates the specified endpoint of the given clipped edge and returns the
  /// resulting clipped edge.
  func updateBound(edge: ClippedEdge, uEnd: Int, u: Double, vEnd: Int, v: Double) -> ClippedEdge {
    let x = R1Interval(lo: uEnd == 0 ? u : edge.bound.x.lo, hi: uEnd == 0 ? edge.bound.x.hi : u)
    let y = R1Interval(lo: vEnd == 0 ? v : edge.bound.y.lo, hi: vEnd == 0 ? edge.bound.y.hi : v)
    let bound = R2Rect(x: x, y: y)
    return ClippedEdge(faceEdge: edge.faceEdge, bound: bound)
  }

  /// Clips the given endpoint (lo=0, hi=1) of the u-axis so that
  /// it does not extend past the given value of the given edge.
  func clipUBound(edge: ClippedEdge, uEnd: Int, u: Double) -> ClippedEdge {
    // First check whether the edge actually requires any clipping. (Sometimes
    // this method is called when clipping is not necessary, e.g. when one edge
    // endpoint is in the overlap area between two padded child cells.)
    if uEnd == 0 {
      if edge.bound.x.lo >= u {
        return edge
      }
    } else {
      if edge.bound.x.hi <= u {
        return edge
      }
    }
    // We interpolate the new v-value from the endpoints of the original edge.
    // This has two advantages: (1) we don't need to store the clipped endpoints
    // at all, just their bounding box; and (2) it avoids the accumulation of
    // roundoff errors due to repeated interpolations. The result needs to be
    // clamped to ensure that it is in the appropriate range.
    let e = edge.faceEdge
    let value = interpolateFloat64(x: u, a: e.a.x, b: e.b.x, a1: e.a.y, b1: e.b.y)
    let v = edge.bound.y.clamp(value)
    // Determine which endpoint of the v-axis bound to update. If the edge
    // slope is positive we update the same endpoint, otherwise we update the
    // opposite endpoint.
    var vEnd = 0
    let positiveSlope = (e.a.x > e.b.x) == (e.a.y > e.b.y)
    if (uEnd == 1) == positiveSlope {
      vEnd = 1
    }
    return updateBound(edge: edge, uEnd: uEnd, u: u, vEnd: vEnd, v: v)
  }

  /// Clips the given endpoint (lo=0, hi=1) of the v-axis so that
  /// it does not extend past the given value of the given edge.
  func clipVBound(edge: ClippedEdge, vEnd: Int, v: Double) -> ClippedEdge {
    if vEnd == 0 {
      if edge.bound.y.lo >= v {
        return edge
      }
    } else {
      if edge.bound.y.hi <= v {
        return edge
      }
    }
    // We interpolate the new v-value from the endpoints of the original edge.
    // This has two advantages: (1) we don't need to store the clipped endpoints
    // at all, just their bounding box; and (2) it avoids the accumulation of
    // roundoff errors due to repeated interpolations. The result needs to be
    // clamped to ensure that it is in the appropriate range.
    let e = edge.faceEdge
    let value = interpolateFloat64(x: v, a: e.a.y, b: e.b.y, a1: e.a.x, b1: e.b.x)
    let u = edge.bound.x.clamp(value)
    // Determine which endpoint of the v-axis bound to update. If the edge
    // slope is positive we update the same endpoint, otherwise we update the
    // opposite endpoint.
    var uEnd = 0
    let positiveSlope = (e.a.x > e.b.x) == (e.a.y > e.b.y)
    if (vEnd == 1) == positiveSlope {
      uEnd = 1
    }
    return updateBound(edge: edge, uEnd: uEnd, u: u, vEnd: vEnd, v: v)
  }
  
  /// Returns the given edge clipped to within the boundaries of the middle
  /// interval along the v-axis, and adds the result to its children.
  func clipVAxis(edge: ClippedEdge, middle: R1Interval) -> (a: ClippedEdge?, b: ClippedEdge?) {
    if edge.bound.y.hi <= middle.lo {
      // Edge is entirely contained in the lower child.
      return (edge, nil)
    } else if edge.bound.y.lo >= middle.hi {
      // Edge is entirely contained in the upper child.
      return (nil, edge)
    }
    // The edge bound spans both children.
    return (clipVBound(edge: edge, vEnd: 1, v: middle.hi), clipVBound(edge: edge, vEnd: 0, v: middle.lo))
  }
  
  /// Absorbs an index cell by transferring its contents to edges
  /// and/or "tracker", and then delete this cell from the index. If edges includes
  /// any edges that are being removed, this method also updates their
  /// InteriorTracker state to correspond to the exit vertex of this cell.
  func absorbIndexCell(p: PaddedCell, iter: ShapeIndexIterator, edges: inout [ClippedEdge], t: Tracker) {
    // When we absorb a cell, we erase all the edges that are being removed.
    // However when we are finished with this cell, we want to restore the state
    // of those edges (since that is how we find all the index cells that need
    // to be updated).  The edges themselves are restored automatically when
    // UpdateEdges returns from its recursive call, but the InteriorTracker
    // state needs to be restored explicitly.
    //
    // Here we first update the InteriorTracker state for removed edges to
    // correspond to the exit vertex of this cell, and then save the
    // InteriorTracker state.  This state will be restored by UpdateEdges when
    // it is finished processing the contents of this cell.
    if t.isActive && edges.count != 0 && isShapeBeingRemoved(shapeId: edges[0].faceEdge.shapeId) {
      // We probably need to update the tracker. ("Probably" because
      // it's possible that all shapes being removed do not have interiors.)
      if !t.isAt(cellId: p.id) {
        t.move(to: p.entryVertex())
      }
      t.draw(to: p.exitVertex())
      t.setNextCellId(nextCellId: p.id.next())
      for edge in edges {
        let fe = edge.faceEdge
        if !isShapeBeingRemoved(shapeId: fe.shapeId) {
          break // All shapes being removed come first.
        }
        if fe.hasInterior {
          t.testEdge(shapeId: fe.shapeId, edge: fe.edge)
        }
      }
    }
    // Save the state of the edges being removed, so that it can be restored
    // when we are finished processing this cell and its children.  We don't
    // need to save the state of the edges being added because they aren't being
    // removed from "edges" and will therefore be updated normally as we visit
    // this cell and its children.
    t.saveAndClearStateBefore(limitShapeId: pendingAdditionsPos)
    // Create a faceEdge for each edge in this cell that isn't being removed.
    var faceEdges: [FaceEdge] = []
    var trackerMoved = false
    let cell = iter.indexCell()!
    for clipped in cell.shapes {
      let shapeId = clipped.shapeId
      guard let shape = self.shape(id: shapeId) else {
        continue // This shape is being removed.
      }
      let numClipped = clipped.numEdges()
      // If this shape has an interior, start tracking whether we are inside the
      // shape. updateEdges wants to know whether the entry vertex of this
      // cell is inside the shape, but we only know whether the center of the
      // cell is inside the shape, so we need to test all the edges against the
      // line segment from the cell center to the entry vertex.
      let hasInterior = shape.hasInterior()
      if hasInterior {
        t.addShape(shapeId: shapeId, containsFocus: clipped.containsCenter)
        // There might not be any edges in this entire cell (i.e., it might be
        // in the interior of all shapes), so we delay updating the tracker
        // until we see the first edge.
        if !trackerMoved && numClipped > 0 {
          t.move(to: p.center())
          t.draw(to: p.entryVertex())
          t.setNextCellId(nextCellId: p.id)
          trackerMoved = true
        }
      }
      for i in 0..<numClipped {
        let edgeId = clipped.edges[i]
        let edge = shape.edge(edgeId)
        let maxLevel = ShapeIndex.maxLevel(edge: edge)
        if hasInterior {
          t.testEdge(shapeId: shapeId, edge: edge)
        }
        let (a, b, intersects) = clipToPaddedFace(a: edge.v0, b: edge.v1, f: p.id.face(), padding: cellPadding)
        if !intersects {
          fatalError("invariant failure in ShapeIndex")
        }
        let fe = FaceEdge(shapeId: shapeId, edgeId: edgeId, maxLevel: maxLevel, hasInterior: hasInterior, edge: edge, a: a, b: b)
        faceEdges.append(fe)
      }
    }
    // Now create a clippedEdge for each faceEdge, and put them in "new_edges".
    var newEdges: [ClippedEdge] = []
    for faceEdge in faceEdges {
      let clipped = ClippedEdge(faceEdge: faceEdge, bound: clippedEdgeBound(a: faceEdge.a, b: faceEdge.b, clip: p.bound))
      newEdges.append(clipped)
    }
    // Discard any edges from "edges" that are being removed, and append the
    // remainder to "newEdges"  (This keeps the edges sorted by shape id.)
    for (i, clipped) in edges.enumerated() {
      if !isShapeBeingRemoved(shapeId: clipped.faceEdge.shapeId) {
        newEdges += edges[i...]
        break
      }
    }
    // Update the edge list and delete this cell from the index.
    (edges, newEdges) = (newEdges, edges)
    cellMap.removeValue(forKey: p.id)
    // TODO(roberts): delete from s.Cells
  }

  /// Calls the trackers testEdge on all edges from shapes that have interiors.
  func testAllEdges(edges: [ClippedEdge], t: Tracker) {
    for edge in edges {
      if edge.faceEdge.hasInterior {
        t.testEdge(shapeId: edge.faceEdge.shapeId, edge: edge.faceEdge.edge)
      }
    }
  }

  /// Reports the number of distinct shapes that are either associated with the
  /// given edges, or that are currently stored in the InteriorTracker.
  func countShapes(edges: [ClippedEdge], shapeIds: [Int32]) -> Int {
    var count = 0
    var lastShapeId = Int32(-1)
    var cNext = Int32(0)
    for edge in edges {
      if edge.faceEdge.shapeId == lastShapeId {
        continue
      }
      count += 1
      lastShapeId = edge.faceEdge.shapeId
      // Skip over any containing shapes up to and including this one,
      // updating count as appropriate.
      while cNext < Int32(shapeIds.count) {
        if cNext > lastShapeId {
          break
        }
        if cNext < lastShapeId {
          count += 1
        }
        cNext += 1
      }
    }
    // Count any remaining containing shapes.
    count += Int(shapeIds.count - Int(cNext))
    return count
  }

  /// Reports the maximum level for a given edge.
  static func maxLevel(edge: Edge) -> Int {
    // Compute the maximum cell size for which this edge is considered long.
    // The calculation does not need to be perfectly accurate, so we use Norm
    // rather than Angle for speed.
    let cellSize = edge.v0.sub(edge.v1).norm * cellSizeToLongEdgeRatio
    // Now return the first level encountered during subdivision where the
    // average cell size is at most cellSize.
    return S2CellMetric.avgEdge.minLevel(cellSize)
  }
  
}
