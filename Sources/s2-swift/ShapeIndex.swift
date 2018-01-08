//
//  S2ShapeIndex.swift
//  s2-swift
//

import Foundation


/// Defines the types of geometry dimensions that a Shape supports.
enum ShapeDimension: Int {
  case pointGeometry
  case polylineGeometry
  case polygonGeometry
}

//Describes the possible relationships between a target cell
/// and the cells of the ShapeIndex. If the target is an index cell or is
/// contained by an index cell, it is Indexed. If the target is subdivided
/// into one or more index cells, it is Subdivided. Otherwise it is Disjoint.
enum CellRelation: Int {
  case indexed = 0
  case subdivided = 1
  case disjoint = 2
}

// ShapeIndex indexes a set of Shapes, where a Shape is some collection of
// edges. A shape can be as simple as a single edge, or as complex as a set of loops.
// For Shapes that have interiors, the index makes it very fast to determine which
// Shape(s) that contain a given point or region.
public class ShapeIndex {
  // shapes contains all the shapes in this index, accessible by their shape id.
  // Removed shapes are replaced by nil.
  //
  // TODO: Is there a better storage structure to use? C++ uses a btree
  // deep down for the index. There do appear to be a number of Go BTree
  // implementations available that may be suitable. Further investigation
  // is needed before selecting an appropriate option.
  //
  // The slice is an interim storage solution to get the index up and usable.
  var shapes: [S2ShapeType] = []
  
  let maxEdgesPerCell: Int
  
  init() {
    maxEdgesPerCell = 10
  }

  // Add adds the given shape to the index and assign a unique id to the shape.
  // Shape ids are assigned sequentially starting from 0 in the order shapes are added.
  func add(_ shape: S2ShapeType) {
    shapes.append(shape)
  }

  // Len reports the number of Shapes in this index.
  var count: Int {
    return shapes.count
  }

  // At returns the shape with the given index. If the given index is not valid, nil is returned.
  func at(_ i: Int) -> S2ShapeType {
    // TODO: This blindly assumes that no Shapes have been removed and
    // that the slice has no holes in it. As this gets implemented, change this
    // to be smarter and safer about verifying existence before returning it.
    return shapes[i]
  }

  // Reset clears the contents of the index and resets it to its original state.
  // Any options specified via Init are preserved.
  func reset() {
    shapes = []
  }

}
