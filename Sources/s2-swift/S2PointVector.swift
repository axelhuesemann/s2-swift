//
//  S2PointVector.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/29/18.
//

import Foundation


/// A Shape representing a set of Points. Each point
/// is represented as a degenerate point with the same starting and ending
/// vertices.
/// This type is useful for adding a collection of points to an ShapeIndex.
struct S2PointVector {
  let points: [S2Point]
}

extension S2PointVector: S2Shape {
  func numEdges() -> Int { return points.count }
  func edge(_ i: Int) -> Edge { return Edge(v0: points[i], v1: points[i]) }
  func hasInterior() -> Bool { return false }
  func containsOrigin() -> Bool { return false }
  func referencePoint() -> ReferencePoint { return ReferencePoint(origin: true, contained: false) }
  func numChains() -> Int { return points.count }
  func chain(_ i: Int) -> Chain { return Chain(start: i, length: 1) }
  func chainEdge(chainId i: Int, offset j: Int) -> Edge { return Edge(v0: points[i], v1: points[j]) }
  func chainPosition(_ e: Int) -> ChainPosition { return ChainPosition(chainId: e, offset: 0) }
  func dimension() -> ShapeDimension { return .pointGeometry }
}
