//
//  WedgeRelations.swift
//  s2-swift
//
//  Created by Axel Huesemann on 1/11/18.
//

import Foundation


// Enumerates the possible relation between two wedges A and B.
// Define the different possible relationships between two wedges.
//
// Given an edge chain (x0, x1, x2), the wedge at x1 is the region to the
// left of the edges. More precisely, it is the set of all rays from x1x0
// (inclusive) to x1x2 (exclusive) in the *clockwise* direction.
enum WedgeRelation {
  case equals // A and B are equal.
  case properlyContains // A is a strict superset of B.
  case isProperlyContained // A is a strict subset of B.
  case properlyOverlaps // A-B, B-A, and A intersect B are non-empty.
  case isDisjoint // A and B are disjoint.
}

// WedgeRelation reports the relation between two non-empty wedges
// A=(a0, ab1, a2) and B=(b0, ab1, b2).
func wedgeRelation(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> WedgeRelation {
  // There are 6 possible edge orderings at a shared vertex (all
  // of these orderings are circular, i.e. abcd == bcda):
  //  (1) a2 b2 b0 a0: A contains B
  //  (2) a2 a0 b0 b2: B contains A
  //  (3) a2 a0 b2 b0: A and B are disjoint
  //  (4) a2 b0 a0 b2: A and B intersect in one wedge
  //  (5) a2 b2 a0 b0: A and B intersect in one wedge
  //  (6) a2 b0 b2 a0: A and B intersect in two wedges
  // We do not distinguish between 4, 5, and 6.
  // We pay extra attention when some of the edges overlap.  When edges
  // overlap, several of these orderings can be satisfied, and we take
  // the most specific.
  if a0 == b0 && a2 == b2 {
    return .equals
  }
  // Cases 1, 2, 5, and 6
  if S2Point.orderedCCW(a0, a2, b2, ab1) {
    // The cases with this vertex ordering are 1, 5, and 6,
    if S2Point.orderedCCW(b2, b0, a0, ab1) {
      return .properlyContains
    }
    // We are in case 5 or 6, or case 2 if a2 == b2.
    if a2 == b2 {
      return .isProperlyContained
    }
    return .properlyOverlaps
  }
  // We are in case 2, 3, or 4.
  if S2Point.orderedCCW(a0, b0, b2, ab1) {
    return .isProperlyContained
  }
  if S2Point.orderedCCW(a0, b0, a2, ab1) {
    return .isDisjoint
  }
  return .properlyOverlaps
}

// WedgeContains reports whether non-empty wedge A=(a0, ab1, a2) contains B=(b0, ab1, b2).
// Equivalent to WedgeRelation == WedgeProperlyContains || WedgeEquals.
func wedgeContains(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
  // For A to contain B (where each loop interior is defined to be its left
  // side), the CCW edge order around ab1 must be a2 b2 b0 a0.  We split
  // this test into two parts that test three vertices each.
  return S2Point.orderedCCW(a2, b2, b0, ab1) && S2Point.orderedCCW(b0, a0, a2, ab1)
}

// WedgeIntersects reports whether non-empty wedge A=(a0, ab1, a2) intersects B=(b0, ab1, b2).
// Equivalent but faster than WedgeRelation != WedgeIsDisjoint
func wedgeIntersects(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
  // For A not to intersect B (where each loop interior is defined to be
  // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
  // Note that it's important to write these conditions as negatives
  // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
  // results when two vertices are the same.
  return !(S2Point.orderedCCW(a0, b2, b0, ab1) && S2Point.orderedCCW(b0, a2, a0, ab1))
}
