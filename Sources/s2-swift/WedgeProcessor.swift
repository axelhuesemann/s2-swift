//
//  WedgeProcessor.swift
//  s2-swiftGo
//

import Foundation

enum WedgeRelation {
  case WEDGE_EQUALS
  case WEDGE_PROPERLY_CONTAINS  // A is a strict superset of B.
  case WEDGE_IS_PROPERLY_CONTAINED  // A is a strict subset of B.
  case WEDGE_PROPERLY_OVERLAPS  // All of A intsect B, A-B and B-A are non-empty.
  case WEDGE_IS_DISJOINT  // A is disjoint from B
}

class S2EdgeUtil {
  static func wedgeContains(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool { return false }
  static func wedgeIntersects(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool { return false }
  static func getWedgeRelation(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> WedgeRelation { return .WEDGE_EQUALS }
}
// This is a helper class for the AreBoundariesCrossing function.
// Each time there is a point in common between the two loops passed
// as parameters, the two associated wedges centered at this point are
// passed to the ProcessWedge function of this processor. The function
// updates an internal state based on the wedges, and returns true to
// signal that no further processing is needed.
//
// To use AreBoundariesCrossing, subclass this class and keep an
// internal state that you update each time ProcessWedge is called,
// then query this internal state in the function that called
// AreBoundariesCrossing.
class WedgeProcessor {
  
  func processWedge(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool { return false }
}

// WedgeProcessor to be used to check if loop A contains loop B.
// DoesntContain() then returns true if there is a wedge of B not
// contained in the associated wedge of A (and hence loop B is not
// contained in loop A).
class ContainsWedgeProcessor: WedgeProcessor {
  var doesnt_contain = false
  
  override func processWedge(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
    doesnt_contain = !S2EdgeUtil.wedgeContains(a0: a0, ab1: ab1, a2: a2, b0: b0, b2: b2)
    return doesnt_contain
  }
}

// WedgeProcessor to be used to check if loop A intersects loop B.
// Intersects() then returns true when A and B have at least one pair
// of associated wedges that intersect.
class IntersectsWedgeProcessor: WedgeProcessor {
  var intersects = false
  
  override func processWedge(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
    intersects = S2EdgeUtil.wedgeIntersects(a0: a0, ab1: ab1, a2: a2, b0: b0, b2: b2)
    return intersects
  }
}

// WedgeProcessor to be used to check if the interior of loop A
// contains the interior of loop B, or their boundaries cross each
// other (therefore they have a proper intersection).
// CrossesOrMayContain() then returns -1 if A crossed B, 0 if it is
// not possible for A to contain B, and 1 otherwise.
class ContainsOrCrossesProcessor: WedgeProcessor {
  // True if any crossing on the boundary is discovered.
  var has_boundary_crossing = false
  // True if A (B) has a strictly superwedge on a pair of wedges that
  // share a common center point.
  var a_has_strictly_super_wedge = false
  var b_has_strictly_super_wedge = false
  // True if there is a pair of disjoint wedges with common center
  // point.
  var has_disjoint_wedge = false
  
  func crossesOrMayContain() -> Int {
    if has_boundary_crossing { return -1 }
    if has_disjoint_wedge || b_has_strictly_super_wedge { return 0 }
    return 1
  }
  
  override func processWedge(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point) -> Bool {
    let wedge_relation = S2EdgeUtil.getWedgeRelation(a0: a0, ab1: ab1, a2: a2, b0: b0, b2: b2)
    if wedge_relation == .WEDGE_PROPERLY_OVERLAPS {
      has_boundary_crossing = true
      return true
    }
    a_has_strictly_super_wedge = a_has_strictly_super_wedge || (wedge_relation == .WEDGE_PROPERLY_CONTAINS)
    b_has_strictly_super_wedge = b_has_strictly_super_wedge || (wedge_relation == .WEDGE_IS_PROPERLY_CONTAINED)
    if a_has_strictly_super_wedge && b_has_strictly_super_wedge {
      has_boundary_crossing = true
      return true
    }
    has_disjoint_wedge = has_disjoint_wedge || (wedge_relation == .WEDGE_IS_DISJOINT)
    return false
  }
}
