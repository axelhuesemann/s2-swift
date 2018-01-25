//
//  S2Metric.swift
//  s2-swift
//

import Foundation

/// A measure for cells.
public struct S2CellMetric {
  
  // Dim is either 1 or 2, for a 1D or 2D metric respectively.
  let dim: Int
  // Deriv is the scaling factor for the metric.
  let deriv: Double

  // MARK: inits / factory 

  // default constructur is automatic
  
  // Defined metrics.
  // We only support the quadratic projection.
  static let minWidth = S2CellMetric(dim: 1, deriv: 2 * sqrt(2.0) / 3)
  static let avgWidth = S2CellMetric(dim: 1, deriv: 1.434523672886099389)
  static let maxWidth = S2CellMetric(dim: 1, deriv: 1.704897179199218452)

  static let minArea = S2CellMetric(dim: 2, deriv: 8 * sqrt(2.0) / 9)
  static let avgArea = S2CellMetric(dim: 2, deriv: 4 * .pi / 6)
  static let maxArea = S2CellMetric(dim: 2, deriv: 2.635799256963161491)

  // Each cell is bounded by four planes passing through its four edges and
  // the center of the sphere. These metrics relate to the angle between each
  // pair of opposite bounding planes, or equivalently, between the planes
  // corresponding to two different s-values or two different t-values.
  static let minAngleSpan = S2CellMetric(dim: 1, deriv: 4.0 / 3)
  static let avgAngleSpan = S2CellMetric(dim: 1, deriv: .pi / 2)
  static let maxAngleSpan = S2CellMetric(dim: 1, deriv: 1.704897179199218452)

  // The edge length metrics can be used to bound the minimum, maximum,
  // or average distance from the center of one cell to the center of one of
  // its edge neighbors. In particular, it can be used to bound the distance
  // between adjacent cell centers along the space-filling Hilbert curve for
  // cells at any given level.
  static let minEdge = S2CellMetric(dim: 1, deriv: 2 * sqrt(2.0) / 3.0)
  static let avgEdge = S2CellMetric(dim: 1, deriv: 1.459213746386106062)
  static let maxEdge = S2CellMetric(dim: 1, deriv: 1.704897179199218452)

  static let maxEdgeAspect = 1.442615274452682920

  // The maximum diagonal is also the maximum diameter of any cell,
  // and also the maximum geometric width (see the comment for widths). For
  // example, the distance from an arbitrary point to the closest cell center
  // at a given level is at most half the maximum diagonal length.
  static let minDiag = S2CellMetric(dim: 1, deriv: 8 * sqrt(2.0) / 9.0)
  static let avgDiag = S2CellMetric(dim: 1, deriv: 2.060422738998471683)
  static let maxDiag = S2CellMetric(dim: 1, deriv: 2.438654594434021032)

  static let maxDiagAspect = sqrt(3)

  // TODO: port GetValue, GetClosestLevel

  // Value returns the value of the metric at the given level.
  func value(_ level: Int) -> Double {
    return ldexp(deriv, -dim * level)
  }

  // MinLevel returns the minimum level such that the metric is at most
  // the given value, or maxLevel (30) if there is no such level.
  func minLevel(_ val: Double) -> Int {
    if val <= 0 {
      return CellId.maxLevel
    }
    let level = -(Int(logb(val / deriv)) >> (dim - 1))
    if level > CellId.maxLevel {
      return CellId.maxLevel
    }
    if level < 0 {
      return 0
    }
    return level
  }

  // MaxLevel returns the maximum level such that the metric is at least
  // the given value, or zero if there is no such level.
  func maxLevel(_ val: Double) -> Int {
    if val <= 0 {
      return CellId.maxLevel
    }
    let level = Int(logb(deriv / val)) >> (dim - 1)
    if level > CellId.maxLevel {
      return CellId.maxLevel
    }
    if level < 0 {
      return 0
    }
    return level
  }

}
