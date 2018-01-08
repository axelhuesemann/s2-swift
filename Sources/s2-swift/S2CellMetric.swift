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
  static let maxWidth = S2CellMetric(dim: 1, deriv: 1.704897179199218452)

  static let minArea = S2CellMetric(dim: 2, deriv: 8 * sqrt(2.0) / 9)
  static let avgArea = S2CellMetric(dim: 2, deriv: 4 * .pi / 6)
  static let maxArea = S2CellMetric(dim: 2, deriv: 2.635799256963161491)

  // TODO: more metrics, as needed
  // TODO: port GetValue, GetClosestLevel

  // Value returns the value of the metric at the given level.
  func value(_ level: Int) -> Double {
    return ldexp(deriv, -dim * level)
  }

  // MinLevel returns the minimum level such that the metric is at most
  // the given value, or maxLevel (30) if there is no such level.
  func minLevel(_ val: Double) -> Int {
    if val < 0 {
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
