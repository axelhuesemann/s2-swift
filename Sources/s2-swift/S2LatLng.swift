//
//  S2LatLng.swift
//  s2-swift
//

import Foundation

/// Represents a point on the unit sphere as a pair of angles.
/// LatLng representation is good for distances.
public struct LatLng {
  
  private static let northPoleLat = .pi / 2.0
  private static let southPoleLat = -.pi / 2.0
  
  let lat: S1Angle
  let lng: S1Angle

  // MARK: inits / factory
  
  public init(lat: S1Angle, lng: S1Angle) {
    self.lat = lat
    self.lng = lng
  }
  
  public init(vector: R3Vector) {
    lat = atan2(vector.z, sqrt(vector.x * vector.x + vector.y * vector.y))
    lng = atan2(vector.y, vector.x)
  }
    
  // MARK: computed members
  
  /// Returns the normalized version of the LatLng,
  /// with Lat clamped to [-π/2,π/2] and Lng wrapped in [-π,π].
  func normalized() -> LatLng {
    let lat2 = clamp(lat, min: LatLng.southPoleLat, max: LatLng.northPoleLat)
    let lng2 = remainder(lng, 2.0 * .pi)
    return LatLng(lat: lat2, lng: lng2)
  }
  
  // MARK: tests
  
  /// Returns true iff normalized, with Lat ∈ [-π/2,π/2] and Lng ∈ [-π,π].
  var isValid: Bool {
    return abs(lat) <= .pi/2 && abs(lng) <= .pi
  }

  // MARK: arithmetic
  
  /// Returns the angle between two LatLngs.
  func distance(_ latLng: LatLng) -> S1Angle {
    // Haversine formula
    let lat1 = lat
    let lat2 = latLng.lat
    let lng1 = lng
    let lng2 = latLng.lng
    let dlat = sin(0.5 * (lat2 - lat1))
    let dlng = sin(0.5 * (lng2 - lng1))
    let x = dlat * dlat + dlng * dlng * cos(lat1) * cos(lat2)
    return 2 * atan2(sqrt(x), sqrt(max(0, 1 - x)))
  }

  // MARK: conversion
  
  /// Returns a point for the given LatLng.
  func toPoint() -> S2Point {
    let phi = lat
    let theta = lng
    let cosphi = cos(phi)
    let vector = R3Vector(x: cos(theta) * cosphi, y: sin(theta) * cosphi, z: sin(phi))
    return S2Point(raw: vector)
  }
  
}

extension LatLng: CustomStringConvertible {
  
  public var description: String {
    let lat2 = String(format: "%.7f", lat * toDegrees)
    let lng2 = String(format: "%.7f", lng * toDegrees)
    return "[\(lat2), \(lng2)]"
  }
  
}
