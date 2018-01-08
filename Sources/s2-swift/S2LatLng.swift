//
//  S2LatLng.swift
//  s2-swift
//

import Foundation


fileprivate let northPoleLat = .pi / 2.0
fileprivate let southPoleLat = -.pi / 2.0


/// Represents a point on the unit sphere as a pair of angles.
public struct LatLng {
  
  let lat: S1Angle
  let lng: S1Angle

  // MARK: inits / factory
  
  public init(lat: S1Angle, lng: S1Angle) {
    self.lat = lat
    self.lng = lng
  }
  
  public init(point: S2Point) {
    self.init(lat: point.lat, lng: point.lng)
  }
  
  // MARK: computed members
  
  /// Returns the normalized version of the LatLng,
  /// with Lat clamped to [-π/2,π/2] and Lng wrapped in [-π,π].
  func normalized() -> LatLng {
    var lat2 = lat
    if lat2 > northPoleLat {
      lat2 = northPoleLat
    } else if lat2 < southPoleLat {
      lat2 = southPoleLat
    }
    let lng2 = remainder(lng, 2.0 * .pi)
    return LatLng(lat: lat2, lng: lng2)
  }
  
  // MARK: tests
  
  /// Returns true iff the LatLng is normalized, with Lat ∈ [-π/2,π/2] and Lng ∈ [-π,π].
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
    return 2 * atan2(sqrt(x), sqrt(max(0, 1-x)))
  }

  static func lat(_ vector: R3Vector) -> S1Angle {
    return atan2(vector.z, sqrt(vector.x * vector.x + vector.y * vector.y))
  }
  
  static func lng(_ vector: R3Vector) -> S1Angle {
    return atan2(vector.y, vector.x)
  }
  
  // MARK: conversion
  
  /// Returns a point for the given LatLng.
  func toPoint() -> S2Point {
    return S2Point(latLng: self)
  }
  
}

extension LatLng: CustomStringConvertible {
  
  public var description: String {
    let lat2 = String(format: "%.7f", lat * toDegrees)
    let lng2 = String(format: "%.7f", lng * toDegrees)
    return "[\(lat2), \(lng2)]"
  }
  
}
