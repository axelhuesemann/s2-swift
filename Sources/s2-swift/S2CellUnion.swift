//
//  S2CellUnion.swift
//  s2-swift
//

import Foundation


/// A collection of CellIDs.
/// It is normalized if it is sorted and does not contain redundancy.
/// Specifically, it may not contain the same CellId twice, nor a CellId that is contained by another,
/// nor the four sibling CellIds that are children of a single higher level CellId.
public class CellUnion: S2RegionType {

  var cellIds: [CellId]
  
  // MARK: inits 
  
  public init(cellIds: [CellId]) {
    self.cellIds = cellIds
  }
  
  public init(ids: [UInt64]) {
    cellIds = ids.map { CellId(id: $0) }
  }
  
  // MARK: protocols
  
  public subscript(i: Int) -> CellId {
    get {
      return cellIds[i]
    }
    set(newValue) {
      cellIds[i] = newValue
    }
  }
  
  // MARK: arithmetic
  
  public func add(_ cellId: CellId) {
    cellIds.append(cellId)
  }
  
  public func truncate(_ length: Int) -> CellUnion {
    let s = Array(cellIds[0..<length])
    return CellUnion(cellIds: s)
  }
  
  /// Normalizes the CellUnion.
  public func normalize() {
    cellIds.sort { $0.id < $1.id }
    var output = [CellId]() // the list of accepted cells
    // Loop invariant: output is a sorted list of cells with no redundancy.
    for ci in cellIds {
      // The first two passes here either ignore this new candidate,
      // or remove previously accepted cells that are covered by this candidate.
      // Ignore this cell if it is contained by the previous one.
      // We only need to check the last accepted cell. The ordering of the
      // cells implies containment (but not the converse), and output has no redundancy,
      // so if this candidate is not contained by the last accepted cell
      // then it cannot be contained by any previously accepted cell.
      if output.count > 0 && output[output.count - 1].contains(ci) {
        continue
      }
      // Discard any previously accepted cells contained by this one.
      // This could be any contiguous trailing subsequence, but it can't be
      // a discontiguous subsequence because of the containment property of
      // sorted S2 cells mentioned above.
      var j = output.count - 1 // last index to keep
      while j >= 0 {
        if !ci.contains(output[j]) {
          break
        }
        j -= 1
      }
      output = Array(output.prefix(j+1))
      // See if the last three cells plus this one can be collapsed.
      // We loop because collapsing three accepted cells and adding a higher level cell
      // could cascade into previously accepted cells.
      var ci_ = ci
      while output.count >= 3 {
        // last three
        let fin = Array(output.suffix(3))
        // fast XOR test; a necessary but not sufficient condition
        if fin[0].id^fin[1].id^fin[2].id^ci_.id != 0 {
          break
        }
        // more expensive test; exact.
        // Compute the two bit mask for the encoded child position,
        // then see if they all agree.
        var mask = ci_.lsb() << 1
        mask = ~(mask + mask<<1)
        let should = ci_.id & mask
        if fin[0].id & mask != should || fin[1].id & mask != should || fin[2].id & mask != should || ci_.isFace() {
          break
        }
        output = Array(output.prefix(output.count-3))
        ci_ = ci_.immediateParent() // checked !ci.isFace above
      }
      output.append(ci_)
    }
    cellIds = output
  }

  /// Replaces this CellUnion with an expanded version of the
  /// CellUnion where any cell whose level is less than minLevel or where
  /// (level - minLevel) is not a multiple of levelMod is replaced by its
  /// children, until either both of these conditions are satisfied or the
  /// maximum level is reached.
  public func denormalize(minLevel: Int, levelMod: Int) {
    var denorm = [CellId]()
    for id in cellIds {
      let level = id.level()
      var newLevel = level
      if newLevel < minLevel {
        newLevel = minLevel
      }
      if levelMod > 1 {
        newLevel += (CellId.maxLevel - (newLevel - minLevel)) % levelMod
        if newLevel > CellId.maxLevel {
          newLevel = CellId.maxLevel
        }
      }
      if newLevel == level {
        denorm.append(id)
      } else {
        let end = id.childEnd(newLevel)
        var ci = id.childBegin(newLevel)
        while ci.id != end.id {
          denorm.append(ci)
          ci = ci.next()
        }
      }
    }
    cellIds = denorm
  }

  // MARK: computed members 
  
  var count: Int {
    return cellIds.count
  }
  
  /// Returns an S2Rect that bounds this entity.
  public func rectBound() -> S2Rect {
    var bound = S2Rect.empty
    for cellId in cellIds {
      let c = Cell(id: cellId)
      bound = bound.union(c.rectBound())
    }
    return bound
  }

  public func capBound() -> S2Cap {
    return rectBound().capBound()
  }
  
  // MARK: tests
  
  /// Reports whether this cell union intersects the given cell ID.
  /// This method assumes that the CellUnion has been normalized.
  func intersects(_ cellId: CellId) -> Bool {
    // Find index of array item that occurs directly after our probe cell:
    var i = cellIds.count
    for (index, ci) in cellIds.enumerated() {
      if cellId.id < ci.id {
        i = index
        break
      }
    }
    if i != cellIds.count && cellIds[i].rangeMin().id <= cellId.rangeMax().id {
      return true
    }
    return i != 0 && cellIds[i-1].rangeMax().id >= cellId.rangeMin().id
  }
  
  /// Reports whether the cell union contains the given cell ID.
  /// Containment is defined with respect to regions, e.g. a cell contains its 4 children.
  /// This method assumes that the CellUnion has been normalized.
  func contains(_ cellId: CellId) -> Bool {
    // Find index of array item that occurs directly after our probe cell:
    var i = cellIds.count
    for (index, ci) in cellIds.enumerated() {
      if cellId.id < ci.id {
        i = index
        break
      }
    }
    if i != cellIds.count && cellIds[i].rangeMin().id <= cellId.id {
      return true
    }
    return i != 0 && cellIds[i-1].rangeMax().id >= cellId.id
  }
  
  /// Reports whether this cell union contains the given cell.
  public func contains(_ cell: Cell) -> Bool {
    return contains(cell.id)
  }
  
  /// Reports whether this cell union intersects the given cell.
  public func intersects(_ cell: Cell) -> Bool {
    return intersects(cell.id)
  }
  
}

extension CellUnion: Equatable {

  public static func ==(lhs: CellUnion, rhs: CellUnion) -> Bool {
    return lhs.cellIds == rhs.cellIds
  }

}
