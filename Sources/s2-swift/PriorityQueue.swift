//
//  PriorityQueue.swift
//  s2-swift
//

/// A generic priority queue.
/// This code was inspired by Section 2.4 of Algorithms by Sedgewick & Wayne, 4th Edition
public struct PriorityQueue<T: Comparable> {
  
  fileprivate var heap = [T]()
  private let ordered: (T, T) -> Bool
  
  /// Initialize with optiuoal starting values. 
  /// You can choose the order to create min and max queues.
  public init(ascending: Bool = false, startingValues: [T] = []) {
    if ascending {
      ordered = { $0 > $1 }
    } else {
      ordered = { $0 < $1 }
    }
    for value in startingValues { push(value) }
  }
  
  // MARK: the usual suspects
  
  public var count: Int { return heap.count }
  
  public var isEmpty: Bool { return heap.isEmpty }
  
  /// Insert new element into queue O(lg n)
  public mutating func push(_ element: T) {
    heap.append(element)
    swim(heap.count - 1)
  }
  
  /// Remove and return the element with the highest priority (or lowest if ascending) O(lg n)
  public mutating func pop() -> T? {
    if heap.isEmpty { return nil }
    if heap.count == 1 { return heap.removeFirst() }
    heap.swapAt(0, heap.count - 1)
    let temp = heap.removeLast()
    sink(0)
    return temp
  }
  
  /// Like pop() without removing it. O(1)
  public func peek() -> T? {
    return heap.first
  }
  
  /// Remove all elements
  public mutating func clear() {
    heap.removeAll(keepingCapacity: false)
  }

  // MARK: private sink and swim
  // based on Sedgewick p. 316
  
  private mutating func sink(_ index: Int) {
    var index = index
    while 2 * index + 1 < heap.count {
      var j = 2 * index + 1
      if j < (heap.count - 1) && ordered(heap[j], heap[j + 1]) { j += 1 }
      if !ordered(heap[index], heap[j]) { break }
      heap.swapAt(index, j)
      index = j
    }
  }
  
  private mutating func swim(_ index: Int) {
    var index = index
    while index > 0 && ordered(heap[(index - 1) / 2], heap[index]) {
      heap.swapAt((index - 1) / 2, index)
      index = (index - 1) / 2
    }
  }
  
}

// MARK: IteratorProtocol

extension PriorityQueue: IteratorProtocol {
  public typealias Element = T
  mutating public func next() -> Element? { return pop() }
}

// MARK: Sequence

extension PriorityQueue: Sequence {
  public typealias Generator = PriorityQueue
  public func generate() -> Generator { return self }
}

// MARK: Collection

extension PriorityQueue: Collection {
  public typealias Index = Int
  public var startIndex: Int { return heap.startIndex }
  public var endIndex: Int { return heap.endIndex }
  public subscript(i: Int) -> T { return heap[i] }
  public func index(after: Int) -> Int { return after + 1 }
}

// MARK: CustomStringConvertible, CustomDebugStringConvertible

extension PriorityQueue: CustomStringConvertible, CustomDebugStringConvertible {
  public var description: String { return heap.description }
  public var debugDescription: String { return heap.debugDescription }
}
