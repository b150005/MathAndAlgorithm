//
//  main.swift
//  MathAndAlgorithm
//
//  Created by 伊藤 直輝 on 2022/05/22.
//

import Foundation

/// 001 - Print 5+N
func main_001() -> Void {
  let numOrange: Int = Int(readLine()!)!
  let sum: Int = numOrange + 5
  
  print(sum)
}

/// 10進法→2進法
func decToBin() -> Void {
  guard var dec: Int = Int(readLine()!) else { return }
  
  if (dec == 0) {
    print(dec.description)
    return
  }
  
  var ans: String = ""
  while (dec >= 1) {
    ans.insert( ((dec % 2 == 0) ? "0" : "1"), at: ans.startIndex)
    dec /= 2
  }
  
  print(ans)
}

/// 002 - Sum of 3 Integers
func main_002() -> Void {
  let nums: [Int] = readLine()!.split(separator: " ").map{Int($0)!}
  let sum: Int = nums.reduce(0, +)
  print(sum)
}

/// 003 - Sum of N Integers
func main_003() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let sum: Int = nums.reduce(0, +)
  print(sum)
}

/// 004 - Product of 3 Integers
func main_004() -> Void {
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let p: Int = nums.reduce(1, *)
  print(p)
}

/// 005 - Modulo 100
func main_005() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let sum: Int = nums.reduce(0, +)
  print(sum % 100)
}

/// 006 - Print 2N+3
func main_006() -> Void {
  let num: Int = Int(readLine()!)!
  print(2 * num + 3)
}

/// 007 - Number of Multiples 1
func main_007() -> Void {
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  var ans: Int = 0
  
  for n in 1 ... nums[0] {
    if (n % nums[1] == 0 || n % nums[2] == 0) { ans += 1 }
  }
  
  print(ans)
}

/// 008 - Brute Force 1
func main_008() -> Void {
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  var count: Int = 0
  if (nums[0] < nums[1]) {
    for n in 1 ... nums[0] {
      count += (nums[1] - n >= nums[0]) ? nums[0] : nums[1] - n
    }
  }
  else {
    for n in 1 ... nums[1] - 1 {
      count += nums[1] - n
    }
  }
  
  print(count)
}

/// 009 - Brute Force 2
func main_009() -> Void {
  let s: Int = readLine()!.split(separator: " ").map{ Int($0)! }[1]
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }.filter{ $0 <= s }
  
  var isMakeable: [Bool] = [Bool](repeating: false, count: s + 1)
  // 総和が0の場合は作成可能
  isMakeable[0] = true
  
  for num in nums {
    // 1つ前の要素の時点で作成可能な総和 + 現在の要素 = 現在の要素を使って作成可能な総和
    for makeableNum in (1 ... s).reversed()
    where (makeableNum - num >= 0 && isMakeable[makeableNum - num]) {
      isMakeable[makeableNum] = true
    }
  }
  
  print(isMakeable[s] == true ? "Yes" : "No")
}

/// 010 - Factorial
func main_010() -> Void {
  let num: Int = Int(readLine()!)!
  
  var answer: Int = 1
  for i in 1 ... num {
    answer *= i
  }
  
  print(answer)
}

/// エラトステネスの篩を用いた素数判定(計算量: O(*√N*))
/// - parameter n: 素数判定を行う整数
/// - returns: nが素数である場合は`true`、素数でない場合は`false`
func isPrime(_ n: Int) -> Bool {
  if (n <= 1) { return false }
  
  var i: Int = 2
  while (i * i <= n) {
    if (n % i == 0) { return false }
    
    i += 1
  }
  
  return true
}

/// 011 - Print Prime Numbers
func main_011() -> Void {
  let num: Int = Int(readLine()!)!
  let answer: String = (1 ... num).filter{ isPrime($0) }.map{ String($0) }.joined(separator: " ")
  print(answer)
}

/// 012 - Primality Test
func main_012() -> Void {
  let num: Int = Int(readLine()!)!
  print(isPrime(num) ? "Yes" : "No")
}

/// エラトステネスの篩を用いて昇順の正の約数を取得する(計算量: O(*N* log *N*))
/// - parameter n: 正の約数を求める整数
/// - returns: 昇順に並べ替えられた正の約数の`[Int]`型配列
func getSortedDivisors(_ n: Int) -> [Int] {
  // 約数列挙のアルゴリズム自体の計算量はO(√N)
  var i: Int = 1
  var divisors: [Int] = [Int]()
  while (i * i <= n) {
    if (n % i == 0) {
      divisors.append(i)
      
      if (i * i != n) {
        divisors.append(n / i)
      }
    }
    
    i += 1
  }
  
  // Array#sorted()の計算量がO(N log N)
  return divisors.sorted()
}

/// 013 - Divisor Enumeration
func main_013() -> Void {
  let num: Int = Int(readLine()!)!
  let divisors: [Int] = getSortedDivisors(num)
  divisors.forEach{ print($0) }
}

/// 重複のない素因数を取得する(計算量: O(*N* log *N*))
/// - parameter n: 素因数を求める整数
/// - returns: 昇順に並べ替えられた重複のない素因数の`[Int]`型配列
func getUniquePrimeFactors(_ n: Int) -> [Int] {
  let divisors: [Int] = getDivisors(n)
  
  return divisors.filter{ isPrime($0) }
}

/// 全ての素因数を取得する(計算量: O(*N* log *N*))
/// - parameter n: 素因数を求める整数
/// - returns: 昇順に並べ替えられた全ての素因数の`[Int]`型配列
func getAllPrimeFactors(_ n: Int) -> [Int] {
  let uniquePrimeFactors: [Int] = getUniquePrimeFactors(n)
  
  var num: Int = n
  var primeFactors: [Int] = [Int]()
  while (num != 1) {
    for p in uniquePrimeFactors where (num % p == 0) {
      primeFactors.append(p)
      num /= p
      break
    }
  }
  
  return primeFactors
}

/// 014 - Factorization
func main_014() -> Void {
  let num: Int = Int(readLine()!)!
  let primeFactors: [Int] = getAllPrimeFactors(num)
  print(primeFactors.map{ String($0) }.joined(separator: " "))
}

/// ユークリッドの互除法を用いて2つの整数の最大公約数を取得する(計算量: O(log *N*))
/// - parameter num1: 最大公約数を求める整数
/// - parameter num2: 最大公約数を求める整数
/// - returns: `num1, num2`の最大公約数
func getGcd(_ num1: Int, _ num2: Int) -> Int {
  var (p, q): (Int, Int) = (num1, num2)
  while (p >= 1 && q >= 1) {
    if (p >= q) { p %= q }
    else { q %= p }
  }
  
  if (p >= q) { return p }
  else { return q }
}

/// 再帰的にユークリッドの互除法を用いて2つの整数の最大公約数を取得する(計算量: O(log *N*))
/// - parameter num1: 最大公約数を求める整数
/// - parameter num2: 最大公約数を求める整数
/// - returns: `num1, num2`の最大公約数
func getGcdRecursively(_ num1: Int, _ num2: Int) -> Int {
  // ベースケース
  if (num2 == 0) { return num1 }
  return getGcdRecursively(num2, num1 % num2)
}

/// 015 - Greatest Common Divisor
func main_015() -> Void {
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  print(getGcd(nums[0], nums[1]))
}

/// 016 - Greatest Common Divisor of N Integers
func main_016() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  // 初期値を最初の要素にする
  let gcd: Int = nums.reduce(nums.first!, getGcd(_:_:))
  print(gcd)
}

/// ユークリッドの互除法を用いて2つの整数の最小公倍数を取得する(計算量: O(log *N*))
/// - parameter num1: 最小公倍数を求める整数
/// - parameter num2: 最小公倍数を求める整数
/// - returns: `num1, num2`の最小公倍数
func getLcm(_ num1: Int, _ num2: Int) -> Int {
  let gcd: Int = getGcd(num1, num2)
  // 除算を先に行うことでInt型のオーバーフローを避ける
  return num1 / gcd * num2
}

/// 017 - Least Common Multiple of N Integers
func main_017() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  // 初期値を1にする
  let lcm: Int = nums.reduce(1, getLcm(_:_:))
  print(lcm)
}

/// 018 - Convenience Store 1
func main_018() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  // 値段ごとの商品数を求める
  var counts: [Int] = [Int](repeating: 0, count: 4)
  nums.forEach{ counts[$0 / 100 - 1] += 1 }
  print(counts[0] * counts[3] + counts[1] * counts[2])
}

/// 019 - Choose Cards 1
func main_019() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  // 色ごとの枚数を求める
  var counts: [Int] = [Int](repeating: 0, count: 3)
  nums.forEach{ counts[$0 - 1] += 1 }
  print(
    counts[0] * (counts[0] - 1) / 2
    + counts[1] * (counts[1] - 1) / 2
    + counts[2] * (counts[2] - 1) / 2
  )
}

/// 020 - Choose Cards 2
func main_020() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  var count: Int = 0
  for i in 0 ..< nums.count - 4 {
    for j in i + 1 ..< nums.count - 3 {
      for k in j + 1 ..< nums.count - 2 {
        for l in k + 1 ..< nums.count - 1 {
          for m in l + 1 ..< nums.count
          where (nums[i] + nums[j] + nums[k] + nums[l] + nums[m] == 1000) {
            count += 1
          }
        }
      }
    }
  }
  
  print(count)
}



/// 順列の総数を取得する(計算量: O(*N*))
/// - parameter m: 全体の要素数
/// - parameter n: 抽出する要素数
/// - returns: *P(m, n)* の順列の総数
func getPermutation(_ m: Int, _ n: Int) -> Int {
  if (m < n) { return 0 }
  if (m == 0 || n == 0) { return 1 }
  
  var p: Int = 1
  for i in stride(from: m, through: (m - n + 1), by: -1) {
    p *= i
  }
  
  return p
}

/// 組合せの総数を取得する(計算量: O(*N*))
/// - parameter m: 全体の要素数
/// - parameter n: 抽出する要素数
/// - returns: *C(m, n)* の組合せの総数
func getCombination(_ m: Int, _ n: Int) -> Int {
  if (m < n) { return 0 }
  if (m == 0 || n == 0) { return 1 }
  
  var p: Int = getPermutation(m, n)
  for i in 1 ... n {
    p /= i
  }
  
  return p
}

/// 021 - Combination Easy
func main_021() -> Void {
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  print(getCombination(nums[0], nums[1]))
}

/// 022 - Choose Cards 3
func main_022() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  // 値ごとの枚数を求める
  let sum: Int = 100000
  var counts: [Int] = [Int](repeating: 0, count: sum)
  nums.forEach{ counts[$0] += 1 }
  
  let halfCount: Int = counts[sum / 2]
  var answer: Int = (halfCount % 2 == 0)
  ? halfCount / 2 * (halfCount - 1)
  : (halfCount - 1) / 2 * halfCount
  
  for lth in 1 ..< (sum / 2) {
    answer += counts[lth] * counts[sum - lth]
  }
  
  print(answer)
}

/// 023 - Dice Expectation
func main_023() -> Void {
  let n: Double = Double(readLine()!)!
  let blues: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let reds: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  let sum: Int = blues.reduce(0, +) + reds.reduce(0, +)
  print(Double(sum) / n)
}

/// 024 - Answer Exam Randomly
func main_024() -> Void {
  let n: Int = Int(readLine()!)!
  var questions: [[Double]] = [[Double]]()
  for _ in 0 ..< n {
    questions.append(readLine()!.split(separator: " ").map{ Double($0)! })
  }
  
  var expectation: Double = 0
  for question in questions {
    expectation += question[1] / question[0]
  }
  
  print(expectation)
}

/// 025 - Jiro's Vacation
func main_025() -> Void {
  _ = readLine()
  let sum1_2: Int = readLine()!.split(separator: " ").map{ Int($0)! }.reduce(0, +)
  let sum3_6: Int = readLine()!.split(separator: " ").map{ Int($0)! }.reduce(0, +)
  
  print(Double(sum1_2 + sum3_6 * 2) / 3)
}

/// 026 - Coin Gacha
func main_026() -> Void {
  let n: Int = Int(readLine()!)!
  
  var probability: Double = 0
  for i in 1 ... n { probability += 1 / Double(i) }
  probability *= Double(n)
  
  print(probability)
}

/// 2つの配列をマージする
/// - parameter left: 左半分の配列
/// - parameter right: 右半分の配列
/// - returns: `left, right`がマージされた`[Int]`型のソート済配列
func merge(_ left: [Int], _ right: [Int]) -> [Int] {
  var result: [Int] = [Int]()
  
  // Merge操作
  var (l, r): (Int, Int) = (0, 0)
  while (l < left.count || r < right.count) {
    // 左側の配列の値全てが走査済である場合は右側の配列から値を追加
    if (l == left.count) {
      result.append(right[r])
      r += 1
    }
    // 右側の配列の値全てが走査済である場合は左側の配列から値を追加
    else if (r == right.count) {
      result.append(left[l])
      l += 1
    }
    // 両方の配列の走査を終えていない場合は両方の配列の値同士を比較
    else {
      // 左側の配列の値の方が小さい場合
      if (left[l] <= right[r]) {
        result.append(left[l])
        l += 1
      }
      // 右側の配列の値の方が小さい場合
      else {
        result.append(right[r])
        r += 1
      }
    }
  }
  
  return result
}

/// 配列に対してマージソートを行う
/// - parameter list: マージソートを行う配列
/// - returns: マージソートされた`[Int]`型配列
func mergeSort(_ list: [Int]) -> [Int] {
  // 要素数が1以下の場合はソートの必要がない
  if (list.count <= 1) { return list }
  
  // 配列を分割(分割統治法)
  let (l, r): ([Int], [Int]) = (
    [Int](list[0 ..< list.count / 2]),
    [Int](list[list.count / 2 ..< list.count])
  )
  
  // Merge操作
  return merge(mergeSort(l), mergeSort(r))
}

/// 027 - Sorting
func main_027() -> Void {
  _ = readLine()
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let sorted: [Int] = mergeSort(nums)
  print(sorted.map{String($0)}.joined(separator: " "))
}

/// 028 - Frog 1
func main_028() -> Void {
  let n: Int = Int(readLine()!)!
  let heights: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  var costs: [Int] = [Int]()
  var (p, q): (Int, Int) = (0, 0)
  for i in 0 ..< n {
    if (i == 0) { costs.append(0) }
    if (i == 1) { costs.append(abs(heights[i] - heights[i - 1])) }
    if (i >= 2) {
      // 1つ前の足場
      p = costs[i - 1] + abs(heights[i] - heights[i - 1])
      // 2つ前の足場
      q = costs[i - 2] + abs(heights[i] - heights[i - 2])
      costs.append(min(p, q))
    }
  }
  
  print(costs[n - 1])
}

/// 029 - Climb Stairs
func main_029() -> Void {
  let n: Int = Int(readLine()!)!
  var counts: [Int] = [Int](repeating: 0, count: n + 1)
  for i in 0 ... n {
    if (i <= 1) { counts[i] = 1 }
    else { counts[i] += counts[i - 1] + counts[i - 2] }
  }
  
  print(counts[n])
}

/// 030 - Knapsack 1
func main_030() -> Void {
  let input1: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, maxWeight): (Int, Int) = (input1[0], input1[1])
  var values: [Int] = [Int](repeating: 0, count: maxWeight + 1)
  
  // 商品を1つずつ走査
  for _ in 0 ..< n {
    let input2: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    let (weight, value): (Int, Int) = (input2[0], input2[1])
    
    // (新商品の重量weight) ≦ (設定重量w)であれば、重量の合計が設定重量wとなる
    // 「新商品を取らない、現時点での最適な取り方」と「新商品を取る、新たな取り方」を比較
    // ※ 値が再代入されたvaluesのインデックスwがループ処理の中で(w - weight)として再登場しないよう
    //   設定重量wは「重い順」で考える
    for w in (0 ... maxWeight).reversed() where (w - weight >= 0) {
      values[w] = max(values[w], values[w - weight] + value)
    }
  }
  
  print(values[maxWeight])
}

/// 031 - Taro's Vacation
func main_031() -> Void {
  let n: Int = Int(readLine()!)!
  let scores: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  var maxScores: [Int] = [Int](repeating: 0, count: n + 1)
  for d in 1 ... n {
    if (d == 1) { maxScores[d] = scores[d - 1] }
    else { maxScores[d] = max(maxScores[d - 1], maxScores[d - 2] + scores[d - 1]) }
  }
  
  print(maxScores[n])
}

/// 032 - Binary Search
func main_032() -> Void {
  let x: Int = readLine()!.split(separator: " ").map { Int($0)! }[1]
  let nums: [Int] = readLine()!.split(separator: " ").map { Int($0)! }.sorted()
  
  var (left, mid, right): (Int, Int, Int) = (1, 0, nums.count)
  while (left <= right) {
    mid = (left + right) / 2
    
    if (nums[mid] == x) {
      print("Yes")
      return
    }
    else if (nums[mid] > x) { right = mid - 1 }
    else if (nums[mid] < x) { left = mid + 1 }
  }
  
  print("No")
}

/// 033 - Distance
func main_033() -> Void {
  let aInput: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let bInput: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let cInput: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  // 点A, B, Cの座標
  let (ax, ay): (Int, Int) = (aInput[0], aInput[1])
  let (bx, by): (Int, Int) = (bInput[0], bInput[1])
  let (cx, cy): (Int, Int) = (cInput[0], cInput[1])
  
  // ベクトルの成分
  let (BAx, BAy): (Int, Int) = (ax - bx, ay - by)
  let (BCx, BCy): (Int, Int) = (cx - bx, cy - by)
  let (CAx, CAy): (Int, Int) = (ax - cx, ay - cy)
  let (CBx, CBy): (Int, Int) = (bx - cx, by - cy)
  
  // ∠Bが90°より大きい ⇔ 点Aの最短は点B
  var answer: Double = 0
  if (BAx * BCx + BAy * BCy < 0) {
    answer = sqrt(Double(BAx * BAx + BAy * BAy))
  }
  // ∠Cが90°より大きい ⇔ 点Aの最短は点C
  else if (CAx * CBx + CAy * CBy < 0) {
    answer = sqrt(Double(CAx * CAx + CAy * CAy))
  }
  // ∠B,Cが90°以下 ⇔ 点Aの最短は線分BC上の点
  else {
    answer = Double(abs(CAx * CBy - CAy * CBx)) / sqrt(Double(BCx * BCx + BCy * BCy))
  }
  
  print(answer)
}

/// x座標でソート済の点群に対して、最近接の2点間の距離の2乗を取得する(計算量: O(*N* log *N*))
/// - parameter xSortedPoints: x座標の値でソートされた`[[Int]]`型の点群
/// - returns: `xSortedPoints`内での最近接の2点間の距離の2乗
func getNearestSquaredDistance(_ xSortedPoints: [[Int]]) -> Int {
  /*
   点が1つしかない場合:
   (最短距離の2乗) = Int.max
   */
  if (xSortedPoints.count == 1) { return Int.max }
  
  /*
   点が2つ以上ある場合:
   (最短距離の2乗) = min(同領域内のx,y座標が近接する2点間の最短距離の2乗,
   他領域同士のx, y座標が近接する2点間の最短距離の2乗)
   */
  let mid: Int = xSortedPoints.count / 2
  let midX: Int = xSortedPoints[mid][0]
  
  // x座標を基準とした点の分割 → 少なくともx座標が近接した点同士でしか比較しない
  var nearerSquaredDistance: Int = min(
    getNearestSquaredDistance([[Int]](xSortedPoints[0 ..< mid])),
    getNearestSquaredDistance([[Int]](xSortedPoints[mid ..< xSortedPoints.count]))
  )
  
  // x座標が近接した点同士のうち、y座標も近接している点同士の距離(の2乗)を求める
  let ySortedPoints: [[Int]] = xSortedPoints.sorted{ $0[1] < $1[1] }
  for i in 0 ..< (ySortedPoints.count - 1) {
    // 基準点Aは、x座標が中央値である点Bとのx座標の距離が最短距離未満である必要がある
    let midDx: Int = ySortedPoints[i][0] - midX
    if (midDx * midDx >= nearerSquaredDistance) { continue }
    
    for j in (i + 1) ..< ySortedPoints.count {
      // 基準点Aとのy座標の距離dyが最短距離以上となった時点で、
      // 以降の要素と基準点Aのy座標の距離はdy以上であることが確定しているため次の基準点の走査に移る
      let dy: Int = ySortedPoints[i][1] - ySortedPoints[j][1]
      let squaredDy: Int = dy * dy
      if (squaredDy >= nearerSquaredDistance) { break }
      
      let dx: Int = ySortedPoints[i][0] - ySortedPoints[j][0]
      let squaredDx: Int = dx * dx
      nearerSquaredDistance = min(nearerSquaredDistance, squaredDx + squaredDy)
    }
  }
  
  return nearerSquaredDistance
}

/// 034 - Nearest Points
func main_034() -> Void {
  let n: Int = Int(readLine()!)!
  var points: [[Int]] = [[Int]]()
  for _ in 0 ..< n {
    points.append(readLine()!.split(separator: " ").map{ Int($0)! })
  }
  
  // x座標順にソート
  points.sort{ $0[0] < $1[0] }
  
  let nearestSquaredDistance: Int = getNearestSquaredDistance(points)
  print(sqrt(Double(nearestSquaredDistance)))
}

/// 035 - Two Circles
func main_035() -> Void {
  let circle1: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (x1, y1, r1): (Int, Int, Int) = (circle1[0], circle1[1], circle1[2])
  let circle2: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (x2, y2, r2): (Int, Int, Int) = (circle2[0], circle2[1], circle2[2])
  
  let squaredDx: Int = (x1 - x2) * (x1 - x2)
  let squaredDy: Int = (y1 - y2) * (y1 - y2)
  let d: Double = sqrt(Double(squaredDx + squaredDy))
  let rSum: Double = Double(r1 + r2)
  let dr: Double = Double(abs(r1 - r2))
  switch d {
      // 一方の円が他方の円を完全に含み、2つの円は接していない
    case let d where (d < dr): print("1")
      // 一方の円が他方の円を完全に含み、2つの円は接している
    case let d where (d == dr): print("2")
      // 2つの円が互いに交差する
    case let d where (d < rSum && d > dr): print("3")
      // 2つの円の内部に共通部分は存在しないが、2つの円は接している
    case let d where (d == rSum): print("4")
      // 2つの円の内部に共通部分は存在せず、2つの円は接していない
    case let d where (d > rSum): print("5")
    default: break
  }
}

/// 時・分・秒を表す列挙型
enum Hand {
  /// 時
  case hour
  /// 分
  case minute
  /// 秒
  case second
}

/// 時刻をもとに短針・長針・秒針のベクトルの成分を取得する
/// - note: 針の長さを1, 正午を基準に針が進んだ角度をθとすると、
///         針の先端のx座標はcos(90° - θ) = sinθ, y座標はsin(90° - θ) = cosθ で表される
/// - parameters:
///   - hour: 時
///   - min: 分
///   - sec: 秒
///   - len: 針の長さ
///   - type: 短針・長針・秒針を示す`Hand`型の列挙ケース
/// - returns: `[Double]`型で表される短針・長針・秒針のベクトルのx, y成分
func getHandVector(_ hour: Int = 0, _ min: Int = 0, _ sec: Int = 0, _ len: Int = 1, _ type: Hand) -> [Double] {
  var degree: Double = 0
  switch type {
    case .hour:
      degree = Double((hour % 12) * 60 * 60 + min * 60 + sec) / 120
    case .minute:
      degree = Double(min * 60 + sec) / 10
    case .second:
      degree = Double(sec) * 6
  }
  
  // 度数法[°] → 弧度法[rad] への変換
  let radian: Double = degree * Double.pi / 180
  
  return [Double(len) * sin(radian), Double(len) * cos(radian)]
}

/// 036 - : (Colon)
func main_036() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (a, b, h, m): (Int, Int, Int, Int) = (input[0], input[1], input[2], input[3])
  
  // 短針・長針のベクトルの成分(=先端の平面座標)を取得
  let vecHour: [Double] = getHandVector(h, m, 0, a, .hour)
  let vecMin: [Double] = getHandVector(h, m, 0, b, .minute)
  
  // 短針・長針の2つの端点の差分を算出
  let dx: Double = vecHour[0] - vecMin[0]
  let dy: Double = vecHour[1] - vecMin[1]
  
  print(sqrt(dx * dx + dy * dy))
}

/// 直線で2つに分割される領域に対して、与えられた2点がそれぞれ(直線上を含む)異なる領域に存在するかどうか
/// - parameters:
///  - endPoints: 領域を分割する直線を形成する2つの点
///  - comparedPoints: 異なる領域に存在するかどうかを調べる2つの点
/// - returns: 2点がそれぞれ異なる領域に存在する場合は`true`, 同じ領域に存在する場合は`false`
func existsSeparately(_ endPoints: [[Int]], _ comparedPoints: [[Int]]) -> Bool {
  let (start, end): ([Int], [Int]) = (endPoints[0], endPoints[1])
  let (comp1, comp2): ([Int], [Int]) = (comparedPoints[0], comparedPoints[1])
  
  let vecBase: [Int] = [end[0] - start[0], end[1] - start[1]]
  let vecCompared1: [Int] = [comp1[0] - start[0], comp1[1] - start[1]]
  let vecCompared2: [Int] = [comp2[0] - start[0], comp2[1] - start[1]]
  
  // 基準ベクトルと残りのベクトルの外積同士の積が正の場合は同じ領域に存在
  // → オーバーフロー対策のため、外積の符号(正・負・0)だけ区別できるようにする
  let (op1, op2): (Int, Int)
  switch (vecBase[0] * vecCompared1[1] - vecBase[1] * vecCompared1[0]) {
    case let op where (op < 0): op1 = -1
    case let op where (op > 0): op1 = 1
    default: op1 = 0
  }
  switch (vecBase[0] * vecCompared2[1] - vecBase[1] * vecCompared2[0]) {
    case let op where (op < 0): op2 = -1
    case let op where (op > 0): op2 = 1
    default: op2 = 0
  }
  
  /*
   op1 = op2 = 0となる(⇔全ての点が一直線上に存在する)場合は2つの線分が重なっているかどうかで交差を判定
   直線を形成する2点のうちxまたはy座標が小さい順にそれぞれ点A, Bとし、
   残りの2点のうちxまたはy座標が小さい順にそれぞれ点C, Dとすると、
   A < B < C (< D)または C < D < A (< B) となる場合は交差しない
   */
  if (op1 == 0 && op2 == 0) {
    // 全ての点がy軸と平行に並んでいなければx座標を基準にソート
    if (start[0] != end[0]) {
      let xSortedEndPoints: [[Int]] = endPoints.sorted{ $0[0] < $1[0] }
      let (ax, bx): (Int, Int) = (xSortedEndPoints[0][0], xSortedEndPoints[1][0])
      let xSortedCompPoints: [[Int]] = comparedPoints.sorted{ $0[0] < $1[0] }
      let (cx, dx): (Int, Int) = (xSortedCompPoints[0][0], xSortedCompPoints[1][0])
      
      if (ax < cx && bx < cx) || (cx < ax && dx < ax) { return false }
    }
    // 全ての点がy軸と平行に並んでいればy座標を基準にソート
    else {
      let ySortedEndPoints: [[Int]] = endPoints.sorted{ $0[1] < $1[1] }
      let (ay, by): (Int, Int) = (ySortedEndPoints[0][1], ySortedEndPoints[1][1])
      let ySortedCompPoints: [[Int]] = comparedPoints.sorted{ $0[1] < $1[1] }
      let (cy, dy): (Int, Int) = (ySortedCompPoints[0][1], ySortedCompPoints[1][1])
      
      if (ay < cy && by < cy) || (cy < ay && dy < ay) { return false }
    }
  }
  
  if (op1 * op2 <= 0) { return true }
  else { return false }
}

/// 037 - Intersection
func main_037() -> Void {
  let n: Int = 4
  var points: [[Int]] = [[Int]](repeating: [], count: n)
  for i in 0 ..< n {
    points[i] = readLine()!.split(separator: " ").map{ Int($0)! }
  }
  
  if (existsSeparately([points[0], points[1]], [points[2], points[3]]) == true &&
      existsSeparately([points[2], points[3]], [points[0], points[1]]) == true) {
    print("Yes")
  }
  else {
    print("No")
  }
}

/// 038 - How Many Guests?
func main_038() -> Void {
  let q: Int = readLine()!.split(separator: " ").map{ Int($0)! }[1]
  // 0〜N日目の累積和を算出
  var cumSums: [Int] = [0]
  readLine()!.split(separator: " ").map{ Int($0)! }.forEach{
    cumSums.append(cumSums.last! + $0)
  }
  
  var questions: [[Int]] = [[Int]]()
  for _ in 0 ..< q {
    questions.append(readLine()!.split(separator: " ").map{ Int($0)! })
  }
  
  for question in questions {
    let (l, r): (Int, Int) = (question[0], question[1])
    print(cumSums[r] - cumSums[l - 1])
  }
}

/// 039 - Snowy Days
func main_039() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, q): (Int, Int) = (input[0], input[1])
  /// 1〜N日目の積雪量の階差
  var differences: [Int] = [Int](repeating: 0, count: n)
  for _ in 0 ..< q {
    let news: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    let (l, r, snow): (Int, Int, Int) = (news[0], news[1], news[2])
    differences[l - 1] += snow
    if (r != n) {
      differences[r] -= snow
    }
  }
  
  var output: String = ""
  for i in 1 ..< n {
    // i(=次の日)の積雪量の階差の符号で(i - 1)日目とi日目の積雪量を比較
    switch (differences[i]) {
      case let d where (d < 0): output += ">"
      case let d where (d > 0): output += "<"
      default: output += "="
    }
  }
  
  print(output)
}

/// 040 - Travel
func main_040() -> Void {
  _ = readLine()!
  var sums: [Int] = [0, 0]
  readLine()!.split(separator: " ").map{ Int($0)! }.forEach{ sums.append(sums.last! + $0) }
  
  let m: Int = Int(readLine()!)!
  var stations: [Int] = [Int]()
  for _ in 0 ..< m {
    stations.append(Int(readLine()!)!)
  }
  
  var answer: Int = 0
  for i in 1 ..< stations.count {
    answer += abs(sums[stations[i]] - sums[stations[i - 1]])
  }
  
  print(answer)
}

/// 041 - Convenience Store 2
func main_041() -> Void {
  let (t, n): (Int, Int) = (Int(readLine()!)!, Int(readLine()!)!)
  var changes: [Int] = [Int](repeating: 0, count: t + 1)
  for _ in 0 ..< n {
    let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    let (start, end): (Int, Int) = (input[0], input[1])
    
    changes[start] += 1
    changes[end] -= 1
  }
  
  var sum: Int = 0
  for i in 0 ..< t {
    sum += changes[i]
    print(sum)
  }
}

/// ニュートン法を利用した平方根の取得
func getSqrt(_ d: Double, _ precision: Int) -> Double {
  var a: Double = d
  for _ in 0 ..< precision {
    let new: Double = (d + a * a) / (2 * a)
    a = new
  }
  return a
}

/// 素因数とその指数を取得する
/// - parameter n: 素因数分解を行う正の整数
/// - returns: [素因数, 指数]が列挙された`[[Int]]`型配列
/// - note: `n = 1`の場合は例外的に`[[1, 0]]`を返却
func getPrimeFactorization(_ n: Int) -> [[Int]] {
  // n = 1の場合は例外的に1^0として返却
  if (n == 1) { return [[1, 0]] }
  
  // 2 ≦ i ≦ [√n] となる整数iでnum(=n)を除算し、割り切れなくなった時点で素因数iとその指数を配列に追加
  var result: [[Int]] = [[Int]]()
  var (num, i): (Int, Int) = (n, 2)
  while (i * i <= num) {
    if (num % i == 0) {
      var exponent: Int = 0
      while (num % i == 0) {
        exponent += 1
        num /= i
      }
      result.append([i, exponent])
    }
    i += 1
  }
  
  // この時点でnumが1でなければnumは素数
  if (num != 1) { result.append([num, 1]) }
  return result
}

/// エラトステネスの篩を用いて正の約数を取得する(計算量: O(√*N*))
/// - parameter n: 正の約数を求める整数
/// - returns: 正の約数の`[Int]`型配列
func getDivisors(_ n: Int) -> [Int] {
  var i: Int = 1
  var divisors: [Int] = [Int]()
  while (i * i <= n) {
    if (n % i == 0) {
      divisors.append(i)
      
      if (i * i != n) {
        divisors.append(n / i)
      }
    }
    
    i += 1
  }
  
  return divisors
}

/// 042 - Sum of Divisors
func main_042() -> Void {
  let n: Int = Int(readLine()!)!
  
  // iの倍数全てに1を加算
  var divisorCounts: [Int] = [Int](repeating: 0, count: n + 1)
  for i in 1 ... n {
    for j in stride(from: i, through: n, by: i) {
      divisorCounts[j] += 1
    }
  }
  
  var answer: Int = 0
  for i in 1 ... n {
    answer += i * divisorCounts[i]
  }
  print(answer)
}

/// 深さ優先探索を用いてグラフが連結であるかどうかを判定(計算量: O(*M + N*); *M*はエッジ数、*N*はノード数)
/// - parameter nodes: `[[Int]]`型の隣接リスト表現(最初のインデックスは判定に使用しない)
/// - returns: 連結である場合は`true`、そうでない場合は`false`
func isConnected(_ nodes: [[Int]]) -> Bool {
  var isVisited: [Bool] = [Bool](repeating: false, count: nodes.count)
  
  /// 深さ優先探索
  /// - parameter node: 現在のノード番号
  func dfs(_ node: Int) -> Void {
    // 現在のノードを到達済ノードとしてセット
    isVisited[node] = true
    
    // 隣接ノードのうち、未到達ノードに対して再帰処理を実行
    for next in nodes[node] {
      if (isVisited[next] == false) { dfs(next) }
    }
  }
  
  // ノード1から深さ優先探索
  dfs(1)
  
  // 全ての頂点が到達済であれば連結であり、そうでなければ連結でない
  for i in 1 ..< isVisited.count {
    if (isVisited[i] == false) { return false }
  }
  
  return true
}

/// 043 - Is It Connected?
func main_043() -> Void {
  let nm: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, m): (Int, Int) = (nm[0], nm[1])
  
  var nodes: [[Int]] = [[Int]](repeating: [], count: n + 1)
  for _ in 0 ..< m {
    let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    nodes[input[0]].append(input[1])
    nodes[input[1]].append(input[0])
  }
  
  if (isConnected(nodes) == true) {
    print("The graph is connected.")
  }
  else {
    print("The graph is not connected.")
  }
}

/// 幅優先探索を用いてノード群の最短経路長を取得(計算量: O(*M + N*); *M*はエッジ数、*N*はノード数)
/// - parameter nodes: `[[Int]]`型の隣接リスト表現(最初のインデックスはダミー)
/// - returns: 最短経路長を格納した`[Int]`型配列(最初のインデックスはダミー)
func getShortestDistances(_ nodes: [[Int]], _ start: Int) -> [Int] {
  var distances: [Int] = [Int](repeating: -1, count: nodes.count)
  
  /// 幅優先探索
  /// - parameter start: 始点のノード番号
  func bfs(_ start: Int) -> Void {
    var queue: [Int] = [Int]()
    // 始点の最短経路長を0にセット、キューに追加
    distances[start] = 0
    queue.append(start)
    
    while (queue.isEmpty == false) {
      // キューから先頭要素(ノード番号)を抽出
      let pop: Int = queue.removeFirst()
      
      // 抽出したノードと隣接するノードを走査
      for i in 0 ..< nodes[pop].count {
        let neighbor: Int = nodes[pop][i]
        
        // 隣接ノードが未到達ノードであれば最短経路長を更新、キューに追加
        if (distances[neighbor] == -1) {
          distances[neighbor] = distances[pop] + 1
          queue.append(neighbor)
        }
      }
    }
  }
  
  // ノード1を始点とした幅優先探索
  bfs(start)
  
  return distances
}

/// 044 - Shortest Path Problem
func main_044() -> Void {
  let nm: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, m): (Int, Int) = (nm[0], nm[1])
  
  var nodes: [[Int]] = [[Int]](repeating: [], count: n + 1)
  for _ in 0 ..< m {
    let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    nodes[input[0]].append(input[1])
    nodes[input[1]].append(input[0])
  }
  
  let distances: [Int] = getShortestDistances(nodes, 1)
  for i in 1 ... n {
    print(distances[i])
  }
}

/// 045 - Easy Graph Problem
func main_045() -> Void {
  let nm: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, m): (Int, Int) = (nm[0], nm[1])
  
  var nodes: [[Int]] = [[Int]](repeating: [], count: n + 1)
  for _ in 0 ..< m {
    let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    nodes[input[0]].append(input[1])
    nodes[input[1]].append(input[0])
  }
  
  var answer: Int = 0
  for i in 1 ..< nodes.count {
    if (nodes[i].filter{ $0 < i }.count == 1) {
      answer += 1
    }
  }
  
  print(answer)
}

/// 046 - 幅優先探索
func main_046() -> Void {
  let rc: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (row, col): (Int, Int) = (rc[0], rc[1])
  
  let start: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let goal: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  var squares: [[Character]] = [[Character]](repeating: [Character](repeating: Character("-"), count: col + 1), count: row + 1)
  var distances: [[Int]] = [[Int]](repeating: [Int](repeating: -1, count: col + 1), count: row + 1)
  
  for r in 1 ... row {
    let input: String = readLine()!
    
    for c in 1 ... col {
      squares[r][c] = input[input.index(input.startIndex, offsetBy: c - 1)]
    }
  }
  
  // 幅優先探索
  var queue: [[Int]] = [[Int]]()
  // 始点の最短経路長を0にセット、キューに追加
  distances[start[0]][start[1]] = 0
  queue.append([start[0], start[1]])
  
  while (queue.isEmpty == false) {
    // キューから先頭要素(行・列番号)を抽出
    let pop: [Int] = queue.removeFirst()
    let (pR, pC): (Int, Int) = (pop[0], pop[1])
    
    // 隣接しうるマス
    let neighbors: [[Int]] = [
      // 上
      [pR - 1, pC],
      // 下
      [pR + 1, pC],
      // 左
      [pR, pC - 1],
      // 右
      [pR, pC + 1]
    ]
    
    // 隣接マスが移動可能マスかつ未到達マスである場合は最短経路長を更新、マスを追加
    for neighbor in neighbors {
      let (nR, nC): (Int, Int) = (neighbor[0], neighbor[1])
      if (squares[nR][nC] == "." && distances[nR][nC] == -1) {
        distances[nR][nC] = distances[pR][pC] + 1
        queue.append([nR, nC])
      }
    }
  }
  
  print(distances[goal[0]][goal[1]])
}

/// 深さ優先探索を用いて二部グラフであるかどうかを判定(計算量: O(*M + N*); *M*はエッジ数、*N*はノード数)
/// - parameter nodes: `[[Int]]`型の隣接リスト表現(最初のインデックスは判定に使用しない)
/// - returns: 二部グラフである場合は`true`、そうでない場合は`false`
func isBipartiteGraph(_ nodes: [[Int]]) -> Bool {
  var values: [Int] = [Int](repeating: -1, count: nodes.count)
  values[0] = -2
  
  /// 深さ優先探索
  /// - parameters:
  ///   - node: 現在のノード番号
  ///   - value: 現在のノードにセットする値
  func dfs(_ node: Int, _ value: Int) -> Void {
    // 現在のノードに値をセット
    values[node] = value
    
    // 隣接ノードのうち、未到達ノードに対して再帰処理を実行
    for next in nodes[node] {
      if (values[next] == -1) {
        if (value == 0) { dfs(next, 1) }
        if (value == 1) { dfs(next, 0) }
      }
    }
  }
  
  // ノード1を始点とした深さ優先探索
  dfs(1, 0)
  
  // 全てのノードが走査済になるまで深さ優先探索を実行
  while let start = values.firstIndex(of: -1) {
    dfs(start, 0)
  }
  
  // 全ての隣接ノードの値を走査し、0と1で交互になっていなければ二部グラフでない
  for node in 1 ..< nodes.count {
    for next in nodes[node] {
      if (values[node] == 0 && values[next] != 1) { return false }
      else if (values[node] == 1 && values[next] != 0) { return false }
    }
  }
  
  return true
}

/// 047 - Bipartite Graph
func main_047() -> Void {
  let nm: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, m): (Int, Int) = (nm[0], nm[1])
  
  var nodes: [[Int]] = [[Int]](repeating: [], count: n + 1)
  for _ in 0 ..< m {
    let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
    nodes[input[0]].append(input[1])
    nodes[input[1]].append(input[0])
  }
  
  if (isBipartiteGraph(nodes) == true) { print("Yes") }
  else { print("No") }
}

/// 優先度付きキュー
/// - authors: [koher](https://github.com/koher)
/// - seealso: [swift-atcoder-support/Sources/AtCoderSupport/PriorityQueue.swift](https://github.com/koher/swift-atcoder-support/blob/main/Sources/AtCoderSupport/PriorityQueue.swift)
struct PriorityQueue<Element> {
  /// 優先度付きキューの要素
  private var elements: [Element] = []
  /// 順序付け方式
  private let areInIncreasingOrder: (Element, Element) -> Bool
  
  init<S>(_ elements: S, by areInIncreasingOrder: @escaping (Element, Element) -> Bool) where S: Sequence, S.Element == Element {
    self.areInIncreasingOrder = areInIncreasingOrder
    for element in elements {
      append(element)
    }
  }
  
  init(by areInIncreasingOrder: @escaping (Element, Element) -> Bool) {
    self.init(EmptyCollection(), by: areInIncreasingOrder)
  }
  
  var isEmpty: Bool { elements.isEmpty }
  var count: Int { elements.count }
  var first: Element? { elements.first }
  
  /// ヒープを用いたキューへの値の追加
  mutating func append(_ element: Element) {
    var i = elements.count
    elements.append(element)
    elements.withUnsafeMutableBufferPointer { elements in
      while i > 0 {
        let parentIndex = (i - 1) >> 1
        let parent = elements[parentIndex]
        guard areInIncreasingOrder(element, parent) else { break }
        elements[parentIndex] = element
        elements[i] = parent
        i = parentIndex
      }
    }
  }
  
  /// ヒープを用いたキューからの値の取り出し
  mutating func popFirst() -> Element? {
    guard let element = elements.popLast() else { return nil }
    guard let first = elements.first else { return element }
    
    elements.withUnsafeMutableBufferPointer { elements in
      elements[0] = element
      
      var i = 0
      while true {
        var childIndex: Int = (i << 1) + 1
        guard childIndex < elements.count else { break }
        var child: Element = elements[childIndex]
        let rightIndex = childIndex + 1
        if rightIndex < elements.count {
          let right = elements[rightIndex]
          if areInIncreasingOrder(right, child) {
            childIndex = rightIndex
            child = right
          }
        }
        if areInIncreasingOrder(element, child) { break }
        elements[childIndex] = element
        elements[i] = child
        i = childIndex
      }
    }
    
    return first
  }
}

extension PriorityQueue where Element: Comparable {
  init<S>(_ elements: S) where S: Sequence, S.Element == Element {
    self.init(elements, by: <)
  }
  
  init() {
    self.init(by: <)
  }
}

/// Dijkstra法を用いてノード0を始点とした各ノードへの重み付き有向グラフの最短経路長を取得する(計算量: O((*M + N*) log *N*))
/// - parameter nodes: `[[[Int]]]`型の隣接リスト表現
/// - returns: 各ノードへの最短経路長を格納した`[Int]`型配列
func getShortestPath(_ nodes: [[[Int]]]) -> [Int] {
  /// ノード0からインデックスに対応するノードまでの暫定累計コスト
  var distances: [Int] = [Int](repeating: Int.max, count: nodes.count)
  /// 始点ノードとして設定済かどうか
  var isDetermined: [Bool] = [Bool](repeating: false, count: nodes.count)
  /// 累計コストを優先度とする優先度付きキュー([隣接ノード, 累計コスト])
  var pQueue: PriorityQueue<(Int, Int)> = PriorityQueue<(Int, Int)>{ $0.1 < $1.1 }
  
  // 始点となるノード0を累計コスト0としてキューに追加
  pQueue.append((0, 0))
  
  while (pQueue.isEmpty == false) {
    // ノード0からの累計コストが最小となる始点ノード
    let from: Int = pQueue.popFirst()!.0
    
    // すでに始点ノードとして走査済の場合はスキップ
    if (isDetermined[from] == true) { continue }
    
    // 始点ノードまでの累計コスト(始点ノードが0の場合は0とする)
    let distanceToFrom: Int
    if (from == 0) { distanceToFrom = 0 }
    else { distanceToFrom = distances[from] }
    
    // 最短距離が確定していない隣接ノードを走査
    for node in nodes[from] where (isDetermined[node[0]] == false) {
      let (to, distance): (Int, Int) = (node[0], distanceToFrom + node[1])
      
      // 暫定累計コストの更新
      distances[to] = min(distances[to], distance)
      // 累計コストが小さい順にソートしながらキューに挿入
      pQueue.append((to, distances[to]))
    }
    
    // 始点ノードを探索済にする
    isDetermined[from] = true
  }
  
  return distances
}

/// 048 - Small Multiple
func main_048() -> Void {
  let n: Int = Int(readLine()!)!
  /// `nodes[p] = [next, distance]`はノード`p` → ノード`next` の移動コストが`distance`であることを示す
  var nodes: [[[Int]]] = [[[Int]]](repeating: [[Int]](), count: n)
  
  // 各ノードp → 隣接ノードq への最小移動コストを算出
  for p in 0 ..< n {
    for d in 0 ... 9 {
      let q: Int = (10 * p + d) % n
      
      // 自己ループは最短経路を導出する上で不要経路のため除外
      if (p == q) { continue }
      
      // 隣接ノード・移動コストがすでにnodesに含まれる場合は現在の移動コストと比較
      if let next: Int = nodes[p].indices.map({ nodes[p][$0][0] }).firstIndex(of: q) {
        if (d < nodes[p][next][1]) { nodes[p][next][1] = d }
      }
      // 含まれない場合はnodesに隣接ノード・移動コストの情報がないため追加
      else {
        nodes[p].append([q, d])
      }
    }
  }
  
  print(getShortestPath(nodes)[0])
}

/// 049 - Fibonacci Easy (mod 1000000007)
func main_049() -> Void {
  let n: Int = Int(readLine()!)!
  let division: Int = 1000000007
  
  // フィボナッチ数列の生成
  var a: [Int] = [Int](repeating: 1, count: n + 1)
  for i in 3 ... n {
    a[i] = (a[i - 2] + a[i - 1]) % division
  }
  
  print(a[n] % division)
}

/// 繰り返し二乗法を用いて累乗の剰余を取得する(計算量: O(log *exp*))
/// - parameters:
///  - base: `Int`型の底(base)
///  - exp: `Int`型の冪指数(exponent)
///  - div: `Int`型の除数(division)
/// - returns: 累乗`base`^`exp`を除数`div`で割ったときの剰余
func getModPower(_ base: Int, _ exp: Int, _ div: Int) -> Int {
  var (mod, answer): (Int, Int) = (base, 1)
  
  // 冪指数expについて、2^i(0 ≦ i ≦ 30)のフラグが立っている場合は剰余を更新
  // → 10^9 < 2^30であるため、iの最大値を30とする
  for i in 0 ..< 30 {
    if (exp & (1 << i) != 0) {
      answer *= mod
      answer %= div
    }
    mod *= mod
    mod %= div
  }
  
  return answer
}

/// 050 - Power
func main_050() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (a, b): (Int, Int) = (input[0], input[1])
  let division: Int = 1000000007
  
  print(getModPower(a, b, division))
}

/// 繰り返し二乗法を用いて組合せの総数の剰余を取得する(計算量: O(*N*))
/// - parameters:
///  - n: `Int`型の要素数
///  - r: `Int`型の組合せ数
///  - div: `Int`型の除数
/// - returns: nCrを`div`で割ったときの剰余
func getCombinationMod(_ n: Int, _ r: Int, _ div: Int) -> Int {
  var factMods: [Int] = [Int](repeating: 1, count: n + 1)
  
  // 階乗の剰余を取得(計算量: O(N))
  for i in 1 ... n {
    factMods[i] = (factMods[i - 1] * i) % div
  }
  
  // nCrの分母(=除数)となるr!,(n-r)!のモジュラ逆数を求める(計算量: O(log div))
  let (rFactMMI, nrFactMMI): (Int, Int) = (
    getModPower(factMods[r], div - 2, div) % div,
    getModPower(factMods[n - r], div - 2, div) % div
  )
  
  return factMods[n] * rFactMMI % div * nrFactMMI % div
}

/// 051 - Combination Hard
func main_051() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (x, y): (Int, Int) = (input[0], input[1])
  let division: Int = 1000000007
  
  print(getCombinationMod(x + y, x, division))
}

/// 052 - Knight
func main_052() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (x, y): (Int, Int) = (input[0], input[1])
  let division: Int = 1000000007
  
  // (+1, +2)の移動回数をa, (+2, +1)の移動回数をbとすると、
  // (3a, 3b) = (2x - y, 2y - x)で表され、共に非負整数の3の倍数とならなければならない
  var (a, b): (Int, Int) = (2 * x - y, 2 * y - x)
  if (a % 3 != 0 || b % 3 != 0 || a < 0 || b < 0) {
    print(0)
    return
  }
  (a, b) = (a / 3, b / 3)
  
  print(getCombinationMod(a + b, a, division))
}

/// 繰り返し二乗法を用いて累乗の剰余を取得する(計算量: O(log *exp*))
/// - parameters:
///  - base: `Int`型の底(base)
///  - exp: `Int`型の冪指数(exponent)
///  - div: `Int`型の除数(division)
/// - returns: 累乗`base`^`exp`を除数`div`で割ったときの剰余
func getGreaterModPower(_ base: Int, _ exp: Int, _ div: Int) -> Int {
  var (mod, answer): (Int, Int) = (base, 1)
  
  // 冪指数expについて、2^i(0 ≦ i ≦ 60)のフラグが立っている場合は剰余を更新
  // → 10^18 < 2^60であるため、iの最大値を60とする
  for i in 0 ..< 60 {
    if (exp & (1 << i) != 0) {
      answer *= mod
      answer %= div
    }
    mod *= mod
    mod %= div
  }
  
  return answer
}

/// 053 - Sum of 4^N
func main_053() -> Void {
  let n: Int = Int(readLine()!)!
  let division: Int = 1000000007
  
  // 等比数列の総和(4^{n + 1} - 1) / 3の分子をdivで割った剰余を求める
  let top: Int = getGreaterModPower(4, n + 1, division) - 1
  
  // 分母(除数)のモジュラ逆数を求める
  let bottomMMI: Int = getModPower(3, division - 2, division) % division
  
  print(top * bottomMMI % division)
}

/// 行列とその乗算を定義する構造体
struct Matrix {
  /// 行列の成分群
  var values: [[Int]] = [[1, 1], [1, 0]]
  
  init(_ matrix: [[Int]]) {
    self.values = [[Int]]()
    self.values += matrix
  }
  
  /// 行列 A, Bの積の各成分を任意の除数の剰余で表した行列を求める
  /// - parameters:
  ///  - a: `Matrix`型の行列
  ///  - b: `Matrix`型の行列
  ///  - div: `Int`型の除数
  /// - returns: `Matrix`型の行列`a`と行列`b`の積
  static func getProductMod(_ a: Matrix, _ b: Matrix, _ div: Int) -> Matrix {
    // 行列aの行数と行列bの列数が一致しない場合は実行時エラーを発生させる
    let (row, col, common): (Int, Int, Int) = (
      a.values.count,
      b.values[0].count,
      a.values[0].count == b.values.count ? a.values.count : -1
    )
    var answer: [[Int]] = [[Int]](
      repeating: [Int](repeating: 0, count: row),
      count: col
    )
    
    for r in 0 ..< row {
      for com in 0 ..< common {
        for c in 0 ..< col {
          answer[r][c] += a.values[r][com] * b.values[com][c]
          answer[r][c] %= div
        }
      }
    }
    
    return Matrix(answer)
  }
  
  /// 行列の累乗について各要素を任意の除数の剰余で表した行列を求める(計算量: O(log *exp*))
  /// - parameters:
  ///  - base: 底となる`Matrix`型の行列
  ///  - exp: `Int`型の冪指数
  ///  - div: `Int`型の除数
  /// - returns: `base`^`exp`
  static func getPowerMod(_ base: Matrix, _ exp: Int, _ div: Int) -> Matrix {
    var (baseMatrix, answer): (Matrix, Matrix) = (base, Matrix([[1, 1], [1, 0]]))
    
    // 冪指数expについて、2^i(0 ≦ i ≦ 60)のフラグが立っている場合は剰余を更新
    // → 10^18 < 2^60であるため、iの最大値を60とする
    var isFirst: Bool = true
    for i in 0 ... 60 {
      if (exp & (1 << i) != 0) {
        // 右から数えて初めてフラグが立っている部分では、その直前の走査でbaseMatrixがすでに解を持っているため流用
        if (isFirst == true) {
          answer = baseMatrix
          isFirst = false
        }
        else {
          answer = Matrix.getProductMod(answer, baseMatrix, div)
        }
      }
      baseMatrix = Matrix.getProductMod(baseMatrix, baseMatrix, div)
    }
    
    return answer
  }
}

/// 054 - Fibonacci Hard (mod 1000000000)
func main_054() -> Void {
  let n: Int = Int(readLine()!)!
  let division: Int = 1000000000
  
  let answer: Matrix = Matrix.getPowerMod(Matrix([[1, 1], [1, 0]]), n - 1, division)
  print((answer.values[1][0] + answer.values[1][1]) % division)
}

/// 055 - Recurrence Formula 1
func main_055() -> Void {
  let n: Int = Int(readLine()!)!
  let division: Int = 1000000007
  
  let answer: Matrix = Matrix.getPowerMod(Matrix([[2, 1], [1, 0]]), n - 1, division)
  print((answer.values[1][0] + answer.values[1][1]) % division)
}

/// 056 - Recurrence Formula 2
func main_056() -> Void {
  let n: Int = Int(readLine()!)!
  let division: Int = 1000000007
  
  let answer: Matrix = Matrix.getPowerMod(Matrix([[1, 1, 1], [1, 0, 0], [0, 1, 0]]), n - 1, division)
  print(
    (2 * answer.values[2][0] + answer.values[2][1] + answer.values[2][2]) % division
  )
}

/// 057 - Domino Tiling
func main_057() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (k, n): (Int, Int) = (input[0], input[1])
  let division: Int = 1000000007
  
  let matrix: Matrix
  switch k {
    case 2: matrix = Matrix([
      //  1  2  3
      [0, 0, 0, 1],
      [0, 0, 1, 0],
      [0, 1, 0, 0],
      [1, 0, 0, 1]
    ])
    case 3: matrix = Matrix([
      //  1  2  3  4  5  6  7
      [0, 0, 0, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0, 0, 0],
      [0, 0, 0, 1, 0, 0, 0, 1],
      [0, 0, 1, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0, 0, 1],
      [1, 0, 0, 0, 1, 0, 1, 0]
    ])
    case 4: matrix = Matrix([
      //  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
      [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
      [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
      [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
      [1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1]
    ])
    default: return
  }
  
  let answer: Matrix = Matrix.getPowerMod(matrix, n, division)
  print(answer.values[answer.values.count - 1][answer.values[0].count - 1])
}

/// 058 - Move on Squares 1
func main_058() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, x, y): (Int, Int, Int) = (input[0], input[1], input[2])
  
  // |x| + |y| ≦ n かつ |x + y|とnの偶奇が一致 であれば到達可能
  let isReachable: Bool = (abs(x) + abs(y) <= n && abs(x + y) % 2 == n % 2) ? true : false
  if (isReachable == true) { print("Yes") }
  else { print("No") }
}

/// Power of Two
func main_059() -> Void {
  let n: Int = Int(readLine()!)!
  
  // 2^nの一の位は{2, 4, 8, 6}の繰り返し
  let units: [Int] = [6, 2, 4, 8]
  print(units[n % 4])
}

/// 060 - Stones Game 1
func main_060() -> Void {
  let n: Int = Int(readLine()!)!
  
  if (n % 4 == 0) {
    print("Second")
  }
  else {
    print("First")
  }
}

/// 061 - Stones Game 2
func main_061() -> Void {
  let n: Int = Int(readLine()!)!
  
  var (isWinnable, isOne): (Bool, Bool) = (false, true)
  // 2進数表記の2^m - 1は2^0の位から順に1が続き、2^kの位で0が初めて登場すると以降の位で1は登場しない
  for i in 0 ... 60 {
    if (isOne == true) {
      if (n & (1 << i) == 0) { isOne = false }
    }
    else {
      // 2^m - 1でない場合は先手必勝
      if (n & (1 << i) != 0) { isWinnable = true }
    }
  }
  
  if (isWinnable == true) {
    print("First")
  }
  else {
    print("Second")
  }
}

/// 062 - Teleporter
func main_062() -> Void {
  let k: Int = readLine()!.split(separator: " ").map{ Int($0)! }[1]
  var destinations: [Int] = [0]
  readLine()!.split(separator: " ").forEach{ destinations.append(Int($0)!) }
  
  var visitedCount: [Int] = [Int](repeating: 0, count: destinations.count)
  
  // 最初の街1は到達済とする
  visitedCount[1] += 1
  var (notLoops, loops): ([Int], [Int]) = ([1] ,[Int]())
  
  var next: Int = destinations[1]
  while (true) {
    // 未到達の街であればループ対象外
    if (visitedCount[next] == 0) {
      notLoops.append(next)
      visitedCount[next] += 1
      next = destinations[next]
    }
    // 到達済の街であればループ対象なので、初回のループにのみ注目
    else if (visitedCount[next] == 1) {
      loops.append(next)
      visitedCount[next] += 1
      next = destinations[next]
    }
    // 2回目以降のループは注目しない
    else { break }
  }
  
  // kがループになる前の移動回数である場合
  if (k < notLoops.count) {
    print(notLoops[k])
  }
  // kがループになった後の移動回数である場合
  else {
    print(loops[(k - notLoops.count) % loops.count])
  }
}

/// 063 - Move on Squares 2
func main_063() -> Void {
  let n: Int = Int(readLine()!)!
  
  if (n % 2 == 0) {
    print("Yes")
  }
  else {
    print("No")
  }
}

/// 064 - All Zero
func main_064() -> Void {
  var sum: Int = 0
  let k: Int = readLine()!.split(separator: " ").map{ Int($0)! }[1]
  readLine()!.split(separator: " ").forEach{ sum += Int($0)! }
  
  if (k < sum || k % 2 != sum % 2) {
    print("No")
  }
  else {
    print("Yes")
  }
}

/// 065 - Bishop
func main_065() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (h, w): (Int, Int) = (input[0], input[1])
  
  // 1行または1列の場合は移動不可
  if (h == 1 || w == 1) { print(1) }
  // 2x2以上のマスの場合
  else {
    // 奇数x奇数である場合は移動可能マス数は「マス数の半分 + 1」
    if (h % 2 == 1 && w % 2 == 1) { print((h * w / 2) + 1) }
    // 縦横のどちらかが偶数の場合は移動可能マス数は「マス数の半分」
    else { print(h * w / 2) }
  }
}

/// 066 - Three Cards
func main_066() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, k): (Int, Int) = (input[0], input[1])
  
  /// 各カードの差の絶対値がk未満(= k - 1以下)になる余事象の場合の数
  var complement: Int = 0
  /*
   黒・白・灰のカードの値をそれぞれa, b, cとすると、aの値を固定することで
   b, cの走査範囲を以下のように絞ることができる
   max(1, a - (k - 1)) ≦ b, c ≦ min(a + (k - 1), n)
   */
  for a in 1 ... n {
    for b in max(1, a - (k - 1)) ... min(a + (k - 1), n) {
      for c in max(1, a - (k - 1)) ... min(a + (k - 1), n) {
        if (abs(b - c) <= k - 1) { complement += 1 }
      }
    }
  }
  
  print(n * n * n - complement)
}

/// 067 - Cross Sum（★2）
func main_067() -> Void {
  let hw: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (h, w): (Int, Int) = (hw[0], hw[1])
  var nums: [[Int]] = [[Int]](repeating: [0], count: h + 1)
  for i in 1 ... h {
    nums[i] += readLine()!.split(separator: " ").map{ Int($0)! }
  }
  
  var (rowSums, colSums): ([Int], [Int]) = ([0], [0])
  var (colSum, answer): (Int, Int)
  var output: String = ""
  // 行単位での合計
  for r in 1 ... h {
    // 行合計
    rowSums.append(nums[r].reduce(0, +))
  }
  
  // 列単位での合計
  for c in 1 ... w {
    colSum = 0
    nums[1 ... h].forEach{ colSum += $0[c] }
    colSums.append(colSum)
  }
  
  for r in 1 ... h {
    output = ""
    for c in 1 ... w {
      answer = rowSums[r] + colSums[c] - nums[r][c]
      output += (c != w) ? "\(answer) " : "\(answer)"
    }
    print(output)
  }
}

/// 068 - Number of Multiples 2
func main_068() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, k): (Int, Int) = (input[0], input[1])
  let nums: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  
  /// `nums`から選択した値の数
  var choices: Int
  /// `nums`から選択した値の最小公倍数
  var lcm: Int
  
  var answer: Int = 0
  // ビット全探索(計算量: O(2^k))
  // k桁の2進数を用意し、フラグが立っていれば桁に対応するnumsの値を選択することとする
  for bit in 1 ..< (1 << k) {
    (choices, lcm) = (0, 1)
    
    for index in 0 ..< k {
      if (bit & (1 << index) != 0) {
        choices += 1
        lcm = getLcm(lcm, nums[index])
      }
    }
    
    // numsから値を奇数個選択した場合は要素数を加算
    if (choices % 2 == 1) {
      answer += n / lcm
    }
    // numsから値を偶数個選択した場合は要素数を減算
    else {
      answer -= n / lcm
    }
  }
  
  print(answer)
}

/// 069 - Product Max
func main_069() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (a, b, c, d): (Int, Int, Int, Int) = (input[0], input[1], input[2], input[3])
  
  print(max(a * c, a * d, b * c, b * d))
}

/// 070 - Axis-Parallel Rectangle
func main_070() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, k): (Int, Int) = (input[0], input[1])
  
  var nodes: [[Int]] = [[Int]]()
  for _ in 0 ..< n {
    nodes.append(readLine()!.split(separator: " ").map{ Int($0)! })
  }
  
  var (count, area, answer): (Int, Int, Int) = (0, Int.max, Int.max)
  var (left, right, top, bottom): (Int, Int, Int, Int)
  for l in 0 ..< n {
    for r in 0 ..< n {
      for t in 0 ..< n {
        for b in 0 ..< n {
          // x, y座標の範囲を決める
          (left, right, top, bottom) = (
            nodes[l][0],
            nodes[r][0],
            nodes[t][1],
            nodes[b][1]
          )
          
          // 範囲内に含まれる点の個数をカウント
          count = 0
          for i in 0 ..< n {
            if (left <= nodes[i][0] && nodes[i][0] <= right &&
                nodes[i][1] <= top && bottom <= nodes[i][1]) {
              count += 1
            }
          }
          
          // 最小面積を求める
          if (count >= k) {
            area = (right - left) * (top - bottom)
            if (area < answer) { answer = area }
          }
        }
      }
    }
  }
  
  print(answer)
}

/// 071 - Linear Programming
func main_071() -> Void {
  let n: Int = Int(readLine()!)!
  var nums: [[Double]] = [[Double]]()
  for _ in 0 ..< n {
    nums.append(readLine()!.split(separator: " ").map{ Double($0)! })
  }
  
  var intersection: [Double] = [0, 0]
  var answer: Double = Double(Int.min)
  var isSatisfied: Bool
  for i in 0 ..< n - 1 {
    for j in i + 1 ..< n {
      // 2直線の交点の座標を求める
      (intersection[0], intersection[1]) = (
        (nums[i][2] * nums[j][1] - nums[j][2] * nums[i][1]) / (nums[i][0] * nums[j][1] - nums[j][0] * nums[i][1]),
        (nums[i][2] * nums[j][0] - nums[j][2] * nums[i][0]) / (nums[i][1] * nums[j][0] - nums[j][1] * nums[i][0])
      )
      
      // 全ての条件式を満たすかどうかを判定
      isSatisfied = true
      for k in 0 ..< n where (nums[k][0] * intersection[0] + nums[k][1] * intersection[1] > nums[k][2]) {
        isSatisfied = false
        break
      }
      
      if (isSatisfied == true) { answer = max(answer, intersection[0] + intersection[1]) }
    }
  }
  
  print(answer)
}

/// 072 - Max GCD 2
func main_072() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (a, b): (Int, Int) = (input[0], input[1])
  
  // iの倍数がa以上b未満の間に2つ以上存在する場合は最大公約数をiにセット
  var answer: Int = 1
  for gcd in 1 ..< b where (b / gcd > (a + gcd - 1) / gcd) {
    answer = max(answer, gcd)
  }
  
  print(answer)
}

/// 073 - Sum of Maximum
func main_073() -> Void {
  let n: Int = Int(readLine()!)!
  let nums: [Int] = [0] + readLine()!.split(separator: " ").map{ Int($0)! }
  let div: Int = 1000000007
  
  var answer: Int = 0
  // index枚選ぶ場合の数(2^index)
  var cases: [Int] = [1]
  for i in 1 ... n {
    cases.append(2 * cases[i - 1] % div)
  }
  
  for i in 1 ... n {
    // nums[i]が最大値となる場合の数 = nums[i]を選ぶ場合の数(=1) × (i - 1)枚の中から選ぶ場合の数
    answer += nums[i] * cases[i - 1]
    answer %= div
  }
  
  print(answer)
}

/// 074 - Sum of difference Easy
func main_074() -> Void {
  let n: Int = Int(readLine()!)!
  let nums: [Int] = [0] + readLine()!.split(separator: " ").map{ Int($0)! }
  
  var (initial, answer): (Int, Int) = (1 - n, 0)
  // nums[i]が加算される回数は初項initial, 公差2の等差数列の一般項initial + (i - 1) * 2で表される
  for i in 1 ... n {
    answer += nums[i] * (initial + (i - 1) * 2)
  }
  
  print(answer)
}

/// 繰り返し二乗法を用いて組合せの総数の剰余を取得する(計算量: O(log *div*))
/// - parameters:
///  - n: `Int`型の要素数
///  - r: `Int`型の組合せ数
///  - factMods: `n!`を`div`で割った剰余を格納する`[Int]`型配列
///  - div: `Int`型の除数
/// - returns: nCrを`div`で割ったときの剰余
func getCombinationMod2(_ n: Int, _ r: Int, _ factMods: [Int], _ div: Int) -> Int {
  // nCrの分母(=除数)となるr!,(n-r)!のモジュラ逆数を求める(計算量: O(log div))
  let (rFactMMI, nrFactMMI): (Int, Int) = (
    getModPower(factMods[r], div - 2, div) % div,
    getModPower(factMods[n - r], div - 2, div) % div
  )
  
  return factMods[n] * rFactMMI % div * nrFactMMI % div
}

/// 075 - Pyramid
func main_075() -> Void {
  let n: Int = Int(readLine()!)!
  let nums: [Int] = [0] + readLine()!.split(separator: " ").map{ Int($0)! }
  let div: Int = 1000000007
  
  // 階乗の剰余を取得(計算量: O(N))
  var factMods: [Int] = [Int](repeating: 1, count: n + 1)
  for i in 1 ... n {
    factMods[i] = (factMods[i - 1] * i) % div
  }
  
  // nums[i]は最上段に到達するまでの移動回数n - 1のうち、右方向に(i - 1)回進む場合の数だけ加算
  var answer: Int = 0
  for i in 1 ... n {
    answer += nums[i] * getCombinationMod2(n - 1, i - 1, factMods, div)
    answer %= div
  }
  
  print(answer)
}

/// 076 - Sum of difference
func main_076() -> Void {
  let n: Int = Int(readLine()!)!
  let nums: [Int] = [0] + readLine()!.split(separator: " ").map{ Int($0)! }.sorted()
  
  var (initial, answer): (Int, Int) = (1 - n, 0)
  // nums[i]が加算される回数は初項initial, 公差2の等差数列の一般項initial + (i - 1) * 2で表される
  for i in 1 ... n {
    answer += nums[i] * (initial + (i - 1) * 2)
  }
  
  print(answer)
}

/// 077 - Distance Sum
func main_077() -> Void {
  let n: Int = Int(readLine()!)!
  var (x, y): ([Int], [Int]) = ([Int.min], [Int.min])
  var input: [Int]
  for _ in 1 ... n {
    input = readLine()!.split(separator: " ").map{ Int($0)! }
    x.append(input[0])
    y.append(input[1])
  }
  x.sort()
  y.sort()
  
  var (initial, answer) = (1 - n, 0)
  // nums[i]が加算される回数は初項initial, 公差2の等差数列の一般項initial + (i - 1) * 2で表される
  for i in 1 ... n {
    answer += (x[i] + y[i]) * (initial + (i - 1) * 2)
  }
  
  print(answer)
}

/// 078 - Difference Optimization 1
func main_078() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, m): (Int, Int) = (input[0], input[1])
  
  // 関係性をグラフとして処理
  var nodes: [[Int]] = [[Int]](repeating: [], count: n + 1)
  nodes[0] = [0]
  var members: [Int]
  for _ in 1 ... m {
    members = readLine()!.split(separator: " ").map{ Int($0)! }
    nodes[members[0]].append(members[1])
    nodes[members[1]].append(members[0])
  }
  
  // 幅優先探索を用いた最短経路長の取得
  let answers: [Int] = getShortestDistances(nodes, 1)
  
  // 最短経路長が120を超える または 存在しない(=-1)場合は最大値(=120)を出力
  for answer in answers[1 ... n] {
    if (answer > 120 || answer == -1) {
      print(120)
    }
    else {
      print(answer)
    }
  }
}

/// 079 - ModSum
func main_079() -> Void {
  let n: Int = Int(readLine()!)!
  
  /*
   N = {1, 2, 3, .. , n - 1, n}に対して、
   P = {1 + 1, 2 + 1, 3 + 1, .. , (n - 1) + 1, 1}と設定した場合が
   N_{i} % P_{i} の総和が最大となる
   → 1からn-1までの総和が最大値
   */
  print((n - 1) * n / 2)
}

/// 080 - Difference Optimization 2
func main_080() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int(String($0))! }
  let (n, m): (Int, Int) = (input[0], input[1])
  
  // nodes[i]はノードiの[隣接ノード - 1, 移動コスト]を格納する
  var nodes: [[[Int]]] = [[[Int]]](repeating: [[Int]](), count: n)
  for _ in 0 ..< m {
    let input2: [Int] = readLine()!.split(separator: " ").map{ Int(String($0))! }
    nodes[input2[0] - 1].append([input2[1] - 1, input2[2]])
    nodes[input2[1] - 1].append([input2[0] - 1, input2[2]])
  }
  
  // Dijkstra法を用いて各ノードへの最短経路を取得する
  let shortestPath: [Int] = getShortestPath(nodes)
  
  // ノードnまでの最短経路長が初期値Int.maxである場合は-1を出力
  if (shortestPath[n - 1] == Int.max) {
    print(-1)
  }
  else {
    print(shortestPath[n - 1])
  }
}

/// 081 - Bill Changing Problem
func main_081() -> Void {
  let n: Int = Int(readLine()!)!
  print((n / 10000) + (n % 10000 / 5000) + (n % 10000 % 5000 / 1000))
}

/// 082 - Interval Scheduling Problem
func main_082() -> Void {
  let n: Int = Int(readLine()!)!
  
  var intervals: [[Int]] = [[Int.min, Int.min]]
  for _ in 1 ... n {
    intervals.append(readLine()!.split(separator: " ").map { Int(String($0))! })
  }
  // 終了時刻の早い順にソート
  intervals.sort{ $0[1] < $1[1] }
  
  var (now, answer): (Int, Int) = (intervals[0][1], 0)
  for i in 1 ... n {
    // 現在の時刻が次の映画の開始時刻より早ければ次の映画は見れる
    if (now <= intervals[i][0]) {
      // 現在時刻を次の映画の終了時刻にセット
      now = intervals[i][1]
      answer += 1
    }
  }
  
  print(answer)
}

/// 083 - We Used to Sing a Song Together（★3）
func main_083() -> Void {
  let n: Int = Int(readLine()!)!
  let kids: [Int] = readLine()!.split(separator: " ").map{ Int(String($0))! }.sorted()
  let schools: [Int] = readLine()!.split(separator: " ").map{ Int(String($0))! }.sorted()
  
  var answer: Int = 0
  for i in 0 ..< n {
    answer += abs(kids[i] - schools[i])
  }
  print(answer)
}

/// 084 - Sqrt Inequality
func main_084() -> Void {
  /*
   √a + √b < √c の両辺を2乗すると
   2√(ab) < c - a - b
   c - a - b ≦ 0 ⇔ 条件を満たさない
   c - a - b > 0 のとき、両辺を2乗すると
   4 * ab < (c - a - b)^2
   */
  let input: [Int] = readLine()!.split(separator: " ").map { Int($0)! }
  let (a, b, c): (Int, Int, Int) = (input[0], input[1], input[2])
  let d: Int = c - a - b
  
  if (d > 0 && 4 * a * b < d * d) { print("Yes") }
  else { print("No") }
}

/// 085 - Two Conditions
func main_085() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, x, y): (Int, Int, Int) = (input[0], input[1], input[2])
  for a in 1 ... n {
    for b in a ... n {
      for c in b ... n {
        for d in c ... n {
          let sum: Int = a + b + c + d
          let product: Int = a * b * c * d
          if (sum == x && product == y) {
            print("Yes")
            return
          }
        }
      }
    }
  }
  print("No")
}

/// 086 - Parentheses Check
func main_086() -> Void {
  let n: Int = Int(readLine()!)!
  let s: [String] = Array(readLine()!).map{ String($0) }

  var lp: Int = 0
  for i in 0 ..< n {
    if (lp >= 0) {
      lp += (s[i] == "(") ? 1 : -1
    }
    else {
      print("No")
      return
    }
  }
  
  print(lp == 0 ? "Yes" : "No")
}

/// 087 - Simple Math Easy
func main_087() -> Void {
  let n: Int = Int(readLine()!)!
  let div: Int = 1000000007
  
  let sum1ToN: Int = n * (n + 1) / 2 % div
  print(sum1ToN * sum1ToN % div)
}

/// 088 - Simple Math
func main_088() -> Void {
  /// 1から`n`までの総和を`div`で割った剰余を取得する
  /// - parameters:
  ///  - n: `Int`型の最大値
  ///  - div: `Int`型の除数
  /// - returns: `1 ... n`の総和を`div`で割った剰余
  func sumToNum(_ n: Int, _ div:Int = 998244353) -> Int {
    return n * (n + 1) / 2 % div
  }
  
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (a, b, c): (Int, Int, Int) = (input[0], input[1], input[2])
  let div: Int = 998244353
  
  let sumA: Int = sumToNum(a)
  let sumB: Int = sumToNum(b)
  let sumC: Int = sumToNum(c)
  
  print(sumA * sumB % div * sumC % div)
}

/// 089 - Log Inequality 2
func main_089() -> Void {
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (a, b, c): (Int, Int, Int) = (input[0], input[1], input[2])
  
  // a ≧ 1の制約上、a ≧ 1^bとなる
  if (c == 1) {
    print("No")
    return
  }
  
  /*
   log_2{a} < b * log_2{c}
   ⇔ a < c^b
   */
  var fact: Int = 1
  // ここで、c ≧ 2であることが保証されているが、2^60 > 10^18よりループ回数は最大でも60回程度に収まる
  for _ in 1 ... b {
    // オーバーフローを避けるためa < fact * cを変形
    if (a / c < fact) {
      print("Yes")
      return
    }
    fact *= c
  }
  print("No")
}

/// 090 - Digit Product Equation（★7）
func main_090() -> Void {
  /// 各桁の値の積を求める
  /// - parameter n: `Int`型の整数値
  /// - returns: `n`の各桁の値の積
  func productEachDigit(_ n: Int) -> Int {
    let digits: [Int] = Array(String(n)).map{ String($0) }.map{ Int($0)! }
    return digits.reduce(1, *)
  }
  
  /// `f(x)`の候補を取得
  /// - parameters:
  ///  - digit: `num`の桁数
  ///  - num: `元の数`
  func enumerateFx(_ digit: Int, _ num: Int) -> Void {
    // numが11桁に達した場合
    // → numをnumsに追加し、再帰処理を終了する
    if (digit == 11) {
      nums.insert(productEachDigit(num))
      return
    }
    
    // numが11桁に達していない場合
    // → 現在の1の位の値以上の値を1の位に追加
    let min: Int = num % 10
    for i in min ... 9 {
      enumerateFx(digit + 1, num * 10 + i)
    }
  }
  
  let input: [Int] = readLine()!.split(separator: " ").map{ Int($0)! }
  let (n, b): (Int, Int) = (input[0], input[1])
  
  /// `f(x)`として考えられる値
  var nums: Set<Int> = Set<Int>()
  
  // f(x)として考えられる値を列挙
  enumerateFx(0, 0)
  
  var answer: Int = 0
  for num in nums {
    // f(x)として考えられる値をもとにmを算出(m ≧ 1がここで保証される)
    let m: Int = num + b
    // mをもとにf(m)を算出
    let prodM: Int = productEachDigit(m)
    if (m - prodM == b && m <= n) {
      answer += 1
    }
  }
  
  print(answer)
}
