\ Lattice triangle with point 0 at the origin and edge lengths 2 k_ij + 1.
Minimize
 obj: 0
Subject To
 dist_01: - 4 k_01 + [ x1^2 + y1^2 - 4 k_01^2 ] = 1
 dist_02: - 4 k_02 + [ x2^2 + y2^2 - 4 k_02^2 ] = 1
 dist_12: - 4 k_12 + [ x1^2 - 2 x1 * x2 + x2^2 + y1^2 - 2 y1 * y2 + y2^2 - 4 k_12^2 ] = 1
Bounds
 0 <= x1
 0 <= y1
 0 <= x2
 0 <= y2
 0 <= k_01
 0 <= k_02
 0 <= k_12
Generals
 x1
 y1
 x2
 y2
 k_01
 k_02
 k_12
End
