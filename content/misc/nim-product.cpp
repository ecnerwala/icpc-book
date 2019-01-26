/**
 * Author: Andrew He
 * Description: Nim Product.
 */
using ull = uint64_t;
ull _nimProd2[64][64];
ull nimProd2(int i, int j) {
  if (_nimProd2[i][j]) return _nimProd2[i][j];
  if ((i & j) == 0) return _nimProd2[i][j] = 1ull << (i|j);
  int a = (i&j) & -(i&j);
  return _nimProd2[i][j] = nimProd2(i ^ a, j) ^ nimProd2((i ^ a) | (a-1), (j ^ a) | (i & (a-1)));
}
ull nimProd(ull x, ull y) {
  ull res = 0;
  for (int i = 0; x >> i; i++)
    if ((x >> i) & 1)
      for (int j = 0; y >> j; j++)
        if ((y >> j) & 1)
          res ^= nimProd2(i, j);
  return res;
}
