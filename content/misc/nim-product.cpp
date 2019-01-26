/**
 * Author: 
 * Description: Nim Product.
 */
ll f[64][64]; // initialize to -1
ll _nim(int x, int y) {
  if (!x || !y)  return 1 << (x + y);
  if (f[x][y] != -1)  return f[x][y];
  ll ret = 1, e = 1;
  for (int i = 0; (x >> i) || (y >> i); i++)
    if (((x ^ y) >> i) & 1)  e <<= (1 << i);
    else if ((x >> i) & 1)  ret = nim(ret, 3 * (1ll << (1 << i)) / 2);
  f[x][y] = nim(ret, e);
  return f[x][y];
}
ll nim(ll x, ll y) {
  ll ret = 0;
  for (int i = 0; x >> i; i++)
    if ((x >> i) & 1)
      for (int j = 0; y >> j; j++)
        if ((y >> j) & 1)  ret ^= _nim(i, j);
  return ret;
}
