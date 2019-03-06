/**
 * Author: Simon Lindholm
 * Date: 2018-07-06
 * License: CC0
 * Description: Permutation -> integer conversion. (Not order preserving.)
 * Time: O(n)
 */
#pragma once

int permToInt(vi& v) {
	int use = 0, i = 0, r = 0;
	trav(x,v)r=r * ++i + __builtin_popcount(use & -(1 << x)),
		use |= 1 << x;                // (note: minus, not ~!)
	return r;
}
