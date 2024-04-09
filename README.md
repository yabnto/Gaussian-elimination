# 高斯消元
## 板题 [SDOI2006] 线性方程组（别问我为什么不是【模板】高斯消元法，这个太**了）
### 思路
首先需要引入一个定义[增广矩阵](https://baike.baidu.com/item/%E5%A2%9E%E5%B9%BF%E7%9F%A9%E9%98%B5/7254773)。

所以一个 $n$ 元线性方程组就可以抽象成一个矩阵，$a$ 为系数，$b$ 为方程的常数项:

$\begin{bmatrix}a_{11}&\cdots&a_{1n}&b_{1}\\\vdots&\ddots&\vdots\\a_{n1}&\cdots&a_{nn}&b{n}\end{bmatrix}$

我们知道对于任意一个一元一次方程，我们都需要将其解成一个类似于 $ax=b$ 的形式，然后变成 $x=\frac{b}{a}$，所以对于一个方程组，其实是一样的也需要解成一堆七七八八的东西，然后将 $a$ 变成一个[单位矩阵](https://baike.baidu.com/item/%E5%8D%95%E4%BD%8D%E7%9F%A9%E9%98%B5/8540268#:~:text=%E5%9C%A8%E7%9F%A9%E9%98%B5%E7%9A%84%E4%B9%98%E6%B3%95%E4%B8%AD%EF%BC%8C%E6%9C%89%E4%B8%80%E7%A7%8D%E7%9F%A9%E9%98%B5%E8%B5%B7%E7%9D%80%E7%89%B9%E6%AE%8A%E7%9A%84%E4%BD%9C%E7%94%A8%EF%BC%8C%E5%A6%82%E5%90%8C%E6%95%B0%E7%9A%84%E4%B9%98%E6%B3%95%E4%B8%AD%E7%9A%841%EF%BC%8C%E8%BF%99%E7%A7%8D%E7%9F%A9%E9%98%B5%E8%A2%AB%E7%A7%B0%E4%B8%BA%E5%8D%95%E4%BD%8D%E7%9F%A9%E9%98%B5%E3%80%82,%E5%AE%83%E6%98%AF%E4%B8%AA%E6%96%B9%E9%98%B5%EF%BC%8C%E4%BB%8E%E5%B7%A6%E4%B8%8A%E8%A7%92%E5%88%B0%E5%8F%B3%E4%B8%8B%E8%A7%92%E7%9A%84%E5%AF%B9%E8%A7%92%E7%BA%BF%EF%BC%88%E7%A7%B0%E4%B8%BA%E4%B8%BB%E5%AF%B9%E8%A7%92%E7%BA%BF%EF%BC%89%E4%B8%8A%E7%9A%84%E5%85%83%E7%B4%A0%E5%9D%87%E4%B8%BA1%E3%80%82%20%E9%99%A4%E6%AD%A4%E4%BB%A5%E5%A4%96%E5%85%A8%E9%83%BD%E4%B8%BA0%E3%80%82)，此时右边的 $b$ 就是对应行的[主元](https://baike.baidu.com/item/%E4%B8%BB%E5%85%83/19061118)的解了。

为了达成这个目标，我们有三种行变换操作，目的是为了在保证主元解不变的情况下进行变化：

1. 交换两行，这个很显然不会变答案
2. 将一行整体缩放一个数，这个很显然可以用等式的性质二来证明（欸，要证明吗，好像本来就满足）
3. 将一行叠加到另一行去，这个可以用等式的性质一来证明（嘿，好像不用证明）

好，此时我们就可以快速口胡一份代码了，我们可以对于每一行都将其下面的全部都变成 $0$，因为动了上面的，那么上面的也得行变换，下面的也都得行变换了，注意此时是要用行变换的！！！做完以后就发现最后一行刚好是一个 $a_{nn}=b_n$ 的东西，所以可以将 $x_n$ 求出，然后倒着带入求解其它的元即可。

到这里，我们就解决完了所有等式有唯一解的情况了，接下来就要说一些特殊的情况：

1. 等式当前行主元系数为 $0$，这个又分两种情况，一个呢是下面的当前主元系数不为 $0$，那么此时只要用行变换把它弄上来即可，至于为什么不动上面的，理由和上面一样，二个呢是下面的当前主元系数也为 $0$，那么这个主元，就变成了无论是什么，整个式子都成立，所以就成了无唯一解，可此时我们还不能下定结论，因为下面还有更离谱的事情，不过这就告诉我们，我们不能单纯的只维护当前行，应该把当前列也维护了，因为如果当前行的当前主元无唯一解时，我们就需要让当前列往后走，找到一个有唯一解的主元给当前行，那么跑到最后时，当前行和当前列不一定都跑完了所有的方程。
2. 等式接完后有剩余等式，这个时候我们会发现一件事情，就是这些等式与上面的等式是[线性相关](https://baike.baidu.com/item/%E7%BA%BF%E6%80%A7%E7%9B%B8%E5%85%B3/6416511#:~:text=1.%20%E5%AF%B9%E4%BA%8E%E4%BB%BB%E4%B8%80%E5%90%91%E9%87%8F%E7%BB%84%E8%80%8C%E8%A8%80%2C%EF%BC%8C%E4%B8%8D%E6%98%AF%E7%BA%BF%E6%80%A7%E6%97%A0%E5%85%B3%E7%9A%84%E5%B0%B1%E6%98%AF%E7%BA%BF%E6%80%A7%E7%9B%B8%E5%85%B3%E7%9A%84%E3%80%82%202.%20%E5%90%91%E9%87%8F%20%E7%BB%84%E5%8F%AA%E5%8C%85%E5%90%AB%E4%B8%80%E4%B8%AA%E5%90%91%E9%87%8Fa%E6%97%B6%EF%BC%8Ca%E4%B8%BA0%E5%90%91%E9%87%8F%EF%BC%8C%E5%88%99%E8%AF%B4A%E7%BA%BF%E6%80%A7%E7%9B%B8%E5%85%B3%3B,%E8%8B%A5a%E2%89%A00%2C%20%E5%88%99%E8%AF%B4A%E7%BA%BF%E6%80%A7%E6%97%A0%E5%85%B3%E3%80%82%203.%20%E5%8C%85%E5%90%AB%20%E9%9B%B6%E5%90%91%E9%87%8F%20%E7%9A%84%E4%BB%BB%E4%BD%95%E5%90%91%E9%87%8F%E7%BB%84%E6%98%AF%E7%BA%BF%E6%80%A7%E7%9B%B8%E5%85%B3%E7%9A%84%E3%80%82)的，所以解不解都没关系，但是，这里因为已经被上面削成了 $0$，所以有可能会出现 $0\ne 0$ 的情况，所以要判无解，然后再判无唯一解。

好，到这里就全算完了。

### code
```cpp
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const int MaxN = 110;
const double eps = 1e-4;

double a[MaxN][MaxN], ans[MaxN];
int n;

int main() {
  cin >> n;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n + 1; j++) {
      cin >> a[i][j];
    }
  }
  int i, j;
  for (i = 1, j = 1; i <= n && j <= n;) {
    int u = 0;
    for (int k = n; k >= i; k--) {
      if (fabs(a[k][j]) > eps) {
        u = k;
      }
    }
    if (!u) {
      j++;
      continue;
    } 
    for (int k = 1; k <= n + 1; k++) {
      swap(a[i][k], a[u][k]);
    }
    for (int k = i + 1; k <= n; k++) {
      double d = a[k][j] / a[i][j];
      for (int l = j; l <= n + 1; l++) {
        a[k][l] -= a[i][l] * d;
      }
    }
    i++, j++;
  }
  for (int k = i; k <= n; k++) {
    if (fabs(a[k][n + 1]) > eps) {
      cout << "-1" << endl;
      return 0;
    }
  }
  if (j < n || i < n) {
    cout << "0" << endl;
    return 0;
  }
  for (int i = n; i >= 1; i--) {
    double sum = 0;
    for (int j = 1; j <= n; j++) {
      sum += ans[j] * a[i][j];
    }
    if (fabs(a[i][i]) < eps) {
      cout << "0" << endl;
      return 0;
    }
    ans[i] = (a[i][n + 1] - sum) / a[i][i];
  }
  for (int i = 1; i <= n; i++) {
    cout << fixed << setprecision(2) << "x" << i << "=" << ans[i] << endl;
  }
  return 0;
}
```
