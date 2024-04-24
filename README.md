# 47-ICPC-World-Final-2023

<img src=".../assets/images/B.png" alt="Paris" class="center">
<img src=".../assets/images/B.png" alt="Paris" class="center">

## Solution

```cpp
#include <iomanip>
#include <iostream>
using namespace std;

int main() {
  int n, m;
  while (cin >> n >> m) {
    int a, tot = 0;
    for (int i = 0; i < n*m; i++) {
      cin >> a;
      tot += a;
    }
    cout << fixed << setprecision(9) << (long double)(tot) / (n*m) << endl;
  }
  return 0;
}
```
