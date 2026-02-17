# 連鎖行列積DP

$N$ 個の行列 $A_0, A_1, \ldots, A_{N-1}$ のそれぞれの行および列の組を $(r_0, c_0), (r_1, c_1),\ldots, (r_{N-1},c_{N-1})$ とする．
このとき， $c_i=r_{i+1}$ のときのみ行列積を計算できるので，これを省略した数列 $\{x_n\}$ を以下のように定める．
```math
\begin{align*}
  \{x_n\}=\{r_0, r_1,\ldots, r_{N-1}, c_{N-1}\}
\end{align*}
```

$(A_0,\ldots,A_{N-1})$
の最適な計算順序を決定するためには，
$(A_0)(A_1,\ldots,A_{N-1})$
から
$(A_0,\ldots,A_{N-2})(A_{N-1})$
までの最も低コストな木を選べばよい．  
またこの時長さ3以上の行列の列が有れば，再帰的に同じ操作を繰り返せばよい．  
このとき内側では再帰される部分が同じ計算をしていることが有るので，DPによって重複した計算を除去する．

DPしない場合は以下の計算式で計算量が見積もられ，カタラン数になる．
```math
\begin{align}
  P(n)=\begin{cases}
    1 & n=1\textrm{のとき}\\
    \displaystyle\sum_{k=1}^{n-1} P(k)P(n-k) & n\ge 2\textrm{のとき}
  \end{cases}
\end{align}
```

DPのやり方は以下の通り
```math
\begin{align}
  \mathrm{dp}[i][j]=\begin{cases}
    0 & i=j\textrm{のとき}\\
    \min\left\{
      \rm{dp}[i][k]+dp[k+1][j] + x_{i-1}x_{k}x_{j}
      \mid i\le k < j
    \right\} & j<j \textrm{のとき}
  \end{cases}
\end{align}
```