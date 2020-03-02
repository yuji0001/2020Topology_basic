Numpyで整数行列をSmith標準化してみた
===

最近，「タンパク質構造とトポロジー（平岡）」という本を友人と読んでいます．
その本の中でSmith標準化という整数行列の分解手法があったので，Numpyで実装してみました．

## 基本行列
Smith標準化の説明をする前に，整数行列の基本行列の説明をする必要があります．

ベクトル空間上の，基本行列は行列の基本変形を行うための行列になります
(大学1年生ぐらいで習っているはず．)．
行列の基本変形は以下の3つがあります．
1. ある行か列を定数倍（0を除く）する．
2. 2つの行か列を入れ替える．
3. ある行か列に，他のある行か列の定数倍を加える．

それぞれの変形は左右のどちらからか，**正則な**行列（つまり逆行列が存在する行列）を作用させることによって達成できます．
この正則な行列は基本行列と呼ばれています．

ここで，1,2,3それぞれに対応する基本行列を$P_i(c), Q_{ij}, R_{ij}(c)$と表記すると，以下のようになります．

$$
P_i(c) =\begin{pmatrix}
1&&&&&&\\
&\ddots&&&&&\\
&&1&&&&\\
&&&c&&&\\
&&&&1&&\\
&&&&&\ddots&\\
&&&&&&1\\
\end{pmatrix}
$$

$$
Q_{ij} =\begin{pmatrix}
1&&&&&&\\
&\ddots&&&&&\\
&&0&&1&&\\
&&&\ddots&&&\\
&&1&&0&&\\
&&&&&\ddots&\\
&&&&&&1\\
\end{pmatrix}
$$

$$
R_{ij}(c) =\begin{pmatrix}
1&&&&&&\\
&\ddots&&&&&\\
&&1&&c&&\\
&&&\ddots&&&\\
&&&&1&&\\
&&&&&\ddots&\\
&&&&&&1\\
\end{pmatrix}
$$

また，それぞれの基本行列の逆行列は
1. $P_i(c)^{-1} = P_i(1/c)$
2. $Q_{ij}^{-1} = Q_{ij}$
3. $R_{ij}(c)^{-1} = R_{ij}(-c)$

となるため，基本行列の正則性が導かれます.
掃き出し法などの線形代数の基礎理論は，この基本行列を用いて導かれています． 

さて，整数行列に対する基本行列を考えて見ましょう．

ベクトル空間上の基本行列$Q_{ij}$と$R_{ij}(c)$は整数$c$に対して，逆行列が整数行列となるため，整数行列に対する基本行列として採用できます．
一方で，$P_i(c)$は定数$c$が 1 か -1 でない限り，逆行列は整数行列にはなりません．
$P_i(1)$は単位行列となるため，$P_i(-1)$は整数行列の基本行列となります．

よって， 整数$c$に対して，$P_i(-1), Q_{ij}, R_{ij}(c)$が整数行列の基本行列となります．

## Smith標準形
Smith標準化は整数行列の対角化手法のです．
任意の整数行列Aは行列$P$と$Q$を用いて，
$$
B = Q^{-1} A P = 
\left(
    \begin{array}{c|c}
     \begin{matrix}
     c_1 & &0 \\
      & \ddots &\\
     0 & & c_k
     \end{matrix}
     & 0 \\
\hline
     0 &0
\end{array}\right)
$$
と対角化する事ができます．
ここで，$P$と$Q$は基本行列の積であり，$c_{i+1}$は$c_i$で割り切ることができます．
行列$P$と$Q$は基本行列の積であるため，整数行列の逆行列が存在します．
この基本行列の積によって対角化された形式をSmith標準形と呼びます．

例：

$$
\begin{pmatrix}
1&0&0\\
7&-7&-4\\
0&2&1
\end{pmatrix}^{-1}
\begin{pmatrix}
7&3&2&1\\
7&6&7&7\\
4&8&2&0
\end{pmatrix}
\begin{pmatrix}
0&0&-6&13\\
0&0&-13&28\\
0&1&65&-138\\
1&-2&-49&101
\end{pmatrix}
=
\begin{pmatrix}
1&0&0&0\\
0&1&0&0\\
0&0&2&0
\end{pmatrix}
$$

Smith標準形は，次の4つの操作を繰り返すことによって達成されます．

1. $A_{t+1} = {\rm moveMN}(A_t)$ ： (1,1)要素に最小の値を移動する変換．

非ゼロの絶対値が最小の要素$A[i,j]$を(1,1)成分に移動する．
つまり，
$$
    {\rm moveMN}(A_t) = 
    \begin{cases}
    Q_{i,0} A_tQ_{0,j},& A[i,j]>0\\
    Q_{i,0} A_tQ_{0,j}P_1(-1),& A[i,j]<0    
    \end{cases}
$$
と定義する．

例：
$$
\begin{pmatrix}
7&9&5&8\\
5&1&0&2\\
8&1&8&5
\end{pmatrix}
\rightarrow
\begin{pmatrix}
1&5&0&2\\
9&7&5&8\\
1&8&8&5
\end{pmatrix}
$$
   
2. $A_{t+1}= {\rm rowR}(A_t)$：$A[1,1]$ を除く左端の列を小さくする変換.

任意の$i$に対して，$q_i = A[i,0]\div A[0,0]$(余りは切り捨て)とすると，
$$
{\rm rowR}(A_t) := \prod_{i>1} R_{i,0}(-q_i) A_t
$$
と定義する．

例：
$$
\begin{pmatrix}
1&5&0&2\\
9&7&5&8\\
1&8&8&5
\end{pmatrix}
\rightarrow
\begin{pmatrix}
1&5&0&2\\
0&-38&5&10\\
0&3&8&3
\end{pmatrix}
$$

3. $A_{t+1}= {\rm colR}(A_t)$：$A[1,1]$ を除く上端の行を小さくする変換.

任意の$i$に対して，$q_i = A[0,i]\div A[0,0]$(余りは切り捨て)とすると，
$$
{\rm colR}(A_t) := \prod_{i>1} A_t R_{0,i}(-q_i) 
$$
と定義する．

例：
$$
\begin{pmatrix}
1&5&0&2\\
0&-38&5&10\\
0&3&8&3
\end{pmatrix}
\rightarrow
\begin{pmatrix}
1&0&0&0\\
0&-38&5&10\\
0&3&8&3
\end{pmatrix}
$$

1. $A_{t+1}= {\rm remR}(A_t)$：$A[2:,2:]$において，$A[1,1]$で割り切れない要素を除外する変換．

$A[i,j]$を$A[1,1]$で割り切れない要素とし，$q = A[i,j]\div A[1,1]$（余りは切り捨て）とする．
この時，
$$
{\rm remR}(A_t) :=  R_{0,j}(q) A_t R_{i,0}(-1)
$$
と定義する．

例：
$$
\begin{pmatrix}
2&0&0&0\\
0&2&3&4\\
0&2&2&6
\end{pmatrix}
\rightarrow
\begin{pmatrix}
2&0&-2&0\\
2&2&1&4\\
0&2&2&6
\end{pmatrix}
$$

以上の4つの変換を以下のアルゴリズムにしたがって繰り返すことにより，Smith標準形を生成できる．

1. A = moveMN(A)
2. A = rowR(A)
3. もし, A[2:,1] != 0ならば，1.にもどる
4. A = colR(A)
5. もし，A[1,2:] != 0ならば，1.にもどる
6. もし,A[2:,2:]のすべての要素がA[1,1]で割り切れないならば，A = remR(A)とし，1．に戻る．
7. A = A[2:,2:]とし，1.に戻る．

以下にこのpython codeを紹介する．

## Python コード
```python
import numpy as np

# 整数行列上の基本行列
def Eij(i,j,n):
    E = np.eye(n)
    E[i,i] = 0
    E[j,j] = 0
    E[i,j] = 1
    E[j,i] = 1
    return E

def Ei(i,n):
    E = np.eye(n)
    E[i,i] = -1
    return E

def Ec(i,j,c,n):
    E = np.eye(n)
    E[i,j] = c
    return E    

# A[k:,k:]上の0以外の最小の絶対値をA[k,k]に移動する変換
def moveMN(A,k):
    tmp_A = A[k:,k:]
    a = np.abs(tmp_A[tmp_A != 0]).min()
    i = np.where(np.abs(tmp_A) == a)[0][0] + k
    j = np.where(np.abs(tmp_A) == a)[1][0]+ k
    P = Eij(k,j,A.shape[1])
    invQ = Eij(i,k,A.shape[0])
    B = invQ.dot(A).dot(P)
    if B[k,k]<0:
        Pi =Ei(k,A.shape[1])
        B = B.dot(Pi)
        P = P.dot(Pi)
    return invQ.astype(int),B.astype(int),P.astype(int)

#A[k,k]を使って，A[k+1:,k]を0に整える変換．(整いきれなかったら要素はA[k,k]より小さくなる)
def rowR(A,k):
    B = A.copy()
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    for i in range(k+1,A.shape[0]):
        q = A[i,k]//A[k,k]
        #残渣
        r = A[i,k]%A[k,k]
        invQi = Ec(i,k,-q,A.shape[0])
        B = invQi.dot(B)
        invQ = invQi.dot(invQ)
    return invQ.astype(int),B.astype(int),P.astype(int)

#A[k,k]を使って，A[k,k+1]を0に整える変換．(整いきれなかったら要素はA[k,k]より小さくなる)
def colR(A,k):
    B = A.copy()
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    for i in range(k+1,A.shape[1]):
        q = A[k,i]//A[k,k]
        #残渣
        r = A[k,i]%A[k,k]
        Pi = Ec(k,i,-q,A.shape[1])
        B = B.dot(Pi)
        P = P.dot(Pi)
    return invQ.astype(int),B.astype(int),P.astype(int)

# A[k+1:,k+1:]においてA[k,k]で割り切れない要素A[i,j]をA[k,k]の残差に変換する変換
def remR(A,k):
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    #  Find i,j
    i = np.where(A[k+1:,k+1:]%A[k,k] !=0)[0][0] +k+1
    j = np.where(A[k+1:,k+1:]%A[k,k] !=0)[1][0] +k+1
    q = A[i,j]//A[k,k]
    r = A[i,j]%A[k,k]
    invQi = Ec(i,k,q,A.shape[0])
    Pi = Ec(k,j,-1,A.shape[1])
    B = invQi.dot(A).dot(Pi)
    P = P.dot(Pi)
    invQ = invQi.dot(invQ)
    return invQ.astype(int),B.astype(int),P.astype(int)


# Main Function
def Smith_Normalization(A):
    invQ = np.eye(A.shape[0])
    P = np.eye(A.shape[1])
    A0 = A.copy()
    # limit of optimization
    N = 1000
    for k in range(min(A0.shape)):
        # If A0[k:,k:] is zero matrix, then stop calculation
        if np.sum(np.abs(A0[k:,k:]))==0:
            break
        for i in range(N):
            if i == N-1 : 
                print("Error: Time Out")
            # minimize A[k,k]
            invQi,A1,Pi = moveMN(A0,k)
            invQ = invQi.dot(invQ)
            P = P.dot(Pi)
            # make zero row A[k+1:,k]
            invQi,A2,Pi = rowR(A1,k)
            invQ = invQi.dot(invQ)
            P = P.dot(Pi)
            # if row A2[k+1:,k] is zero vector ?
            if np.abs(A2[k+1:,k]).sum() ==0:
                # make zero col A[k,k+1:]
                invQi,A3,Pi = colR(A2,k)
                invQ = invQi.dot(invQ)
                P = P.dot(Pi)
                # if col A3[k+1:,k] is zero vector ?
                if np.abs(A3[k,k+1:]).sum() ==0:
                    # A[k,k]|A[k+1:,k+1:]?
                    if np.sum(A3[k+1:,k+1:]%A3[k,k]) == 0:                    
                        A0 = A3.copy()            
                        break;
                    else:
                        # reduce A[k+1:,k+1:]
                        invQi,A0,Pi = remR(A3,k)
                        invQ = invQi.dot(invQ)
                        P = P.dot(Pi)
                else:
                    A0 = A3.copy()            
            else:
                A0 = A2.copy()

    B = A0.copy().astype(int)
    P = P.astype(int)
    invQ = invQ.astype(int)
    return invQ,B,P



# Check result
A = np.random.randint(0,10,(4,5))
invQ,B,P = Smith_Normalization(A)
print(A)
print(invQ)
print(B)
print(P)

```

## まとめ
- 整数行列の基本行列を定義した．
- Smith標準形とその変換方法を紹介した．
- Smith標準形へ変換するpythonコードを作成した．

pythonコードに比べ，まとめを書くのがとても大変でした．

## Code詳細
https://github.com/yuji0001/2020Topology_basic

## Reference

平岡裕章，「タンパク質構造とトポロジー,パーシステントホモロジー群入門」，共立出版,
2013 
## Author 
Yuji Okamoto yuji.0001@gmail.com

