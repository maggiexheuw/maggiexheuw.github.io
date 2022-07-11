---
layout:     post
title:      Kondo model
subtitle:   magnetic
date:       2022-07-10
author:     Maggie
header-img: img/desktop.jpg
catalog: true
tags:
    - Kondo model
---




## Kondo model

---

​        现在我们讲述第一部分，如何从Anderson 模型到s-d模型。

​        根据安德森所给的哈密顿量                                
$$
\hat{H}=\sum_{k,d}\varepsilon_{k}\hat{c}^{\dagger}_kc_k+\sum_{\sigma}\varepsilon_{d}\hat{c}_{d\sigma}^{\dagger}\hat{c}_{d\sigma}+\sum_{k,d}V_{k}\left(\hat{c}_{k\sigma}^{\dagger}\hat{c}_{d\sigma}+h.c\right)
+U\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}
$$
​       安德森模型哈密顿量的希尔伯特空间可以分为三个部分，零占据、单占据、双占据，因此我们定义三套投影算符算符。                                                                      
$$
\begin{align*}
\begin{cases}
\hat{P}_0=(1-\hat{n}_{d\uparrow})(1-\hat{n}_{d\downarrow})\\
\hat{P}_1=\hat{n}_{d\downarrow}(1-\hat{n}_{d\uparrow})+\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
\hat{P}_2=\hat{n}_{d\downarrow}\hat{n}_{d\uparrow}
\end{cases}
\end{align*}
$$
​       作为投影算符，这三个算符满足于投影算符的一切性质（细节留在后面）                                                                          
$$
\sum\limits_{n=0}^{2}\hat{P}_n=\mathbf{I}\qquad \hat{P}_{i}\hat{P}_j=\delta_{ij}\hat{P}_i
$$
​       定义$\hat{P}=(\hat{P}_0,\hat{P}_1,\hat{P}_2)^{\mathbf{T}}$,我们对整个薛定谔方程进行投影                                                                         
$$
\hat{P}\hat{H}\mid \psi\rangle=\hat{P}E\mid\psi\rangle
$$
​      本征态分成三部分$\mid \psi\rangle =(\mid \psi_1\rangle,\mid \psi_2\rangle,\mid \psi_3\rangle )^{\mathbf{T}}$,那么我们再次消去$\mid \psi \rangle $前面的$\hat{P}$,可以得到                                                            
$$
\hat{P}\hat{H}\hat{P}^{\dagger}\mid\psi\rangle=E\mid \psi\rangle
$$
​       我们将上面这个本征方程展开成矩阵形式                                                                    
$$
\begin{bmatrix}
\hat{P}_0\hat{H}\hat{P}_0&\hat{P}_0\hat{H}\hat{P}_1&\hat{P}_0\hat{H}\hat{P}_2\\
\hat{P}_1\hat{H}\hat{P}_0&\hat{P}_1\hat{H}\hat{P}_1&\hat{P}_1\hat{H}\hat{P}_2\\
\hat{P}_2\hat{H}\hat{P}_0&\hat{P}_2\hat{H}\hat{P}_1&\hat{P}_2\hat{H}\hat{P}_2
\end{bmatrix}
\begin{bmatrix}
\mid \psi_0\rangle\\
\mid \psi_1\rangle \\
\mid \psi_2\rangle
\end{bmatrix}=
E
\begin{bmatrix}
\mid \psi_0\rangle\\
\mid \psi_1\rangle \\
\mid \psi_2\rangle
\end{bmatrix}
$$
​     为了书写方便我们定义$\hat{P}_m\hat{H}\hat{P}_n=\hat{H}_{mn}$,接下来我们考虑具体计算矩阵元$\hat{H}_{mn}$,先来计算较为复杂的非对角项目

​      先来看一下$\hat{H}_{01}$，我们首先可以分析出,此时只有对hopping项才有作用，hopping项中$c_{k},c_{k}^{\dagger}$与$\hat{P}_0$和$\hat{P}_1$是相互对易可以互换的，因此我们需要先来研究$\hat{P}_0$和$\hat{P}_1$,对$c_{d\sigma}$ 的作用，我们来看                                                              
$$
\begin{align*}
\hat{P}_0c_{d\uparrow}\hat{P}_1&=(1-\hat{n}_{d\uparrow})(1-\hat{n}_{d\downarrow})c_{d\uparrow}[\hat{n}_{d\downarrow}(1-\hat{n}_{d\uparrow})+\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})]\\
&=(1-\hat{n}_{d\uparrow})(1-\hat{n}_{d\downarrow})c_{d\uparrow}\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
&=(1-\hat{n}_{d\uparrow})c_{d\uparrow}\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
&=c_{d\uparrow}((1-\hat{n}_{d\uparrow})+1)\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
&=c_{d\uparrow}\hat{n}_{d\uparrow}(1-\hat{n}_{d \downarrow})
\end{align*}
$$

$$
\begin{align*}
\hat{P}_0c_{d\uparrow}^{\dagger}\hat{P}_1&=(1-\hat{n}_{d\uparrow})(1-\hat{n}_{d\downarrow})c_{d\uparrow}^{\dagger}[\hat{n}_{d\downarrow}(1-\hat{n}_{d\uparrow})+\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})]\\
&=(1-\hat{n}_{d\uparrow})(1-\hat{n}_{d\downarrow})c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
&=(1-\hat{n}_{d\uparrow})c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
&=-c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})\\
&=-c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}(1-\hat{n}_{d \downarrow})
\end{align*}
$$

​       我们来分析一下$c_{d\uparrow}\hat{n}_{d\uparrow}(1-\hat{n}_{d \downarrow})$,自旋向上的局域电子占据数为$1$的时候，$c_{d\uparrow}$作用之后才不为0，于是这一项可以写成$c_{d\uparrow}\hat{n}_{d\uparrow}(1-\hat{n}_{d \downarrow})$,另一方面$c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}(1-\hat{n}_{d \downarrow})$这一项为0，因为$c^{\dagger}_{d\uparrow}$作用之后不为0，说明这里本来没有局域电子占据， $c^{\dagger}_{d\uparrow}$作用产生了一个局域电子，然而$\hat{n}_{d\uparrow}$在这之前作用，这说明这一项为0。那么                                                                
$$
\hat{H}_{01}=\sum_{k,\sigma}V_{k}\hat{c}_{k\sigma}^{\dagger}\hat{c}_{d\sigma}(1-\hat{n}_{d,-\sigma})
$$
​       现在我们计算$\hat{H}_{12}$,同样我们先来考察对$\hat{P}_1$,对$c_{d\sigma}$ 的作用                                              
$$
\begin{align*}
\hat{P}_0\hat{c}_{d\uparrow}\hat{P}_2&=[\hat{n}_{d\downarrow}(1-\hat{n}_{d\uparrow})+\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})]c_{d\uparrow}\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}\\
&=(1-\hat{n}_{d\uparrow})c_{d\uparrow}\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}\\
&=\hat{n}_{d\downarrow}\hat{c}_{d\uparrow}\hat{n}_{d\uparrow}\\
&=\hat{n}_{d\downarrow}(1+\hat{n}_{d\uparrow})\hat{c}_{d\uparrow}
\end{align*}
$$

$$
\begin{align*}
\hat{P}_0\hat{c}_{d\uparrow}^{\dagger}\hat{P}_2&=[\hat{n}_{d\downarrow}(1-\hat{n}_{d\uparrow})+\hat{n}_{d\uparrow}(1-\hat{n}_{d\downarrow})]c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}\\
&=(1-\hat{n}_{d\uparrow})c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}\\
&=-c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}\\
\end{align*}
$$

​       同理，我们知道$-c_{d\uparrow}^{\dagger}\hat{n}_{d\uparrow}\hat{n}_{d\downarrow}$为0，那么$\hat{H}_{12}$为                                                                        
$$
\hat{H}_{12}=\sum_{k,\sigma}V_{k}\hat{c}_{k\sigma}^{\dagger}\hat{c}_{d\sigma}\hat{n}_{d,-\sigma}
$$
​       $\hat{H}_{02}$ 和$\hat{H}_{20}$为0，因为hopping项涉及的是单粒子的hopping。非对角项计算完成之后，我们来计算对角项，我们这里直接给出结果                                                                         

$$
\begin{align*}
\displaystyle
\begin{cases}
\hat{H}_{00}=\sum\limits_{k,\sigma}\varepsilon_k \hat{c}^{\dagger}_{k,\sigma}c_{k,\sigma}\\
\hat{H}_{11}=\sum\limits_{k,\sigma}\varepsilon_k \hat{c}^{\dagger}_{k,\sigma}c_{k,\sigma}+\varepsilon_d\\
\hat{H}_{22}=\sum\limits_{k,\sigma}\varepsilon_k \hat{c}^{\dagger}_{k,\sigma}c_{k,\sigma}+2\varepsilon_d+U
\end{cases}
\end{align*}
$$
​          因为我们只考虑单粒子情况，我们要研究$\mid \psi_1\rangle$                                                                            
$$
\begin{align*}
\begin{cases}
\hat{H}_{00}\mid \psi_0\rangle+\hat{H}_{01}\mid \psi_1\rangle=E\mid \psi_0\rangle
&\implies \mid \psi_0\rangle=\frac{\hat{H}_{01}\mid \psi_1\rangle}{E-\hat{H}_{00}}
\\
\hat{H}_{21}\mid \psi_1\rangle+\hat{H}_{22}\mid \psi_2\rangle=E\mid \psi_2\rangle
&\implies \mid \psi_1\rangle=\frac{\hat{H}_{21}\mid \psi_1\rangle}{E-\hat{H}_{22}}
\\
\end{cases}
\end{align*}
$$
​        上面是因为$\hat{H}_{02}=\hat{H}_{20}=0$,因此少了一项，我们现在考虑把它代到另一个方程中                                             
$$
\begin{align*}
(\hat{H}_{10}\mid \psi_0\rangle+\hat{H}_{11}\mid\psi_1\rangle+\hat{H}_{12}\mid\psi_2\rangle)=E\mid\psi_1\rangle
\end{align*}
$$
​       我们先来考察                                                   
$$
\begin{align*}
&\frac{\hat{H}_{12}\hat{H}_{21}}{E-\hat{H}_{22}}=\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{n}_{k-\sigma}\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma}}^{\dagger}\hat{n}_{k^{\prime}-\sigma^{\prime}}\hat{c}_{k^{\prime}\sigma}}{E-\hat{H}_{22}}\\
&
\frac{\hat{H}_{10}\hat{H}_{01}}{E-\hat{H}_{00}}=\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}(1-\hat{n}_{k\sigma})\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma}}^{\dagger}(1-\hat{n}_{k^{\prime}\sigma^{\prime}})\hat{c}_{k^{\prime}\sigma}}{E-\hat{H}_{00}}
\end{align*}
$$
​       现在我们来看                                                           
$$
\begin{align*}
&\frac{1}{E-\hat{H}_{22}}=\frac{1}{E-(\hat{H}_{00}-\varepsilon_k+2\varepsilon_{d}+U)}
=-\frac{1}{\varepsilon_k-\varepsilon_{d}+U}\frac{1}{1-\frac{E-\varepsilon_d-\hat{H}_{00}}{\varepsilon_k-\varepsilon_{d}+U}}\\
&\frac{1}{E-\hat{H}_{00}}=-\frac{1}{\varepsilon_k-\varepsilon_d}\frac{1}{1-\frac{E-\varepsilon_d-\hat{H}_{00}}{\varepsilon_k-\varepsilon_{d}}}
\end{align*}
$$
​       我们取一阶近似                                          
$$
\frac{\hat{H}_{12}\hat{H}_{21}}{E-\hat{H}_{22}}+\frac{\hat{H}_{10}\hat{H}_{01}}{E-\hat{H}_{00}}=-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{n}_{k-\sigma}\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{n}_{k^{\prime}-\sigma^{\prime}}\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d+U}-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}(1-\hat{n}_{k\sigma})\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}(1-\hat{n}_{k^{\prime}\sigma^{\prime}})\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d}
$$
​        这里我们知道占据数算符$\hat{n}$对$\mid\psi_1\rangle$作用的结果只能是$0,1$因此我们可以考虑把占据数算符拿掉                                                                  
$$
\begin{align*}
\frac{\hat{H}_{12}\hat{H}_{21}}{E-\hat{H}_{22}}+\frac{\hat{H}_{10}\hat{H}_{01}}{E-\hat{H}_{00}}&=-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d+U}-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d}\\
&=-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d+U}-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c_{k^{\prime}\sigma}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{{d}{\sigma}}\hat{c}^{\dagger}_{k\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d}
\end{align*}
$$
​         最后一个等号利用指标互换，现在我们进一步考虑利用泡利算符的性质进行化简                                                         
$$
\hat{{\sigma}}_{\alpha\beta}\cdot\hat{\sigma}_{\gamma\delta}=2\delta_{\alpha\delta}\delta_{\beta\gamma}-\delta_{\alpha\beta}\delta_{\gamma\delta}
$$
​        现在我们定义$\hat{S}_{kk^{\prime}}=\frac{1}{2}\sum\limits_{\alpha,\beta}c_{k\alpha}^{\dagger}\hat{\sigma}_{\alpha\beta}\hat{c}_{k^{\prime}\beta}$,$\hat{s}_{d}=\frac{1}{2}\sum\limits_{\gamma,\delta}\hat{c}_{d\gamma}^{\dagger}\sigma_{\gamma\delta}\hat{c}_{d\delta}$,现在我们来看                                                         
$$
\begin{align*}
2\hat{S}_{kk^{\prime}}\cdot\hat{s}_{d}&=\frac{1}{2}\sum_{\alpha,\beta,\gamma,\delta}
\hat{c}_{k\alpha}^{\dagger}\hat{c}_{k^{\prime}\beta}\hat{c}_{d\gamma}^{\dagger}\hat{c}_{d\delta}
(2\delta_{\alpha\delta}\delta_{\beta\gamma}-\delta_{\alpha\beta}\delta_{\gamma\delta})\\
&=\sum_{\alpha,\beta}
\hat{c}_{k\alpha}^{\dagger}\hat{c}_{k^{\prime}\beta}\hat{c}_{d\beta}^{\dagger}\hat{c}_{d\alpha}
-\frac{1}{2}\sum_{\alpha,\gamma}
\hat{c}_{k\alpha}^{\dagger}\hat{c}_{k^{\prime}\alpha}\hat{c}_{d\gamma}^{\dagger}\hat{c}_{d\gamma}
\end{align*}
$$
​         我们可以得到                                                

$$
\begin{align*}
\frac{\hat{H}_{12}\hat{H}_{21}}{E-\hat{H}_{22}}+\frac{\hat{H}_{10}\hat{H}_{01}}{E-\hat{H}_{00}}
&=-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{c}_{{d}{\sigma}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d+U}-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c_{k^{\prime}\sigma^{\prime}}\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{{d}{\sigma}}\hat{c}^{\dagger}_{k\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d}\\
&=
-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}(\delta_{\sigma\sigma^{\prime}}-\hat{c}_{{d}{\sigma^{\prime}}}^{\dagger}\hat{c}_{{d}{\sigma}})\hat{c}_{k^{\prime}\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d+U}-\sum\limits_{k,k^{\prime},\sigma,\sigma^{\prime}}\frac{V_kV_{k^{\prime}}^{*}(\delta_{\sigma\sigma^{\prime}}-\hat{c}_{k\sigma}^{\dagger}\hat{c}_{{k}^{\prime}\sigma^{\prime}})c^{\dagger}_{d\sigma}\hat{c}_{d\sigma^{\prime}}}{\varepsilon_k-\varepsilon_d}\\
&=\sum_{k,k^{\prime},\sigma,\sigma^{\prime}}\left(
\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{c}_{k^{\prime}\sigma^{\prime}}\hat{c}_{{d}{\sigma}^{\prime}}^{\dagger}\hat{c}_{{d}{\sigma}}}{\varepsilon_k-\varepsilon_d+U}
+\frac{V_kV_{k^{\prime}}^{*}c^{\dagger}_{k\sigma}\hat{c}_{k^{\prime}\sigma^{\prime}}\hat{c}_{{d}{\sigma}^{\prime}}^{\dagger}\hat{c}_{{d}{\sigma}}}{\varepsilon_k-\varepsilon_d}\right)
-\sum_{k,k^{\prime},\sigma}\frac{\hat{c}_{k\sigma}^{\dagger}\hat{c}_{k^{\prime}\sigma}}{\varepsilon_k-\varepsilon_d+U}\\
&=\sum_{k,k^{\prime}}\left[\left(
\frac{1}{\varepsilon_k-\varepsilon_d+U}
+\frac{1}{\varepsilon_k-\varepsilon_d}\right)
\left(2\hat{S}_{kk^{\prime}}\cdot\hat{s}_{d}+{\color{red}{\frac{1}{2}\sum\limits_{\sigma,\sigma^{\prime}}
\hat{c}_{k\sigma}^{\dagger}\hat{c}_{k^{\prime}\sigma}\hat{c}_{d\sigma^{\prime}}^{\dagger}\hat{c}_{d\sigma^{\prime}}}}\right)-\sum\limits_{\sigma}\frac{\hat{c}_{k\sigma}^{\dagger}\hat{c}_{k^{\prime}\sigma}}{\varepsilon_k-\varepsilon_d+U}\right]
\end{align*}
$$
​                   

​          这里我们利用$\sum\limits_{\sigma}\hat{c}_{d\sigma^{\prime}}^{\dagger}\hat{c}_{d\sigma^{\prime}}\equiv1$,可以扔掉红色部分有关项，红色部分去掉一个求和指标变为                                               
$$
\frac{1}{2}\sum\limits_{\sigma}
\hat{c}_{k\sigma}^{\dagger}\hat{c}_{k^{\prime}\sigma}
$$
​           现在我们可以考虑引入                                                      
$$
\begin{align*}
&2J_{kk^{\prime}}=V_kV_{k^{\prime}}^{*}\left(
\frac{1}{\varepsilon_k-\varepsilon_d+U}
+\frac{1}{\varepsilon_k-\varepsilon_d}\right)\\
&K_{kk^{\prime}}=\frac{V_kV_{k^{\prime}}^{*}}{2}\left(
\frac{1}{\varepsilon_k-\varepsilon_d}-\frac{1}{\varepsilon_k-\varepsilon_d+U}\right)
\end{align*}
$$
​           最后有效哈密顿量可以写为                                     
$$
\hat{H}_{\text{eff}}=\sum\limits_{k,\sigma}\varepsilon_k \hat{c}^{\dagger}_{k,\sigma}c_{k,\sigma}+\varepsilon_d+\sum\limits_{k,k^{\prime}}J_{kk^{\prime}}
\hat{S}_{kk^{\prime}}\cdot\hat{s}_{d}+K_{kk^{\prime}}\sum\limits_{\sigma}
\hat{c}_{k\sigma}^{\dagger}\hat{c}_{k^{\prime}\sigma}
$$
