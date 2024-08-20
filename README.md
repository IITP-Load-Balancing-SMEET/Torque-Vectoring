# Torque Vectoring (IJAT)

<div align='center'>
 
![image](https://github.com/IITP-Load-Balancing-SMEET/Torque-Vectoring/assets/86957779/dc042e82-fc9f-4b85-af4a-72224ab99e6f)

</div>


## Todo

change the ioniq model
- 4-wheel model 

320hp maximum Torque 605[Nm]

tire - Dugoff Tire Model

- search the gear and inwheel motor



## Kalman Filter Model
$$
X = \begin{bmatrix}
V_x \\
V_y \\
\gamma \\
F_{yfl} \\
F_{yfr} \\
F_{yrl} \\
F_{yrr} \\
\end{bmatrix} \in{R^{7}}
$$ 

$$
u = \begin{bmatrix}
\frac{\delta_1 + \delta_2}{2} \\
T_{fl} = T_{dfl}-T_{bfl}\\
T_{fr} = T_{dfr}-T_{bfr}\\
T_{rl} = T_{drl}-T_{brl}\\
T_{rl} = T_{drl}-T_{brl}\\
\end{bmatrix}\in{R^{9}}
$$

$$
z = h(X)+\textbf{V}=\begin{bmatrix}
V_x \\
V_y \\
\gamma \\
a_x \\
a_y \\
\end{bmatrix}\in{R^{5}}
$$

---

### System model

$$
f(X) = \begin{bmatrix}
\frac{1}{m} \left( \frac{1}{R}(u_2 + u_3) \cos(u_1) - (x_4 + x_5) \sin(u_1) + \frac{1}{R}(u_4 + u_5) - C_{av} x_1^2 \right) + x_2 x_3 \\
\frac{1}{m} \left( \frac{1}{R}(u_2 + u_3) \sin(u_1) + (x_4 + x_5) \cos(u_1) + (x_6 + x_7) \right) - x_1 x_3 \\
\frac{1}{m} \left( l_f \left( \frac{1}{R}(u_2 + u_3) \sin(u_1) + (x_4 + x_5) \cos(u_1) \right) + t \left( \frac{1}{R}(u_2 - u_3) \cos(u_1) + (-x_4 + x_5) \cos(u_1) + \frac{1}{R}(u_4 - u_5) \right) - l_r (x_6 + x_7) \right) \\
\frac{x_1}{\sigma} \left( -x_i + \overline{F}_{yj} \right) \\
\end{bmatrix}
$$

where $i$=[4, 5, 6, 7], $j$=[fl, fr, rl, rr]
 
---

### Mesurement model

$$
h(X) = \begin{bmatrix}
x_1 \ (=V_x) \\
x_2 \ (=V_y) \\
x_3 \ (=\gamma) \\
\frac{1}{m}[\frac{(u_2+u_3)}{R}cos(u_1) - (x_4+x_5)sin(u_1)+\frac{(u_4+u_5)}{R}-C_{av}x_1^2] \\
\frac{1}{m}[\frac{(u_2+u_3)}{R}sin(u_1) + (x_4+x_5)cos(u_1)+x_6+x_7]
\end{bmatrix}, where \ R \ is \ radius \ of \ tire
$$

---
