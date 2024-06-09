# Torque_Vectoring_for_IJAT

# Todo

change the ioniq model
- 4-wheel model 

320hp maximum Torque 605[Nm]

tire - Magic Formula Model <- just use tire?


- search the gear and inwheel motor



Kalman Filter Model

$$\
X = \begin{bmatrix}
V_x \\
V_y \\
\gamma \\
F_{yfl} \\
F_{yfr} \\
F_{yrl} \\
F_{yrr} \\
\end{bmatrix}
\$$

$$\
u = \begin{bmatrix}
\frac{\delta_1 + \delta_2}{2} \\
T_{di}-T_{bi}\\
\end{bmatrix}
\$$

$$\
z = \begin{bmatrix}
V_x \\
V_y \\
\gamma \\
a_x \\
a_y \\
\end{bmatrix}
\$$


$$
f(x) = \begin{bmatrix}
\frac{1}{m} \left( \frac{1}{R}(u_2 + u_3) \cos(u_1) - (x_4 + x_5) \sin(u_1) + \frac{1}{R}(u_4 + u_5) - C_{av} x_1^2 \right) + x_2 x_3 \\
\frac{1}{m} \left( \frac{1}{R}(u_2 + u_3) \sin(u_1) + (x_4 + x_5) \cos(u_1) + (x_6 + x_7) \right) + x_1 x_3 \\
\frac{1}{m} \left( l_f \left( \frac{1}{R}(u_2 + u_3) \sin(u_1) + (x_4 + x_5) \cos(u_1) \right) + t \left( \frac{1}{R}(u_2 - u_3) \cos(u_1) + (-x_4 + x_5) \cos(u_1) + \frac{1}{R}(u_4 - u_5) \right) - l_r (x_6 + x_7) \right) \\
\frac{x_1}{\sigma} \left( -x_4 + \overline{F}_{yfl} \right) \\
\end{bmatrix}
$$

