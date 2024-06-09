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


$$\
f(x) = \begin{bmatrix}
\frac{1}{m}{\frac{1}{R}(u_2 + u_3)cos(u_1)-(x_4+x_5)sin(u_1)+\frac{1}{R}(u_4+u_5)-C_{av}x_1^2}+x_2x_3 \\
\frac{1}{m}{\frac{1}{R}(u_2 + u_3)sin(u_1)+(x_4+x_5)cos(u_1)+(x_6+x_7)}+x_1x_3 \\
\frac{1}{m}[l_f{\frac{1}{R}(u_2 + u_3)sin(u_1)+(x_4+x_5)cos(u_1)}+ t{\frac{1}{R}(u2-u3)cos(u_1) + (-x_4+x_5)cos(u_1)+\frac{1}{R}(u4-u5)} -l_r(x_6+x_7)] \\
\frac{x_1}{\sigma}- (x_4+\overline{F}_{yfl}) \\

\frac{x_1}{\sigma}- (x_5+\overline{F}_{yfr}) \\ 

\frac{x_1}{\sigma}- (x_6+\overline{F}_{yrl}) \\ 

\frac{x_1}{\sigma}- (x_7+\overline{F}_{yrr}) \\ 
\end{bmatrix}
\$$
