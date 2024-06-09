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
T_{di}-T{bi}\\
\end{bmatrix}
\$$

$$\
z = \begin{bmatrix}
V_x \\
V_y \\
\gamma \\
a_x \\
a_y \\
F_{xfl} \\
F_{xfr} \\
F_{xrl} \\
F_{xrr} \\
\end{bmatrix}
\$$
