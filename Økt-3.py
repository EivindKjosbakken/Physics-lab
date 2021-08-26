import matplotlib.pyplot as plt
import numpy as np

#read txt file
f = open("VideoLab\Vid1_table.txt", "r")
f.readline()
f.readline()

t = []
x = []
y = []

for lines in f:
    val = lines.split("\t")
    t.append(float(val[0]))
    x.append(float(val[1]))
    y.append(float(val[2]))

f.close()


#define constants
g = 9.81        #m/s2
c = 2/5         #
m = 0.031       #kg
r = 0.011       #m
#define lists
dt = []
for i in range(1,len(t)):
    dt_val = t[i]-t[i-1]
    dt.append(dt_val)

#calculate polynoms
y_p_val = np.polyfit(x,y,9)
y_p = np.poly1d(y_p_val)
d1y_p = y_p.deriv()
d2y_p = d1y_p.deriv()


#velocity
v_x = []
for i in x:
    v_val = np.sqrt((2*g*(y_p(x[0])-y_p(i)))/(1+c))
    v_x.append(v_val)

v_x_max = 0
for i in range(len(v_x)):
    if v_x[i] > v_x_max:
        v_x_max = v_x[i]
        i_max = i
print("v(x)_max =",v_x_max, x[i_max])

#curvature
κ = []
for i in x:
    κ_val = (d2y_p(i))/((1+(d1y_p(i))**(2))**(3/2))
    κ.append(κ_val)
κ_max = 0
for i in κ:
    if i > κ_max:
        κ_max = i
print("κ_max =",κ_max)

#acceleration
a = []
for i in range(0,len(x)):
    a_val = (v_x[i]**2)*κ[i]
    a.append(a_val)


#angle
β = []
for i in x:
    β_val = np.arctan(d1y_p(i))
    β.append(β_val)

#Normal force
N = []
for i in range(0,len(x)):
    N_val = m*(g*np.cos(β[i]) + a[i])
    N.append(N_val)

#friction
f = []
for i in β:
    f_val = (c*m*g*np.sin(i))/(1+c)
    f.append(f_val)

#friction/normal force
fN = []
for i in range(0,len(x)):
    fN_val = f[i]/N[i]
    fN.append(abs(fN_val))

fN_max = 0
for i in range(len(fN)):
    if fN[i] > fN_max:
        fN_max = fN[i]
        i_max = i
print("fN_max =",fN_max, x[i_max])

#time dependent x,y,v_x(t)
x_t = []
y_t = []
v_t = []
x_t.append(0), v_t.append(0), y_t.append(y_p(0))
for i in range(0,len(t)-1):
    a = -((5*g)/7)*np.sin(β[i])
    #next values
    x_n = x_t[i] + v_t[i]*dt[i]
    y_n = y_p(x_n)
    v_n = v_t[i] + a*dt[i]*np.cos(β[i])
    #append
    x_t.append(x_n), y_t.append(y_n), v_t.append(v_n)

#reverserer v_x til v siden vi kun bruker x_komponenten for å finne x(t)
for i in range(len(v_t)):
    v_t[i] = v_t[i]/(np.cos(β[i]))

v_t_max = 0
for i in range(len(v_t)):
    if v_t[i] > v_t_max:
        v_t_max = v_t[i]
        i_max = i
print("v(t)_max =",v_t_max, t[i_max])


#plots
plt.plot(x,y_p(x))
plt.legend(["y(x)"])
plt.xlabel("x[m]")
plt.ylabel("y[m]")
plt.show()
plt.plot(x,β)
plt.legend(["β(x)"])
plt.xlabel("x[m]")
plt.ylabel("vinkel[rad]")
plt.show()
plt.plot(x,κ)
plt.legend(["κ(x)"])
plt.xlabel("x[m]")
plt.ylabel("krummning[1/m]")
plt.show()
plt.plot(x,v_x)
plt.legend(["v(x)"])
plt.xlabel("x[m]")
plt.ylabel("v[m/s]")
plt.show()
plt.plot(x,N)
plt.legend(["N(x)"])
plt.xlabel("x[m]")
plt.ylabel("F[N]")
plt.show()
plt.plot(x,fN)
plt.legend(["n(x)"])
plt.xlabel("x[m]")
plt.ylabel("n[]")
plt.show()
plt.plot(t,x_t)
plt.legend(["x_num(t)"])
plt.xlabel("t[s]")
plt.ylabel("x[m]")
plt.show()
plt.plot(t,x)
plt.legend(["x_exp(t)"])
plt.xlabel("t[s]")
plt.ylabel("x[m]")
plt.show()
plt.plot(t,x_t)
plt.plot(t,x)
plt.legend(["x_num(t)","x_exp(t)"])
plt.xlabel("t[s]")
plt.ylabel("x[m]")
plt.show()
plt.plot(t,v_t)
plt.legend(["v(t)"])
plt.xlabel("t[s]")
plt.ylabel("v[m/s]")
plt.show()

print(v_x[-1])
