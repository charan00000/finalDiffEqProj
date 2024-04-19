import matplotlib.pyplot as plot

"""
The following code implements Euler's method to solve a system of Lotka-Volterra equations, modified to account for cannibalistic behavior within each species.
x0 and y0 are initial populations of the prey and predator species, ALPHA is the prey reproduction rate, BETA is the prey death rate, GAMMA is the 
predator death rate, DELTA is the predator growth rate, A is the predator's cannibalism conversion constant, B is the cannibalism rate amongst predators, C is the 
half-saturation constant, size is the number of iterations (steps) simulated, and step is the change in x for each step of the Euler method, and show_graph is a 
boolean that toggles the graphing of a phase portrait involving the predator population as a function of the prey population.

method changes the type of plot formed (predator vs prey or predator & prey vs time). Lower size values are recommended when viewing phase plot for ease of viewing.
Phase plots arent supported for graphs with 
"""
method = "time"


#edit parameters, initial conditions here
"""
#generic values, sample A
file_name = "Simulation_A_Time_Plot.png"
ALPHA = 0.2
BETA = 0.1
GAMMA = 0.2
DELTA = 0.1
A = 0.05
B = 0.1
C = 0.5
x0 = 2
y0 = 2
step = .1
size = 1*10**6
show_graph = True
real_values_X = []  #leave empty if unknown or not applicable
real_values_y = []  #leaver empty if unknown or not applicable
method = "phase" #'time' or 'phase', phase plot not supported if real x, y values are provided. change file name based on method
"""

#polar bears v ringed seals, sample B, generated from polar_bears_bearded_seals.py
file_name = "Polar_Bear_V_Ringed_Seals_Time_Plot.png"
ALPHA = 0.6791379276083253
BETA = 0.024129879215543246
GAMMA = 1.081362947845921
DELTA = 0.07035954120293572
A = 0.022698993545745047
B = 1.0841016229644673
C = 0.12454929999882855
x0 = 32
y0 = 23
step = .1
size = 70
show_graph = True
real_values_X = [32, 45, 13, 9, 35, 43, 40]
real_values_y = [23, 43, 34, 46, 12, 31, 10]


def x_(x, y, ALPHA, BETA):
    return x*(ALPHA-BETA*y)

def y_(x, y, GAMMA, DELTA, A, B, C):
    return y*(-GAMMA+DELTA*x)-(B*y*y)/(C+y)+A*y

def euler_phase(ALPHA, BETA, GAMMA, DELTA, A, B, C, size, step, x0, y0, graph = False):
    points = [[0, 0, 0, 0, 0] for _ in range(size)]
    points[0][0] = x0
    points[0][1] = y0
    for i in range(1, size):
        x = points[i-1][0]
        y = points[i-1][1]
        points[i][0] = x + step * x_(x, y, ALPHA, BETA)
        points[i][1] = y + step * y_(x, y, GAMMA, DELTA, A, B, C)
        points[i][2] = x_(x, y, ALPHA, BETA)
        points[i][3] = y_(x, y, GAMMA, DELTA, A, B, C)
        points[i][4] = abs(points[i][2]) + abs(points[i][3])
    points = sorted(points, key = lambda x: x[4])
    for i in range(len(points) - 1, -1, -1):
        if points[i][0] < 0 or points[i][1] < 0:
            del points[i]
            size -= 1
    print(points)
    if graph:
        plot.scatter([points[i][0] for i in range(size)], [points[i][1] for i in range(size)])
        plot.xlabel("x (prey) population")
        plot.ylabel("y (predator) population")
        plot.savefig("Simulation_A_Phase_Plot.png")
        plot.show()
    return points

def euler_time(ALPHA, BETA, GAMMA, DELTA, A, B, C, size, step, x0, y0, graph = False):
    points = [[0, 0, 0, 0, 0, 0] for _ in range(size)]
    points[0][0] = 0
    points[0][1] = x0
    points[0][2] = y0
    for i in range(1, size):
        t = points[i-1][0]
        x = points[i-1][1]
        y = points[i-1][2]
        points[i][0] = t + step
        points[i][1] = x + step * x_(x, y, ALPHA, BETA)
        points[i][2] = y + step * y_(x, y, GAMMA, DELTA, A, B, C)
        points[i][3] = x_(x, y, ALPHA, BETA)
        points[i][4] = y_(x, y, GAMMA, DELTA, A, B, C)
        points[i][5] = abs(points[i][3]) + abs(points[i][4])
    for i in range(len(points) - 1, -1, -1):
        if points[i][1] < 0 or points[i][2] < 0:
            del points[i]
    points = sorted(points, key=lambda x: x[5])
    print(points)
    size = len(points)
    if graph:
        plot.xlabel("time (arbitrary units)")
        plot.ylabel("population (arbitrary units)")
        plot.scatter([points[i][0] for i in range(size)], [points[i][1] for i in range(size)])
        plot.scatter([points[i][0] for i in range(size)], [points[i][2] for i in range(size)])
        if len(real_values_X) > 0:
            plot.scatter([t for t in range(len(real_values_X))], real_values_X)
            plot.scatter([t for t in range(len(real_values_y))], real_values_y)
            plot.plot([t for t in range(len(real_values_X))], real_values_X)
            plot.plot([t for t in range(len(real_values_y))], real_values_y)  
            plot.xlabel("time (months)")     
            plot.legend(["simulated x (prey) population", "simulated y (predator) population", "true x (prey) population", "true y (predator) population"])
        else:
            plot.legend(["simulated x (prey) population", "simulated y (predator) population"])
        plot.savefig(file_name)
        plot.show()
    return points

if method == "phase":
    euler_phase(ALPHA, BETA, GAMMA, DELTA, A, B, C, size, step, x0, y0, graph = show_graph)
else:
    euler_time(ALPHA, BETA, GAMMA, DELTA, A, B, C, size, step, x0, y0, graph = show_graph)
