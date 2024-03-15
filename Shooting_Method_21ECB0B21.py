""" 
    Write a program to solve the Boundary value problem (BVP) using Shooting Method
    y" - (1/2)(y^2  - 1)y' + y = 0
    y(0) = 0 y(2) = 1

    (1) Regula Falsi method
    (2) Newton Raphson method

    h = 0.01

    use (1) Euler's Method (2) RK 4th order Method
    
"""

# Inputs
h = 0.01                    # Step size
iter_length = int(2 / h)    # iter_length = 200

x0, y0, z0 = 0, 0, 1
x = x0

omega = 0.5     # omega = mu (from the Question)
alpha0 = 1      # Alpha is Guess value for y' (y derivative)

tolerance = 1e-4

eta = 0
mu = 1

# Regula Falsi Inputs
# First Guess Alpha 1
alpha1 = 1.2
# Second Guess Alpha 2
alpha2 = 1.5


# Selecting Root Finding Method

print("1. Regula Falsi Method")
print("2. Newton Raphson Method")

root_method_list = ['Regula Falsi Method', 'Newton Raphson Method']
root_method = int(input("\nEnter (1 / 2) for selecting Root Finding Method: "))

print("\n1. Euler Method")
print("2. RK 4th Order Method")

single_method_list = ['Euler\'s Method', 'RK 4th Order Method']
single_step_method = int(input("\nEnter (1 / 2) for selecting Single Step Method: "))

print("\nRoot Finding Method - {}".format(root_method_list[root_method - 1]))
print("Single-Step Method - {}\n".format(single_method_list[single_step_method - 1]))


if (root_method == 1):
    if (single_step_method == 1):
        output = open('Regula_Falsi_Euler.txt', 'w')
        output1 = open('Regula_Falsi_Euler_Single_Step_Calculations.txt', 'w')

        output.write("Root Finding Method - Regula Falsi Method\n")
        output.write("Single-Step Method - Euler's Method\n")
        output.write("\nShooting Method Calculations\n")

        output1.write("Root Finding Method - Regula Falsi Method\n")
        output1.write("Single-Step Method - Euler's Method\n")
        output1.write("\nSingle Step Calculations\n")

    elif (single_step_method == 2):
        output = open('Regula_Falsi_RK.txt', 'w')
        output1 = open('Regula_Falsi_RK_Single_Step_Calculations.txt', 'w')

        output.write("Root Finding Method - Regula Falsi Method\n")
        output.write("Single-Step Method - RK 4th Order Method\n")
        output.write("\nShooting Method Calculations\n")

        output1.write("Root Finding Method - Regula Falsi Method\n")
        output1.write("Single-Step Method - RK 4th Order Method\n")
        output1.write("\nSingle Step Calculations\n")


elif (root_method == 2):
    if (single_step_method == 1):
        output = open('Newton_Raphson_Euler.txt', 'w')
        output1 = open('Newton_Raphson_Euler_Single_Step_Calculations.txt', 'w')
        output2 = open('Newton_Raphson_Euler_Second_Order_Derivative_Eta_Calculations.txt', 'w')

        output.write("Root Finding Method - Newton Raphson Method\n")
        output.write("Single-Step Method - Euler's Method\n")
        output.write("\nShooting Method Calculations\n")

        output1.write("Root Finding Method - Newton Raphson Method\n")
        output1.write("Single-Step Method - Euler's Method\n")
        output1.write("\nSingle Step Calculations\n")

        output2.write("Root Finding Method - Newton Raphson Method\n")
        output2.write("Single-Step Method - Euler's Method\n")
        output2.write("\nSecond Order Derivative of Eta Calculation\n")

    elif (single_step_method == 2):
        output = open('Newton_Raphson_RK.txt', 'w')
        output1 = open('Newton_Raphson_RK_Single_Step_Calculations.txt', 'w')
        output2 = open('Newton_Raphson_RK_Second_Order_Derivative_Eta_Calculations.txt', 'w')

        output.write("Root Finding Method - Newton Raphson Method\n")
        output.write("Single-Step Method - RK 4th Order Method\n")
        output.write("\nShooting Method Calculations\n")

        output1.write("Root Finding Method - Newton Raphson Method\n")
        output1.write("Single-Step Method - RK 4th Order Method\n")
        output1.write("\nSingle Step Calculations\n")

        output2.write("Root Finding Method - Newton Raphson Method\n")
        output2.write("Single-Step Method - RK 4th Order Method\n")
        output2.write("\nSecond Order Derivative of Eta Calculation\n")

# Functions

# Van der Pol equation
def f(y, z, omega):
    return round(((omega * ((y ** 2) - 1) * z) - y), 8)
        
# Euler's method
def euler(y, f, h):
    return round((y + (h * f)), 8)

# RK Method
def rk(y, z, h, omega):
    # fx = f(y, z, omega)

    # k1 = h * z
    # l1 = h * f(z, fx, omega)

    # k2 = h * (z + (l1 * 0.5))
    # l2 = h * f(z + (0.5 * k1), fx + (l1 * 0.5), omega)

    # k3 = h * (z + (l2 * 0.5))
    # l3 = h * f(z + (0.5 * k2), fx + (l2 * 0.5), omega)

    # k4 = h * (z + (l3))
    # l4 = h * f(z + (k3), fx + (l3), omega)

    # k = (k1 + (2 * (k2 + k3)) + k4) / 6
    # l = (l1 + (2 * (l2 + l3)) + l4) / 6

    k1 = h * z
    l1 = h * f(y, z, omega)

    k2 = h * (z + (l1 * 0.5))
    l2 = h * f(y + (0.5 * k1), z + (l1 * 0.5), omega)

    k3 = h * (z + (l2 * 0.5))
    l3 = h * f(y + (0.5 * k2), z + (l2 * 0.5), omega)

    k4 = h * (z + (l3))
    l4 = h * f(y + (k3), z + (l3), omega)

    k = (k1 + (2 * (k2 + k3)) + k4) / 6
    l = (l1 + (2 * (l2 + l3)) + l4) / 6

    return (round((y + k), 8), round((z + l), 8))

# Derivative of the Van der Pol Equation wrt z
def df_dz(y, z, omega):
    return round((omega * ((y ** 2) - 1)), 8)

# Derivative of the Van der Pol Equation wrt y
def df_dy(y, z, omega):
    return round(((2 * omega * y * z) - 1), 8)


# Calculating ETA

# ETA_d2 = eta'' = ((df / dy) * eta) + ((df / dz) * mu(eta'))

def calc_eta_d2(z, iter_length, single_step_method):
    y = y0
    eta = 0
    mu = 1

    index = 0
    print("Calculating Second Order Derivative of Eta")
    print("Iteration: {}".format(index))
    print("y{} = {}".format(index, y))
    print("z = y'{} = {}".format(index, z))
    print("eta = {}".format(eta))
    print("mu = {}".format(mu))

    output2.write("\nCalculating Second Order Derivative of Eta\n")
    output2.write("\nIteration: {}\n".format(index))
    output2.write("y{} = {}\n".format(index, y))
    output2.write("z = y'{} = {}\n".format(index, z))
    output2.write("eta = {}\n".format(eta))
    output2.write("mu = {}\n".format(mu))

    
    for _ in range(iter_length):
        index = index + 1

        fx= f(y, z, omega)
        dfx_dy = df_dy(y, z, omega)
        dfx_dz = df_dz(y, z, omega)
                
        temp_mu = (dfx_dy * eta) + (dfx_dz * mu)

        # Euler Method
        if (single_step_method == 1):
            eta = euler(eta, mu, h)
            mu = euler(mu, temp_mu, h)
            
            y = euler(y, z, h)
            z = euler(z, fx, h)

        # RK 4th Order Method
        elif (single_step_method == 2):
            eta, mu = rk(eta, mu, h, omega)

            y, z = rk(y, z, h, omega)

        print("Iteration: {}".format(index))
        print("y{} = {}".format(index, y))
        print("z = y'{} = {}".format(index, z))
        print("eta = {}".format(eta))
        print("mu = {}".format(mu))

        output2.write("\nIteration: {}\n".format(index))
        output2.write("y{} = {}\n".format(index, y))
        output2.write("z = y'{} = {}\n".format(index, z))
        output2.write("eta = {}\n".format(eta))
        output2.write("mu = {}\n".format(mu))

        
    return round(eta, 8)


# Newton-Raphson Method

def Calc_Newton_Raphson(x0, y0, z0, omega, single_step_method):
    x = x0
    y = y0
    z = z0
    
    index = 0

    print("Interpolating to find y(2)")
    print("Iteration: {}\n".format(index))
    print("x = {}".format(round(x, 2)))
    print("y({}) = {}".format(round(x, 2), round(y, 8)))
    print("z = y'({}) = {}".format(round(x, 2), round(z, 8)))

    output1.write("\nInterpolating to find y(2)\n")
    output1.write("\nIteration: {}\n".format(index))
    output1.write("x = {}\n".format(round(x, 2)))
    output1.write("y({}) = {}\n".format(round(x, 2), round(y, 8)))
    output1.write("z = y'({}) = {}\n".format(round(x, 2), round(y, 8)))

    for _ in range(iter_length):
        index = index + 1
        x = x + h

        fx = f(y, z, omega)

        # Euler Method
        if (single_step_method == 1):
            y = euler(y, z, h)                 # Updating y using Euler Method
            z = euler(z, fx, h)                # Updating z using Euler Method

        # RK 4th Order Method
        elif(single_step_method == 2):
            y, z = rk(y, z, h, omega)

        print("Iteration: {}\n".format(index))
        print("x = {}".format(round(x, 2)))
        print("y({}) = {}".format(round(x, 2), round(y, 8)))
        print("z = y'({}) = {}".format(round(x, 2), round(z, 8)))

        output1.write("\nIteration: {}\n".format(index))
        output1.write("x = {}\n".format(round(x, 2)))
        output1.write("y({}) = {}\n".format(round(x, 2), round(y, 8)))
        output1.write("z = y'({}) = {}\n".format(round(x, 2), round(z, 8)))

    return round(y, 8)

# Regula Falsi

def Calc_Reg_Falsi(y0, z0, omega, single_step_method, id):
    x = x0
    y = y0
    z = z0
    
    index = 0

    if (id == 1):
        print("Shooting Method 1")
        output1.write("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        output1.write("\nShooting Method 1\n")
    elif (id == 2):
        print("Shooting Method 2")
        output1.write("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        output1.write("\nShooting Method 2\n")

    print("Interpolating to find y(2)")
    print("Iteration: {}".format(index))
    print("x = {}".format(x))
    print("y({}) = {}".format(x, y))
    print("z = y'({}) = {}".format(x, z))

    output1.write("\nInterpolating to find y(2)\n")
    output1.write("\nIteration: {}\n".format(index))
    output1.write("x = {}\n".format(x))
    output1.write("y({}) = {}\n".format(x, y))
    output1.write("z = y'({}) = {}\n".format(x, z))
    
    for _ in range(iter_length):
        index = index + 1
        x = x + h

        # Evaluate the functions and their derivatives at the current point
        fx = f(y, z, omega)

        if (single_step_method == 1):
            # Update y using Newton-Raphson iteration
            y = euler(y, z, h)
            # print('y(2) =  ',y)
            
            # Update z using Newton-Raphson iteration
            z = euler(z, fx, h)
            # print('z : ',z)

        elif(single_step_method == 2):
            y, z = rk(y, z, h, omega)

        print("Iteration: {}\n".format(index))
        print("x = {}".format(round(x, 2)))
        print("y({}) = {}".format(round(x, 2), round(y, 8)))
        print("z = y'({}) = {}".format(round(x, 8), round(z, 8)))

        output1.write("\nIteration: {}\n".format(index))
        output1.write("x = {}\n".format(round(x, 2)))
        output1.write("y({}) = {}\n".format(round(x, 8), round(y, 8)))
        output1.write("z = y'({}) = {}\n".format(round(x, 8), round(z, 8)))

    return round(y, 8)
    
# Newton Raphson Error

def Error_Newton_Raphson(y_shooting, single_step_method):
    y_actual = 1
    tolerance = 1e-4
    alpha = 1

    error = y_shooting - y_actual
    eta = calc_eta_d2(z0, iter_length, single_step_method)

    index = 0

    output.write("\nAll the Iterations for finding Alpha\n")
    print("All the Iterations for finding Alpha\n")

    while(abs(error) >= tolerance):
        print("Error Iteration: {}".format(index))
        output1.write("\n_________________________________________________________________________________________\n")
        output1.write("\nError Iteration: {}\n".format(index))

        print("Error Iteration: {}".format(index))
        output2.write("\n_________________________________________________________________________________________\n")
        output2.write("\nError Iteration: {}\n".format(index))

        alpha = alpha - (error / eta)

        eta = calc_eta_d2(alpha, iter_length, single_step_method)
        y_shooting = Calc_Newton_Raphson(x0, y0, alpha, omega, single_step_method)

        error = y_shooting - y_actual  # Error Calculation

        print("\nIteration - {}".format(index))
        print("y(2) = {}".format(round(y_shooting, 8)))
        print("error = {}".format(error, 8))
        print("Alpha = {}".format(round(alpha, 8)))

        output.write("\nIteration - {}\n".format(index))
        output.write("y(2) = {}\n".format(round(y_shooting, 8)))
        output.write("Error = {}\n".format(error, 8))
        output.write("Alpha = {}\n".format(round(alpha, 8)))

        index = index + 1

    print("\nFinal Solution")
    print("y(2) = {}".format(round(y_shooting, 8)))
    print("Error = {}".format(error, 8))
    print("Alpha = {}".format(round(alpha, 8)))

    output.write("\nFinal Solution")
    output.write("\ny(2) = {}\n".format(round(y_shooting, 8)))
    output.write("Error = {}\n".format(error, 8))
    output.write("Alpha = {}\n".format(round(alpha, 8)))

    return round(y_shooting, 8)

# Regula Falsi Error

def Error_Regula_Falsi(y_shooting1, y_shooting2, alpha1, alpha2):
    y_actual = 1
    tolerance = 1e-4

    # alpha1 = 1              # Guess 1 for y'
    # alpha2 = 1.5            # Guess 2 for y'


    error1 = abs(y_shooting1 - y_actual)  # Error = Diff btw Actual Value and Calculated Value
    error2 = abs(y_shooting2 - y_actual)  # Error

    index = 0

    output.write("\nAll the Iterations for finding Alpha\n")
    print("All the Iterations for finding Alpha\n")

    while(abs(error1) >= tolerance):
        print("Error Iteration: {}".format(index))
        output1.write("\n_________________________________________________________________________________________\n")
        output1.write("\nError Iteration: {}\n".format(index))

        alpha1 = ((alpha1 * error2) - (alpha2 * error1)) / (error2 - error1)

        y_shooting1 = Calc_Reg_Falsi(y0, alpha1, omega, single_step_method, 0)

        error1 = y_shooting1 - y_actual

        print("Iteration: {}\n".format(index))
        print("y(2) = {}".format(round(y_shooting1, 8)))
        print("Error = {}".format(error1, 8))
        print("Alpha = {}".format(round(alpha1, 8)))

        output.write("\nIteration - {}\n".format(index))
        output.write("y(2) = {}\n".format(round(y_shooting1, 8)))
        output.write("Error = {}\n".format(error1, 8))
        output.write("Alpha = {}\n".format(round(alpha1, 8)))

        index = index + 1

    print("\nFinal Solution")
    print("y(2) = {}".format(round(y_shooting1, 8)))
    print("Error = {}".format(error1, 8))
    print("Alpha = {}".format(round(alpha1, 8)))

    output.write("\nFinal Solution")
    output.write("\ny(2) = {}\n".format(round(y_shooting1, 8)))
    output.write("Error = {}\n".format(error1, 8))
    output.write("Alpha = {}\n".format(round(alpha1, 8)))

    return round(y_shooting1, 8)


if (root_method == 1):

    # Regula Falsi Method
    y_shooting1 = Calc_Reg_Falsi(y0, alpha1, omega, single_step_method, 1)
    y_shooting2 = Calc_Reg_Falsi(y0, alpha2, omega, single_step_method, 2)

    y_final = Error_Regula_Falsi(y_shooting1, y_shooting2, alpha1, alpha2)

elif (root_method == 2):

    # Newton Raphson
    y_shooting = Calc_Newton_Raphson(x0, y0, z0, omega, single_step_method)
    y_final = Error_Newton_Raphson(y_shooting, single_step_method)

output.close()