import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ! In the following code, the reaction kinetics is considered to follow the power law

# Calculating the gradient of the linear correlation
def gradient(conc,point): 
    return (conc[point+1]-conc[point])/(time[point+1]-time[point])

# Checking if the user wants to select the standard concentration data or wants to insert theirs
type_of_data = int(input("Type 1 if you want to use the standard concentration data, 2 if you want to insert yours and 3 if you want to import from a csv file "))

# Standard concentration data
if type_of_data == 1:
    time = [0,10,20,30,40,50] # in minutes
    ca = [1.,0.8,0.65,0.52,0.42,0.35] # mol/lt

# Concentration data in case they are added by the user
elif type_of_data == 2:
    n_of_elements = int(input("How many data points do you want to insert? "))
    time = []
    ca = []
    print("\nAdd each data point separately")
    for i in range(n_of_elements):
        time.append(float(input("\nEnter the time at which you have data: ")))
        ca.append(float(input("Enter the concentration of the reactant in that time: ")))

# Concentration data in case they are imported from a csv file
elif type_of_data == 3:
    file_path = input("Enter the file path ")
    df = pd.read_csv(file_path, delimiter=';')
    time = df['time'].tolist()
    ca = df['ca'].tolist()

# In case the user inserts a wrong input
else:
    print("Wrong input")
    exit()

# Initialization of some variables
mean_k = 0
best_corr= 0

# Plotting concentration of reactant as a function of time
plt.plot(time,ca,'o')
plt.title("Concentration of the reactant function with time")
plt.xlabel("Time (min)")
plt.ylabel("Concentration of the reactant (mol/lt)")
plt.grid(True)
plt.savefig('plot_conc.png')
plt.show()

# ! If it is a first grade reaction, the logarithm of the reactant's concentration has a linear correlation with time, with a gradient of -k
# Checking if reaction is first grade in respect of the reactant
ca_new = np.log(ca) 
data = {'time': time, 'log_ca': ca_new}
df = pd.DataFrame(data)
correlation = df.corr().loc['time','log_ca'] # Finding the correlation between time and the logarithm of the reactant's concentration

if abs(correlation) > 0.999: # Checking if the correlation is linear
    # Calculating the kinetic's constant for each neighbor pair of data and then its mean value
    for i in range(len(time)-1):
        k = - gradient(ca_new,i)
        mean_k = mean_k + k
    mean_k = mean_k/(len(time)-1)
    print("\nThe rate of reaction is equal to r = -","%.4f" % mean_k,"*ca\n")

# ! If it is an n grade reaction, then the reactant's consentration raised to (1-n) has a linear correlation with time, with a gradient of -(1-n)k
else: # In case the correlation is NOT linear
    for i in range(0,61):
        n = i/20 # Checking for different grades of reaction. It checks for a grade between 0 and 3 with a step of 0.05
        if n == 1: # Excluding the first grade reaction
            continue
        ca_n = [i**(1-n) for i in ca] 
        data = {'time': time, 'ca_n': ca_n}
        df = pd.DataFrame(data)
        correlation = df.corr().loc['time','ca_n'] # Finding the correlation between time and the reactant's concentration raised to 1-n
        if abs(correlation) > best_corr and abs(correlation) > 0.999: # We need the best correlation, which simultaneously is greater than 0.995 
            best_corr = abs(correlation)
            n_fin = n
            # Calculating the kinetic's constant for each neighbor pair of data and then its mean value
            mean_k = 0
            for i in range(len(time)-1):
                k = gradient(ca_n,i)/(n-1)
                mean_k = mean_k + k
            mean_k = mean_k/(len(time)-1)

    # Checking if there is a correlation that is linear
    if best_corr != 0:
        print("\nthe rate of reaction is equal to r = -","%.4f" % mean_k,"*ca^",n_fin,"\n")
    else:
        print("\nthe rate of reaction doesn't follow the power law\n")
