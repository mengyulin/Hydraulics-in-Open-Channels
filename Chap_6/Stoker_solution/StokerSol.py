"""
**************************************************************************
      Stoker's (1957) analytical solution for ideal dam-break
                 code converted to Python from MATLAB
          Original MATLAB code by: Payam Sarkhosh
           Research assistant at Prof. Yee-Chung Jin's Lab
             University of Regina, Saskatchewan, Canada
                            Fall 2021
                    Python conversion: 2025
**************************************************************************
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def stoker_dam_break():
    """
    Stoker's analytical solution for ideal dam-break problem
    """
    
    # Clear any existing plots
    plt.clf()
    
    # ******************************* inputs ***********************************
    print('*******************************************************************')
    print('******* Please enter the following inputs, and press ENTER ********')
    print('*******************************************************************')
    print('                                                                   ')
    
    h0 = float(input('   initial upstream water depth:   h0 (m) >> '))
    hd = float(input('   initial downstream water depth: hd (m) >> '))
    Lc = float(input('   channel length:                 Lc (m) >> '))
    Lr = float(input('   reservoir length:               Lr (m) >> '))
    T = float(input('   total simulation time            T (s) >> '))
    
    n = 2000            # Number of space intervals
    m = int(T * 200)    # Number of time intervals
    
    # **************************************************************************
    dx = Lc / n         # space step
    dt = T / m          # Time step
    g = 9.81            # Gravitational acceleration
    x = np.zeros(n)     # X vector
    u = np.zeros(n)     # Flow velocity vector
    h = np.zeros(n)     # Flow depth vector
    
    x_UBC = -Lr
    x_DBC = Lc - Lr
    x0 = 0
    
    # ************************ Plotting initial condition plot *****************
    x1 = -Lr
    x2 = 0
    xDam2 = np.linspace(x1, x2, n)
    Dam2 = h0 + 1e-50 * xDam2
    
    x1 = x0 + 1e-10
    xDam1 = np.linspace(x0, x1, n)
    Dam1 = h0 * (xDam1 - x0) / 1e-10
    
    # *************************** Mesh generation ******************************
    x[0] = x_UBC
    i_0 = 0
    for i in range(n-1):
        x[i+1] = x_UBC + (i+1) * dx
        if abs(x[i] - x0) < 0.5 * dx:
            i_0 = i
    
    x[n-1] = Lc - Lr
    h[0] = h0
    
    # ****************** definition of constant values**************************
    x2_end = -1e-20
    hA = 1e+20
    
    Cup = np.sqrt(g * h0)    # wave speed at upstream
    Cdown = np.sqrt(g * hd)  # wave speed at downstream
    
    # Setup interactive plotting
    plt.ion()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    for k in range(1, m+1):
        Time = k * dt
        
        # **************** Newton-Raphson iteration  ***************************
        f_CB = 1
        df_CB = 1
        CB = 10 * h0
        
        while abs(f_CB / CB) > 1e-10:
            sqrt_term = np.sqrt((8 * CB**2) / Cdown**2 + 1)
            bracket_term = (sqrt_term - 1)
            sqrt_hd_term = np.sqrt(g * hd * bracket_term / 2)
            
            f_CB = (CB * hd - hd * bracket_term * 
                   (CB/2 - Cup + sqrt_hd_term))
            
            df_CB = (hd - hd * ((2 * CB * g * hd) / 
                               (Cdown**2 * sqrt_term * sqrt_hd_term) + 1/2) * 
                    bracket_term - (8 * CB * hd * (CB/2 - Cup + sqrt_hd_term)) / 
                    (Cdown**2 * sqrt_term))
            
            CB = CB - f_CB / df_CB
        
        # *************** Newton-Raphson iteration end *************************
        
        hA = 0.5 * hd * (np.sqrt(1 + 8 * CB**2 / Cdown**2) - 1)
        if hd == 0:
            CB = 0
            hA = 0
        
        X2_start = (2 * np.sqrt(g * h0) - 3 * np.sqrt(g * hA)) * Time
        X2_end = CB * Time
        uA = 2 * Cup - 2 * np.sqrt(g * hA)
        
        for i in range(1, n):
            h[i] = (2 * np.sqrt(g * h0) - x[i] / Time)**2 / (9 * g)
            u[i] = 2 * (x[i] / Time + np.sqrt(g * h0)) / 3
            
            h[0] = h[1]
            u[0] = u[1]
            
            # ******************************************************************
            if h[i] >= h0:
                i_A = i
                h[i] = h0
                u[i] = 0
            
            if hA == 0 and i > 0 and h[i] > h[i-1]:
                h[i] = 0
                u[i] = 0
            
            if hA > 0:
                if x[i] <= X2_end and h[i] <= hA:
                    i_B = i
                    h[i] = hA
                    u[i] = uA
                elif x[i] > X2_end:
                    h[i] = hd
                    u[i] = 0
        
        # Plotting at regular intervals
        if (Time * m) % 1 == 0 or k == m:
            # Clear previous plots
            ax1.clear()
            ax2.clear()
            
            # Plot water depth
            ax1.plot(xDam2, Dam2, '--k', linewidth=1, label='Dam')
            ax1.plot(xDam1, Dam1, '--k', linewidth=1)
            ax1.plot(x, h, 'b', linewidth=3, label='Water depth')
            ax1.set_xlim([x_UBC, x_DBC])
            ax1.set_ylim([0, 1.1 * h0])
            ax1.set_ylabel('water depth (m)', fontsize=14)
            ax1.tick_params(axis='both', which='major', labelsize=14)
            ax1.set_title(f"Stoker's (1957) solution for ideal dam-break problem\nt = {Time:.2f} s", 
                         fontsize=15)
            
            # Plot velocity
            brown = [0.5, 0, 0]
            ax2.plot(x, u, color=brown, linewidth=3, label='Velocity')
            ax2.set_xlim([x_UBC, x_DBC])
            ax2.set_ylim([0, np.sqrt(g * h0) * 2.2])
            ax2.set_xlabel('x (m)', fontsize=14)
            ax2.set_ylabel('velocity (m/s)', fontsize=14)
            ax2.tick_params(axis='both', which='major', labelsize=14)
            
            plt.tight_layout()
            plt.pause(0.01)  # Small pause for animation effect
            
            if Time == T:
                print('                                                       ')
                print(f'******* Outputs at t = {Time:.2f} s **********')
                
                # Create output table
                df = pd.DataFrame({
                    'x': x,
                    'h': h,
                    'u': u
                })
                
                print(df.to_string(index=False, float_format='%.6f'))
                print(f'******* Outputs at t = {Time:.2f} s **********')
    
    plt.ioff()
    plt.show()
    
    return x, h, u, Time

def stoker_dam_break_batch(h0, hd, Lc, Lr, T, plot_animation=True):
    """
    Batch version of Stoker's dam-break solution without user input
    
    Parameters:
    h0: initial upstream water depth (m)
    hd: initial downstream water depth (m)  
    Lc: channel length (m)
    Lr: reservoir length (m)
    T: total simulation time (s)
    plot_animation: whether to show animated plots
    
    Returns:
    x, h, u, Time: position, depth, velocity arrays and final time
    """
    
    n = 2000            # Number of space intervals
    m = int(T * 200)    # Number of time intervals
    
    # **************************************************************************
    dx = Lc / n         # space step
    dt = T / m          # Time step
    g = 9.81            # Gravitational acceleration
    x = np.zeros(n)     # X vector
    u = np.zeros(n)     # Flow velocity vector
    h = np.zeros(n)     # Flow depth vector
    
    x_UBC = -Lr
    x_DBC = Lc - Lr
    x0 = 0
    
    # ************************ Plotting initial condition plot *****************
    x1 = -Lr
    x2 = 0
    xDam2 = np.linspace(x1, x2, n)
    Dam2 = h0 + 1e-50 * xDam2
    
    x1 = x0 + 1e-10
    xDam1 = np.linspace(x0, x1, n)
    Dam1 = h0 * (xDam1 - x0) / 1e-10
    
    # *************************** Mesh generation ******************************
    x[0] = x_UBC
    for i in range(n-1):
        x[i+1] = x_UBC + (i+1) * dx
    
    x[n-1] = Lc - Lr
    h[0] = h0
    
    # ****************** definition of constant values**************************
    Cup = np.sqrt(g * h0)    # wave speed at upstream
    Cdown = np.sqrt(g * hd)  # wave speed at downstream
    
    if plot_animation:
        plt.ion()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    for k in range(1, m+1):
        Time = k * dt
        
        # **************** Newton-Raphson iteration  ***************************
        f_CB = 1
        CB = 10 * h0
        
        while abs(f_CB / CB) > 1e-10:
            sqrt_term = np.sqrt((8 * CB**2) / Cdown**2 + 1)
            bracket_term = (sqrt_term - 1)
            sqrt_hd_term = np.sqrt(g * hd * bracket_term / 2)
            
            f_CB = (CB * hd - hd * bracket_term * 
                   (CB/2 - Cup + sqrt_hd_term))
            
            df_CB = (hd - hd * ((2 * CB * g * hd) / 
                               (Cdown**2 * sqrt_term * sqrt_hd_term) + 1/2) * 
                    bracket_term - (8 * CB * hd * (CB/2 - Cup + sqrt_hd_term)) / 
                    (Cdown**2 * sqrt_term))
            
            CB = CB - f_CB / df_CB
        
        # *************** Newton-Raphson iteration end *************************
        
        hA = 0.5 * hd * (np.sqrt(1 + 8 * CB**2 / Cdown**2) - 1)
        if hd == 0:
            CB = 0
            hA = 0
        
        X2_end = CB * Time
        uA = 2 * Cup - 2 * np.sqrt(g * hA)
        
        for i in range(1, n):
            h[i] = (2 * np.sqrt(g * h0) - x[i] / Time)**2 / (9 * g)
            u[i] = 2 * (x[i] / Time + np.sqrt(g * h0)) / 3
            
            h[0] = h[1]
            u[0] = u[1]
            
            # ******************************************************************
            if h[i] >= h0:
                h[i] = h0
                u[i] = 0
            
            if hA == 0 and i > 0 and h[i] > h[i-1]:
                h[i] = 0
                u[i] = 0
            
            if hA > 0:
                if x[i] <= X2_end and h[i] <= hA:
                    h[i] = hA
                    u[i] = uA
                elif x[i] > X2_end:
                    h[i] = hd
                    u[i] = 0
        
        # Plotting at regular intervals
        if plot_animation and ((Time * m) % 1 == 0 or k == m):
            ax1.clear()
            ax2.clear()
            
            # Plot water depth
            ax1.plot(xDam2, Dam2, '--k', linewidth=1)
            ax1.plot(xDam1, Dam1, '--k', linewidth=1)
            ax1.plot(x, h, 'b', linewidth=3)
            ax1.set_xlim([x_UBC, x_DBC])
            ax1.set_ylim([0, 1.1 * h0])
            ax1.set_ylabel('water depth (m)', fontsize=14)
            ax1.tick_params(axis='both', which='major', labelsize=14)
            ax1.set_title(f"Stoker's (1957) solution for ideal dam-break problem\nt = {Time:.2f} s", 
                         fontsize=15)
            
            # Plot velocity
            brown = [0.5, 0, 0]
            ax2.plot(x, u, color=brown, linewidth=3)
            ax2.set_xlim([x_UBC, x_DBC])
            ax2.set_ylim([0, np.sqrt(g * h0) * 2.2])
            ax2.set_xlabel('x (m)', fontsize=14)
            ax2.set_ylabel('velocity (m/s)', fontsize=14)
            ax2.tick_params(axis='both', which='major', labelsize=14)
            
            plt.tight_layout()
            plt.pause(0.01)
    
    if plot_animation:
        plt.ioff()
        plt.show()
    
    return x, h, u, Time

# Example usage
if __name__ == "__main__":
    print("Choose mode:")
    print("1. Interactive mode (with user input)")
    print("2. Batch mode (predefined parameters)")
    
    mode = input("Enter mode (1 or 2): ")
    
    if mode == "1":
        # Interactive mode - equivalent to original MATLAB code
        x, h, u, final_time = stoker_dam_break()
    else:
        # Batch mode with example parameters
        print("Running with example parameters:")
        print("h0 = 2.0 m, hd = 0.1 m, Lc = 200 m, Lr = 100 m, T = 10 s")
        x, h, u, final_time = stoker_dam_break_batch(h0=2.0, hd=0.1, Lc=200, Lr=100, T=10)
        
        # Display final results
        print('                                                       ')
        print(f'******* Final Outputs at t = {final_time:.2f} s **********')
        df = pd.DataFrame({'x': x, 'h': h, 'u': u})
        print(df.head(20).to_string(index=False, float_format='%.6f'))
        print('... (showing first 20 rows)')
        print(f'******* Final Outputs at t = {final_time:.2f} s **********')