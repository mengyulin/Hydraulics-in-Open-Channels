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
from matplotlib.animation import FuncAnimation

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

def plot_characteristics(h0, hd, Lc, Lr, T, x_dam=0):
    """
    Plot characteristic lines and show wave propagation
    
    Parameters:
    h0: initial upstream water depth (m)
    hd: initial downstream water depth (m)
    Lc: channel length (m)
    Lr: reservoir length (m)
    T: total time (s)
    x_dam: dam position (m, default 0)
    """
    
    g = 9.81
    Cup = np.sqrt(g * h0)    # upstream wave speed
    Cdown = np.sqrt(g * hd)  # downstream wave speed
    
    # Time array
    t = np.linspace(0, T, 1000)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Upper plot: Characteristic diagram (x-t plot)
    ax1.set_xlabel('Position x (m)', fontsize=12)
    ax1.set_ylabel('Time t (s)', fontsize=12)
    ax1.set_title('Characteristic Lines - Dam-break Wave Propagation', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([-Lr, Lc-Lr])
    ax1.set_ylim([0, T])
    
    # Plot dam location
    ax1.axvline(x=x_dam, color='black', linewidth=3, label='Dam location')
    
    # Characteristic line 1: Rightward wave front (positive characteristic)
    # Rarefaction wave front: x = 2*sqrt(g*h0)*t
    x_front = 2 * Cup * t
    valid_front = x_front <= (Lc - Lr)
    ax1.plot(x_front[valid_front], t[valid_front], 'r-', linewidth=2, 
             label=f'Positive wave front (speed = 2√(gh₀) = {2*Cup:.2f} m/s)')
    
    # Characteristic line 2: Rightward wave tail
    # Rarefaction wave tail: x = -sqrt(g*h0)*t  
    x_tail = -Cup * t
    valid_tail = x_tail >= -Lr
    ax1.plot(x_tail[valid_tail], t[valid_tail], 'b-', linewidth=2,
             label=f'Negative wave tail (speed = -√(gh₀) = {-Cup:.2f} m/s)')
    
    # If downstream depth exists, calculate shock characteristic
    if hd > 0:
        # Use Newton-Raphson method to calculate shock velocity
        CB_values = []
        for time_val in t:
            if time_val > 0:
                f_CB = 1
                CB = 10 * h0
                
                # Newton-Raphson iteration
                iteration_count = 0
                while abs(f_CB / CB) > 1e-10 and iteration_count < 100:
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
                    iteration_count += 1
                
                CB_values.append(CB)
            else:
                CB_values.append(0)
        
        CB_values = np.array(CB_values)
        
        # Shock position
        x_shock = CB_values * t
        valid_shock = (x_shock <= (Lc - Lr)) & (t > 0)
        ax1.plot(x_shock[valid_shock], t[valid_shock], 'g-', linewidth=2,
                 label=f'Shock front (avg speed ≈ {np.mean(CB_values[CB_values>0]):.2f} m/s)')
        
        # Intersection of rarefaction wave and shock
        # Time and position where rarefaction front meets shock
        for i in range(len(t)-1):
            if t[i] > 0 and abs(x_front[i] - x_shock[i]) < 0.1:
                ax1.plot(x_front[i], t[i], 'ko', markersize=8, 
                         label='Wave intersection' if i == np.argmin(np.abs(x_front - x_shock)) else "")
                break
    
    # Internal characteristics of rarefaction wave
    n_chars = 8  # Number of characteristics
    h_levels = np.linspace(h0*0.1, h0*0.9, n_chars)
    colors = plt.cm.viridis(np.linspace(0, 1, n_chars))
    
    for i, h_level in enumerate(h_levels):
        # Internal rarefaction characteristics: x = (2*sqrt(g*h0) - 3*sqrt(g*h))*t
        C_level = np.sqrt(g * h_level)
        x_char = (2 * Cup - 3 * C_level) * t
        
        # Only plot within rarefaction wave range
        valid_char = (x_char >= x_tail) & (x_char <= x_front) & (t > 0)
        if np.any(valid_char):
            ax1.plot(x_char[valid_char], t[valid_char], '--', 
                     color=colors[i], alpha=0.6, linewidth=1)
    
    ax1.legend(loc='upper left', fontsize=10)
    
    # Lower plot: Water depth profiles at different times
    ax2.set_xlabel('Position x (m)', fontsize=12)
    ax2.set_ylabel('Water depth h (m)', fontsize=12)
    ax2.set_title('Water Depth Profiles at Different Times', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([-Lr, Lc-Lr])
    ax2.set_ylim([0, h0*1.2])
    
    # Plot initial condition
    x_plot = np.linspace(-Lr, Lc-Lr, 1000)
    h_initial = np.where(x_plot < x_dam, h0, hd)
    ax2.plot(x_plot, h_initial, 'k--', linewidth=2, alpha=0.5, label='Initial depth')
    
    # Plot water depth distribution at different times
    time_snapshots = [T*0.2, T*0.4, T*0.6, T*0.8, T*1.0]
    colors_time = ['red', 'orange', 'green', 'blue', 'purple']
    
    for i, time_snap in enumerate(time_snapshots):
        if time_snap > 0:
            x_positions = np.linspace(-Lr, Lc-Lr, 500)
            h_profile = np.zeros_like(x_positions)
            
            # Calculate wave speeds and positions at this time
            x_front_snap = 2 * Cup * time_snap
            x_tail_snap = -Cup * time_snap
            
            # Calculate shock velocity (if downstream depth exists)
            if hd > 0:
                f_CB = 1
                CB = 10 * h0
                iteration_count = 0
                while abs(f_CB / CB) > 1e-10 and iteration_count < 100:
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
                    iteration_count += 1
                
                x_shock_snap = CB * time_snap
                hA = 0.5 * hd * (np.sqrt(1 + 8 * CB**2 / Cdown**2) - 1)
            else:
                x_shock_snap = float('inf')
                hA = 0
            
            # Build depth profile
            for j, x_pos in enumerate(x_positions):
                if x_pos < x_tail_snap:
                    # Undisturbed upstream region
                    h_profile[j] = h0
                elif x_pos <= x_front_snap and hd == 0:
                    # Rarefaction wave region (dry bed case)
                    h_profile[j] = (2 * Cup - x_pos / time_snap)**2 / (9 * g)
                    if h_profile[j] < 0:
                        h_profile[j] = 0
                elif x_pos <= x_shock_snap and hd > 0:
                    # Rarefaction wave region (wet bed case)
                    h_profile[j] = (2 * Cup - x_pos / time_snap)**2 / (9 * g)
                    if h_profile[j] <= hA:
                        h_profile[j] = hA
                else:
                    # Downstream undisturbed region
                    h_profile[j] = hd
            
            ax2.plot(x_positions, h_profile, color=colors_time[i], 
                     linewidth=2, label=f't = {time_snap:.1f} s')
    
    ax2.axvline(x=x_dam, color='black', linewidth=2, alpha=0.5, label='Dam location')
    ax2.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def animate_dam_break_with_characteristics(h0, hd, Lc, Lr, T, save_animation=False):
    """
    Create animation of dam-break wave propagation with characteristic lines
    
    Parameters:
    h0: initial upstream water depth (m)
    hd: initial downstream water depth (m)
    Lc: channel length (m)
    Lr: reservoir length (m)
    T: total time (s)
    save_animation: whether to save animation as gif file
    """
    
    g = 9.81
    Cup = np.sqrt(g * h0)
    Cdown = np.sqrt(g * hd) if hd > 0 else 0
    
    # Animation parameters
    n_frames = 200
    time_array = np.linspace(0.01, T, n_frames)  # Avoid singularity at t=0
    
    # Create figure
    fig = plt.figure(figsize=(15, 10))
    
    # Subplot 1: Characteristic diagram
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_xlabel('Position x (m)')
    ax1.set_ylabel('Time t (s)')
    ax1.set_title('Characteristic Lines - Real-time Wave Propagation')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([-Lr, Lc-Lr])
    ax1.set_ylim([0, T])
    
    # Subplot 2: Current water depth distribution
    ax2 = plt.subplot(2, 2, 2)
    ax2.set_xlabel('Position x (m)')
    ax2.set_ylabel('Water depth h (m)')
    ax2.set_title('Current Water Depth Distribution')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([-Lr, Lc-Lr])
    ax2.set_ylim([0, h0*1.2])
    
    # Subplot 3: Current velocity distribution
    ax3 = plt.subplot(2, 2, 3)
    ax3.set_xlabel('Position x (m)')
    ax3.set_ylabel('Velocity u (m/s)')
    ax3.set_title('Current Velocity Distribution')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([-Lr, Lc-Lr])
    ax3.set_ylim([0, 2*Cup*1.1])
    
    # Subplot 4: Wave position tracking
    ax4 = plt.subplot(2, 2, 4)
    ax4.set_xlabel('Time t (s)')
    ax4.set_ylabel('Wave position x (m)')
    ax4.set_title('Wave Front Position vs Time')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, T])
    ax4.set_ylim([-Lr, Lc-Lr])
    
    # Pre-calculate all time step data
    all_data = []
    CB_history = []
    
    for t_current in time_array:
        # Calculate shock velocity
        if hd > 0:
            f_CB = 1
            CB = 10 * h0
            iteration_count = 0
            while abs(f_CB / CB) > 1e-10 and iteration_count < 100:
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
                iteration_count += 1
            
            hA = 0.5 * hd * (np.sqrt(1 + 8 * CB**2 / Cdown**2) - 1)
        else:
            CB = 0
            hA = 0
        
        CB_history.append(CB)
        
        # Calculate wave positions
        x_front = 2 * Cup * t_current
        x_tail = -Cup * t_current
        x_shock = CB * t_current if hd > 0 else float('inf')
        
        # Calculate depth and velocity distributions
        n_points = 500
        x_positions = np.linspace(-Lr, Lc-Lr, n_points)
        h_profile = np.zeros(n_points)
        u_profile = np.zeros(n_points)
        
        for i, x_pos in enumerate(x_positions):
            if x_pos < x_tail:
                # Undisturbed upstream region
                h_profile[i] = h0
                u_profile[i] = 0
            elif x_pos <= x_front and hd == 0:
                # Rarefaction wave region (dry bed)
                h_profile[i] = (2 * Cup - x_pos / t_current)**2 / (9 * g)
                u_profile[i] = 2 * (x_pos / t_current + Cup) / 3
                if h_profile[i] < 0:
                    h_profile[i] = 0
                    u_profile[i] = 0
            elif x_pos <= x_shock and hd > 0:
                # Rarefaction wave region (wet bed)
                h_calc = (2 * Cup - x_pos / t_current)**2 / (9 * g)
                u_calc = 2 * (x_pos / t_current + Cup) / 3
                
                if h_calc <= hA:
                    h_profile[i] = hA
                    u_profile[i] = 2 * Cup - 2 * np.sqrt(g * hA)
                else:
                    h_profile[i] = h_calc
                    u_profile[i] = u_calc
            else:
                # Downstream undisturbed region
                h_profile[i] = hd
                u_profile[i] = 0
        
        all_data.append({
            'time': t_current,
            'x_positions': x_positions.copy(),
            'h_profile': h_profile.copy(),
            'u_profile': u_profile.copy(),
            'x_front': x_front,
            'x_tail': x_tail,
            'x_shock': x_shock,
            'CB': CB
        })
    
    CB_history = np.array(CB_history)
    
    def animate(frame):
        # Clear all subplots
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        
        current_data = all_data[frame]
        current_time = current_data['time']
        
        # Reset subplot 1: Characteristic diagram
        ax1.set_xlabel('Position x (m)')
        ax1.set_ylabel('Time t (s)')
        ax1.set_title(f'Characteristic Lines - t = {current_time:.2f} s')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim([-Lr, Lc-Lr])
        ax1.set_ylim([0, T])
        
        # Plot dam location
        ax1.axvline(x=0, color='black', linewidth=3, label='Dam')
        
        # Plot characteristics up to current time
        t_history = time_array[:frame+1]
        
        # Positive wave front characteristic
        x_front_history = 2 * Cup * t_history
        ax1.plot(x_front_history, t_history, 'r-', linewidth=2, label='Positive wave front')
        
        # Negative wave tail characteristic
        x_tail_history = -Cup * t_history
        ax1.plot(x_tail_history, t_history, 'b-', linewidth=2, label='Negative wave tail')
        
        # Shock characteristic
        if hd > 0:
            CB_current = CB_history[:frame+1]
            x_shock_history = CB_current * t_history
            ax1.plot(x_shock_history, t_history, 'g-', linewidth=2, label='Shock')
        
        # Current time horizontal line
        ax1.axhline(y=current_time, color='orange', linestyle='--', alpha=0.7, label='Current time')
        
        # Current wave position points
        ax1.plot(current_data['x_front'], current_time, 'ro', markersize=8)
        ax1.plot(current_data['x_tail'], current_time, 'bo', markersize=8)
        if hd > 0 and current_data['x_shock'] < Lc-Lr:
            ax1.plot(current_data['x_shock'], current_time, 'go', markersize=8)
        
        ax1.legend(fontsize=8)
        
        # Subplot 2: Current water depth distribution
        ax2.set_xlabel('Position x (m)')
        ax2.set_ylabel('Water depth h (m)')
        ax2.set_title(f'Water Depth Distribution - t = {current_time:.2f} s')
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim([-Lr, Lc-Lr])
        ax2.set_ylim([0, h0*1.2])
        
        ax2.plot(current_data['x_positions'], current_data['h_profile'], 'b-', linewidth=2)
        ax2.axvline(x=0, color='black', linewidth=2, alpha=0.5, label='Dam')
        ax2.axvline(x=current_data['x_front'], color='red', linestyle='--', alpha=0.7, label='Wave front')
        ax2.axvline(x=current_data['x_tail'], color='blue', linestyle='--', alpha=0.7, label='Wave tail')
        if hd > 0 and current_data['x_shock'] < Lc-Lr:
            ax2.axvline(x=current_data['x_shock'], color='green', linestyle='--', alpha=0.7, label='Shock')
        ax2.legend(fontsize=8)
        
        # Subplot 3: Current velocity distribution
        ax3.set_xlabel('Position x (m)')
        ax3.set_ylabel('Velocity u (m/s)')
        ax3.set_title(f'Velocity Distribution - t = {current_time:.2f} s')
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim([-Lr, Lc-Lr])
        ax3.set_ylim([0, 2*Cup*1.1])
        
        ax3.plot(current_data['x_positions'], current_data['u_profile'], 'brown', linewidth=2)
        ax3.axvline(x=0, color='black', linewidth=2, alpha=0.5)
        ax3.axvline(x=current_data['x_front'], color='red', linestyle='--', alpha=0.7)
        ax3.axvline(x=current_data['x_tail'], color='blue', linestyle='--', alpha=0.7)
        if hd > 0 and current_data['x_shock'] < Lc-Lr:
            ax3.axvline(x=current_data['x_shock'], color='green', linestyle='--', alpha=0.7)
        
        # Subplot 4: Wave position history
        ax4.set_xlabel('Time t (s)')
        ax4.set_ylabel('Wave position x (m)')
        ax4.set_title('Wave Position Tracking')
        ax4.grid(True, alpha=0.3)
        ax4.set_xlim([0, T])
        ax4.set_ylim([-Lr, Lc-Lr])
        
        # Plot wave position history up to current time
        x_front_hist = [data['x_front'] for data in all_data[:frame+1]]
        x_tail_hist = [data['x_tail'] for data in all_data[:frame+1]]
        
        ax4.plot(t_history, x_front_hist, 'r-', linewidth=2, label='Positive wave front')
        ax4.plot(t_history, x_tail_hist, 'b-', linewidth=2, label='Negative wave tail')
        
        if hd > 0:
            x_shock_hist = [data['x_shock'] for data in all_data[:frame+1]]
            ax4.plot(t_history, x_shock_hist, 'g-', linewidth=2, label='Shock')
        
        # Current time points
        ax4.plot(current_time, current_data['x_front'], 'ro', markersize=8)
        ax4.plot(current_time, current_data['x_tail'], 'bo', markersize=8)
        if hd > 0:
            ax4.plot(current_time, current_data['x_shock'], 'go', markersize=8)
        
        ax4.axvline(x=current_time, color='orange', linestyle='--', alpha=0.7)
        ax4.legend(fontsize=8)
        
        plt.tight_layout()
    
    # Create animation
    anim = FuncAnimation(fig, animate, frames=n_frames, interval=100, blit=False, repeat=True)
    
    if save_animation:
        print("Saving animation...")
        anim.save('dam_break_characteristics.gif', writer='pillow', fps=10)
        print("Animation saved as dam_break_characteristics.gif")
    
    plt.show()
    
    return anim

def compare_dry_wet_scenarios():
    """
    Compare differences between dry bed and wet bed dam-break scenarios
    """
    print("Comparative Analysis: Dry Bed vs Wet Bed Dam-break")
    
    # Common parameters
    h0 = 2.0
    Lc = 200
    Lr = 100
    T = 15
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    scenarios = [
        {'hd': 0.0, 'name': 'Dry Bed Dam-break', 'color': 'blue'},
        {'hd': 0.5, 'name': 'Wet Bed Dam-break', 'color': 'red'}
    ]
    
    g = 9.81
    Cup = np.sqrt(g * h0)
    
    for i, scenario in enumerate(scenarios):
        hd = scenario['hd']
        color = scenario['color']
        name = scenario['name']
        
        # Characteristic diagram
        ax = axes[i, 0]
        t = np.linspace(0.01, T, 1000)
        
        # Positive wave front
        x_front = 2 * Cup * t
        ax.plot(x_front, t, color=color, linewidth=2, label='Positive wave front')
        
        # Negative wave tail
        x_tail = -Cup * t
        ax.plot(x_tail, t, color=color, linewidth=2, linestyle='--', label='Negative wave tail')
        
        if hd > 0:
            # Shock (simplified calculation)
            Cdown = np.sqrt(g * hd)
            CB_approx = np.sqrt(g * (h0 + hd) / 2)
            x_shock = CB_approx * t
            ax.plot(x_shock, t, color='green', linewidth=2, linestyle=':', label='Shock')
        
        ax.set_xlim([-Lr, Lc-Lr])
        ax.set_ylim([0, T])
        ax.set_xlabel('Position x (m)')
        ax.set_ylabel('Time t (s)')
        ax.set_title(f'{name} - Characteristic Lines')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.axvline(x=0, color='black', linewidth=2, alpha=0.5)
        
        # Final water depth distribution
        ax = axes[i, 1]
        x_pos = np.linspace(-Lr, Lc-Lr, 500)
        h_final = np.zeros_like(x_pos)
        
        t_final = T
        x_front_final = 2 * Cup * t_final
        x_tail_final = -Cup * t_final
        
        if hd > 0:
            CB_final = np.sqrt(g * (h0 + hd) / 2)  # Simplified
            x_shock_final = CB_final * t_final
            hA = hd * 0.5  # Simplified
        
        for j, x in enumerate(x_pos):
            if x < x_tail_final:
                h_final[j] = h0
            elif hd == 0:
                if x <= x_front_final:
                    h_final[j] = max(0, (2 * Cup - x / t_final)**2 / (9 * g))
                else:
                    h_final[j] = 0
            else:
                if x <= x_shock_final:
                    h_calc = (2 * Cup - x / t_final)**2 / (9 * g)
                    h_final[j] = max(hA, h_calc)
                else:
                    h_final[j] = hd
        
        ax.plot(x_pos, h_final, color=color, linewidth=2)
        ax.set_xlim([-Lr, Lc-Lr])
        ax.set_ylim([0, h0*1.2])
        ax.set_xlabel('Position x (m)')
        ax.set_ylabel('Water depth h (m)')
        ax.set_title(f'{name} - Final Depth (t={T}s)')
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0, color='black', linewidth=2, alpha=0.5)
        
        # Wave speed comparison
        ax = axes[i, 2]
        wave_speeds = ['Positive\nwave front', 'Negative\nwave tail']
        speeds = [2*Cup, -Cup]
        
        if hd > 0:
            wave_speeds.append('Shock')
            speeds.append(CB_approx)
        
        bars = ax.bar(wave_speeds, speeds, color=[color, color, 'green'][:len(speeds)])
        ax.set_ylabel('Wave speed (m/s)')
        ax.set_title(f'{name} - Wave Speed Comparison')
        ax.grid(True, alpha=0.3)
        
        # Annotate values on bars
        for bar, speed in zip(bars, speeds):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1*abs(height),
                   f'{speed:.2f}', ha='center', va='bottom' if height > 0 else 'top')
    
    plt.tight_layout()
    plt.suptitle('Dry Bed vs Wet Bed Dam-break Comparative Analysis', fontsize=16, y=0.98)
    plt.show()
    
    return fig

def theoretical_analysis(h0, hd, Lc, Lr):
    """
    Theoretical analysis and key time calculations
    """
    g = 9.81
    Cup = np.sqrt(g * h0)
    
    print("="*60)
    print("Dam-break Wave Propagation Theoretical Analysis")
    print("="*60)
    
    print(f"Initial Conditions:")
    print(f"  Upstream depth h₀ = {h0:.2f} m")
    print(f"  Downstream depth h₁ = {hd:.2f} m")
    print(f"  Channel length Lc = {Lc:.2f} m")
    print(f"  Reservoir length Lr = {Lr:.2f} m")
    
    print(f"\nWave Speed Calculations:")
    print(f"  Upstream shallow water speed c₀ = √(gh₀) = {Cup:.3f} m/s")
    print(f"  Positive wave front speed = 2c₀ = {2*Cup:.3f} m/s")
    print(f"  Negative wave tail speed = -c₀ = {-Cup:.3f} m/s")
    
    if hd > 0:
        Cdown = np.sqrt(g * hd)
        print(f"  Downstream shallow water speed c₁ = √(gh₁) = {Cdown:.3f} m/s")
        print(f"  Dimensionless parameter h₁/h₀ = {hd/h0:.3f}")
        print(f"  Froude number Fr = √(h₀/h₁) = {np.sqrt(h0/hd):.3f}")
    
    print(f"\nKey Times:")
    t_front_downstream = (Lc - Lr) / (2 * Cup)
    t_tail_upstream = Lr / Cup
    
    print(f"  Positive wave front reaches downstream boundary: t = {t_front_downstream:.2f} s")
    print(f"  Negative wave tail reaches upstream boundary: t = {t_tail_upstream:.2f} s")
    
    if hd > 0:
        # Approximate shock velocity calculation
        CB_approx = np.sqrt(g * (h0 + hd) / 2)
        t_shock_downstream = (Lc - Lr) / CB_approx
        print(f"  Shock reaches downstream boundary (estimate): t ≈ {t_shock_downstream:.2f} s")
    
    print(f"\nRarefaction Wave Parameters:")
    print(f"  Maximum rarefaction wave width: {3*Cup*t_front_downstream:.2f} m")
    print(f"  Rarefaction wave travel time (to downstream): {t_front_downstream:.2f} s")
    
    if hd == 0:
        print(f"  Dry bed case: dry region behind rarefaction wave")
    else:
        print(f"  Wet bed case: rarefaction wave interacts with shock")
    
    return {
        'Cup': Cup,
        't_front_downstream': t_front_downstream,
        't_tail_upstream': t_tail_upstream,
        't_shock_downstream': t_shock_downstream if hd > 0 else None
    }

# Example usage
if __name__ == "__main__":
    print("Choose execution mode:")
    print("1. Interactive mode (user input parameters)")
    print("2. Batch mode (predefined parameters)")
    print("3. Characteristic analysis (static plot)")
    print("4. Animation mode (characteristics + wave propagation)")
    print("5. Comparative analysis (dry vs wet bed)")
    print("6. Theoretical analysis")
    
    mode = input("Enter mode (1-6): ")
    
    if mode == "1":
        # Interactive mode - equivalent to original MATLAB code
        x, h, u, final_time = stoker_dam_break()
    elif mode == "2":
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
    elif mode == "3":
        # Characteristic lines analysis
        print("Characteristic Analysis Mode")
        print("Please enter parameters:")
        h0 = float(input('Initial upstream depth h0 (m) [default: 2.0]: ') or 2.0)
        hd = float(input('Initial downstream depth hd (m) [default: 0.1]: ') or 0.1)
        Lc = float(input('Channel length Lc (m) [default: 200]: ') or 200)
        Lr = float(input('Reservoir length Lr (m) [default: 100]: ') or 100)
        T = float(input('Total time T (s) [default: 10]: ') or 10)
        
        print(f"Generating characteristic lines plot...")
        print(f"Parameters: h0={h0}m, hd={hd}m, Lc={Lc}m, Lr={Lr}m, T={T}s")
        
        fig = plot_characteristics(h0, hd, Lc, Lr, T)
        
        # Calculate and display key wave speed information
        g = 9.81
        Cup = np.sqrt(g * h0)
        print(f"\nKey Parameters:")
        print(f"Upstream wave speed c₀ = √(gh₀) = {Cup:.3f} m/s")
        print(f"Positive wave front speed = 2c₀ = {2*Cup:.3f} m/s")
        print(f"Negative wave tail speed = -c₀ = {-Cup:.3f} m/s")
        
        if hd > 0:
            Cdown = np.sqrt(g * hd)
            print(f"Downstream wave speed c₁ = √(gh₁) = {Cdown:.3f} m/s")
            print(f"Froude number Fr = √(h₀/h₁) = {np.sqrt(h0/hd):.3f}")
    
    elif mode == "4":
        # Animation mode
        print("Animation Mode - Synchronized display of characteristics and wave propagation")
        print("Please enter parameters:")
        h0 = float(input('Initial upstream depth h0 (m) [default: 2.0]: ') or 2.0)
        hd = float(input('Initial downstream depth hd (m) [default: 0.1]: ') or 0.1)
        Lc = float(input('Channel length Lc (m) [default: 200]: ') or 200)
        Lr = float(input('Reservoir length Lr (m) [default: 100]: ') or 100)
        T = float(input('Total time T (s) [default: 15]: ') or 15)
        
        save_gif = input('Save as GIF animation? (y/n) [default: n]: ').lower().startswith('y')
        
        print(f"Generating animation...")
        print(f"Parameters: h0={h0}m, hd={hd}m, Lc={Lc}m, Lr={Lr}m, T={T}s")
        print("Animation includes:")
        print("- Top left: Characteristic lines (x-t plane)")
        print("- Top right: Current water depth distribution")
        print("- Bottom left: Current velocity distribution") 
        print("- Bottom right: Wave position tracking")
        
        if save_gif:
            print("Note: Saving GIF may take considerable time...")
        
        anim = animate_dam_break_with_characteristics(h0, hd, Lc, Lr, T, save_gif)
        
        # Display theoretical information
        g = 9.81
        Cup = np.sqrt(g * h0)
        print(f"\nTheoretical Analysis:")
        print(f"Positive wave front reaches downstream boundary time: t ≈ {(Lc-Lr)/(2*Cup):.2f} s")
        print(f"Negative wave tail reaches upstream boundary time: t ≈ {Lr/Cup:.2f} s")
        
        if hd > 0:
            # Estimate shock arrival time at downstream boundary
            CB_approx = np.sqrt(g * (h0 + hd) / 2)  # Rough estimate
            print(f"Shock reaches downstream boundary time: t ≈ {(Lc-Lr)/CB_approx:.2f} s")
    
    elif mode == "5":
        # Comparative analysis
        print("Comparative Analysis Mode - Dry vs Wet Bed Dam-break")
        fig = compare_dry_wet_scenarios()
        
    elif mode == "6":
        # Theoretical analysis
        print("Theoretical Analysis Mode")
        print("Please enter parameters:")
        h0 = float(input('Initial upstream depth h0 (m) [default: 2.0]: ') or 2.0)
        hd = float(input('Initial downstream depth hd (m) [default: 0.1]: ') or 0.1)
        Lc = float(input('Channel length Lc (m) [default: 200]: ') or 200)
        Lr = float(input('Reservoir length Lr (m) [default: 100]: ') or 100)
        
        results = theoretical_analysis(h0, hd, Lc, Lr)
    
    else:
        print("Invalid selection, please restart the program")