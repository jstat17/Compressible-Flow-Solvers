## imports
from time import perf_counter_ns
from math import log

## convergence conditions
ε = 1e-8

## air conditions
R = 287.1 # [J/kg·K]
k = 1.4

## General

def sound_speed(T):
    return (k * R * T)**0.5 # [m/s]

def get_mach_number(v, T):
    return v / sound_speed(T) # M

def central_diff_derivative(f, x):
    h = ε * 1e5
    return (f(x+h) - f(x-h)) / (2*h)

## Variable Area Flow

def stagn_temp_ratio(M): # (4.19)
    return 1 + (k-1)/2*M**2 # T0 / T

def stagn_press_ratio(M): # (4.21)
    return stagn_temp_ratio(M)**(k/(k-1)) # P0 / P

def stagn_dens_ratio(M): # (4.22)
    return stagn_temp_ratio(M)**(1/(k-1)) # ρ0 / ρ

def mass_flow_rate(M, P_0, T_0, A): # (4.36)
    return P_0*(k/(R*T_0))**0.5 * A * M * (1 + (k-1)/2*M**2)**( (k+1)/(2-2*k) ) # ṁ [kg/s]

def expansion_ratio(M): # (4.38)
    return 1/M * ( ((k+1)/2) / (stagn_temp_ratio(M)) ) ** ((k+1)/(2-2*k)) # A / A*

def get_mach_from_expansion_ratio(A_Astar, Mi=1.5):
    # TODO:
    # Use bisection search for subsonic inlet Mach number due to odd behavior of the Newton-Raphson method (jumping over to
    # the supersonic regime because of the high gradient changes in 0 < M < 1), and because it is bounded, so bisection search
    # would be the fastest method. Keep the current solver for supersonic Mi as the region is unbounded (M > 1).
    objective_function = lambda M, A_Astar=A_Astar: A_Astar - expansion_ratio(M)
    begin = True
    Mi_next = Mi
    while begin or abs(Mi - Mi_next) > ε:
        begin = False
        Mi = Mi_next
        Mi_next = Mi - objective_function(Mi)/(central_diff_derivative(objective_function, Mi))
        # print(Mi_next)
        
    return Mi_next

## Normal Shock Waves

def shock_mach_number(M1): # (5.10)
    return ((M1**2 + 2/(k-1)) / (2*k/(k-1)*M1**2 - 1))**0.5 # M2

def shock_stagn_press_ratio(M1): # (5.13)
    return ( ((k+1)/2*M1**2) / ((k-1)/2*M1**2 + 1) )**(k/(k-1))  *  (shock_stat_press_ratio(M1))**(-1/(k-1)) # P02 / P01

def shock_stat_press_ratio(M1): # (5.14)
    return 2*k/(k+1)*M1**2 - (k-1)/(k+1) # P2 / P1

def shock_stat_temp_ratio(M1): # (5.15)
    return (1 + (k-1)/2*M1**2) * (2*k/(k-1)*M1**2 - 1) / ((k+1)**2/(2*(k-1))*M1**2) # T2 / T1

def shock_velocity_ratio(M1): # (5.16)
    return (1 + (k-1)/2*M1**2) / ((k+1)/2*M1**2) # v2 / v1

def get_shock_mach_from_velocity_ratio(v_ratio, Mi=2.5): # v2/v1 ratio input
    objective_function = lambda M1, v_ratio=v_ratio: v_ratio - shock_velocity_ratio(M1)
    begin = True
    Mi_next = Mi
    while begin or abs(Mi - Mi_next) > ε:
        begin = False
        Mi = Mi_next
        Mi_next = Mi - objective_function(Mi)/(central_diff_derivative(objective_function, Mi))
        # print(Mi_next)
        
    return Mi_next
    

def shock_stat_dens_ratio(M1): # (5.16)
    return 1/shock_velocity_ratio(M1) # ρ1 / ρ2

def reflected_shock_wave(Ms1, T1, v1=0):
    vs1 = Ms1 * sound_speed(T1)
    v1_prime = vs1 - v1
    M1_prime = get_mach_number(v1_prime, T1)
    
    # after shock jump
    T2 = shock_stat_temp_ratio(M1_prime) * T1
    v2_prime = shock_velocity_ratio(M1_prime) * v1_prime
    v2 = -(vs1 - v2_prime)
    
    # guessing-method for reflected shock wave
    v3 = 0
    a2 = sound_speed(T2)
    
    def get_delta_mach_from_guess(vs2):
        v2_dprime = vs2 - v2
        v3_dprime = vs2 - v3
        v_dprime_ratio = v3_dprime / v2_dprime
        v_dprime_ratio = 0.2 if v_dprime_ratio < 0.2 else v_dprime_ratio
        v_dprime_ratio = 1 if v_dprime_ratio > 1 else v_dprime_ratio
        
        M2_dprime_calc = get_mach_number(v2_dprime, T2)
        M2_dprime_ratio = get_shock_mach_from_velocity_ratio(v_dprime_ratio, Mi=M2_dprime_calc)
        
        return M2_dprime_ratio - M2_dprime_calc
    
    lower_v = 1.01 * a2 + v2
    lower_Δ = get_delta_mach_from_guess(lower_v)
    upper_v = Ms1 * a2 + v2
    upper_Δ = get_delta_mach_from_guess(upper_v)
    
    prev_central_Δ = (upper_Δ + lower_Δ)/2
    while True:
        central_v = (upper_v + lower_v)/2
        central_Δ = get_delta_mach_from_guess(central_v)
        
        if lower_Δ < 0 < central_Δ or lower_Δ > 0 > central_Δ:
            upper_v = central_v
            upper_Δ = central_Δ
        elif central_Δ < 0 < upper_Δ or central_Δ > 0 > upper_Δ:
            lower_v = central_v
            lower_Δ = central_Δ
        
        curr_central_Δ = central_Δ
        if abs(prev_central_Δ - curr_central_Δ) < ε:
            break
        else:
            prev_central_Δ = curr_central_Δ
    
    vs2 = -((upper_v - lower_v) / (upper_Δ - lower_Δ) * upper_Δ - upper_v)
    dict = {'vs2': vs2, 'Ms2': get_mach_number(vs2, T2)}
    return dict
       
def find_shock_position(Mi, n, P01, Pb): # n*Ae = Ai
    # to get the size of the duct at a percentage p along the duct, relative to the inlet area
    q = lambda p, n=n: 1 - (1 - 1/n)*p
    
    # rational function to approximate the supersonic mach number from an expansion ratio (valid for A_Astar ∈ [1, 25])
    approx_A_Astar_to_M = lambda A_Astar: (0.29893446*A_Astar**3 + 2.73075312*A_Astar**2 - 3.87408494*A_Astar + 0.94654509) / (0.0415028*A_Astar**3 + 0.99318356*A_Astar**2 - 0.40238173*A_Astar - 0.53364279)
    
    Ai_Astar1 = expansion_ratio(Mi)
    
    def shock_pos_to_outlet_exp_ratio_approx(p):
        As_Astar1 = q(p) * Ai_Astar1
        M1 = approx_A_Astar_to_M(As_Astar1)
        M2 = shock_mach_number(M1)
        As_Astar2 = expansion_ratio(M2)
        Ae_Astar2 = As_Astar2 / (n*q(p))
        
        return Ae_Astar2
    
    def shock_pos_to_outlet_conditions(p):
        As_Astar1 = q(p) * Ai_Astar1
        approx = approx_A_Astar_to_M(As_Astar1)
  
        M1 = get_mach_from_expansion_ratio(As_Astar1, approx)
        M2 = shock_mach_number(M1)
        P02_P01 = shock_stagn_press_ratio(M1)
        As_Astar2 = expansion_ratio(M2)
        Ae_Astar2 = As_Astar2 / (n*q(p))
        
        if n < 1:
            approx = 0.25
        else:
            approx = M2
        Me = get_mach_from_expansion_ratio(Ae_Astar2, approx)
        Pe = P02_P01 * P01 / stagn_press_ratio(Me)
        Pe_Δ = Pb - Pe
        
        return {'Me': Me,
                'Pe_Δ': Pe_Δ}
        
    
    if n > 1:
        p_lower, p_upper = 0, 1
        outlet_exp_ratio_lower = shock_pos_to_outlet_exp_ratio_approx(p_lower)
        outlet_exp_ratio_upper = shock_pos_to_outlet_exp_ratio_approx(p_upper)
        
        p_mid = (p_upper + p_lower) / 2
        while True:
            outlet_exp_ratio_mid = shock_pos_to_outlet_exp_ratio_approx(p_mid)
            
            if outlet_exp_ratio_lower < 1 < outlet_exp_ratio_mid:
                p_upper = p_mid
                outlet_exp_ratio_upper = outlet_exp_ratio_mid
                
            elif outlet_exp_ratio_mid < 1 < outlet_exp_ratio_upper:
                p_lower = p_mid
                outlet_exp_ratio_lower = outlet_exp_ratio_mid
                
            p_prev = p_mid
            p_mid = (p_upper + p_lower) / 2
            if abs(p_mid - p_prev) < 1e-6:
                p_lower = p_upper
                break
            
    else:
        p_lower = 0
    
    p_upper = 1
    print(f"bounds: [{p_lower:.4f}, {p_upper:.4f}]") # where the solver will look for the shockwave
    Pe_Δ_upper = shock_pos_to_outlet_conditions(p_upper)['Pe_Δ']
    Pe_Δ_lower = shock_pos_to_outlet_conditions(p_lower)['Pe_Δ']
    p_mid = (p_upper + p_lower) / 2
    Pe_Δ_mid, Pe_Δ_mid_prev = 0, 0
    begin = True
    
    while begin or abs(Pe_Δ_mid - Pe_Δ_mid_prev) > 1e-1:
        begin = False
        Pe_Δ_mid_prev = Pe_Δ_mid
        Pe_Δ_mid = shock_pos_to_outlet_conditions(p_mid)['Pe_Δ']
        
        if Pe_Δ_lower < 0 < Pe_Δ_mid or Pe_Δ_lower > 0 > Pe_Δ_mid:
            p_upper = p_mid
            Pe_Δ_upper = Pe_Δ_mid
            
        elif Pe_Δ_mid < 0 < Pe_Δ_upper or Pe_Δ_mid > 0 > Pe_Δ_upper:
            p_lower = p_mid
            Pe_Δ_lower = Pe_Δ_mid
            
        p_mid = (p_upper + p_lower) / 2
        
    p = -((p_lower - p_upper) / (Pe_Δ_lower - Pe_Δ_upper) * Pe_Δ_lower - p_lower)
    dict_exit_cond = shock_pos_to_outlet_conditions(p)
    dict_exit_cond['p'] = p
    
    return dict_exit_cond



## Fanno Flow (adiabatic flow with friction)

def get_hydraulic_diameter(shape_dict): # (7.3)
    '''
    Returns the hydraulic diameter (Dh) for a cross-section.
    
    Parameters:
    shape_dict (dict): Gives the dimensions of the input shape
        shape_dict['shape']: an implemented shape in this function
        shape_dict['side']: the side length of shape 'square'
        shape_dict['diameter']: the diameter of shape 'circle'
        shape_dict['side1']: the length of side 1 for shape 'rectangle'
        shape_dict['side2']: the length of side 2 for shape 'rectangle'
        
    Returns:
    float: the hydraulic diameter of the input shape
    '''

    try:
        shape = shape_dict['shape']
        if shape == 'square':
            return shape_dict['side']
        elif shape == 'circle':
                return shape_dict['diameter']
        elif shape == 'rectangle':
            A = shape_dict['side1']
            B = shape_dict['side2']
            return 2*A*B/(A+B)
        else:
            print(f"\n\nEXCEPTION:\n\nShape '{shape_dict['shape']}' not implemented in the get_hydraulic_diameter function.\nThe program will now terminate.")
            exit()
    except KeyError:
        print("\n\nEXCEPTION:\n\nYou are missing a key in the input dictionary for the selected shape in the 'get_hydraulic_diameter' function.\nThe program will now terminate.")
        exit()
        
 # For Python 3.10 and above
    # try:
    #     match shape_dict['shape']:
    #         case 'square':
    #             return shape_dict['side']
    #         case 'circle':
    #             return shape_dict['diameter']
    #         case 'rectangle':
    #             A = shape_dict['side1']
    #             B = shape_dict['side2']
    #             return 2*A*B/(A+B)
    #         case _:
    #             print(f"\n\nEXCEPTION:\n\nShape '{shape_dict['shape']}' not implemented in the get_hydraulic_diameter function.\nThe program will now terminate.")
    #             exit()
    # except KeyError:
    #     print("\n\nEXCEPTION:\n\nYou are missing a key in the input dictionary for the selected shape in the 'get_hydraulic_diameter' function.\nThe program will now terminate.")
    #     exit()
    
def fanno_sonic_stagn_press_ratio(M): # (7.25)
    return 1/M * ( (1 + (k-1)/2 * M**2) / ((k+1)/2) ) ** ((k+1) / (2*k - 2)) # P0 / P0star

def fanno_sonic_press_ratio(M): # (7.22)
    return 1/M * fanno_sonic_temp_ratio(M) ** 0.5 # P / Pstar

def fanno_sonic_temp_ratio(M): # (7.21)
    return ((k+1)/2) / (1 + (k-1)/2 * M**2) # T / Tstar

def fanno_length_parameter(M): # (7.28)
    return -1/k - (k+1)/(2*k) * log(1/(1 + (k-1)/2)) + 1/(k*M**2) + (k+1)/(2*k) * log(M**2/(1 + (k-1)/2*M**2))
    
def get_mach_from_fanno_length_parameter(fanno_param, Mi):
    if Mi < 1:
        M_lower = 0.001
        M_upper = 1
    
    elif Mi > 1:
        M_lower = 1
        M_upper = 10
        
    else: # Mi = 1 in this case
        return 1
    
    M_mid = (M_lower + M_upper) / 2
    fanno_M_lower = fanno_length_parameter(M_lower)
    fanno_M_upper = fanno_length_parameter(M_upper)
    
    while True:
        fanno_M_mid = fanno_length_parameter(M_mid)
        
        if fanno_M_lower < fanno_param < fanno_M_mid or fanno_M_lower > fanno_param > fanno_M_mid:
            fanno_M_upper = fanno_M_mid
            M_upper = M_mid
            
        elif fanno_M_mid < fanno_param < fanno_M_upper or fanno_M_mid > fanno_param > fanno_M_upper:
            fanno_M_lower = fanno_M_mid
            M_lower = M_mid
            
        M_prev = M_mid
        M_mid = (M_lower + M_upper) / 2
        
        if abs(M_mid - M_prev) < ε:
            return M_mid
    
def solve_fanno_duct(Mi, f, Dh, L):
    if Mi > 1:
        fanno_param_inlet = fanno_length_parameter(Mi)
        Lstar = fanno_param_inlet * Dh / 4 / f
        
        if L > Lstar: # the duct chokes + supersonic inlet, so there will be a shock wave
            L_lower = 0
            L_upper = Lstar
            
            def get_fanno_exit_Δ_from_shock_pos_guess(L1):
                fanno_param_upstream_of_shock = fanno_param_inlet - 4*f*L1/Dh
                M1 = get_mach_from_fanno_length_parameter(fanno_param_upstream_of_shock, 1.5)
                M2 = shock_mach_number(M1)
                fanno_param_downstream_of_shock = fanno_length_parameter(M2)
                
                return fanno_param_downstream_of_shock - 4*f*(L-L1)/Dh
            
            # fanno_exit_Δ_lower = get_fanno_exit_Δ_from_shock_pos_guess(L_lower)
            # fanno_exit_Δ_upper = get_fanno_exit_Δ_from_shock_pos_guess(L_upper)
            L_mid = (L_lower + L_upper) / 2
            fanno_exit_Δ_mid = 0
            while True:
                fanno_prev = fanno_exit_Δ_mid
                fanno_exit_Δ_mid = get_fanno_exit_Δ_from_shock_pos_guess(L_mid)
                
                if fanno_exit_Δ_mid < 0:
                    # fanno_exit_Δ_upper = fanno_exit_Δ_mid
                    L_upper = L_mid
                    
                elif fanno_exit_Δ_mid > 0:
                    # fanno_exit_Δ_lower = fanno_exit_Δ_mid
                    L_lower = L_mid
                    
                L_mid = (L_lower + L_upper) / 2
                
                if abs(fanno_exit_Δ_mid - fanno_prev) < ε:
                    Ls = L_mid # L shock
                    break
                
            print(f"Distance from inlet to shockwave: {Ls} m")
            P_ratio_list = []
            P_ratio_list.append(fanno_sonic_press_ratio(Mi)) # Pi / Pstar,1
            M1 = get_mach_from_fanno_length_parameter( P_ratio_list[0] - 4*f*Ls/Dh , 1.5)
            P_ratio_list.append(1/fanno_sonic_press_ratio(M1)) # Pstar,1 / P1
            P_ratio_list.append(1/shock_stat_press_ratio(M1)) # P1 / P2
            M2 = shock_mach_number(M1)
            P_ratio_list.append(fanno_sonic_press_ratio(M2)) # P2 / Pstar,2 = P2 / Pe
            
            T_ratio_list = []
            T_ratio_list.append(fanno_sonic_temp_ratio(Mi)) # Ti / Tstar,1
            T_ratio_list.append(1/fanno_sonic_temp_ratio(M1)) # Tstar,1 / T1
            T_ratio_list.append(1/shock_stat_temp_ratio(M1)) # T1 / T2
            T_ratio_list.append(fanno_sonic_temp_ratio(M2)) # T2 / Tstar,2 = T2 / Te
            
            total_stat_press_ratio, total_stat_temp_ratio = 1, 1
            for ratio_press, ratio_temp in zip(P_ratio_list, T_ratio_list):
                total_stat_press_ratio *= ratio_press
                total_stat_temp_ratio *= ratio_temp
            
            return {'Ls': Ls, 'Pi/Pe': total_stat_press_ratio, 'Ti/Te': total_stat_temp_ratio}
    
if __name__ == "__main__":
    tic = perf_counter_ns()
    # print(reflected_shock_wave(1.5, 300))
    # print(find_shock_position(1.8, 1/3, 100e3, 60e3))
    # print(find_shock_position(2.8, 3, 100e3, 60e3))
    print(solve_fanno_duct(5, 0.0035, 12.7e-3, 1))
    toc = perf_counter_ns()
    print(f"time: {(toc - tic)/1e6} ms")