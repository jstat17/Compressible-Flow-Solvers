## imports
from time import perf_counter_ns

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
        M1 = get_mach_from_expansion_ratio(As_Astar1, approx_A_Astar_to_M(As_Astar1))
        M2 = shock_mach_number(M1)
        P02_P01 = shock_stagn_press_ratio(M1)
        As_Astar2 = expansion_ratio(M2)
        Ae_Astar2 = As_Astar2 / (n*q(p))
        
        Me = get_mach_from_expansion_ratio(Ae_Astar2, 0.5)
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
    print(f"bounds: [{p_lower:.4f}, {p_upper:.4f}]")
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

if __name__ == "__main__":
    tic = perf_counter_ns()
    # print(reflected_shock_wave(1.5, 300))
    print(find_shock_position(2.8, 3, 100e3, 60e3))
    toc = perf_counter_ns()
    print(f"time: {(toc - tic)/1e6} ms")