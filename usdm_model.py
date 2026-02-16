import numpy as np
import matplotlib.pyplot as plt

# TASK 2 
def user_mod_prototype(material_props, sigo, deps):
    E, nu, c, phi_deg = material_props
    phi_rad = np.radians(phi_deg)
    
    # Matriz de Rigidez Elástica (D) 
    G = E / (2 * (1 + nu))
    fac = (2 * G) / (1 - 2 * nu)
    term1 = fac * (1 - nu)
    term2 = fac * nu
    
    D = np.array([
        [term1, term2, term2, 0, 0, 0],
        [term2, term1, term2, 0, 0, 0],
        [term2, term2, term1, 0, 0, 0],
        [0, 0, 0, G, 0, 0],
        [0, 0, 0, 0, G, 0],
        [0, 0, 0, 0, 0, G]
    ])

    # Tensões
    sig_trial = sigo + np.dot(D, deps)
    
    # Atualização das tensões
    p_trial = -(sig_trial[0] + sig_trial[1] + sig_trial[2]) / 3.0
    s = np.zeros(6)
    s[0:3] = -sig_trial[0:3] - p_trial
    s[3:6] = -sig_trial[3:6]
    
    # Critério de resistência
    q_trial = np.sqrt(0.5 * ((s[0]-s[1])**2 + (s[1]-s[2])**2 + (s[2]-s[0])**2 + 6*(s[3]**2 + s[4]**2 + s[5]**2)))
    M = (6 * np.sin(phi_rad)) / (3 - np.sin(phi_rad))
    q_max = M * p_trial + (6 * c * np.cos(phi_rad)) / (3 - np.sin(phi_rad))
    
    sig_final = np.copy(sig_trial)
    ipl = 0 #Verificador de plastificação (ipl=1 -> plastificou)

    # Verificação da plastificação
    if p_trial < -1e-6:
        sig_final = np.zeros(6) 
        ipl = 1
    elif q_trial > q_max and q_trial > 1e-6:
        scaling = q_max / q_trial
        # Devolve para a convenção PLAXIS (Negativa)
        sig_final[0:3] = -(p_trial + s[0:3] * scaling)
        sig_final[3:6] = -(s[3:6] * scaling)
        ipl = 1
            
    return sig_final, ipl

# Simulação triaxial
def run_triaxial_test(material_props, confinement=-100.0, max_strain=-0.05, steps=500):
    E, nu, c, phi = material_props
    
    # Estado inicial: Isotrópico (Confinamento)
    curr_sig = np.array([confinement, confinement, confinement, 0.0, 0.0, 0.0])
    curr_eps = np.zeros(6)
    
    # Incremento de deformação axial
    de_a = max_strain / steps
    
    history_q = []
    history_eps = []

    for i in range(steps):
        
        # Ajuste do incremento de deformações laterais
        de_l = -nu * de_a 
        deps = np.array([de_a, de_l, de_l, 0, 0, 0])
        
        curr_sig, ipl = user_mod_prototype(material_props, curr_sig, deps)
        
        # Ajuste do confinamento lateral
        curr_sig[1] = confinement
        curr_sig[2] = confinement
        
        # Cálculo de p e q para o gráfico
        q = abs(curr_sig[0] - curr_sig[2])
        curr_eps[0] += de_a
        
        history_q.append(q)
        history_eps.append(abs(curr_eps[0]) * 100) # Em %

    return history_eps, history_q


# Plot dos resultados
props_usp = [20000, 0.3, 5, 30] # E, nu, c, phi

eps_ax, q_dev = run_triaxial_test(props_usp, confinement=-100)

plt.figure(figsize=(10, 6))
plt.plot(eps_ax, q_dev, label='Protótipo UDSM (Python)', linewidth=2, color='darkblue')
plt.title('Validação de Protótipo: Ensaio Triaxial Drenado (Mohr-Coulomb)')
plt.xlabel('Deformação Axial εa [%]')
plt.ylabel('Tensão Desviadora q [kPa]')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.axhline(y=max(q_dev), color='r', linestyle=':', label=f'Resistência de Pico: {max(q_dev):.1f} kPa')
plt.legend()
plt.show()