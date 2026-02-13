import numpy as np
import matplotlib.pyplot as plt

def user_mod(id_task, props, sig0, deps, stvar0):
    """
    Simulação da sub-rotina User_Mod do PLAXIS
    id_task: tarefa atual (1 a 6) [cite: 32]
    props: array de parâmetros do material [cite: 32]
    sig0: tensões anteriores [cite: 32]
    deps: incrementos de deformação [cite: 40]
    stvar0: variáveis de estado iniciais [cite: 32]
    """
    
    # IDTask 4: Retornar número de variáveis de estado [cite: 66, 140]
    if id_task == 4:
        n_stat = 0 # Para um modelo elástico simples, começamos com zero
        return n_stat

        # IDTask 5: Atributos da matriz (Simetria e Dependência) [cite: 146]
    if id_task == 5:
        non_sym = 0   # 0 = Matriz Simétrica (comum em elasticidade) [cite: 155]
        i_strs_dep = 0 # 0 = Não depende da tensão (Linear) [cite: 158]
        i_time_dep = 0 # 0 = Não depende do tempo (Sem creep) [cite: 159]
        i_tang = 0     # 0 = Não é matriz tangente (Newton-Raphson simples) [cite: 160]
        
        return non_sym, i_strs_dep, i_time_dep, i_tang

    # IDTask 3: Criar matriz de rigidez efetiva D [cite: 104]
    if id_task == 3:
        # Extração dos parâmetros do array Props 
        # No manual (Fortran) os índices são 1 e 2. No Python (0-indexed) usamos 0 e 1.
        E = props[0]  # Módulo de Young [cite: 111]
        nu = props[1] # Coeficiente de Poisson [cite: 113]
        
        # Cálculo das constantes elásticas conforme o manual [cite: 116, 117, 118]
        # Atenção: nu deve ser menor que 0.5 para evitar divisão por zero [cite: 117]
        G = 0.5 * E / (1.0 + nu) [cite: 116]
        fac = 2.0 * G / (1.0 - 2.0 * nu) [cite: 117]
        term1 = fac * (1.0 - nu) [cite: 117]
        term2 = fac * nu [cite: 118]
        
        # Inicialização da matriz 6x6 [cite: 39, 135]
        D = np.zeros((6, 6))
        
        # Preenchimento dos componentes da matriz (conforme pág. 6 do manual)
        # Relações entre tensões e deformações normais [cite: 119-130]
        D[0, 0] = term1 # Equivale ao D[1,1] do manual [cite: 119]
        D[0, 1] = term2 # D[1,2] [cite: 120]
        D[0, 2] = term2 # D[1,3] [cite: 121]
        
        D[1, 0] = term2 # D[2,1] [cite: 122]
        D[1, 1] = term1 # D[2,2] [cite: 123]
        D[1, 2] = term2 # D[2,3] [cite: 124]
        
        D[2, 0] = term2 # D[3,1] [cite: 126]
        D[2, 1] = term2 # D[3,2] [cite: 128]
        D[2, 2] = term1 # D[3,3] [cite: 130]
        
        # Componentes de cisalhamento (Módulo de Cisalhamento G) [cite: 131-133]
        D[3, 3] = G # D[4,4] [cite: 131]
        D[4, 4] = G # D[5,5] [cite: 132]
        D[5, 5] = G # D[6,6] [cite: 133]
        
        return D

    # IDTask 2: Calcular tensões constitutivas (Sig) [cite: 19, 92]
    if id_task == 2:
        # 1. Obter parâmetros de rigidez e resistência [cite: 107]
        E = props[0]
        nu = props[1]
        c = props[2]
        phi_rad = np.radians(props[3])
        
        # 2. Preditor Elástico (Cálculo da tensão "tentativa") 
        D = user_mod(3, props, sig0, deps, stvar0)
        sig_trial = np.array(sig0[:6]) + np.dot(D, deps[:6]) # [cite: 95, 101]
        
        # 3. Calcular Invariantes de Tensão do estado "trial"
        p_trial = (sig_trial[0] + sig_trial[1] + sig_trial[2]) / 3.0
        
        # Cálculo simplificado de q (von Mises) para o protótipo
        s_xx, s_yy, s_zz = sig_trial[0]-p_trial, sig_trial[1]-p_trial, sig_trial[2]-p_trial
        q_trial = np.sqrt(0.5*((s_xx-s_yy)**2 + (s_yy-s_zz)**2 + (s_zz-s_xx)**2 + 
                         6*(sig_trial[3]**2 + sig_trial[4]**2 + sig_trial[5]**2)))

        # 4. Critério de Mohr-Coulomb no espaço p-q (simplificado)
        # M é a inclinação da linha de estado crítico baseada em phi
        M = (6.0 * np.sin(phi_rad)) / (3.0 - np.sin(phi_rad))
        q_max = M * p_trial + (6.0 * c * np.cos(phi_rad)) / (3.0 - np.sin(phi_rad))
        
        # Função de falha f
        f = q_trial - q_max
        
        if f <= 0:
            # Caso Elástico: A tensão trial é válida [cite: 39]
            # ipl = 0 (sem plasticidade) [cite: 39]
            return sig_trial
        else:
            # Caso Plástico: Retorno de Tensão (Stress Return)
            # Reduzimos q_trial para q_max mantendo p constante (retorno radial simples)
            scaling_factor = q_max / q_trial
            
            # As tensões desviatórias são escalonadas para voltar ao envelope
            sig_corrected = np.zeros(6)
            sig_corrected[0] = p_trial + (s_xx * scaling_factor)
            sig_corrected[1] = p_trial + (s_yy * scaling_factor)
            sig_corrected[2] = p_trial + (s_zz * scaling_factor)
            sig_corrected[3] = sig_trial[3] * scaling_factor
            sig_corrected[4] = sig_trial[4] * scaling_factor
            sig_corrected[5] = sig_trial[5] * scaling_factor
            
            # No PLAXIS real, você marcaria ipl = 1 (Mohr-Coulomb failure) [cite: 39]
            return sig_corrected

        # IDTask 6: Criar matriz de rigidez elástica De
    if id_task == 6:
        # Para modelos elastoplásticos, aqui deve constar APENAS a parte elástica.
        # Como nosso modelo ainda é puramente elástico, chamamos a Task 3.
        # O manual nota que a variável 'D' é reutilizada para armazenar De aqui[cite: 188].
        De = user_mod(3, props, sig0, deps, stvar0)
        return De

# --- CONFIGURAÇÃO DO TESTE ---
# Propriedades: [E, nu, c, phi]
propriedades = [20000, 0.3, 10, 30] # E=20MPa, nu=0.3, c=10kPa, phi=30°

# Estado inicial: Solo confinado com 100kPa (Tensões são negativas em compressão no PLAXIS)
sig_atual = np.array([-100.0, -100.0, -100.0, 0.0, 0.0, 0.0])
stvar = np.zeros(0) # nStat = 0

# Incremento de deformação vertical (compressão de 1%)
# dEps: [eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx]
deps = np.array([0.0, -0.0001, 0.0, 0.0, 0.0, 0.0]) # Deformação em YY (vertical)

# Listas para armazenar resultados do gráfico
deformacoes = []
tensoes_deviatorias = []

# --- EXECUÇÃO DO PASSO A PASSO (SIMULANDO O PLAXIS) ---
for i in range(200):
    # Chama a Task 2 para calcular a nova tensão
    sig_novo = user_mod(2, propriedades, sig_atual, deps, stvar)
    
    # Invariante q (Tensão desviatória aproximada: |sig_yy - sig_xx|)
    q = abs(sig_novo[1] - sig_novo[0])
    
    # Armazena dados
    deformacoes.append(i * 0.0001 * 100) # em %
    tensoes_deviatorias.append(q)
    
    # Atualiza a tensão para o próximo incremento (como o PLAXIS faz) [cite: 8, 41]
    sig_atual = sig_novo

# --- VISUALIZAÇÃO ---
plt.figure(figsize=(8, 5))
plt.plot(deformacoes, tensoes_deviatorias, label='Modelo UDSM Python', color='blue', linewidth=2)
plt.axhline(y=tensoes_deviatorias[-1], color='red', linestyle='--', label='Limite de Ruptura (MC)')
plt.title('Simulação de Ensaio Triaxial - Protótipo Mohr-Coulomb')
plt.xlabel('Deformação Axial (%)')
plt.ylabel('Tensão Desviatória q (kPa)')
plt.grid(True)
plt.legend()
plt.show()