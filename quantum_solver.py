import numpy as np
from config import QUANTUM_PARAMS, SIMULATION_CONFIG

class QuantumSolver:
    def __init__(self):
        self.steps = SIMULATION_CONFIG["QUANTUM_STEPS"]
        self.decoherence = QUANTUM_PARAMS["DECOHERENCE_RATE"]
        self.tunneling = QUANTUM_PARAMS["TUNNELING_FACTOR"]
        self.dt = 0.1

    def _build_hamiltonian(self, num_nodes, edges):
        H = np.zeros((num_nodes, num_nodes), dtype=complex)
        
        for u, v, w in edges:
            if u < num_nodes and v < num_nodes:
                H[u, v] = -w
                H[v, u] = -w
                
        degrees = np.sum(np.abs(H), axis=1)
        for i in range(num_nodes):
            H[i, i] = degrees[i]
            
        return H

    def _apply_tunneling(self, H, nodes):
        num_nodes = H.shape[0]
        positions = np.array([n["hyp_pos"] for n in nodes])
        
        for i in range(num_nodes):
            dists = np.sum((positions - positions[i])**2, axis=1)
            neighbors = np.where(dists < 0.5)[0] 
            
            for j in neighbors:
                if i != j and H[i, j] == 0:
                    H[i, j] = -self.tunneling
                    H[j, i] = -self.tunneling
        return H

    def solve_wave_function(self, graph_data):
        nodes = graph_data["nodes"]
        edges = graph_data["edges"]
        num_nodes = len(nodes)
        
        if num_nodes == 0:
            return np.zeros(num_nodes)

        H = self._build_hamiltonian(num_nodes, edges)
        H = self._apply_tunneling(H, nodes)
        
        psi = np.zeros(num_nodes, dtype=complex)
        psi[0] = 1.0 + 0j
        
        psi /= np.linalg.norm(psi)
        
        history = []
        
        for _ in range(self.steps):
            d_psi = -1j * np.dot(H, psi) * self.dt
            psi += d_psi
            
            psi *= (1.0 - self.decoherence)
            
            norm = np.linalg.norm(psi)
            if norm > 0:
                psi /= norm
                
            history.append(np.abs(psi)**2)
            
        final_probability = np.mean(history, axis=0)
        return final_probability