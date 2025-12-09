import numpy as np
from config import SIMULATION_CONFIG, TURING_PARAMS, HYPERBOLIC_PARAMS

class StructureEngine:
    def __init__(self):
        self.size = SIMULATION_CONFIG["GRID_SIZE"]
        self.steps = SIMULATION_CONFIG["TURING_STEPS"]
        self.du = TURING_PARAMS["Du"]
        self.dv = TURING_PARAMS["Dv"]
        self.f = TURING_PARAMS["f"]
        self.k = TURING_PARAMS["k"]
        self.radius = HYPERBOLIC_PARAMS["POINCARE_RADIUS"]

    def _laplacian(self, Z):
        return (
            np.roll(Z, 1, axis=0) + np.roll(Z, -1, axis=0) +
            np.roll(Z, 1, axis=1) + np.roll(Z, -1, axis=1) -
            4 * Z
        )

    def generate_turing_pattern(self):
        U = np.ones((self.size, self.size))
        V = np.zeros((self.size, self.size))
        
        r = self.size // 2
        U[r-2:r+2, r-2:r+2] = 0.50
        V[r-2:r+2, r-2:r+2] = 0.25
        
        U += np.random.normal(0, 0.05, (self.size, self.size))
        V += np.random.normal(0, 0.05, (self.size, self.size))

        for _ in range(self.steps):
            Lu = self._laplacian(U)
            Lv = self._laplacian(V)
            
            uvv = U * V * V
            U += (self.du * Lu - uvv + self.f * (1 - U))
            V += (self.dv * Lv + uvv - (self.f + self.k) * V)
            
        return V

    def poincare_transform(self, x, y):
        norm_x = (x / self.size) * 2 - 1
        norm_y = (y / self.size) * 2 - 1
        
        r_euclid = np.sqrt(norm_x**2 + norm_y**2)
        theta = np.arctan2(norm_y, norm_x)
        
        r_hyp = np.tanh(r_euclid) * self.radius
        
        px = r_hyp * np.cos(theta)
        py = r_hyp * np.sin(theta)
        
        return px, py

    def build_hyperbolic_graph(self):
        pattern = self.generate_turing_pattern()
        threshold = np.mean(pattern)
        
        nodes = []
        edges = []
        node_indices = {}
        
        count = 0
        for i in range(self.size):
            for j in range(self.size):
                if pattern[i, j] > threshold:
                    hx, hy = self.poincare_transform(i, j)
                    nodes.append({
                        "id": count,
                        "grid_pos": (i, j),
                        "hyp_pos": (hx, hy),
                        "concentration": float(pattern[i, j])
                    })
                    node_indices[(i, j)] = count
                    count += 1
        
        for n in nodes:
            gx, gy = n["grid_pos"]
            neighbors = [
                (gx+1, gy), (gx-1, gy), (gx, gy+1), (gx, gy-1)
            ]
            
            for nx, ny in neighbors:
                if (nx, ny) in node_indices:
                    target_id = node_indices[(nx, ny)]
                    if n["id"] < target_id:
                        u_pos = np.array(n["hyp_pos"])
                        v_pos = np.array(nodes[target_id]["hyp_pos"])
                        
                        dist_sq = np.sum((u_pos - v_pos)**2)
                        weight = np.exp(-dist_sq)
                        
                        edges.append((n["id"], target_id, float(weight)))
                        
        return {"nodes": nodes, "edges": edges}