import numpy as np
from config import MATERIALS_DB, ENV_CONDITIONS, THRESHOLDS, ACTIVATION_PHYSICS

class ChemicalValidator:
    def __init__(self):
        self.temp_k = ENV_CONDITIONS["TEMPERATURE_K"]
        self.max_gibbs = THRESHOLDS["MAX_GIBBS_ENERGY"]

    def calculate_gibbs_energy(self, recipe_ratios):
        total_H = sum(MATERIALS_DB[m]["Enthalpy_H"] * r for m, r in recipe_ratios.items())
        total_S = sum(MATERIALS_DB[m]["Entropy_S"] * r for m, r in recipe_ratios.items())
        return total_H - (self.temp_k * total_S)

    def calculate_shannon_entropy(self, probability_dist):
        dist = np.array(probability_dist)
        dist = dist[dist > 0]
        entropy = -np.sum(dist * np.log(dist))
        return entropy / np.log(len(probability_dist)) if len(dist) > 1 else 0

    def calculate_lyapunov_exponent(self, time_series):
        if len(time_series) < 2: return 0.0
        diffs = np.diff(time_series)
        diffs = diffs[diffs > 0]
        return np.mean(np.log(diffs)) if len(diffs) > 0 else -1.0

    def predict_performance(self, recipe, quantum_score, temp_c, acid_time_min):
        # 1. Surface Area (The Biochar Engine)
        optimal_temp_area = 580 # Refined optimal temp
        temp_factor = np.exp(-((temp_c - optimal_temp_area)**2) / 3000)
        
        acid_factor = min(1.0, (acid_time_min / 120.0)) * ACTIVATION_PHYSICS["SILICA_REMOVAL_EFFICIENCY"]
        acid_factor = min(1.3, acid_factor) 
        
        generated_surface_area = ACTIVATION_PHYSICS["BASE_SURFACE_AREA"] + \
                                 (ACTIVATION_PHYSICS["MAX_THEORETICAL_AREA"] * temp_factor * acid_factor * 0.55)

        # 2. Surface Chemistry (Biochar Pyrones)
        optimal_temp_chem = 520
        chem_factor = np.exp(-((temp_c - optimal_temp_chem)**2) / 7000)
        
        # 3. THE GELATIN AMINE BOOST (New Logic)
        # Gelatin contains Amino acids (Glycine, Proline) which have Amine groups.
        # Amine groups react with CO2 (Chemisorption).
        gelatin_ratio = recipe.get("GELATIN", 0)
        amine_efficiency_factor = 15.0 # High affinity constant
        gelatin_boost_mg_g = gelatin_ratio * 100 * amine_efficiency_factor
        
        # 4. Final Capacity Calculation
        micropore_efficiency = 1.0 + (quantum_score * 0.6)
        
        # Biochar Contribution
        biochar_capacity_mg_g = (generated_surface_area * 0.035) + \
                                (generated_surface_area * chem_factor * ACTIVATION_PHYSICS["SURFACE_CHEMISTRY_FACTOR"])
        biochar_capacity_mg_g *= micropore_efficiency

        # Composite Total = (Biochar Contribution) + (Gelatin Amine Contribution)
        mass_biochar = recipe.get("BIOCHAR", 0) * 100
        mass_matrix_remainder = recipe.get("STARCH", 0) * 100 # Starch is neutral/blocking
        
        # Starch blocks some pores, but Gelatin helps now
        starch_penalty = 1.0 - (recipe.get("STARCH", 0) * 0.2)
        
        total_uptake_mg = (mass_biochar * biochar_capacity_mg_g * starch_penalty) + \
                          (gelatin_boost_mg_g) # Adding the matrix boost directly
                          
        total_uptake_g = total_uptake_mg / 1000.0
        
        return round(total_uptake_g, 4), int(generated_surface_area)

    def validate_sample(self, recipe, quantum_prob, history, temp_c, acid_time):
        delta_G = self.calculate_gibbs_energy(recipe)
        entropy = self.calculate_shannon_entropy(quantum_prob)
        lyapunov = self.calculate_lyapunov_exponent(history)
        
        quantum_mean = np.mean(quantum_prob)
        
        predicted_grams, surface_area = self.predict_performance(recipe, quantum_mean, temp_c, acid_time)
        
        is_spontaneous = delta_G < self.max_gibbs
        
        # Updated score metric
        score = predicted_grams * 250 
        
        return {
            "valid": is_spontaneous,
            "score": score,
            "metrics": {
                "Gibbs_Energy": delta_G,
                "Shannon_Entropy": entropy,
                "Predicted_CO2_Grams": predicted_grams,
                "Virtual_Surface_Area_m2g": surface_area
            }
        }