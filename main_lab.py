import json
import numpy as np
import sys
import time
from config import SIMULATION_CONFIG
from structure_engine import StructureEngine
from quantum_solver import QuantumSolver
from chemical_validator import ChemicalValidator

class ResearchLab:
    def __init__(self):
        self.structure = StructureEngine()
        self.quantum = QuantumSolver()
        self.validator = ChemicalValidator()
        self.results = []
        self.PYROLYSIS_YIELD = 0.33 

    def _generate_lab_protocol(self, recipe, validation, temp_c, acid_time, input_mass, is_raw):
        mix_temp = 85 
        
        if is_raw:
            available_biochar = input_mass * self.PYROLYSIS_YIELD
            note = f"Start with {input_mass}g Raw Straw -> Yields ~{round(available_biochar, 2)}g Biochar"
        else:
            available_biochar = input_mass
            note = f"Start with {input_mass}g Biochar"

        biochar_ratio = recipe["BIOCHAR"]
        total_batch_mass = available_biochar / biochar_ratio
        
        raw_prediction_100g = validation["metrics"]["Predicted_CO2_Grams"]
        batch_prediction_g = (raw_prediction_100g / 100.0) * total_batch_mass
        efficiency_ratio = batch_prediction_g / total_batch_mass

        biochar_prep = {
            "source_material_note": note,
            "step_1": "Grind Rice Straw to fine powder",
            "step_2_acid_wash": f"Soak in HCL (30%) for {acid_time} minutes",
            "step_3_pyrolysis": f"Heat at {temp_c}°C",
            "activation_goal": f"Target Surface Area: {validation['metrics']['Virtual_Surface_Area_m2g']} m2/g"
        }
        
        protocol = {
            "experiment_meta": {
                "input_type": "Raw Straw" if is_raw else "Biochar",
                "input_mass_g": input_mass,
                "calculated_biochar_g": round(available_biochar, 2),
                "resulting_total_mass_g": round(total_batch_mass, 2),
                "predicted_batch_uptake_g": round(batch_prediction_g, 4),
                "efficiency_g_per_g": round(efficiency_ratio, 5),
                "biochar_concentration_percent": round(biochar_ratio * 100, 1)
            },
            "components_grams": {
                "Starch": round(total_batch_mass * recipe["STARCH"], 2),
                "Gelatin": round(total_batch_mass * recipe["GELATIN"], 2),
                "Glycerol": round(total_batch_mass * recipe["GLYCEROL"], 2),
                "Biochar": round(available_biochar, 2)
            },
            "biochar_prep": biochar_prep,
            "process_steps": {
                "mixing_temp_c": mix_temp,
                "stirring_speed_rpm": 500,
                "stirring_min": 45,
                "dry_hours": 24,
                "curing_temp_c": 25
            },
            "prediction_metrics": validation["metrics"],
            "ai_score": round(validation["score"], 2)
        }
        return protocol

    def _print_progress(self, iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filled_length = int(length * iteration // total)
        bar = fill * filled_length + '-' * (length - filled_length)
        sys.stdout.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
        sys.stdout.flush()
        if iteration == total: print()

    def run_interactive_experiment(self):
        print("\n--- QUANTUM-BIO LAB - SCIENTIFIC OPTIMIZATION MODE ---")
        try:
            batch_input = input("1. Experiments [20]: ")
            batch_size = int(batch_input) if batch_input.strip() else 20
            
            type_input = input("2. Input Type: (R)aw Straw or (B)iochar? [R]: ")
            is_raw = False if type_input.lower().strip() == 'b' else True
            
            label = "Raw Straw" if is_raw else "Biochar"
            mass_input = input(f"3. Mass of {label} (g) [20]: ")
            input_mass = float(mass_input) if mass_input.strip() else 20.0
        except:
            batch_size = 20; input_mass = 20.0; is_raw = True

        print(f"\nOptimizing Process for {input_mass}g {label}...")
        
        for i in range(batch_size):
            graph = self.structure.build_hyperbolic_graph()
            prob_dist = self.quantum.solve_wave_function(graph)
            recipe = self._synthesize_recipe(prob_dist)
            
            proposed_temp = np.random.randint(450, 650) 
            proposed_acid_time = np.random.randint(60, 180)
            
            validation = self.validator.validate_sample(recipe, prob_dist, prob_dist, proposed_temp, proposed_acid_time)
            
            if validation["valid"]:
                full_protocol = self._generate_lab_protocol(recipe, validation, proposed_temp, proposed_acid_time, input_mass, is_raw)
                full_protocol["experiment_id"] = i + 1
                self.results.append(full_protocol)
            
            self._print_progress(i + 1, batch_size, prefix='Simulating:', suffix='Done', length=40)
            time.sleep(0.01) 

        self._save_results()

    def _synthesize_recipe(self, prob_dist):
        peak_prob = np.max(prob_dist)
        noise = np.random.normal(0, 0.1)
        
        biochar_ratio = 0.25 + (peak_prob * 0.4) + noise
        biochar_ratio = max(0.15, min(0.55, biochar_ratio)) 
        
        glycerol_ratio = 0.05 + (biochar_ratio * 0.12)
        matrix_rem = 1.0 - (biochar_ratio + glycerol_ratio)
        if matrix_rem < 0.01: matrix_rem = 0.01
        
        gelatin_ratio = matrix_rem * 0.4
        starch_ratio = matrix_rem * 0.6
        
        return {
            "BIOCHAR": biochar_ratio,
            "GELATIN": gelatin_ratio,
            "STARCH": starch_ratio,
            "GLYCEROL": glycerol_ratio
        }

    def _save_results(self):
        self.results.sort(key=lambda x: x["experiment_meta"]["efficiency_g_per_g"], reverse=True)
        
        with open("detailed_lab_protocol.json", "w") as f: json.dump(self.results, f, indent=4)
        if self.results:
            best = self.results[0]
            with open("Best_Recipe.json", "w") as f: json.dump(best, f, indent=4)

            print(f"\n>>> OPTIMIZATION COMPLETE <<<")
            print(f"Top Result (ID: {best['experiment_id']}):")
            print(f"   > Input: {best['experiment_meta']['input_mass_g']}g {best['experiment_meta']['input_type']}")
            print(f"   > Biochar Yield: {best['experiment_meta']['calculated_biochar_g']}g")
            print(f"   > Optimal Pyrolysis: {best['biochar_prep']['step_3_pyrolysis']}")
            print(f"   > Optimal Acid Wash: {best['biochar_prep']['step_2_acid_wash']}")
            print(f"   > Efficiency: {best['experiment_meta']['efficiency_g_per_g']:.4f} g/g")
        else:
            print("\nNo valid candidates.")

if __name__ == "__main__":
    lab = ResearchLab()
    lab.run_interactive_experiment()