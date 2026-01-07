
## -------------------------------------------------------------------------- ##
##    Author:    A.M.M Elsayed                                                ##
##    Email:     ahmedphysica@outlook.com                                     ##
##    Institute: University of Science and Technology of China                ##
## -------------------------------------------------------------------------- ##
##                                                                            ##
##    Updates On:                                                             ##
##    https://github.com/ammelsayed/LHEAnalysis                               ##
## -------------------------------------------------------------------------- ##


import ROOT
import pylhe
from particle import Particle as pa
import os, math

ROOT.gSystem.CompileMacro("Particle.cxx", "kO") # will run first time for compiling only
ROOT.gSystem.Load("Particle_cxx.so")
ROOT.gInterpreter.Declare('#include "Particle.cxx"')
from ROOT import Particle

# ================== Logging ==================== #
from rich.console import Console
from rich.progress import track
from tqdm import tqdm

class Log:

    def __init__(self, print_option = True):
        self.print_option = print_option
        self.console = Console()

    def title(self, string = ""):
        if self.print_option:
            self.console.print("\n######",string, "######\n", style="green")
        else:
            pass

    def proc_title(self, string = ""):
        if self.print_option:
            self.console.print(">>>>", string, "<<<<", style="blue")
        else:
            pass

    def msg(self, string = ""):
        if self.print_option:
            self.console.print("[INFO]", string, style="white")
        else:
            pass    
        
    def err_msg(self, string = ""):
        self.console.print("[WARN]", string, style="red")

# ================== Logging ==================== #


# ============= Helper Functions ================ #

def UseBranch(branch_name, tree, cls="Particle", size=100):
    arr = ROOT.TClonesArray(cls, size)
    tree.SetBranchAddress(branch_name, arr)
    return arr


# ============= Helper Functions ================ #

class LHEAnalysis:
    def __init__(self, status_list = [-1,2,1]):

        self.verbose_mode = True
        self.log = Log(self.verbose_mode)

        self.status_list = status_list
        self.log.proc_title("LHEAnalysis Initiated")
        self.log.msg(f"Reading {', '.join([self.get_status_name(status) for status in self.status_list[:-1]])} and {self.get_status_name(self.status_list[-1])} state particles in the LHE file.")

        self.lhe_file_path = ""

        self.PDG_MAP = \
        {
        "Electron": { 11: {"U1_Charge": -1}, -11: {"U1_Charge": 1} },
        "Muon":     { 13: {"U1_Charge": -1}, -13: {"U1_Charge": 1} },
        "Tau":      { 15: {"U1_Charge": -1}, -15: {"U1_Charge": 1} },

        "Neutrino": {
            12: {"U1_Charge": 0}, 14: {"U1_Charge": 0}, 16: {"U1_Charge": 0},
            -12: {"U1_Charge": 0}, -14: {"U1_Charge": 0}, -16: {"U1_Charge": 0},
        },

        "Top":     { 6: {"U1_Charge": 2/3},  -6: {"U1_Charge": -2/3} },
        "Bottom":  { 5: {"U1_Charge": -1/3}, -5: {"U1_Charge": 1/3} },
        "Charm":   { 4: {"U1_Charge": 2/3},  -4: {"U1_Charge": -2/3} },
        "Strange": { 3: {"U1_Charge": -1/3}, -3: {"U1_Charge": 1/3} },
        "Up":      { 2: {"U1_Charge": 2/3},  -2: {"U1_Charge": -2/3} },
        "Down":    { 1: {"U1_Charge": -1/3}, -1: {"U1_Charge": 1/3} },

        "Gluon": { 21: {"U1_Charge": 0} },
        "Gamma": { 22: {"U1_Charge": 0} },
        "Z":     { 23: {"U1_Charge": 0} },
        "W":     { 24: {"U1_Charge": 1}, -24: {"U1_Charge": -1} },
        "H":     { 25: {"U1_Charge": 0} }
        }

        self.new_definitions = {}

    def get_status_name(self, status):
        return {-1: "Initial", 2: "Intermediate", 1: "Final" }.get(status)

    def LoadLHE(self, lhe_file_path):
        self.lhe_file_path = lhe_file_path

    def define_particle(self, name, pid, U1_Charge):
        # If particle group doesn't exist, create it
        if name not in self.PDG_MAP:
            self.PDG_MAP[name] = {}
            self.log.msg(f"Particle class {name} added.")

        # Add the new PDG entry with its charge
        self.PDG_MAP[name][pid] = {"U1_Charge": float(U1_Charge)}
        self.log.msg(f"Class {name} : pdgID of {pid} with Q[U(1)] = {float(U1_Charge)} added.")
    
    def get_particle_name(self, pdg_id):
        if isinstance(pdg_id, list):
            return ', '.join([self.get_particle_name(pid) for pid in pdg_id])
        
        if isinstance(pdg_id, int) or isinstance(pdg_id, float):
            try:
                pname = pa.from_pdgid(float(pdg_id)).name
            except Exception:
                try:
                    pname = self.new_definitions[int(pdg_id)]
                except Exception: pname = "Unknown"
            return pname

        else:
            self.log.err_msg("Unknown input type, must be int, float or a list of int/floats.")
    
    def get_particle_info(self, pdgID):
        for name, mapping in self.PDG_MAP.items():
            if pdgID in mapping:
                return name, float(mapping[pdgID]["U1_Charge"])
        return None, None
    
    def get_full_extension(self, file_name):
        # everything after first dot
        base = os.path.basename(file_name)
        parts = base.split(".")
        if len(parts) > 2:
            return ".".join(parts[1:])  
        return parts[-1]
    
    def convert_with_delphes(self):
        import subprocess

        delphes_dir = "~/mg5_v3.5.5/Delphes"
        delphes_card_dir = "/data/ammelsayed/Framework/MC_Samples/my_delphes_card.tcl"
        output_root_file_name = "delphes_output.root"

        command = f"{delphes_dir}/DelphesLHEF {delphes_card_dir} {self.lhe_file_path} {output_root_file_name}"
        subprocess.run(command, check=True)
        

    def SaveAsROOT(self, file_name):

        arrays = {}

        ## Detects if the ROOT file already exisits. If it exsits, then delete everything inside it and write on clean!
        if self.get_full_extension(file_name) != "root":
            self.log.err_msg("Must be a '.root' file. Using default name: 'lhe_analysis.root'")
            outfile = ROOT.TFile("lhe_analysis.root", "RECREATE")
            self.log.msg(f"'lhe_analysis.root' created.")
        else:
            outfile = ROOT.TFile(file_name, "RECREATE")
            self.log.msg(f"'{file_name}' created.")

        tree = ROOT.TTree("LHEAnalysis", "LHE Analysis TTree")
        tree.SetDirectory(outfile)

        def create_array_if_needed(name):
            if name not in arrays:
                arr = ROOT.TClonesArray("Particle", 100)
                arrays[name] = arr
                tree.Branch(name, arrays[name])
            return arrays[name]

        def sort_by_pt(arr):
            n = arr.GetEntriesFast()
            tmp = [arr.At(i) for i in range(n)]
            tmp.sort(key=lambda p: p.PT, reverse=True)
            for i, p in enumerate(tmp):
                arr[i] = p

        def insert_particle_into_clonesarray(arr, index, p):
            slot = arr.ConstructedAt(index)

            slot.Status             = p.Status
            slot.pdgID              = p.pdgID
            slot.pdgID_Mother1      = p.pdgID_Mother1
            slot.pdgID_Mother2      = p.pdgID_Mother2
            slot.pdgID_Daughter1    = p.pdgID_Daughter1
            slot.pdgID_Daughter2    = p.pdgID_Daughter2

            slot.Px                 = p.Px
            slot.Py                 = p.Py
            slot.Pz                 = p.Pz
            slot.Energy             = p.Energy
            slot.Mass               = p.Mass

            slot.PT                 = p.PT
            slot.Eta                = p.Eta
            slot.Phi                = p.Phi

            slot.Charge             = p.Charge   

            slot.Lifetime           = p.Lifetime
            slot.Helicity           = p.Helicity

            slot.SetP4(p.Px, p.Py, p.Pz, p.Energy) 

            slot.Color1             = p.Color1
            slot.Color2             = p.Color2

            slot.Index              = p.Index
            slot.Mother1            = p.Mother1
            slot.Mother2            = p.Mother2
            
            slot.Daughter1          = p.Daughter1
            slot.Daughter2          = p.Daughter2



        self.log.msg("Looping over events, please wait paitently.")

        for event in tqdm(pylhe.read_lhe(self.lhe_file_path)):

            # Clear all arrays for this event
            for arr in arrays.values():
                arr.Clear("C")

            # Temporary Python lists per species
            temp = {}

            def build_daughters_map(event):
                """Return dict: key = particle index (1-based), value = list of daughter indices (1-based)."""
                n = len(event.particles)
                daughters_map = {i: [] for i in range(1, n + 1)}
                for idx, p in enumerate(event.particles, start=1):
                    m1 = int(p.mother1)
                    m2 = int(p.mother2)
                    if 1 <= m1 <= n:
                        daughters_map[m1].append(idx)
                    # append for m2 only if it's different from m1 (prevents double appending when m1==m2)
                    if 1 <= m2 <= n and m2 != m1:
                        daughters_map[m2].append(idx)
                # Optional: ensure uniqueness while preserving order (defensive)
                for k, lst in daughters_map.items():
                    if len(lst) > 1:
                        daughters_map[k] = list(dict.fromkeys(lst))
                return daughters_map


            def get_mothers(particle, pdgIds):
                """Return (m1, m2, pdg_m1, pdg_m2). pdg_* = 0 if mother index invalid."""
                m1 = int(particle.mother1)
                m2 = int(particle.mother2)
                pdgid_m1 = int(pdgIds[m1 - 1]) if (1 <= m1 <= len(pdgIds)) else 0
                pdgid_m2 = int(pdgIds[m2 - 1]) if (1 <= m2 <= len(pdgIds)) else 0
                return m1, m2, pdgid_m1, pdgid_m2

            def get_daughters(index, daughters_map, pdgIds):
                """
                For particle at `index` (1-based) return (d1_idx, d2_idx, pdg_d1, pdg_d2).
                If fewer than two daughters, missing entries return 0. If more than two daughters,
                only the first two are returned (your physics says at most 2).
                """
                ds = daughters_map.get(index, [])
                d1 = ds[0] if len(ds) >= 1 else 0
                d2 = ds[1] if len(ds) >= 2 else 0
                pdgid_d1 = int(pdgIds[d1 - 1]) if (1 <= d1 <= len(pdgIds)) else 0
                pdgid_d2 = int(pdgIds[d2 - 1]) if (1 <= d2 <= len(pdgIds)) else 0
                return d1, d2, pdgid_d1, pdgid_d2

            
            pdgIds = [int(p.id) for p in event.particles]
            daughters_map = build_daughters_map(event)

            for idx, particle in enumerate(event.particles, start=1):
                pdgid = int(particle.id)
                # pdgIds.append(pdgid)
                status = int(particle.status)

                if status not in self.status_list:
                    continue

                name, charge = self.get_particle_info(pdgid)
                if name is None:
                    continue

                if name not in temp:
                    temp[name] = []

                px = float(particle.px)
                py = float(particle.py)
                pz = float(particle.pz)
                E  = float(particle.e)
                m  = float(particle.m)

                pt  = math.sqrt(px*px + py*py)
                phi = math.atan2(py, px)
                # handle zero division
                eta = 0.5 * math.log((E + pz)/(E - pz)) if (E != abs(pz)) else 0

                m1, m2, pdgid_m1, pdgid_m2 = get_mothers(particle, pdgIds)
                d1_idx, d2_idx, pdg_d1, pdg_d2 = get_daughters(idx, daughters_map, pdgIds)

                p = Particle()
                p.pdgID = pdgid
                p.Status = status
                p.Index = int(idx)
                p.Mother1 = m1
                p.Mother2 = m2
                p.pdgID_Mother1 = pdgid_m1
                p.pdgID_Mother2 = pdgid_m2
                p.Daughter1 = int(d1_idx)
                p.Daughter2 = int(d2_idx)
                p.pdgID_Daughter1 = int(pdg_d1)
                p.pdgID_Daughter2 = int(pdg_d2)
                p.Px = px
                p.Py = py
                p.Pz = pz
                p.Energy = E
                p.Mass = m
                p.PT = pt
                p.Eta = eta
                p.Phi = phi
                p.Charge = float(charge)
                p.Color1 = int(particle.color1)
                p.Color2 = int(particle.color2)
                p.Lifetime = float(particle.lifetime)
                p.Helicity = int(particle.spin)
                p.SetP4(px, py, pz, E)

                temp[name].append(p)

            # Now fill ROOT arrays sorted by pt
            for name, plist in temp.items():
                arr = create_array_if_needed(name)

                plist.sort(key=lambda p: p.PT, reverse=True)

                for i, p in enumerate(plist):
                    insert_particle_into_clonesarray(arr, i, p)

            tree.Fill()

        outfile.cd()
        tree.Write("", ROOT.TObject.kOverwrite)
        outfile.Close()
