import os

class MatCat:
    @staticmethod
    def write_pbs_script(modelname):
        pbs_script_name = 'pbs.sh'
        with open(pbs_script_name, 'w') as f:
            f.write(f'''#!/bin/bash
#PBS -V
#PBS -N {modelname}
#PBS -q normal
#PBS -l walltime=18:00:00
#PBS -l select=18:ncpus=18:mpiprocs=18:ompthreads=1
#PBS -A etc

ulimit -l unlimited
ulimit -a unlimited
ulimit -m unlimited
ulimit -n unlimited
ulimit -u unlimited
ulimit -s unlimited
ulimit -v unlimited

module purge
module add craype-network-opa
module add craype-mic-knl
module add intel/19.0.5
module add impi/19.0.5
module add fftw_mpi/3.3.7

cd $PBS_O_WORKDIR

# VASP PART 
time mpirun /home01/x3083a05/bin/vasp_std > vasp-.out

echo "Closing Time is $(date)"

# Anaconda PART 

source /scratch/x3083a05/anaconda3/etc/profile.d/conda.sh
conda activate DFT

# Next caculation 
cd ..
cd bin
time python optimization_to_band.py
''')

    @staticmethod
    def write_incar(directory):
        incar_content = """
ALGO = Normal
EDIFF = 0.002
ENCUT = 520
ISMEAR = 2
ISPIN = 2
LREAL = Auto
LWAVE = False
NELM = 100
NSW = 100
PREC = Normal
SIGMA = 0.05
MAGMOM = 96*0.6
LCHARG = .False.
"""
        incar_path = os.path.join(directory, "INCAR")
        with open(incar_path, "w") as f:
            f.write(incar_content)

    @staticmethod
    def modify_incar_for_magmoms():
        poscar_file = "POSCAR"
        incar_file = "INCAR"

        with open(poscar_file, 'r') as f:
            lines = f.readlines()
            atom_counts = lines[6].split()
            total_atoms = sum(map(int, atom_counts))

        with open(incar_file, 'r') as f:
            incar_lines = f.readlines()

        updated_lines = [line for line in incar_lines if "MAGMOM" not in line]
        updated_lines.append(f"MAGMOM = {total_atoms}*0.6\n")

        with open(incar_file, 'w') as f:
            f.writelines(updated_lines)

        print(f"Modified {incar_file} with MAGMOM = {total_atoms}*0.6")

    @staticmethod
    def write_kpoints(directory):
        kpoints_content = """
kpoints
0
gamma
8 8 8
0 0 0"""
        kpoints_path = os.path.join(directory, "KPOINTS")
        with open(kpoints_path, "w") as f:
            f.write(kpoints_content)

    @staticmethod
    def create_potcar():
        poscar_file = "POSCAR"
        potcar_file = "POTCAR"
        potcar_dir = os.getenv("POTCAR_DIR")
        
        if not potcar_dir:
            print("Environment variable POTCAR_DIR is not set.")
            return

        recommended_potentials = {
            "H": "H", "He": "He", "Li": "Li_sv", "Be": "Be_sv", "B": "B", "C": "C", "N": "N",
            "O": "O", "F": "F", "Ne": "Ne", "Na": "Na_pv", "Mg": "Mg", "Al": "Al", "Si": "Si",
            "P": "P", "S": "S", "Cl": "Cl", "Ar": "Ar", "K": "K_sv", "Ca": "Ca_sv", "Sc": "Sc_sv",
            "Ti": "Ti_sv", "V": "V_sv", "Cr": "Cr_pv", "Mn": "Mn_pv", "Fe": "Fe", "Co": "Co",
            "Ni": "Ni", "Cu": "Cu", "Zn": "Zn", "Ga": "Ga_d", "Ge": "Ge_d", "As": "As",
            "Se": "Se", "Br": "Br", "Kr": "Kr", "Rb": "Rb_sv", "Sr": "Sr_sv", "Y": "Y_sv",
            "Zr": "Zr_sv", "Nb": "Nb_sv", "Mo": "Mo_sv", "Tc": "Tc_pv", "Ru": "Ru_pv", "Rh": "Rh_pv",
            "Pd": "Pd", "Ag": "Ag", "Cd": "Cd", "In": "In_d", "Sn": "Sn_d", "Sb": "Sb",
            "Te": "Te", "I": "I", "Xe": "Xe", "Cs": "Cs_sv", "Ba": "Ba_sv", "La": "La"
        }
        
        with open(poscar_file, 'r') as f:
            lines = f.readlines()
            elements = lines[5].split()
        
        with open(potcar_file, 'w') as potcar:
            for element in elements:
                if element in recommended_potentials:
                    pot_path = os.path.join(potcar_dir, recommended_potentials[element], "POTCAR")
                    if os.path.exists(pot_path):
                        with open(pot_path, 'r') as p:
                            potcar.write(p.read())
                    else:
                        print(f"POTCAR for element {element} not found in {pot_path}")
                        return
                else:
                    print(f"No recommended POTCAR for element {element}")
                    return
        print("POTCAR file successfully created.")

      print("POTCAR file successfully created.")

    @staticmethod
    def convert_cif_to_poscar(cif_folder="./cifs", output_folder="."):
        cif_files = glob.glob(os.path.join(cif_folder, "*.cif"))
        
        if not cif_files:
            print("No CIF files found in the specified folder.")
            return
        
        for cif_file in cif_files:
            file_name = os.path.basename(cif_file).replace(".cif", "")
            folder_path = os.path.join(output_folder, file_name)
            os.makedirs(folder_path, exist_ok=True)
            
            structure = io.read(cif_file)
            poscar_path = os.path.join(folder_path, "POSCAR")
            io.write(poscar_path, structure, format='vasp')
            
            print(f"Processed: {file_name}")
            
            optimization_path = os.path.join(folder_path, "optimization")
            os.makedirs(optimization_path, exist_ok=True)
            os.rename(poscar_path, os.path.join(optimization_path, "POSCAR"))


    @staticmethod
    def copy_bin(destination_folder):
        source_bin = os.path.join(os.path.dirname(__file__), "bin")
        target_bin = os.path.join(destination_folder, "bin")
        
        if os.path.exists(source_bin):
            shutil.copytree(source_bin, target_bin, dirs_exist_ok=True)
            print(f"Copied bin directory to: {target_bin}")
        else:
            print("bin directory not found in MatCat package.")
