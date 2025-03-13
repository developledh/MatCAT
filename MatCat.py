import os
import glob
import shutil
import ase.io as io
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

# Define the path to the POTCAR directory
# Recommended potentials dictionary in a more compact format
declare -A recommended_potentials=(
    ["H"]="H" ["He"]="He" ["Li"]="Li_sv" ["Be"]="Be_sv" ["B"]="B" ["C"]="C" ["N"]="N"
    ["O"]="O" ["F"]="F" ["Ne"]="Ne" ["Na"]="Na_pv" ["Mg"]="Mg" ["Al"]="Al" ["Si"]="Si"
    ["P"]="P" ["S"]="S" ["Cl"]="Cl" ["Ar"]="Ar" ["K"]="K_sv" ["Ca"]="Ca_sv" ["Sc"]="Sc_sv"
    ["Ti"]="Ti_sv" ["V"]="V_sv" ["Cr"]="Cr_pv" ["Mn"]="Mn_pv" ["Fe"]="Fe" ["Co"]="Co"
    ["Ni"]="Ni" ["Cu"]="Cu" ["Zn"]="Zn" ["Ga"]="Ga_d" ["Ge"]="Ge_d" ["As"]="As"
    ["Se"]="Se" ["Br"]="Br" ["Kr"]="Kr" ["Rb"]="Rb_sv" ["Sr"]="Sr_sv" ["Y"]="Y_sv"
    ["Zr"]="Zr_sv" ["Nb"]="Nb_sv" ["Mo"]="Mo_sv" ["Tc"]="Tc_pv" ["Ru"]="Ru_pv" ["Rh"]="Rh_pv"
    ["Pd"]="Pd" ["Ag"]="Ag" ["Cd"]="Cd" ["In"]="In_d" ["Sn"]="Sn_d" ["Sb"]="Sb"
    ["Te"]="Te" ["I"]="I" ["Xe"]="Xe" ["Cs"]="Cs_sv" ["Ba"]="Ba_sv" ["La"]="La"
    ["Ce"]="Ce" ["Pr_3"]="Pr_3"
    ["Nd"]="Nd_3" ["Pm"]="Pm_3" ["Sm"]="Sm_3" ["Eu"]="Eu_2"
    ["Gd"]="Gd_3" ["Tb"]="Tb_3" ["Dy"]="Dy_3" ["Ho"]="Ho_3" ["Er"]="Er_3"
    ["Tm"]="Tm_3" ["Yb"]="Yb_2" ["Lu"]="Lu_3" ["Hf"]="Hf_pv"
    ["Ta"]="Ta_pv" ["W"]="W_sv" ["Re"]="Re" ["Os"]="Os" ["Ir"]="Ir" ["Pt"]="Pt"
    ["Au"]="Au" ["Hg"]="Hg" ["Tl"]="Tl_d" ["Pb"]="Pb_d"
    ["Bi"]="Bi_d" ["Po"]="Po_d" ["At"]="At" ["Rn"]="Rn"
    ["Fr"]="Fr_sv" ["Ra"]="Ra_sv" ["Ac"]="Ac" ["Th"]="Th" ["Pa"]="Pa"
    ["U"]="U" ["Np"]="Np" ["Pu"]="Pu"
    ["Am"]="Am" ["Cm"]="Cm" ["Cf"]="Cf"
)

# Function to create the POTCAR file based on the elements in the POSCAR
function create_potcar {{
    elements=$(awk 'NR==6 {{print $0}}' POSCAR)
    rm -f POTCAR
    for element in $elements; do
        potcar_file=${{recommended_potentials[$element]}}
        if [ -z "$potcar_file" ]; then
            echo "No recommended POTCAR for element $element"
            exit 1
        fi
        if [ -f $POTCAR_DIR/$potcar_file/POTCAR ]; then
            cat $POTCAR_DIR/$potcar_file/POTCAR >> POTCAR
        else
            echo "POTCAR for element $element ($potcar_file) not found in $POTCAR_DIR"
            exit 1
        fi
    done
}}

# Create the POTCAR file
create_potcar

# VASP PART 
time mpirun /home01/x3083a0/bin/vasp_std > vasp-.out

echo "Closing Time is $(date)"

# Anaconda PART 

source /scratch/x3083a05/anaconda3/etc/profile.d/conda.sh
conda activate DFT

# Next caculation 
cd ..
cd bin
time python optimization_to_CHG.py
''')

    @staticmethod
    def write_incar(directory):           ## You write your INCAR option depending on your target material
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
    def write_kpoints(directory):     ## You write your KPOINTS depending on your target material 
        kpoints_content = """
kpoints
0
gamma
4 4 4
0 0 0"""
        kpoints_path = os.path.join(directory, "KPOINTS")
        with open(kpoints_path, "w") as f:
            f.write(kpoints_content)

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
        """ Copies the 'bin' directory from MatCat into the specified folder, avoiding duplicate files. """
        source_bin = os.path.join(os.path.dirname(__file__), "bin")
        target_bin = os.path.join(destination_folder, "bin")

        if not os.path.exists(source_bin):
            print(f"Source bin directory does not exist: {source_bin}")
            return

        if not os.path.exists(target_bin):
            shutil.copytree(source_bin, target_bin)  # If folder doesn't exist, make the folder. 
        else:
            print(f"Bin folder already exists in {destination_folder}. Updating files...")

            for item in os.listdir(source_bin):
                src_item = os.path.join(source_bin, item)
                dest_item = os.path.join(target_bin, item)

                if os.path.isdir(src_item):
                    shutil.copytree(src_item, dest_item, dirs_exist_ok=True)
                
                elif not os.path.exists(dest_item):
                    shutil.copy2(src_item, dest_item)

        print(f"Copied bin directory to: {target_bin}")
