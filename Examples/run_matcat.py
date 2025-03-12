import sys
sys.path.append("/scratch/x3083a05/opt/MatCat")   #Please modify to your path where you downloaded MatCat codes
import os
from MatCat import MatCat

# Current job directory 
current_dir = os.getcwd()
cif_folder = os.path.join(current_dir, "cifs")

# Convert CIF files to folder_name 
MatCat.convert_cif_to_poscar(cif_folder=cif_folder, output_folder=current_dir)

#Repeat in the current directory except structure folder (cifs)
for d in os.listdir(current_dir):
    dir_path = os.path.join(current_dir, d)

    if d == "cifs" or not os.path.isdir(dir_path):
        continue  

    os.chdir(dir_path)
    print(f"Processing directory: {dir_path}")

    # In each materials, copy  bin
    MatCat.copy_bin(dir_path)

    # Generation 'optimization' folder and Execute
    optimization_dir = os.path.join(dir_path, "optimization")
    os.makedirs(optimization_dir, exist_ok=True)
    os.chdir(optimization_dir)

    # Generation INPUT FILES for VASP    ## IF you want use other INCAR or KPOINTS, Please back to MatCat.py 
    MatCat.write_incar(".")
    MatCat.write_kpoints(".")
    MatCat.modify_incar_for_magmoms()
    MatCat.write_pbs_script(d)

    # Running PBS script
    os.system("qsub pbs.sh")

    # Go Back to original directory 
    os.chdir(current_dir)

print("All tasks completed successfully.")

