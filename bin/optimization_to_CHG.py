import os
import shutil

def modify_incar_for_making_chgcar(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    # Remove ISIF and IBRION lines if they exist
    updated_lines = []
    for line in lines:
        if "ISIF" in line or "IBRION" in line:
            continue
        updated_lines.append(line)

    # Add or modify necessary lines
    def add_or_modify_line(lines, keyword, value):
        found = False
        for i, line in enumerate(lines):
            if keyword in line:
                lines[i] = f"{keyword} = {value}\n"
                found = True
                break
        if not found:
            lines.append(f"{keyword} = {value}\n")

    add_or_modify_line(updated_lines, "NSW", 0)
    add_or_modify_line(updated_lines, "ICHARG", 2)
    add_or_modify_line(updated_lines, "LCHARG", ".TRUE.")
    add_or_modify_line(updated_lines, "LORBIT", 11)

    with open(file_path, "w") as f:
        f.writelines(updated_lines)

def process_kpoints_for_chgcar():
    source_file = "optimization/KPOINTS"
    dest_file = "making_CHGCAR/KPOINTS"
    
    with open(source_file, "r") as f:
        lines = f.readlines()
    
    if len(lines) < 4:
        raise ValueError("In KPOINTS files, There are no line 4.")
    
    grid_line = lines[3].strip()
    numbers = grid_line.split()
    try:
        new_numbers = [str(int(num) * 4) for num in numbers]
    except ValueError:
        # if it is not interger, take float
        new_numbers = [str(float(num) * 4) for num in numbers]
    new_grid_line = " ".join(new_numbers) + "\n"
    
    lines[3] = new_grid_line

    with open(dest_file, "w") as f:
        f.writelines(lines)

def copy_and_modify_optimization():
    # making_CHGCAR dir 
    os.makedirs("making_CHGCAR", exist_ok=True)

    shutil.copy("optimization/CONTCAR", "making_CHGCAR/POSCAR")

    # KPOINTS 처리: In optimization/KPOINTS files, After Modify line 4 to quadruple, save to making_CHGCAR/KPOINTS
    process_kpoints_for_chgcar()

    shutil.copy("optimization/INCAR", "making_CHGCAR/")
    shutil.copy("optimization/POTCAR", "making_CHGCAR/")
    shutil.copy("optimization/pbs.sh", "making_CHGCAR/")

    modify_incar_for_making_chgcar("making_CHGCAR/INCAR")
    
def is_optimization_converged():
    out_path = os.path.join("optimization", "vasp-.out")
    if not os.path.exists(out_path):
        print("vasp-.out not found. Aborting.")
        return False
    with open(out_path, "r") as f:
        content = f.read()
    keyword = "reached required accuracy - stopping structural energy minimisation"
    return keyword in content


def submit_qsub_in_making_chgcar(directory):
    os.chdir(directory)
    old_content = "time python optimization_to_CHG.py"
    new_content = "time python CHGCAR_to_dos.py"

    with open("pbs.sh", 'r') as file:
        file_content = file.read()

    new_file_content = file_content.replace(old_content, new_content)

    with open("pbs.sh", 'w') as file:
        file.write(new_file_content)

    os.system("qsub pbs.sh")    


# Call the function
os.chdir("..")
if is_optimization_converged():
    copy_and_modify_optimization()
    submit_qsub_in_making_chgcar("making_CHGCAR")
else:
    sys.exit("Optimization not converged. Skipping making_CHGCAR step.")

