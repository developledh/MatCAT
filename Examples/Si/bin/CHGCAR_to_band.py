import os
import shutil

def modify_incar_for_Band(file_path):
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

     # 내용 추가하고 싶은 것들 사용
    add_or_modify_line(updated_lines, "LORBIT", 11)
    add_or_modify_line(updated_lines, "ICHARG", 11)
    add_or_modify_line(updated_lines, "LCHARG", ".False.")
    add_or_modify_line(updated_lines, "ISPIN", 2 )



    with open(file_path, "w") as f:
        f.writelines(updated_lines)

def copy_and_modify_band():
    # Create the band directory if it doesn't exist
    os.makedirs("band", exist_ok=True)



    # Copy files to dos 
    os.system("cp -r ./making_CHGCAR/{INCAR,POTCAR,POSCAR,CHGCAR,pbs.sh} ./band/")
    os.system("cp -r ./dos/pbs.sh ./band/")
    print("you should write k-path depending on you")
    print("please write KPOINTS file on How you want")


    # Modify INCAR in making_CHGCAR 
    modify_incar_for_Band("band/INCAR")


def submit_qsub_in_band(directory):
    # Change directory to making_CHGCAR
    os.chdir(directory)

    old_content = "time python CHGCAR_to_band.py"
    new_content = ""

    # Open the slurm.sh file for reading
    with open("pbs.sh", 'r') as file:
        file_content = file.read()

    # Replace the old_content with the new_content
    new_file_content = file_content.replace(old_content, new_content)

    # Write the modified content back to the slurm.sh file
    with open("pbs.sh", 'w') as file:
        file.write(new_file_content)
    # Submit qsub command
    # os.system("sbatch pbs.sh")
    
    os.system("sumo-kgen")
    os.system("mv KPOINTS_band KPOINTS")
    os.system("qsub pbs.sh")


# Call the function

os.chdir("..")
copy_and_modify_band()
submit_qsub_in_band("band")

