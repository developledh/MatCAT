import os
import glob
from ase import io


def convert_cif_to_poscar(cif_folder, output_folder):
    # Get CIF file list
    cif_files = glob.glob(os.path.join(cif_folder, "*.cif"))
    
    if not cif_files:
        print("No CIF files found in the specified folder.")
        return
    
    for cif_file in cif_files:
        # Extract filename without extension
        file_name = os.path.basename(cif_file).replace(".cif", "")
        folder_path = os.path.join(output_folder, file_name)
        
        # Create output folder
        os.makedirs(folder_path, exist_ok=True)
        
        # Read CIF and write as POSCAR
        structure = io.read(cif_file)
        structure = io.read(cif_file)
        poscar_path = os.path.join(folder_path, "POSCAR")
        io.write(poscar_path, structure, format='vasp')
        
        print(f"Processed: {file_name}")
 
        # 만든 POSCAR파일 optimization 폴더 만들고 옮기기
        optimization_path = os.path.join(folder_path, "optimization")
        os.makedirs(optimization_path, exist_ok=True)

        # Move POSCAR file into 'optimization' folder
        os.rename(poscar_path, os.path.join(optimization_path, "POSCAR"))

        

# 실행 예제
cif_folder = "./cifs"  # CIF 파일들이 있는 폴더 경로
output_folder = "."   # 결과 폴더 경로
convert_cif_to_poscar(cif_folder, output_folder)


