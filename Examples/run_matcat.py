import sys
sys.path.append("/scratch/x3083a05/opt/MatCat")

import os
from MatCat import MatCat

# 현재 작업 디렉토리
current_dir = os.getcwd()
cif_folder = os.path.join(current_dir, "cifs")

# CIF 파일 변환하여 폴더 생성
MatCat.convert_cif_to_poscar(cif_folder=cif_folder, output_folder=current_dir)

# 현재 디렉토리 내의 모든 폴더를 순회 (cifs 폴더 제외)
for d in os.listdir(current_dir):
    dir_path = os.path.join(current_dir, d)

    if d == "cifs" or not os.path.isdir(dir_path):
        continue  # 'cifs' 폴더는 제외

    os.chdir(dir_path)
    print(f"Processing directory: {dir_path}")

    # 각 구조 폴더 내부에 `bin` 복사
    MatCat.copy_bin(dir_path)

    # 'optimization' 폴더 생성 후 실행
    optimization_dir = os.path.join(dir_path, "optimization")
    os.makedirs(optimization_dir, exist_ok=True)
    os.chdir(optimization_dir)

    # VASP 입력 파일 생성
    MatCat.write_incar(".")
    MatCat.write_kpoints(".")
    MatCat.modify_incar_for_magmoms()
    MatCat.write_pbs_script(d)

    # PBS 스크립트 실행
    os.system("qsub pbs.sh")

    # 원래 디렉토리로 복귀
    os.chdir(current_dir)

print("All tasks completed successfully.")

