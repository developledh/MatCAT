# MaTCat (MaTerials Calculation and Automation Tool) 

## 📌 Overview
**Materials Calculation and Automation Tool (VASP)** is a computational tool developed by the **Innovative Nanomaterials Design Lab**.  
This project aims to **automate material property calculations and modeling based on VASP**, reducing manual efforts and streamlining the computational workflow.  

## 🚀 Key Features  
- **Automation of VASP calculations and modeling**  
- **Linux scripting and Python integration** for seamless execution  
- **Workflow optimization** to handle multi-step calculations automatically

## ✨Installation Guide
To install MatCat in `/opt/MatCat`, first clone the repository:

```bash
git clone https://github.com/developledh/MatCat.git 
```


## 📁 Folder Structure & Required Files  
To use this tool effectively, the following setup is required:  

Setting Up VASP POTCAR Directory
1. POTCAR Directory Setup → Ensure that the necessary POTCAR files are available
To set the POTCAR directory environment variable, add the following line to your shell configuration file (e.g., `.bashrc`, `.bash_profile`, or `.zshrc`):

```sh
export POTCAR_DIR="$HOME/PBE"
source ~/.bashrc  # or source ~/.zshrc
```
2. Modeling Folder Structure  
   - Create an the desired model (.cif format) inside cifs folder
   - For automation, Copy a bin of Python code within the module under each model case folder.
   - The bin folder looks like this : 
     - `optimization_to_CHG.py` → Converts optimized structure to CHGCAR  
     - `CHGCAR_to_dos.py` → Automates Density of States (DOS) calculations  
     - `CHGCAR_to_band.py` → Automates Band Structure calculations
      
3. Run run_matcat.py
   - Utilizing MatCat class to perform automated calculations
   - You will need to specify the download path within Run_matcat.py
     
```sh
cp -r /your_path/MatCat/Example/run_matcat.py /your_usage_path
```
- And then, In run_matcat.py you have to modify your download path

```sh
import sys
sys.path.append("$Your download path$/MatCat")
```

**Example folders structure**
- When you run `run_matcat.py` in your choice a folder, the following directory structure will be automatically created:

```sh
📂 Your_Chosen_Folder/      # The main folder where you execute run_matcat.py
├── 📂 cifs/ # Initial structure files (Folder containing CIF files for conversion)
│
├── 📂 Si/ # Example structure folder
│ ├── 📂 bin/                    # Copied from MatCat/bin (Executable scripts for automatic calculation)
│ ├── 📂 optimization/           # Directory where VASP calculations are performed
│ │ ├── POSCAR                   # Structure file converted from CIF
│ │ ├── INCAR, KPOINTS, pbs.sh       # VASP input files and PBS script
│ │ ├── POTCAR                   # Automatically generated in pbs.sh
│
├── 📂 SiC/ # Another example structure folder
│ ├── 📂 bin/ 
│ ├── 📂 optimization/
│ │ ├── POSCAR
│ │ ├── INCAR, KPOINTS, pbs.sh # Generated VASP input files for optimization
│ │ ├── POTCAR # Automatically generated in pbs.sh
```

## 🛠️ How It Works  
- After **structure optimization**, multiple follow-up calculations are required.  
- This tool **automates** the process using **Linux shell scripts** and **Python integration**.  
- The scripts generate `.sh` files to execute required calculations efficiently.  

## 📝 Usage Instructions  
1. **Set up the required directory and files** as mentioned in the structure above.  
2. **Ensure Python and VASP are properly installed** on your system.  
3. **Run the scripts will perform automated calculations in sequence (Optimization -> Making_CHGCAR -> dos -> band)**  
4. **Monitor and analyze results** from the generated output files.  

## 📌 Future Enhancements  
- Expanding support for additional VASP calculations
- Improving automatic modeling mainly for metal-halide
- Improving automatic modeling mainly for interface 
- Improving error handling and logging
- Improvement of output data collection

## 👥 Authors & Contributors  
This tool was developed by the **Advanced Nano Materials Design Laboratory**.  

### **Main Contributors:**  
- **Dong Ho Lee** – Algorithm & Automation, VASP & Computational MaterialsLinux Scripting & Optimization  

We welcome contributions! If you would like to contribute, please feel free to submit a pull request or contact us.

## 📄 Additional Information  
For detailed explanations, installation steps, and troubleshooting, please refer to the project documentation.  



Let me know if you'd like to customize the **contributors section** with real names or add affiliations! 🚀
