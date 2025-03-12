# MaTCAT (Materials Calculation and Automation Tool) 

## 📌 Overview
**Materials Calculation and Automation Tool (VASP)** is a computational tool developed by the **Innovative Nanomaterials Design Lab**.  
This project aims to **automate material property calculations and modeling based on VASP**, reducing manual efforts and streamlining the computational workflow.  

## 🚀 Key Features  
- **Automation of VASP calculations and modeling**  
- **Linux scripting and Python integration** for seamless execution  
- **Workflow optimization** to handle multi-step calculations automatically  

## 📁 Folder Structure & Required Files  
To use this tool effectively, the following setup is required:  

1. **POTCAR Directory Setup** → Ensure that the necessary POTCAR files are available  
2. **Modeling Folder Structure**  
   - Create an `optimization/` folder inside the desired model directory  
   - Include the following Python scripts for automation:  
     - `optimization_to_CHG.py` → Converts optimized structure to CHGCAR  
     - `CHGCAR_to_dos.py` → Automates Density of States (DOS) calculations  
     - `CHGCAR_to_band.py` → Automates Band Structure calculations  

## 🛠️ How It Works  
- After **structure optimization**, multiple follow-up calculations are required.  
- This tool **automates** the process using **Linux shell scripts** and **Python integration**.  
- The scripts generate `.sh` files to execute required calculations efficiently.  

## 📝 Usage Instructions  
1. **Set up the required directory and files** as mentioned in the structure above.  
2. **Ensure Python and VASP are properly installed** on your system.  
3. **Run the scripts in sequence** to perform automated calculations.  
4. **Monitor and analyze results** from the generated output files.  

## 📌 Future Enhancements  
- Expanding support for additional VASP calculations
- Improving automatic modeling for metal-halide 
- Improving error handling and logging  
- Enhancing parallel computation capabilities  

## 👥 Authors & Contributors  
This tool was developed by the **Advanced Nano Materials Design Laboratory**.  

### **Main Contributors:**  
- **Dong Ho Lee** – Algorithm & Automation, VASP & Computational MaterialsLinux Scripting & Optimization  

We welcome contributions! If you would like to contribute, please feel free to submit a pull request or contact us.

## 📄 Additional Information  
For detailed explanations, installation steps, and troubleshooting, please refer to the project documentation.  



Let me know if you'd like to customize the **contributors section** with real names or add affiliations! 🚀
