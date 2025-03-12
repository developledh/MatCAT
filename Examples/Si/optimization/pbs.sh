#!/bin/bash
#PBS -V
#PBS -N Si
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
function create_potcar {
    elements=$(awk 'NR==6 {print $0}' POSCAR)
    rm -f POTCAR
    for element in $elements; do
        potcar_file=${recommended_potentials[$element]}
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
}

# Create the POTCAR file
create_potcar

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
