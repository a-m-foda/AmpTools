#jlab hpci12k01
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:
export PATH=$PATH:/usr/local/cuda/bin:
export CUDA_INSTALL_PATH=/usr/local/cuda
export AMPTOOLS_HOME=~/AmpTools
#CEDAR
module load cuda
module laod root
export AMPTOOLS_HOME=/scratch/amfoda/AmpTools
#Both
export CUDA_FLAGS="-m64 -arch=sm_75" #k20m=30 TitanRTX=75 P100=60
export AMPPLOTTER=$AMPTOOLS_HOME/AmpPlotter
export AMPTOOLS=$AMPTOOLS_HOME/AmpTools
export DALITZ=$AMPTOOLS_HOME/Tutorials/Dalitz
export omegapi=$AMPTOOLS_HOME/Tutorials/omegapi


#     Fermi (CUDA 3.2 until CUDA 8) (deprecated from CUDA 9):
#         SM20 or SM_20, compute_30 – Older cards such as GeForce 400, 500, 600, GT-630
#     Kepler (CUDA 5 and later):
#         SM30 or SM_30, compute_30 – Kepler architecture (generic – Tesla K40/K80, GeForce 700, GT-730)
#         SM35 or SM_35, compute_35 – More specific Tesla K40
#         SM37 or SM_37, compute_37 – More specific Tesla K80
#     Maxwell (CUDA 6 and later):
#         SM50 or SM_50, compute_50 – Tesla/Quadro M series
#         SM52 or SM_52, compute_52 – Quadro M6000 , GeForce 900, GTX-970, GTX-980, GTX Titan X
#         SM53 or SM_53, compute_53 – Tegra (Jetson) TX1 / Tegra X1, Drive CX, Drive PX, Jetson Nano
#     Pascal (CUDA 8 and later)
#         SM60 or SM_60, compute_60 – Quadro GP100, Tesla P100, DGX-1 (Generic Pascal)
#         SM61 or SM_61, compute_61 – GTX 1080, GTX 1070, GTX 1060, GTX 1050, GTX 1030, Titan Xp, Tesla P40, Tesla P4, Discrete GPU on the NVIDIA Drive PX2
#         SM62 or SM_62, compute_62 – Integrated GPU on the NVIDIA Drive PX2, Tegra (Jetson) TX2
#     Volta (CUDA 9 and later)
#         SM70 or SM_70, compute_70 – DGX-1 with Volta, Tesla V100, GTX 1180 (GV104), Titan V, Quadro GV100
#         SM72 or SM_72, compute_72 – Jetson AGX Xavier, Drive AGX Pegasus, Xavier NX
#     Turing (CUDA 10 and later)
#         SM75 or SM_75, compute_75 – GTX/RTX Turing – GTX 1660 Ti, RTX 2060, RTX 2070, RTX 2080, Titan RTX, Quadro RTX 4000, Quadro RTX 5000, Quadro RTX 6000, Quadro RTX 8000, Quadro T1000/T2000, Tesla T4
#     Ampere (CUDA 11 and later)
#         SM80 or SM_80, compute_80 – RTX Ampere – RTX 3080
