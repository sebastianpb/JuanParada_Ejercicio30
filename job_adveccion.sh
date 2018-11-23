#PBS -l nodes=1:ppn=1,mem=1gb,walltime=00:10:00
#PBS -m abe
#PBS -N ejercicio30


module load anaconda/python3
cd $PBS_O_WORKDIR # este es el directorio desde donde se ejecuto qsub
gcc -fopenmp adveccion.c -o adveccion -lm
./adveccion
python plot.py
