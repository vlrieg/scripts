#!/usr/bin/env python3
# run like:
# ./model_setup.py > run_models.sh

# two models and other parameters to combine
model = ["out-of-africa", "out-of-seasia"]
locus_length = [2000, 20000, 50000] #bp
generation_time = [0.2, 0.5, 1.0] #years
rho = [0.006, 0.00001] # Ben's estimate: 0.006 = ~ 5 * theta of 0.001179 

# set up the combinations
model_param_combos = [[m,l,g,r] for m in model for l in locus_length for g in generation_time for r in rho]

command_list = []
for mcombo in model_param_combos:
    script = './' + mcombo[0] + '.py'
    locus = mcombo[1]
    gen = mcombo[2]
    r = mcombo[3]
    command = [f"sbatch {script} --sim_no 1000 --seq_length {locus} --generation_time {gen} --recomb_rate {r} --print_sim_vcfs 1"]
    command_list.append(command)


print("#!/usr/bin/env bash")
#print("conda activate vivax")
print("source /gpfs/fs1/home/vag15/miniconda3/bin/activate /gpfs/fs1/home/vag15/miniconda3/envs/vivax")
for i in command_list:
    print(i[0])o

