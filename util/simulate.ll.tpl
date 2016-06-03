# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = simulate_${num_procs}_${Nlev}_${Nsub}
#
# @ job_type = parallel
#
# @ node = ${node} 
# @ tasks_per_node = ${tasks_per_node} 
#
# @ wall_clock_limit = ${wall_clock_limit}
#
# Groups to select from: UC, UC_merit, NZ, NZ_merit
# @ group = UC
# Your project number, either bfcs or nesi, followed by 5 digits
# @ account_no = bfcs00299
# @ output = $(job_name).$(schedd_host).$(jobid).out
# @ error = $(job_name).$(schedd_host).$(jobid).err
# To receive an email when yourjob has completed:
# @ notification       = complete
# @ notify_user        = ${email}
 
# @ class = p7linux
#
# To improve performance it is important to specify task and memory affinities 
# @ rset = rset_mcm_affinity
# @ task_affinity = core(1)
# @ network.MPI_LAPI = sn_single,shared,US,,instances=2
# @ mcm_affinity_options = mcm_mem_pref mcm_distribute
# @ environment = COPY_ALL
# @ queue
  
# suggested environment settings:
export MP_EAGER_LIMIT=65536
export MP_SHARED_MEMORY=yes
export MEMORY_AFFINITY=MCM
 
# Launch parallel job  
poe ./simulate ${Nlev} ${Nsub} 

