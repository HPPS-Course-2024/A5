def compute_time_for_input_size(program_cmd, input_size, generator=GEN_NAIVE_CMD_PLACEHOLDER):
    cmd = generator.format(n=input_size)
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    time_cmd = f"bash -c 'time {program_cmd}' 2>&1"
    result = subprocess.run(time_cmd,shell=True,capture_output=True,)
    time = list(map(clean_time_str, result.stdout.decode("utf-8").splitlines()[-3:]))
    subprocess.run(CLEAN_CMD, shell=True)
    return (*time,)