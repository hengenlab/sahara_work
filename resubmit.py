import sahara_work as saw

efiles = []
with open('to_resubmit.txt', 'r') as f:
    for line in f:
        efiles.append(line.strip())

saw.resubmit_jobs(efiles)